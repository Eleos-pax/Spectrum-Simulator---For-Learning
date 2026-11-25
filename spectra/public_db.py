"""More public, free databases (no API key needed).

* OPSIN (EBI) and NCI/CADD CACTUS: name -> structure, also for IUPAC names that PubChem does not list.
* MassBank Europe and MoNA (MassBank of North America): experimental mass spectra (EI, MS/MS).
* Wikidata: Spanish/English names, melting and boiling points.
* PubChem PUG-View: experimental properties with their literature source.
* ChEMBL (EBI): bioactive-molecule record (useful for natural products and extracts).
* Link-outs to databases without an open API (SDBS, SpectraBase, NMRShiftDB, HMDB, ChemSpider, ChEBI,
  NP-MRD, COCONUT, CAS Common Chemistry).
Every function fails quietly (returns None / []) so the app keeps working offline.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from urllib.parse import quote

import numpy as np
import requests

HEADERS = {"User-Agent": "Spectrum-Simulator-For-Learning/1.0 (educational)"}
TIMEOUT = 20


def _json(url, **kw):
    try:
        r = requests.get(url, headers=HEADERS, timeout=TIMEOUT, **kw)
        if r.status_code != 200:
            return None
        return r.json()
    except (requests.RequestException, ValueError):
        return None


# ---------------------------------------------------------------------------
# Name -> structure
# ---------------------------------------------------------------------------
def wikidata_by_label(name: str) -> str | None:
    """SMILES of the Wikidata item whose Spanish or English label is `name` (e.g. "acetofenona", "ácido acético")."""
    variants = {name, name.lower(), name[:1].upper() + name[1:].lower()}
    values = " ".join(f'"{v}"@es "{v}"@en' for v in variants if '"' not in v and "\\" not in v)
    q = f"""SELECT ?smiles WHERE {{ VALUES ?label {{ {values} }} ?item rdfs:label ?label .
      OPTIONAL {{ ?item wdt:P2017 ?iso }} OPTIONAL {{ ?item wdt:P233 ?can }}
      BIND(COALESCE(?iso, ?can) AS ?smiles) FILTER(BOUND(?smiles)) }} LIMIT 1"""
    d = _json("https://query.wikidata.org/sparql", params={"query": q, "format": "json"})
    rows = ((d or {}).get("results") or {}).get("bindings") or []
    return rows[0]["smiles"]["value"] if rows else None


def resolve_name(name: str) -> tuple[str, str] | None:
    """(SMILES, source): OPSIN (systematic/IUPAC names), Wikidata (Spanish or English names) or NCI CACTUS."""
    name = name.strip()
    if not name:
        return None
    d = _json(f"https://www.ebi.ac.uk/opsin/ws/{quote(name, safe='')}.json")
    if d and d.get("status") == "SUCCESS" and d.get("smiles"):
        return d["smiles"], "OPSIN"
    smi = wikidata_by_label(name)
    if smi:
        return smi, "Wikidata"
    try:
        r = requests.get(f"https://cactus.nci.nih.gov/chemical/structure/{quote(name, safe='')}/smiles", headers=HEADERS, timeout=TIMEOUT)
        if r.status_code == 200 and r.text.strip() and "<" not in r.text:
            return r.text.strip().splitlines()[0], "NCI CACTUS"
    except requests.RequestException:
        pass
    return None


# ---------------------------------------------------------------------------
# Experimental mass spectra
# ---------------------------------------------------------------------------
@dataclass
class MassSpectrum:
    source: str
    id: str
    title: str
    x: np.ndarray
    y: np.ndarray
    url: str

    @property
    def label(self):
        return f"{self.source} {self.id} · {self.title}"  # the record id keeps identical titles apart


def massbank_ms(inchikey: str, max_n: int = 8) -> list[MassSpectrum]:
    d = _json("https://massbank.eu/MassBank-api/records/search", params={"inchi_key": inchikey})
    out = []
    for rec in (d or {}).get("data", [])[:max_n]:
        acc = rec.get("accession")
        r = _json(f"https://massbank.eu/MassBank-api/records/{acc}")
        if not r:
            continue
        vals = ((r.get("peak") or {}).get("peak") or {}).get("values") or []
        if not vals:
            continue
        acq = r.get("acquisition") or {}
        ms = acq.get("mass_spectrometry") or {}
        energy = next((t.get("value") for t in ms.get("subtags", []) if "ENERGY" in t.get("subtag", "")), "")
        title = " ".join(x for x in (acq.get("instrument_type", ""), ms.get("ms_type", ""), ms.get("ion_mode", "").lower(), energy) if x)
        x = np.array([v["mz"] for v in vals], float)
        y = np.array([v.get("rel") or v.get("intensity") or 0 for v in vals], float)
        out.append(MassSpectrum("MassBank", acc, title or r.get("title", acc), x, 100 * y / max(y.max(), 1e-9),
                                f"https://massbank.eu/MassBank/RecordDisplay?id={acc}"))
    return out


def mona_ms(inchikey: str, size: int = 10) -> list[MassSpectrum]:
    q = f"exists(compound.metaData.name:'InChIKey' and compound.metaData.value:'{inchikey}')"
    d = _json("https://mona.fiehnlab.ucdavis.edu/rest/spectra/search", params={"query": q, "size": size})
    out = []
    for s in d if isinstance(d, list) else []:
        md = {m.get("name"): m.get("value") for m in s.get("metaData", [])}
        try:
            pairs = [p.split(":") for p in (s.get("spectrum") or "").split()]
            x = np.array([float(a) for a, _ in pairs])
            y = np.array([float(b) for _, b in pairs])
        except ValueError:
            continue
        if not len(x):
            continue
        title = " ".join(str(v) for v in (md.get("instrument type"), md.get("ms level"), md.get("ionization mode"),
                                          md.get("ionization energy") or md.get("collision energy")) if v)
        out.append(MassSpectrum("MoNA", str(s.get("id")), title or "MS", x, 100 * y / max(y.max(), 1e-9),
                                f"https://mona.fiehnlab.ucdavis.edu/spectra/display/{s.get('id')}"))
    return out


# ---------------------------------------------------------------------------
# Names and properties
# ---------------------------------------------------------------------------
_UNITS = {"Q25267": "°C", "Q11579": "K", "Q42289": "°F"}


def _to_c(v: float, unit: str) -> float:
    return v - 273.15 if unit == "K" else (v - 32) * 5 / 9 if unit == "°F" else v


def wikidata(inchikey: str) -> dict | None:
    q = f"""SELECT ?item ?es ?en ?mp ?mpu ?bp ?bpu WHERE {{
      ?item wdt:P235 "{inchikey}" .
      OPTIONAL {{ ?item rdfs:label ?es FILTER(lang(?es)='es') }}
      OPTIONAL {{ ?item rdfs:label ?en FILTER(lang(?en)='en') }}
      OPTIONAL {{ ?item p:P2101/psv:P2101 [ wikibase:quantityAmount ?mp ; wikibase:quantityUnit ?mpu ] }}
      OPTIONAL {{ ?item p:P2102/psv:P2102 [ wikibase:quantityAmount ?bp ; wikibase:quantityUnit ?bpu ] }}
    }} LIMIT 5"""
    d = _json("https://query.wikidata.org/sparql", params={"query": q, "format": "json"})
    rows = ((d or {}).get("results") or {}).get("bindings") or []
    if not rows:
        return None
    r0 = rows[0]
    get = lambda r, k: r[k]["value"] if k in r else None
    out = {"qid": get(r0, "item").rsplit("/", 1)[-1], "es": get(r0, "es"), "en": get(r0, "en"), "mp": None, "bp": None}
    for r in rows:
        for k in ("mp", "bp"):
            if out[k] is None and get(r, k):
                unit = _UNITS.get((get(r, k + "u") or "").rsplit("/", 1)[-1], "")
                if unit:
                    out[k] = round(_to_c(float(get(r, k)), unit), 1)
    out["url"] = f"https://www.wikidata.org/wiki/{out['qid']}"
    return out


PUBCHEM_PROPS = ["Physical Description", "Color/Form", "Odor", "Melting Point", "Boiling Point", "Flash Point", "Solubility",
                 "Density", "Vapor Pressure", "LogP", "Refractive Index", "Dissociation Constants", "Kovats Retention Index"]


def pubchem_experimental(cid: int) -> list[tuple[str, str, str]]:
    """[(property, value, source)] from PubChem's 'Experimental Properties' section."""
    d = _json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON", params={"heading": "Experimental Properties"})
    found = {}

    def walk(sec):
        for s in sec.get("Section", []):
            h = s.get("TOCHeading")
            if h in PUBCHEM_PROPS and h not in found:
                infos = s.get("Information", [])
                # prefer SI / °C values when the property is given in several units
                infos = sorted(infos, key=lambda i: "°C" not in str(i.get("Value", {})))
                for info in infos:
                    v = info.get("Value", {})
                    sw = v.get("StringWithMarkup")
                    text = sw[0].get("String", "") if sw else (f"{v['Number'][0]} {v.get('Unit', '')}".strip() if v.get("Number") else "")
                    if text:
                        if h == "Kovats Retention Index":
                            nums = re.findall(r"\d+(?:\.\d+)?", text)[:6]
                            text = ", ".join(nums) + (" …" if len(nums) == 6 else "")
                        ref = (info.get("Reference") or [""])[0]
                        found[h] = (h, text[:220], ref[:160])
                        break
            walk(s)
    if d:
        walk(d.get("Record", {}))
    return [found[h] for h in PUBCHEM_PROPS if h in found]


def chembl(inchikey: str) -> dict | None:
    d = _json("https://www.ebi.ac.uk/chembl/api/data/molecule", params={"molecule_structures__standard_inchi_key": inchikey, "format": "json"})
    mols = (d or {}).get("molecules") or []
    if not mols:
        return None
    m = mols[0]
    cid = m.get("molecule_chembl_id")
    return {"id": cid, "name": m.get("pref_name"), "max_phase": m.get("max_phase"), "type": m.get("molecule_type"),
            "url": f"https://www.ebi.ac.uk/chembl/explore/compound/{cid}"}


def external_links(inchikey: str, name: str = "", cas: str | None = None, cid: int | None = None) -> list[tuple[str, str, str]]:
    """[(database, url, what it offers)] for databases that have no open API (opened in the browser)."""
    k, n = quote(inchikey or ""), quote(name or inchikey or "")
    links = [
        ("SDBS (AIST)", "https://sdbs.db.aist.go.jp/", "IR, ¹H/¹³C NMR, MS, Raman"),
        ("SpectraBase (Wiley)", f"https://spectrabase.com/search?q={k}", "IR, NMR, MS, UV-Vis, Raman"),
        ("NMRShiftDB2", f"https://nmrshiftdb.nmr.uni-koeln.de/portal/js_pane/P-Results?nmrshiftdbaction=searchbyname&name={n}", "NMR"),
        ("HMDB", f"https://hmdb.ca/unearth/q?query={k}&searcher=metabolites", "NMR, MS (metabolites)"),
        ("NP-MRD", f"https://np-mrd.org/unearth/q?query={k}&searcher=natural_products", "NMR (natural products)"),
        ("COCONUT", f"https://coconut.naturalproducts.net/search?q={k}", "natural products"),
        ("ChemSpider", f"https://www.chemspider.com/Search.aspx?q={k}", "names, properties"),
        ("ChEBI", f"https://www.ebi.ac.uk/chebi/advancedSearchFT.do?searchString={k}", "biological role"),
        ("MassBank Europe", f"https://massbank.eu/MassBank/Result.jsp?inchikey={k}", "MS"),
        ("MoNA", "https://mona.fiehnlab.ucdavis.edu/spectra/search", "MS"),
    ]
    if cas:
        links.append(("CAS Common Chemistry", f"https://commonchemistry.cas.org/detail?cas_rn={quote(cas)}", "CAS data"))
    if cid:
        links.insert(0, ("PubChem", f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Spectral-Information", "spectral information"))
    return links
