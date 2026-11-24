"""Conexión con bases de datos moleculares.

* PubChem (PUG REST): nombre / SMILES / InChIKey -> CID, SMILES, IUPAC, nº CAS.
* NIST Chemistry WebBook: espectros experimentales IR, EI-MS y UV-Vis en formato JCAMP-DX.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from urllib.parse import quote

import numpy as np
import requests

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
NIST = "https://webbook.nist.gov/cgi/cbook.cgi"
HEADERS = {"User-Agent": "Spectrum-Simulator-For-Learning/1.0 (educational)"}
TIMEOUT = 20


@dataclass
class Compound:
    cid: int
    name: str
    smiles: str
    formula: str
    inchikey: str
    cas: str | None
    url: str


def _get(url, **kw):
    r = requests.get(url, headers=HEADERS, timeout=TIMEOUT, **kw)
    r.raise_for_status()
    return r


def pubchem_lookup(query: str, kind: str = "name") -> Compound | None:
    """kind: 'name', 'smiles' o 'inchikey'."""
    props = "SMILES,ConnectivitySMILES,MolecularFormula,InChIKey,IUPACName,Title"
    try:
        if kind == "smiles":
            r = requests.post(f"{PUBCHEM}/compound/smiles/property/{props}/JSON", data={"smiles": query}, headers=HEADERS, timeout=TIMEOUT)
        else:
            r = requests.get(f"{PUBCHEM}/compound/{kind}/{quote(query, safe='')}/property/{props}/JSON", headers=HEADERS, timeout=TIMEOUT)
        if r.status_code == 404:
            return None
        r.raise_for_status()
        p = r.json()["PropertyTable"]["Properties"][0]
    except (requests.RequestException, KeyError, IndexError, ValueError):
        return None
    cid = p.get("CID")
    if not cid:  # PubChem devuelve CID 0 para estructuras que no tiene
        return None
    cas = None
    try:
        syn = _get(f"{PUBCHEM}/compound/cid/{cid}/synonyms/JSON").json()
        for s in syn["InformationList"]["Information"][0].get("Synonym", []):
            if re.fullmatch(r"\d{2,7}-\d{2}-\d", s) and _cas_ok(s):
                cas = s
                break
    except Exception:
        pass
    return Compound(
        cid=cid, name=p.get("Title") or p.get("IUPACName", ""), smiles=p.get("SMILES") or p.get("ConnectivitySMILES", ""),
        formula=p.get("MolecularFormula", ""), inchikey=p.get("InChIKey", ""), cas=cas,
        url=f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
    )


def _cas_ok(cas: str) -> bool:
    digits = cas.replace("-", "")
    body, check = digits[:-1], int(digits[-1])
    return sum((i + 1) * int(d) for i, d in enumerate(reversed(body))) % 10 == check


# ---------------------------------------------------------------------------
# JCAMP-DX
# ---------------------------------------------------------------------------
@dataclass
class Spectrum:
    x: np.ndarray
    y: np.ndarray
    meta: dict
    is_peaks: bool

    @property
    def title(self):
        return self.meta.get("TITLE", "")

    def describe(self) -> str:
        bits = [self.meta.get(k) for k in ("STATE", "ORIGIN", "SOURCEREFERENCE", "$NISTSOURCE") if self.meta.get(k)]
        return " · ".join(bits)


_SQZ = {c: i for i, c in enumerate("@ABCDEFGHI")} | {c: -(i + 1) for i, c in enumerate("abcdefghi")}
_DIF = {c: i for i, c in enumerate("%JKLMNOPQR")} | {c: -(i + 1) for i, c in enumerate("jklmnopqr")}
_DUP = {c: i + 1 for i, c in enumerate("STUVWXYZs")}
_TOKEN = re.compile(r"[@A-Ia-i%J-Rj-rS-Zs]\d*\.?\d*|[+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?")


def _asdf_line(line):
    toks = _TOKEN.findall(line.replace(",", " "))
    if not toks:
        return None, [], False
    vals, last_dif, prev = [], False, None

    def num(sign_digit, rest):
        v = float(f"{abs(sign_digit)}{rest}") if rest else float(abs(sign_digit))
        return -v if sign_digit < 0 else v

    for t in toks[1:]:
        c = t[0]
        if c in _SQZ:
            prev = num(_SQZ[c], t[1:]); vals.append(prev); last_dif = False
        elif c in _DIF:
            prev = prev + num(_DIF[c], t[1:]); vals.append(prev); last_dif = True
        elif c in _DUP:
            n = int(f"{_DUP[c]}{t[1:]}")
            for _ in range(n - 1):
                if last_dif and len(vals) >= 2:
                    prev = prev + (vals[-1] - vals[-2])
                vals.append(prev)
        else:
            prev = float(t); vals.append(prev); last_dif = False
    return float(toks[0]), vals, last_dif


def parse_jcamp(text: str) -> Spectrum | None:
    meta, mode, ys, pairs, prev_dif = {}, None, [], [], False
    for raw in text.splitlines():
        line = raw.split("$$")[0].strip()
        if not line:
            continue
        m = re.match(r"^##([^=]+)=\s*(.*)$", line)
        if m:
            key = re.sub(r"[\s_\-/]", "", m.group(1)).upper()
            meta[key] = m.group(2).strip()
            mode = {"XYDATA": "xy", "PEAKTABLE": "pairs", "XYPOINTS": "pairs"}.get(key, None if key == "END" else (mode if key.startswith("$") else None))
            continue
        if mode == "pairs":
            nums = [float(v) for v in re.split(r"[\s,;]+", line) if v]
            pairs.extend(zip(nums[0::2], nums[1::2]))
        elif mode == "xy":
            x0, vals, ends_dif = _asdf_line(line)
            if x0 is None:
                continue
            if prev_dif and ys:  # en formato DIF el primer Y repite el último (control)
                ys.pop()
            ys.extend(vals)
            prev_dif = ends_dif
    if "TITLE" in meta and meta["TITLE"].lower().startswith("spectrum not found"):
        return None
    f = lambda k, d: float(meta[k]) if k in meta and re.match(r"^[-+\d.Ee]+$", meta[k]) else d
    xf, yf = f("XFACTOR", 1.0), f("YFACTOR", 1.0)
    if ys:
        first, last = f("FIRSTX", 0.0), f("LASTX", float(len(ys) - 1))
        x = np.linspace(first, last, len(ys))
        y = np.array(ys) * yf
        return Spectrum(x, y, meta, False)
    if pairs:
        arr = np.array(pairs)
        return Spectrum(arr[:, 0] * xf, arr[:, 1] * yf, meta, True)
    return None


def nist_spectrum(cas: str, kind: str, index: int = 0) -> Spectrum | None:
    """kind: 'IR', 'Mass' o 'UVVis'."""
    nid = "C" + cas.replace("-", "")
    try:
        r = requests.get(NIST, params={"JCAMP": nid, "Index": index, "Type": kind}, headers=HEADERS, timeout=TIMEOUT)
        if r.status_code != 200 or "##" not in r.text:
            return None
        return parse_jcamp(r.text)
    except requests.RequestException:
        return None


def nist_ir_all(cas: str, max_n: int = 4) -> list[Spectrum]:
    out = []
    for i in range(max_n):
        s = nist_spectrum(cas, "IR", i)
        if s is None:
            break
        out.append(s)
    return out


def nist_page(cas: str) -> str:
    return f"https://webbook.nist.gov/cgi/cbook.cgi?ID=C{cas.replace('-', '')}&Units=SI"


# --- Conversión de unidades para superponer con la simulación --------------
def ir_to_transmittance(s: Spectrum) -> tuple[np.ndarray, np.ndarray]:
    x, y = s.x.copy(), s.y.copy()
    xu = s.meta.get("XUNITS", "").upper()
    if "MICRO" in xu:
        x = 1e4 / x
    yu = s.meta.get("YUNITS", "").upper()
    if "ABSORB" in yu:
        y = 100 * 10 ** (-y / max(1e-9, y.max()) * 1.1)  # normalizado como la simulación
    elif y.max() <= 1.5:
        y = 100 * y
    return x, y


def uv_to_epsilon(s: Spectrum) -> tuple[np.ndarray, np.ndarray, str]:
    x, y = s.x.copy(), s.y.copy()
    yu = s.meta.get("YUNITS", "").upper()
    if "LOG" in yu:
        return x, 10 ** y, "ε"
    return x, y, s.meta.get("YUNITS", "absorbancia")
