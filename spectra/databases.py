"""Connection to molecular databases.

* PubChem (PUG REST): name / SMILES / InChIKey -> CID, SMILES, IUPAC name, CAS number.
* NIST Chemistry WebBook: experimental IR, EI-MS and UV-Vis spectra in JCAMP-DX format.
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
    iupac: str = ""
    synonyms: tuple = ()
    props: tuple = ()  # (name, value) pairs computed by PubChem


def _get(url, **kw):
    r = requests.get(url, headers=HEADERS, timeout=TIMEOUT, **kw)
    r.raise_for_status()
    return r


def pubchem_lookup(query: str, kind: str = "name") -> Compound | None:
    """kind: 'name', 'smiles' or 'inchikey'."""
    props = ("SMILES,ConnectivitySMILES,MolecularFormula,InChIKey,IUPACName,Title,XLogP,TPSA,Complexity,Charge,"
             "HBondDonorCount,HBondAcceptorCount,RotatableBondCount,ExactMass,MolecularWeight")
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
    if not cid:  # PubChem returns CID 0 for structures it does not hold
        return None
    cas, synonyms = None, []
    try:
        syn = _get(f"{PUBCHEM}/compound/cid/{cid}/synonyms/JSON").json()
        for s in syn["InformationList"]["Information"][0].get("Synonym", []):
            if re.fullmatch(r"\d{2,7}-\d{2}-\d", s):
                if cas is None and _cas_ok(s):
                    cas = s
            elif len(synonyms) < 15 and not re.search(r"^(UNII|DTXSID|CHEMBL|NSC|EINECS|SCHEMBL|AKOS|MFCD|ZINC|HSDB|BRN|CCRIS|AI3|WLN|InChI)", s, re.I):
                synonyms.append(s)
    except Exception:
        pass
    keys = [("XLogP", "XLogP3"), ("TPSA", "TPSA (Å²)"), ("Complexity", "Complexity"), ("Charge", "Charge"),
            ("HBondDonorCount", "H-bond donors"), ("HBondAcceptorCount", "H-bond acceptors"), ("RotatableBondCount", "Rotatable bonds")]
    extra = tuple((label, p[k]) for k, label in keys if k in p)
    return Compound(
        cid=cid, name=p.get("Title") or p.get("IUPACName", ""), smiles=p.get("SMILES") or p.get("ConnectivitySMILES", ""),
        formula=p.get("MolecularFormula", ""), inchikey=p.get("InChIKey", ""), cas=cas,
        url=f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
        iupac=p.get("IUPACName", ""), synonyms=tuple(synonyms), props=extra,
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
            if prev_dif and ys:  # in DIF form the first Y repeats the previous line's last Y (check value)
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
    """kind: 'IR', 'Mass' or 'UVVis'."""
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


# --- Unit conversion so experimental data can be overlaid on the simulation ---
def ir_to_transmittance(s: Spectrum) -> tuple[np.ndarray, np.ndarray]:
    x, y = s.x.copy(), s.y.copy()
    xu = s.meta.get("XUNITS", "").upper()
    if "MICRO" in xu:
        x = 1e4 / x
    yu = s.meta.get("YUNITS", "").upper()
    if "ABSORB" in yu:
        y = 100 * 10 ** (-y / max(1e-9, y.max()) * 1.1)  # normalised like the simulation
    elif y.max() <= 1.5:
        y = 100 * y
    return x, y


def uv_to_epsilon(s: Spectrum) -> tuple[np.ndarray, np.ndarray, str]:
    x, y = s.x.copy(), s.y.copy()
    yu = s.meta.get("YUNITS", "").upper()
    if "LOG" in yu:
        return x, 10 ** y, "ε"
    return x, y, s.meta.get("YUNITS", "absorbancia")


# ---------------------------------------------------------------------------
# User files exported by instruments (CSV / TXT / JCAMP-DX)
# ---------------------------------------------------------------------------
@dataclass
class UserTable:
    name: str
    columns: list[str]
    data: np.ndarray  # rows x columns (float)


_NUM = re.compile(r"^[+-]?(\d+([.,]\d*)?|[.,]\d+)([Ee][+-]?\d+)?$")


def _decode(raw: bytes) -> str:
    for enc in ("utf-8-sig", "utf-16", "latin-1"):
        try:
            text = raw.decode(enc)
            if enc == "utf-16" and "\x00" in text:
                continue
            return text
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="ignore")


def read_user_file(raw: bytes, name: str) -> UserTable:
    """Read a spectrum exported by an instrument.

    Accepts JCAMP-DX (.dx/.jdx/.jcamp) and delimited text (.csv/.txt/.tsv/.asc/.prn) with any of
    ``,`` ``;`` tab or spaces as separator, decimal point or decimal comma, and any number of header
    or metadata lines (they are skipped). Several Y columns (e.g. several samples) are kept.
    Raises ValueError when no numeric table is found.
    """
    text = _decode(raw)
    if name.lower().endswith((".dx", ".jdx", ".jcamp")) or text.lstrip().startswith("##"):
        sp = parse_jcamp(text)
        if sp is None or not len(sp.x):
            raise ValueError("empty JCAMP-DX file")
        ylab = sp.meta.get("YUNITS", "Y")
        return UserTable(name, [sp.meta.get("XUNITS", "X"), ylab], np.column_stack([sp.x, sp.y]))

    lines = [l.strip() for l in text.splitlines() if l.strip()]
    best = None
    for sep in (";", "\t", ",", None):
        rows, header, last_text = [], None, None
        for line in lines:
            parts = [p.strip().strip('"') for p in (line.split(sep) if sep else line.split())]
            parts = [p for p in parts if p != ""] if sep is None else parts
            numeric = [p for p in parts if _NUM.match(p)]
            filled = [p for p in parts if p]
            if len(numeric) >= 2 and len(numeric) == len(filled):
                if not rows:
                    header = last_text
                rows.append(parts)
            else:
                if rows and len(rows) > 5:
                    break  # metadata after the table
                rows, last_text = ([], parts) if not rows else ([], parts)
        if len(rows) >= 3 and (best is None or len(rows[0]) * len(rows) > best[0]):
            best = (len(rows[0]) * len(rows), sep, rows, header)
    if best is None:
        raise ValueError("no numeric table found")
    _, sep, rows, header = best
    decimal_comma = sep != "," and any("," in p for r in rows[:20] for p in r)
    ncol = max(len(r) for r in rows)

    def num(p):
        if not _NUM.match(p or ""):
            return np.nan
        return float(p.replace(",", ".")) if decimal_comma or ("," in p and sep != ",") else float(p)

    data = np.array([[num(r[i]) if i < len(r) else np.nan for i in range(ncol)] for r in rows], dtype=float)
    keep = ~np.all(np.isnan(data), axis=0)
    data = data[:, keep]
    if header and len(header) >= ncol:
        cols = [h or f"col {i + 1}" for h, k in zip(header, keep) if k]
    else:
        cols = [f"col {i + 1}" for i in range(data.shape[1])]
    data = data[~np.isnan(data[:, 0])]
    return UserTable(name, cols, data)
