"""Educational simulation of IR, UV-Vis and EI-MS spectra from a structure (RDKit).

The models are empirical:
  * IR  - group-frequency table (SMARTS) + C-H bands by hybridisation, sum of Lorentzians.
  * UV  - conjugated systems (number of pi bonds, fused aromatic rings, auxochromes),
          simplified Woodward-Fieser-type rules, Gaussian bands on an energy scale.
  * MS  - 70 eV electron ionisation: molecular ion with exact isotope pattern and fragments
          from stability rules (alpha-cleavage, benzylic/tropylium, McLafferty, neutral losses).
They are meant for teaching spectral interpretation; they do not replace quantum-chemical
calculations or spectral libraries. Peak labels are generated in Spanish and translated by
`spectra.i18n.label`.
"""
from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# ---------------------------------------------------------------------------
# Isotopes: symbol -> [(nominal mass, exact mass, abundance)]
# ---------------------------------------------------------------------------
ISOTOPES = {
    "H": [(1, 1.007825, 0.99988), (2, 2.014102, 0.00012)],
    "B": [(10, 10.012937, 0.199), (11, 11.009305, 0.801)],
    "C": [(12, 12.0, 0.9893), (13, 13.003355, 0.0107)],
    "N": [(14, 14.003074, 0.99636), (15, 15.000109, 0.00364)],
    "O": [(16, 15.994915, 0.99757), (17, 16.999132, 0.00038), (18, 17.99916, 0.00205)],
    "F": [(19, 18.998403, 1.0)],
    "Si": [(28, 27.976927, 0.92223), (29, 28.976495, 0.04685), (30, 29.97377, 0.03092)],
    "P": [(31, 30.973762, 1.0)],
    "S": [(32, 31.972071, 0.9499), (33, 32.971459, 0.0075), (34, 33.967867, 0.0425)],
    "Cl": [(35, 34.968853, 0.7576), (37, 36.965903, 0.2424)],
    "Br": [(79, 78.918338, 0.5069), (81, 80.916291, 0.4931)],
    "I": [(127, 126.904473, 1.0)],
}
NOMINAL = {s: v[0][0] for s, v in ISOTOPES.items()}


def formula_counts(mol: Chem.Mol, atoms=None) -> Counter:
    c = Counter()
    for a in mol.GetAtoms():
        if atoms is not None and a.GetIdx() not in atoms:
            continue
        if a.GetSymbol() == "H":
            c["H"] += 1
            continue
        c[a.GetSymbol()] += 1
        c["H"] += a.GetTotalNumHs()
    return c


def nominal_mass(c: Counter) -> int:
    return int(sum(NOMINAL.get(s, 0) * n for s, n in c.items()))


def hill(c: Counter, html=False) -> str:
    keys = sorted(c, key=lambda s: (0 if s == "C" else 1 if s == "H" else 2, s))
    sub = (lambda n: f"<sub>{n}</sub>") if html else (lambda n: str(n))
    return "".join(s + (sub(c[s]) if c[s] > 1 else "") for s in keys if c[s])


def isotope_pattern(c: Counter, max_peaks: int = 6) -> list[tuple[int, float]]:
    dist = {0: 1.0}
    for s, n in c.items():
        iso = ISOTOPES.get(s)
        if not iso or len(iso) == 1:
            continue
        base = iso[0][0]
        for _ in range(n):
            new = defaultdict(float)
            for m, p in dist.items():
                for mm, _, ab in iso:
                    new[m + mm - base] += p * ab
            dist = {k: v for k, v in sorted(new.items()) if v > 1e-7}
            dist = dict(list(dist.items())[:16])
    M = nominal_mass(c)
    top = max(dist.values())
    return [(M + k, 100 * v / top) for k, v in list(dist.items())[:max_peaks] if 100 * v / top > 0.05]


# ---------------------------------------------------------------------------
# IR
# ---------------------------------------------------------------------------
# (SMARTS, group, [(centre cm-1, half width, intensity 0-1, assignment)])
IR_RULES = [
    ("[OX2H][CX4,c]", "Alcohol / fenol", [(3350, 120, .85, "ν O–H (puente de H)"), (1060, 25, .75, "ν C–O")]),
    ("[CX3](=O)[OX2H1]", "Ácido carboxílico", [(3000, 300, .75, "ν O–H ácido (muy ancha)"), (1710, 18, 1, "ν C=O ácido"), (1290, 25, .6, "ν C–O"), (930, 30, .35, "δ O–H f.p.")]),
    ("[CX3;!R](=O)[OX2][#6]", "Éster", [(1740, 14, 1, "ν C=O éster"), (1240, 25, .85, "ν C–O–C asim."), (1100, 25, .6, "ν C–O–C sim.")]),
    ("[CX3;R](=O)[OX2;R][#6]", "Lactona", [(1765, 15, 1, "ν C=O lactona"), (1180, 25, .8, "ν C–O")]),
    ("[CX3H1](=O)[#6]", "Aldehído", [(1725, 14, 1, "ν C=O aldehído"), (2820, 15, .35, "ν C–H aldehído (Fermi)"), (2720, 15, .35, "ν C–H aldehído (Fermi)")]),
    ("[#6][CX3](=O)[#6]", "Cetona", [(1715, 14, 1, "ν C=O cetona")]),
    ("[CX3](=O)[NX3H2]", "Amida primaria", [(3350, 40, .55, "ν N–H asim."), (3180, 40, .5, "ν N–H sim."), (1665, 18, 1, "ν C=O (amida I)"), (1620, 15, .6, "δ N–H (amida II)")]),
    ("[CX3](=O)[NX3H1]", "Amida secundaria", [(3300, 40, .55, "ν N–H"), (1650, 18, 1, "ν C=O (amida I)"), (1550, 15, .7, "amida II")]),
    ("[CX3](=O)[NX3H0]", "Amida terciaria", [(1645, 18, 1, "ν C=O amida")]),
    ("[CX3](=O)[Cl,Br]", "Haluro de acilo", [(1800, 12, 1, "ν C=O haluro de acilo")]),
    ("[CX3](=O)[OX2][CX3](=O)", "Anhídrido", [(1820, 12, .9, "ν C=O asim."), (1760, 12, .8, "ν C=O sim.")]),
    ("[NX3;H2;!$(NC=O)][#6]", "Amina primaria", [(3420, 30, .35, "ν N–H asim."), (3340, 30, .3, "ν N–H sim."), (1610, 20, .4, "δ N–H tijera"), (1150, 30, .3, "ν C–N")]),
    ("[NX3;H1;!$(NC=O)]([#6])[#6]", "Amina secundaria", [(3330, 30, .3, "ν N–H"), (1130, 30, .3, "ν C–N")]),
    ("[CX2]#[NX1]", "Nitrilo", [(2250, 10, .6, "ν C≡N")]),
    ("[CX2H1]#[CX2]", "Alquino terminal", [(3300, 15, .7, "ν ≡C–H"), (2120, 10, .3, "ν C≡C"), (630, 25, .6, "δ ≡C–H")]),
    ("[CX2H0;!$(C#N)]#[CX2H0]", "Alquino interno", [(2220, 10, .08, "ν C≡C (débil)")]),
    ("[CX3;!a]=[CX3;!a]", "Alqueno", [(1645, 12, .35, "ν C=C"), (990, 15, .5, "δ =C–H f.p.")]),
    ("[CX3H2]=[CX3]", "Vinilo terminal", [(910, 12, .65, "δ =CH₂ f.p.")]),
    ("c1ccccc1", "Anillo aromático", [(1600, 10, .45, "ν C=C aromático"), (1500, 10, .5, "ν C=C aromático"), (1450, 10, .4, "ν C=C aromático"), (2000, 60, .06, "sobretonos aromáticos"), (1850, 60, .05, "sobretonos aromáticos")]),
    ("[cH1]1[cH1][cH1][cH1][cH1]c1", "Benceno monosustituido", [(750, 12, .85, "δ C–H f.p. (monosust.)"), (695, 12, .85, "deformación del anillo (monosust.)")]),
    ("[$([NX3](=O)=O),$([NX3+](=O)[O-])]", "Nitro", [(1530, 15, 1, "ν NO₂ asim."), (1350, 15, .85, "ν NO₂ sim.")]),
    ("[OD2;!$(OC=O)]([#6])[#6]", "Éter", [(1110, 30, .85, "ν C–O–C")]),
    ("[#6]F", "Fluoruro", [(1100, 40, .9, "ν C–F")]),
    ("[#6]Cl", "Cloruro", [(740, 30, .7, "ν C–Cl")]),
    ("[#6]Br", "Bromuro", [(600, 30, .65, "ν C–Br")]),
    ("[#6]I", "Yoduro", [(520, 30, .6, "ν C–I")]),
    ("[SX2H]", "Tiol", [(2570, 15, .25, "ν S–H")]),
    ("[SX4](=O)(=O)", "Sulfona / sulfonamida", [(1330, 20, .9, "ν SO₂ asim."), (1150, 20, .85, "ν SO₂ sim.")]),
    ("[SX3]=O", "Sulfóxido", [(1050, 20, .9, "ν S=O")]),
    ("N=C=O", "Isocianato", [(2270, 20, 1, "ν N=C=O")]),
]
_PATT = {}


def _q(smarts):
    if smarts not in _PATT:
        _PATT[smarts] = Chem.MolFromSmarts(smarts)
    return _PATT[smarts]


def matches(mol, smarts):
    p = _q(smarts)
    return mol.GetSubstructMatches(p) if p is not None else ()


@dataclass
class IRResult:
    x: np.ndarray
    y: np.ndarray  # %T
    peaks: list = field(default_factory=list)
    groups: list = field(default_factory=list)


def simulate_ir(mol: Chem.Mol) -> IRResult:
    bands, groups = [], []
    add = lambda c, w, i, lab: bands.append((c, w, i, lab))

    sp3 = sp2 = sp = 0
    nCH3 = nCH2 = 0
    for a in mol.GetAtoms():
        h = a.GetTotalNumHs()
        if a.GetSymbol() != "C" or not h:
            continue
        hyb = a.GetHybridization()
        if a.GetIsAromatic() or hyb == Chem.HybridizationType.SP2:
            sp2 += h
        elif hyb == Chem.HybridizationType.SP:
            sp += h
        else:
            sp3 += h
            nCH3 += h == 3
            nCH2 += h == 2
    tot = max(1, sp3 + sp2 + sp)
    if sp3:
        s = min(1, .35 + .6 * sp3 / tot)
        add(2960, 12, s * (.9 if nCH3 else .4), "ν C–H sp³ asim.")
        add(2925, 14, s * (1 if nCH2 else .5), "ν C–H sp³")
        add(2870, 12, s * .6, "ν C–H sp³ sim.")
        if nCH2:
            add(1465, 12, .45, "δ CH₂ tijera")
        if nCH3:
            add(1375, 10, .4, "δ CH₃ paraguas")
        if nCH2 >= 4:
            add(722, 8, .3, "balanceo (CH₂)ₙ, n ≥ 4")
        groups.append("C–H sp³")
    if sp2:
        add(3050, 18, .25 + .2 * sp2 / tot, "ν C–H sp²")
        groups.append("C–H sp²")
    if any(a.GetIsAromatic() and a.GetTotalNumHs() for a in mol.GetAtoms()) and not matches(mol, "[cH1]1[cH1][cH1][cH1][cH1]c1"):
        add(800, 20, .6, "δ C–H aromático f.p.")

    conj_co = bool(matches(mol, "[CX3]=[CX3][CX3]=O") or matches(mol, "a[CX3]=O"))
    for smarts, name, bs in IR_RULES:
        hits = matches(mol, smarts)
        if not hits:
            continue
        groups.append(name + (f" ×{len(hits)}" if len(hits) > 1 else ""))
        for c, w, i, lab in bs:
            cc = c - 25 if conj_co and lab.startswith("ν C=O") and "amida" not in lab else c
            add(cc, w, i * min(1.25, 1 + .1 * (len(hits) - 1)), lab + (" (conjugado)" if cc != c else ""))

    # Fingerprint region: deterministic skeletal bands (seed = canonical SMILES)
    seed = int(hashlib.md5(Chem.MolToSmiles(mol).encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    for _ in range(min(10, 2 + mol.GetNumHeavyAtoms() // 2)):
        add(650 + rng.random() * 700, 6 + rng.random() * 10, .08 + rng.random() * .25, "esqueleto (huella)")

    x = np.arange(4000, 399, -2.0)
    A = np.zeros_like(x)
    for c, w, i, _ in bands:
        A += i / (1 + ((x - c) / w) ** 2)
    y = 100 * 10 ** (-1.1 * A)
    peaks = sorted([b for b in bands if b[2] >= .2 and "huella" not in b[3]], key=lambda b: -b[0])
    return IRResult(x, y, peaks, groups)


# ---------------------------------------------------------------------------
# UV-Vis
# ---------------------------------------------------------------------------
@dataclass
class UVResult:
    x: np.ndarray  # nm
    y: np.ndarray  # epsilon (L mol-1 cm-1)
    bands: list
    lambda_max: float
    color: str


def _conjugated_systems(mol):
    """Connected components of atoms carrying pi bonds (aromatic atoms included)."""
    pi = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic() or any(b.GetBondTypeAsDouble() >= 2 for b in a.GetBonds())}
    seen, systems = set(), []
    ri = mol.GetRingInfo()
    for start in pi:
        if start in seen:
            continue
        comp, stack = set(), [start]
        seen.add(start)
        while stack:
            i = stack.pop()
            comp.add(i)
            for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                j = nb.GetIdx()
                if j in pi and j not in seen:
                    seen.add(j)
                    stack.append(j)
        sub = [b for b in mol.GetBonds() if b.GetBeginAtomIdx() in comp and b.GetEndAtomIdx() in comp]
        arom_bonds = [b for b in sub if b.GetIsAromatic()]
        npi = sum(1 for b in sub if b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)) + len(arom_bonds) / 2
        pairs = [frozenset((b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol())) for b in sub if b.GetBondType() == Chem.BondType.DOUBLE]
        has_co = frozenset(("C", "O")) in pairs
        has_no2 = frozenset(("N", "O")) in pairs or any(mol.GetAtomWithIdx(i).GetFormalCharge() == 1 and mol.GetAtomWithIdx(i).GetSymbol() == "N" for i in comp)
        azo = frozenset(("N",)) in pairs
        aro_rings = sum(1 for r in ri.AtomRings() if all(i in comp and mol.GetAtomWithIdx(i).GetIsAromatic() for i in r))
        aux = 0.0
        for i in comp:
            for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                if nb.GetIdx() in comp:
                    continue
                s = nb.GetSymbol()
                if s == "N":
                    aux += 1.0
                elif s in ("O", "S"):
                    aux += .6
                elif s in ("Cl", "Br", "I"):
                    aux += .3
                elif s == "C":
                    aux += .08
        systems.append(dict(npi=npi, has_co=has_co, has_no2=has_no2, azo=azo, aro_rings=aro_rings, aux=aux, size=len(comp)))
    return systems


def _wl_rgb(l):
    if l < 440: return (-(l - 440) / 60, 0, 1)
    if l < 490: return (0, (l - 440) / 50, 1)
    if l < 510: return (0, 1, -(l - 510) / 20)
    if l < 580: return ((l - 510) / 70, 1, 0)
    if l < 645: return (1, -(l - 645) / 65, 0)
    return (1, 0, 0)


def simulate_uv(mol: Chem.Mol) -> UVResult:
    bands = []
    polyene = lambda n: 171 + 46 * max(0.0, n - 1) ** .85
    for s in _conjugated_systems(mol):
        n_eff = s["npi"] + s["aux"]
        if s["aro_rings"]:
            f = min(s["aro_rings"], 4)
            extra = max(0.0, s["npi"] - (3 * f - (f - 1)))  # pi bonds outside the fused rings
            lB = [254, 286, 375, 440][f - 1] + 22 * extra + 12 * s["aux"]
            lE = [204, 221, 252, 275][f - 1] + 25 * extra + 8 * s["aux"]
            epsB = (250 if f == 1 else 5000) * (1 + 3 * min(1.0, extra + s["aux"] / 2))
            bands.append((lE, 8000 * f * (1 + .5 * extra), "π→π* (banda E₂/K)"))
            bands.append((lB, epsB, "π→π* banda B (estructura fina)" if f == 1 else "π→π* banda p/B"))
            if f == 1 and extra + s["aux"] < .5:
                bands.append((184, 60000, "π→π* (banda E₁)"))
        elif s["npi"] >= 1:
            if s["has_co"] and s["npi"] == 1:
                bands.append((188, 900, "π→π* C=O"))
            else:
                n = int(round(s["npi"]))
                bands.append((polyene(n_eff), min(150000, 10000 * n ** 1.3), f"π→π* ({n} enlaces π conjugados)" if n > 1 else "π→π* aislado"))
        if s["has_co"]:
            bands.append((280 + 20 * max(0, s["npi"] - 1 - 2 * s["aro_rings"]), 30, "n→π* C=O (prohibida, débil)"))
        if s["has_no2"]:
            bands.append((270 if s["aro_rings"] else 275, 8000 if s["aro_rings"] else 20, "transferencia de carga Ar–NO₂" if s["aro_rings"] else "n→π* NO₂"))
        if s["azo"]:
            bands.append((440 if s["aro_rings"] else 350, 500 if s["aro_rings"] else 15, "n→π* N=N (azo)"))
    if not bands:
        bands.append((180, 300, "σ→σ* / n→σ* (< 200 nm, no observable en UV convencional)"))

    x = np.arange(180, 801, 1.0)
    nu = 1e7 / x
    y = np.zeros_like(x)
    for l, e, _ in bands:
        y += e * np.exp(-.5 * ((nu - 1e7 / l) / 2600) ** 2)
    candidates = [b for b in bands if b[0] >= 200]
    lmax = max(candidates, key=lambda b: b[1])[0] if candidates else bands[0][0]
    # perceived colour (transmitted white light)
    vis = (x >= 400) & (x <= 700)
    T = 10 ** (-2 * y[vis] / max(y.max(), 3000))
    base = np.array([_wl_rgb(l) for l in x[vis]]).T
    rgb = (base @ T) / base.sum(axis=1)
    color = "#%02x%02x%02x" % tuple(int(255 * min(1, v)) for v in rgb)
    return UVResult(x, y, sorted(bands, key=lambda b: b[0]), lmax, color)


# ---------------------------------------------------------------------------
# EI-MS
# ---------------------------------------------------------------------------
@dataclass
class MSResult:
    peaks: list  # (m/z, relative intensity, label)
    M: int
    molecular_ion: list


def _side(mol, start, cut_bond_idx):
    seen, stack = {start}, [start]
    while stack:
        a = mol.GetAtomWithIdx(stack.pop())
        for b in a.GetBonds():
            if b.GetIdx() == cut_bond_idx:
                continue
            j = b.GetOtherAtomIdx(a.GetIdx())
            if j not in seen:
                seen.add(j)
                stack.append(j)
    return seen


def _is_carbonyl_c(a):
    return a.GetSymbol() == "C" and any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetSymbol() in ("O", "S") for b in a.GetBonds())


def simulate_ms(mol: Chem.Mol) -> MSResult:
    cnt = formula_counts(mol)
    M = nominal_mass(cnt)
    ions = defaultdict(lambda: {"score": 0.0, "labels": [], "cnt": None})

    def add(mz, score, label, c=None):
        if mz < 12 or mz > M:
            return
        ion = ions[mz]
        ion["score"] += score
        if label:
            ion["labels"].append((score, label))
        if c is not None and ion["cnt"] is None:
            ion["cnt"] = c

    n_arom = sum(a.GetIsAromatic() for a in mol.GetAtoms())
    ring_frac = sum(b.IsInRing() for b in mol.GetBonds()) / max(1, mol.GetNumBonds())
    m_stab = .25 + .6 * (n_arom > 0) + .3 * ring_frac
    if matches(mol, "[OX2H][CX4]"):
        m_stab *= .25
    if matches(mol, "[CX4H0]([#6])([#6])([#6])[#6]"):
        m_stab *= .4
    if mol.GetNumHeavyAtoms() <= 4:
        m_stab = 1.0

    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE or b.IsInRing():
            continue
        A, B = b.GetBeginAtom(), b.GetEndAtom()
        for X, Y in ((A, B), (B, A)):  # the charge stays on X's fragment
            side = _side(mol, X.GetIdx(), b.GetIdx())
            fc = formula_counts(mol, side)
            mz = nominal_mass(fc)
            lost = M - mz
            s, lab = .15, ""
            if _is_carbonyl_c(X):
                s, lab = s + 2.2, "ion acilio (ruptura α)"
            elif X.GetSymbol() == "C" and Y.GetSymbol() == "C" and any(n.GetSymbol() in "NOS" and n.GetIdx() != Y.GetIdx() and not _is_carbonyl_c(n) for n in X.GetNeighbors()):
                het = next(n for n in X.GetNeighbors() if n.GetSymbol() in "NOS" and n.GetIdx() != Y.GetIdx())
                s, lab = (s + 3.0, "ion iminio (ruptura α a N)") if het.GetSymbol() == "N" else (s + 1.6, "ion oxonio (ruptura α a O)")
            elif X.GetSymbol() == "C" and not X.GetIsAromatic() and any(n.GetIsAromatic() for n in X.GetNeighbors()) and not Y.GetIsAromatic():
                s, lab = s + 2.5, "catión bencílico / tropilio"
            elif Y.GetSymbol() in ("Cl", "Br", "I") and len(side) > 1:
                s, lab = s + {"I": 2.0, "Br": 1.5, "Cl": .6}[Y.GetSymbol()], f"pérdida de {Y.GetSymbol()}·"
            elif X.GetSymbol() == "C":
                deg = sum(1 for n in X.GetNeighbors() if n.GetSymbol() == "C" and n.GetIdx() != Y.GetIdx())
                s += [.05, .3, .9, 1.5][min(deg, 3)]
                allyl = any(n.GetIdx() != Y.GetIdx() and not n.GetIsAromatic() and any(bb.GetBondType() == Chem.BondType.DOUBLE and bb.GetOtherAtom(n).GetSymbol() == "C" for bb in n.GetBonds()) for n in X.GetNeighbors())
                if allyl:
                    s, lab = s + 1.0, "catión alílico"
                else:
                    lab = ["", "carbocatión primario", "carbocatión secundario", "carbocatión terciario"][min(deg, 3)]
            s *= .6 + .4 * min(1, lost / 43)  # Stevenson's rule: the larger radical is lost
            if len(side) == 1 and X.GetSymbol() == "C" and X.GetTotalNumHs() == 3:
                s, lab = s * .25, "CH₃⁺"
            add(mz, s, f"{lab} (M − {lost})" if lab else f"M − {lost}", fc)

    # McLafferty
    for m in matches(mol, "[CX3](=[O,S])[CX4][CX4][CX4;!H0]"):
        c, _, a, bt = m[:4]
        bond = mol.GetBondBetweenAtoms(a, bt)
        if bond is None or bond.IsInRing():
            continue
        fc = formula_counts(mol, _side(mol, c, bond.GetIdx()))
        fc["H"] += 1
        add(nominal_mass(fc), 2.4, "reordenamiento de McLafferty", fc)

    for smarts, d, s, lab in [
        ("[OX2H][CX4][CX4;!H0]", 18, 1.6, "M − H₂O (deshidratación)"),
        ("[CX3](=O)[OX2H]", 17, .8, "M − ·OH"),
        ("[CX3](=O)[OX2H]", 45, 1.0, "M − ·COOH"),
        ("[CX3](=O)[OX2][CH3]", 31, 1.2, "M − ·OCH₃"),
        ("[CX3](=O)[OX2][CH2][CH3]", 45, 1.2, "M − ·OC₂H₅"),
        ("[CX3H1]=O", 1, 1.0, "M − H· (aldehído)"),
        ("[CX3H1]=O", 29, .9, "M − ·CHO"),
        ("c[OH]", 28, 1.2, "M − CO (fenol)"),
        ("c[$([N+](=O)[O-]),$(N(=O)=O)]", 46, 1.5, "M − ·NO₂"),
        ("c[$([N+](=O)[O-]),$(N(=O)=O)]", 30, .6, "M − NO"),
        ("[CX4]Cl", 36, .7, "M − HCl"),
        ("[CX4][NX3H2]", 1, .3, "M − H·"),
    ]:
        if matches(mol, smarts):
            add(M - d, s, lab)

    if matches(mol, "c1ccccc1"):
        add(77, .9, "catión fenilo C₆H₅⁺", Counter(C=6, H=5))
        add(51, .5, "C₄H₃⁺ (fenilo − C₂H₂)", Counter(C=4, H=3))
        add(65, .25, "C₅H₅⁺ (tropilio − C₂H₂)", Counter(C=5, H=5))
        if matches(mol, "[CX4;!H0]c1ccccc1"):
            add(91, 1.2, "ion tropilio C₇H₇⁺", Counter(C=7, H=7))
        if matches(mol, "[CX3](=O)c1ccccc1"):
            add(105, 1.4, "ion benzoílo C₆H₅CO⁺", Counter(C=7, H=5, O=1))
    if matches(mol, "[CX4][CX4][CX4]"):
        for n in range(2, min(cnt["C"], 8)):
            add(14 * n + 1, .25 / (1 + abs(n - 3.5) * .4), f"serie alquílica C{n}H{2 * n + 1}⁺", Counter(C=n, H=2 * n + 1))

    best = max([v["score"] for v in ions.values()] + [1.0])
    ions[M]["score"] = ions[M]["score"] * .3 + m_stab * best
    ions[M]["labels"].append((1e9, "ion molecular M⁺·"))
    ions[M]["cnt"] = cnt
    top = max(v["score"] for v in ions.values())

    spec = defaultdict(float)
    labels = {}
    for mz, v in ions.items():
        rel = 100 * v["score"] / top
        if rel < 1.5 and mz != M:
            continue
        lab = max(v["labels"])[1] if v["labels"] else ""
        pat = isotope_pattern(v["cnt"], 5) if v["cnt"] else [(mz, 100.0)]
        off = mz - pat[0][0]
        for k, (pm, pr) in enumerate(pat):
            spec[pm + off] += rel * pr / 100
            if k == 0:
                labels[pm + off] = lab
            elif mz == M:
                labels.setdefault(pm + off, f"M+{k} (isótopos)")
    t2 = max(spec.values())
    peaks = sorted((mz, min(100.0, 100 * r / t2), labels.get(mz, "")) for mz, r in spec.items() if 100 * r / t2 >= .5)
    return MSResult(peaks, M, isotope_pattern(cnt))


# ---------------------------------------------------------------------------
@dataclass
class Analysis:
    mol: Chem.Mol
    smiles: str
    formula: str
    formula_html: str
    mw: float
    exact_mass: float
    dbe: float
    inchi: str
    inchikey: str
    ir: IRResult
    uv: UVResult
    ms: MSResult


def parse_structure(text: str) -> Chem.Mol | None:
    text = (text or "").strip()
    if not text:
        return None
    mol = Chem.MolFromMolBlock(text) if "\n" in text and "M  END" in text else Chem.MolFromSmiles(text)
    if mol is None and text.startswith("InChI="):
        mol = Chem.MolFromInchi(text)
    return mol


def analyze(text: str) -> Analysis:
    mol = parse_structure(text)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Could not parse the structure (invalid SMILES, Molfile or InChI).")
    cnt = formula_counts(mol)
    C = cnt["C"] + cnt["Si"]
    X = cnt["H"] + cnt["F"] + cnt["Cl"] + cnt["Br"] + cnt["I"]
    N = cnt["N"] + cnt["P"] + cnt["B"]
    try:
        inchi = Chem.MolToInchi(mol)
        key = Chem.InchiToInchiKey(inchi)
    except Exception:
        inchi, key = "", ""
    return Analysis(
        mol=mol, smiles=Chem.MolToSmiles(mol), formula=hill(cnt), formula_html=hill(cnt, html=True),
        mw=Descriptors.MolWt(mol), exact_mass=rdMolDescriptors.CalcExactMolWt(mol), dbe=C - X / 2 + N / 2 + 1,
        inchi=inchi, inchikey=key, ir=simulate_ir(mol), uv=simulate_uv(mol), ms=simulate_ms(mol),
    )
