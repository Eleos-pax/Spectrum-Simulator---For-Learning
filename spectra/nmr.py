"""Educational ¹H and ¹³C NMR prediction (additive increment rules).

* Chemically equivalent atoms are grouped with RDKit symmetry classes.
* ¹H shifts: base values by carbon type + alpha/beta substituent increments (Shoolery-type),
  aromatic protons from benzene (7.27 ppm) with ortho/meta/para substituent increments.
* ¹³C shifts: alkyl base values + alpha/beta substituent increments, aromatic carbons from
  128.5 ppm with ipso/ortho/meta/para increments, typical values for C=O, C≡N, C=C, C≡C.
* Multiplicity from the n+1 rule (first-order, a single average J), DEPT-like H count for ¹³C.
Accuracy is typically ±0.3 ppm (¹H) and ±5 ppm (¹³C): good enough to teach interpretation.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import comb

import numpy as np
from rdkit import Chem

# Aromatic substituent increments: (¹H ortho, meta, para), (¹³C ipso, ortho, meta, para)
ARO = {
    "alkyl": ((-0.18, -0.11, -0.21), (9.3, 0.7, -0.1, -2.9)),
    "vinyl": ((0.06, -0.03, -0.10), (9.5, -2.0, 0.2, -0.5)),
    "phenyl": ((0.30, 0.12, 0.10), (13.0, -1.1, 0.5, -1.0)),
    "OH": ((-0.56, -0.12, -0.45), (26.9, -12.7, 1.4, -7.3)),
    "OR": ((-0.48, -0.09, -0.44), (31.4, -14.4, 1.0, -7.7)),
    "OCOR": ((-0.25, 0.03, -0.13), (22.4, -7.1, 0.4, -3.2)),
    "NH2": ((-0.75, -0.25, -0.65), (18.2, -13.4, 0.8, -10.0)),
    "NHCOR": ((0.12, -0.07, -0.28), (9.7, -8.1, 0.2, -4.4)),
    "NO2": ((0.95, 0.26, 0.38), (20.0, -4.8, 0.9, 5.8)),
    "CHO": ((0.56, 0.22, 0.29), (8.6, 1.3, 0.6, 5.5)),
    "COR": ((0.62, 0.14, 0.21), (9.1, 0.1, 0.0, 4.2)),
    "COOH": ((0.85, 0.18, 0.27), (2.1, 1.5, 0.0, 5.1)),
    "COOR": ((0.71, 0.11, 0.21), (2.0, 1.2, -0.1, 4.8)),
    "CONH2": ((0.61, 0.10, 0.17), (5.0, -1.2, 0.1, 3.4)),
    "CN": ((0.36, 0.18, 0.28), (-15.4, 3.6, 0.6, 3.9)),
    "F": ((-0.26, 0.0, -0.20), (34.8, -13.0, 1.6, -4.4)),
    "Cl": ((0.03, -0.02, -0.09), (6.2, 0.4, 1.3, -1.9)),
    "Br": ((0.18, -0.08, -0.04), (-5.5, 3.4, 1.7, -1.6)),
    "I": ((0.39, -0.21, 0.0), (-34.1, 8.7, 1.4, -0.9)),
    "SR": ((-0.08, -0.10, -0.24), (10.2, -1.8, 0.4, -3.6)),
}


def _carbonyl_kind(c: Chem.Atom) -> str | None:
    """Classify a carbonyl carbon: CHO, COOH, COOR, CONH2, COCl or COR (None if not C=O)."""
    if c.GetSymbol() != "C" or not any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(c).GetSymbol() == "O" for b in c.GetBonds()):
        return None
    single = [b.GetOtherAtom(c) for b in c.GetBonds() if b.GetBondType() == Chem.BondType.SINGLE]
    if any(a.GetSymbol() == "O" and a.GetTotalNumHs() for a in single):
        return "COOH"
    if any(a.GetSymbol() == "O" for a in single):
        return "COOR"
    if any(a.GetSymbol() == "N" for a in single):
        return "CONH2"
    if any(a.GetSymbol() in ("Cl", "Br") for a in single):
        return "COCl"
    if c.GetTotalNumHs():
        return "CHO"
    return "COR"


def _substituent(ring_atom: Chem.Atom, s: Chem.Atom) -> str:
    sym = s.GetSymbol()
    if sym in ("F", "Cl", "Br", "I"):
        return sym
    if sym == "O":
        if s.GetTotalNumHs():
            return "OH"
        other = [n for n in s.GetNeighbors() if n.GetIdx() != ring_atom.GetIdx()]
        return "OCOR" if other and _carbonyl_kind(other[0]) else "OR"
    if sym == "N":
        if s.GetFormalCharge() == 1 or sum(1 for n in s.GetNeighbors() if n.GetSymbol() == "O") >= 2:
            return "NO2"
        return "NHCOR" if any(_carbonyl_kind(n) for n in s.GetNeighbors() if n.GetIdx() != ring_atom.GetIdx()) else "NH2"
    if sym == "S":
        return "SR"
    if sym == "C":
        k = _carbonyl_kind(s)
        if k:
            return {"COCl": "COR"}.get(k, k)
        if any(b.GetBondType() == Chem.BondType.TRIPLE and b.GetOtherAtom(s).GetSymbol() == "N" for b in s.GetBonds()):
            return "CN"
        if s.GetIsAromatic():
            return "phenyl"
        if any(b.GetBondType() == Chem.BondType.DOUBLE for b in s.GetBonds()):
            return "vinyl"
    return "alkyl"


def _aromatic_increments(mol, idx, dmat, which):
    """Sum of substituent increments acting on aromatic atom idx (which=0 for ¹H, 1 for ¹³C)."""
    ring_atoms = set()
    for r in mol.GetRingInfo().AtomRings():
        if idx in r and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r):
            ring_atoms |= set(r)
    total = 0.0
    for j in ring_atoms:
        aj = mol.GetAtomWithIdx(j)
        for s in aj.GetNeighbors():
            if s.GetIdx() in ring_atoms or s.GetSymbol() == "H":
                continue
            d = int(dmat[idx][j])
            kind = _substituent(aj, s)
            inc = ARO.get(kind, ARO["alkyl"])[which]
            if which == 0 and 1 <= d <= 3:
                total += inc[d - 1]
            elif which == 1 and d <= 3:
                total += inc[d]
    # heteroaromatic rings: position relative to a ring N (alpha, beta, gamma) or O/S
    for j in ring_atoms:
        h = mol.GetAtomWithIdx(j)
        d = int(dmat[idx][j])
        if h.GetSymbol() == "N" and h.GetTotalNumHs() == 0:
            total += {1: (1.30, 21.3), 2: (-0.02, -4.7), 3: (0.37, 7.5)}.get(d, (0, 0))[which]
        elif h.GetSymbol() in ("O", "S", "N"):
            total += {1: (0.10, 14.0), 2: (-0.9, -18.0)}.get(d, (0, 0))[which]
    return total


# ---------------------------------------------------------------------------
# ¹H
# ---------------------------------------------------------------------------
H_ALPHA = {  # increment for a CHn group bonded to ...
    "C": 0.0, "C=C": 0.75, "Ar": 1.40, "C=O": 1.20, "C#X": 1.10, "COOH": 1.15, "COOR": 1.10,
    "OH": 2.45, "OR": 2.35, "OCOR": 3.00, "OAr": 2.95, "N": 1.45, "NCO": 1.90, "NAr": 2.05, "NO2": 3.40,
    "F": 3.30, "Cl": 2.15, "Br": 1.80, "I": 1.30, "S": 1.25,
}


def _neighbor_class(a: Chem.Atom, n: Chem.Atom) -> str:
    sym = n.GetSymbol()
    if sym in ("F", "Cl", "Br", "I"):
        return sym
    if sym == "S":
        return "S"
    if sym == "O":
        other = [x for x in n.GetNeighbors() if x.GetIdx() != a.GetIdx()]
        if n.GetTotalNumHs():
            return "OH"
        if other and _carbonyl_kind(other[0]):
            return "OCOR"
        if other and other[0].GetIsAromatic():
            return "OAr"
        return "OR"
    if sym == "N":
        if n.GetFormalCharge() == 1:
            return "NO2"
        others = [x for x in n.GetNeighbors() if x.GetIdx() != a.GetIdx()]
        if any(_carbonyl_kind(x) for x in others):
            return "NCO"
        if any(x.GetIsAromatic() for x in others):
            return "NAr"
        return "N"
    if sym == "C":
        if n.GetIsAromatic():
            return "Ar"
        k = _carbonyl_kind(n)
        if k in ("COOH", "COOR"):
            return k
        if k:
            return "C=O"
        if any(b.GetBondType() == Chem.BondType.TRIPLE for b in n.GetBonds()):
            return "C#X"
        if any(b.GetBondType() == Chem.BondType.DOUBLE for b in n.GetBonds()):
            return "C=C"
    return "C"


def _h_shift(mol, a: Chem.Atom, dmat) -> tuple[float, str, bool]:
    """Return (δ ppm, description, exchangeable) for the hydrogens on atom a."""
    sym, nH = a.GetSymbol(), a.GetTotalNumHs()
    if sym == "O":
        nb = a.GetNeighbors()[0] if a.GetDegree() else None
        if nb is not None and _carbonyl_kind(nb) == "COOH":
            return 11.5, "COOH", True
        if nb is not None and nb.GetIsAromatic():
            return 5.5, "Ar–OH", True
        return 2.2, "O–H", True
    if sym == "N":
        if any(_carbonyl_kind(n) for n in a.GetNeighbors()):
            return 6.5 if nH == 1 else 5.8, "N–H amida", True
        if any(n.GetIsAromatic() for n in a.GetNeighbors()):
            return 3.6, "Ar–NH", True
        return 1.3, "N–H", True
    if sym == "S":
        return 1.5, "S–H", True
    if sym != "C":
        return 1.0, sym + "–H", True
    if a.GetIsAromatic():
        return 7.27 + _aromatic_increments(mol, a.GetIdx(), dmat, 0), "Ar–H", False
    k = _carbonyl_kind(a)
    if k == "CHO":
        aryl = any(n.GetIsAromatic() for n in a.GetNeighbors())
        return (10.0 if aryl else 9.75), "CHO", False
    if any(b.GetBondType() == Chem.BondType.TRIPLE for b in a.GetBonds()):
        return 2.4, "≡C–H", False
    dbl = [b for b in a.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
    if dbl:
        partner = dbl[0].GetOtherAtom(a)
        d = 4.95 if nH == 2 else 5.35
        for n in a.GetNeighbors():
            if n.GetIdx() == partner.GetIdx():
                continue
            c = _neighbor_class(a, n)
            d += {"C=O": .6, "COOR": .55, "COOH": .6, "Ar": 1.0, "OR": 1.2, "OCOR": 1.8, "Cl": .8, "Br": .9, "C#X": .3, "C=C": .4}.get(c, 0.0)
        for n in partner.GetNeighbors():  # beta position across the double bond
            if n.GetIdx() == a.GetIdx():
                continue
            c = _neighbor_class(partner, n)
            d += {"C=O": .95, "COOR": .9, "COOH": 1.0, "Ar": .3, "OR": -1.0, "OCOR": -.5, "C#X": .7}.get(c, -0.2 if c == "C" else 0.0)
        return d, "=C–H", False
    base = {3: 0.90, 2: 1.30, 1: 1.55}.get(nH, 1.3)
    inc = 0.0
    for n in a.GetNeighbors():
        cls = _neighbor_class(a, n)
        inc += H_ALPHA.get(cls, 0.0)
        if cls == "C=O" and any(x.GetIsAromatic() for x in n.GetNeighbors()):
            inc += 0.4  # aryl ketones: CH3 at ~2.6 ppm
        if cls == "C":  # beta substituents
            for m in n.GetNeighbors():
                if m.GetIdx() == a.GetIdx():
                    continue
                cm = _neighbor_class(n, m)
                inc += {"OH": .25, "OR": .25, "OCOR": .4, "Cl": .35, "Br": .45, "I": .45, "F": .3, "N": .15, "NCO": .2, "NO2": .7,
                        "C=O": .2, "COOR": .2, "COOH": .2, "Ar": .25, "C=C": .1}.get(cm, 0.0)
    if inc > 2.8:  # two strong groups do not add fully
        inc = 2.8 + .5 * (inc - 2.8)
    return base + inc, {3: "CH₃", 2: "CH₂", 1: "CH"}.get(nH, "CH"), False


# ---------------------------------------------------------------------------
# ¹³C
# ---------------------------------------------------------------------------
# increments relative to H (Grant-Paul scheme); "C" (sp3) uses 9.1 (alpha) / 9.4 (beta)
C_ALPHA = {"C=C": 20.0, "Ar": 23.0, "C=O": 30.0, "COOH": 21.0, "COOR": 21.0, "C#X": 3.5, "OH": 51.0, "OR": 58.0,
           "OCOR": 53.0, "OAr": 58.0, "N": 28.0, "NCO": 26.0, "NAr": 30.0, "NO2": 63.0, "F": 70.0, "Cl": 31.0, "Br": 20.0, "I": -6.0, "S": 11.0}
C_BETA = {"OH": 10.0, "OR": 8.0, "OCOR": 7.0, "OAr": 8.0, "N": 11.0, "NCO": 8.0, "NAr": 9.0, "NO2": 3.0, "F": 8.0, "Cl": 11.0, "Br": 11.0,
          "I": 11.0, "S": 11.0, "C=O": -2.5, "Ar": 9.0, "C=C": 7.0, "COOR": 2.0, "COOH": 2.0, "C#X": 2.5}


def _c_shift(mol, a: Chem.Atom, dmat) -> tuple[float, str]:
    nH = a.GetTotalNumHs()
    if a.GetIsAromatic():
        return 128.5 + _aromatic_increments(mol, a.GetIdx(), dmat, 1), "C aromático"
    k = _carbonyl_kind(a)
    if k:
        conj = any(n.GetIsAromatic() or any(b.GetBondType() == Chem.BondType.DOUBLE for b in n.GetBonds() if b.GetOtherAtom(n).GetIdx() != a.GetIdx())
                   for n in a.GetNeighbors())
        base = {"CHO": 201.0, "COR": 207.0, "COOH": 179.0, "COOR": 171.0, "CONH2": 172.0, "COCl": 170.0}[k]
        return base - (8.0 if conj and k in ("CHO", "COR") else 6.0 if conj else 0.0), "C=O"
    if any(b.GetBondType() == Chem.BondType.TRIPLE and b.GetOtherAtom(a).GetSymbol() == "N" for b in a.GetBonds()):
        return 118.0, "C≡N"
    if any(b.GetBondType() == Chem.BondType.TRIPLE for b in a.GetBonds()):
        return (68.5 if nH else 80.0), "C≡C"
    dbl = [b for b in a.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetSymbol() == "C"]
    if dbl:
        partner = dbl[0].GetOtherAtom(a)
        d = {2: 115.0, 1: 133.0, 0: 142.0}.get(nH, 128.0)
        for n in a.GetNeighbors():
            if n.GetIdx() != partner.GetIdx():
                d += {"C=O": 4.0, "COOR": 1.0, "OR": 25.0, "OCOR": 18.0, "Ar": 4.0, "Cl": 3.0, "C#X": -15.0}.get(_neighbor_class(a, n), 0.0)
        for n in partner.GetNeighbors():
            if n.GetIdx() != a.GetIdx():
                d += {"C=O": 12.0, "COOR": 8.0, "COOH": 9.0, "OR": -35.0, "OCOR": -25.0, "C#X": 15.0}.get(_neighbor_class(partner, n), 0.0)
        return d, "C=C"
    if a.GetSymbol() != "C":
        return 0.0, ""
    # sp3 carbon: Grant-Paul (-2.3 + 9.1 nα + 9.4 nβ - 2.5 nγ) with substituent increments
    d = -2.3
    alpha_ids = {n.GetIdx() for n in a.GetNeighbors()}
    for n in a.GetNeighbors():
        cls = _neighbor_class(a, n)
        d += 9.1 if cls == "C" else C_ALPHA.get(cls, 9.1)
        if cls != "C":
            continue
        for m in n.GetNeighbors():  # beta
            if m.GetIdx() == a.GetIdx():
                continue
            cm = _neighbor_class(n, m)
            d += 9.4 if cm == "C" else C_BETA.get(cm, 0.0)
            if cm == "C":
                d -= 2.5 * sum(1 for g in m.GetNeighbors() if g.GetSymbol() == "C" and g.GetIdx() not in alpha_ids and g.GetIdx() != n.GetIdx())
    return d, {3: "CH₃", 2: "CH₂", 1: "CH", 0: "C"}.get(nH, "C")


# ---------------------------------------------------------------------------
@dataclass
class Signal:
    shift: float
    n: int            # integral (¹H) or number of equivalent carbons (¹³C)
    mult: str         # s, d, t, q, quint, sext, sept, m, br s
    J: float
    atoms: list
    label: str        # letter shown on the structure
    desc: str
    nH: int = 0       # H on the carbon (¹³C, DEPT)
    split: int = 0    # neighbouring H count (n+1 rule)


MULT = {0: "s", 1: "d", 2: "t", 3: "q", 4: "quint", 5: "sext", 6: "sept"}


def _groups(mol):
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    groups = {}
    for a in mol.GetAtoms():
        groups.setdefault(ranks[a.GetIdx()], []).append(a.GetIdx())
    return list(groups.values())


def predict_h(mol) -> list[Signal]:
    dmat = Chem.GetDistanceMatrix(mol)
    out = []
    for g in _groups(mol):
        a = mol.GetAtomWithIdx(g[0])
        nH = a.GetTotalNumHs()
        if not nH or a.GetSymbol() == "H":
            continue
        d, desc, exch = _h_shift(mol, a, dmat)
        if exch:
            mult, n_split, J = "br s", 0, 0.0
        else:
            gset = set(g)
            # n+1 rule: H on neighbouring carbons that are not equivalent to this group
            n_split = sum(nb.GetTotalNumHs() for nb in a.GetNeighbors() if nb.GetSymbol() == "C" and nb.GetIdx() not in gset)
            aro = a.GetIsAromatic()
            J = 8.0 if aro else 7.0
            mult = MULT.get(n_split, "m") if not (aro and n_split > 2) else "m"
            if any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()) and not aro:
                mult = "m" if n_split > 1 else mult
                J = 10.0
        out.append(Signal(round(d, 2), nH * len(g), mult, J, g, "", desc, nH, n_split))
    out.sort(key=lambda s: -s.shift)
    for i, s in enumerate(out):
        s.label = chr(ord("a") + i) if i < 26 else f"z{i}"
    return out


def predict_c(mol) -> list[Signal]:
    dmat = Chem.GetDistanceMatrix(mol)
    out = []
    for g in _groups(mol):
        a = mol.GetAtomWithIdx(g[0])
        if a.GetSymbol() != "C":
            continue
        d, desc = _c_shift(mol, a, dmat)
        nH = a.GetTotalNumHs()
        out.append(Signal(round(d, 1), len(g), {0: "C", 1: "CH", 2: "CH₂", 3: "CH₃", 4: "CH₄"}[nH], 0.0, g, "", desc, nH))
    out.sort(key=lambda s: -s.shift)
    for i, s in enumerate(out):
        s.label = str(i + 1)
    return out


def spectrum_h(signals, mhz=400.0, width_hz=1.2):
    """Build a ¹H spectrum (first-order multiplets, binomial intensities)."""
    x = np.linspace(-0.5, 13.0, 27000)
    y = np.zeros_like(x)
    hw = width_hz / mhz
    for s in signals:
        n = s.split if s.mult not in ("m", "br s") else 0
        lines = [((k - n / 2) * s.J / mhz, comb(n, k)) for k in range(n + 1)]
        if s.mult == "m":  # unresolved: a few overlapping lines
            lines = [(o * s.J / mhz, w) for o, w in ((-1, .6), (-.5, .9), (0, 1), (.5, .9), (1, .6))]
        w = hw * (8 if s.mult == "br s" else 1)
        tot = sum(v for _, v in lines)
        for off, inten in lines:
            y += s.n * inten / tot / (1 + ((x - s.shift - off) / w) ** 2)
    return x, y / max(y.max(), 1e-9)


def spectrum_c(signals, width=0.15):
    x = np.linspace(-5, 230, 23500)
    y = np.zeros_like(x)
    for s in signals:
        inten = s.n * (0.45 if s.nH == 0 else 1.0)  # quaternary carbons are weaker (no NOE)
        y += inten / (1 + ((x - s.shift) / width) ** 2)
    return x, y / max(y.max(), 1e-9)
