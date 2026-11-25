"""Spectrum Simulator – For Learning

Educational simulator of UV-Vis, IR and EI-MS spectra that compares the simulation with
experimental data (NIST Chemistry WebBook) after identifying the compound in PubChem.
Bilingual interface (Spanish / English).

Run:  streamlit run app.py

State handling: Streamlit deletes the state of widgets that are not drawn in a run (for
example the inputs of another entry mode) and may reset a widget whose label changes
(language switch). Every user choice is therefore mirrored in `st.session_state["_keep"]`,
and the structure being studied lives in `st.session_state["_keep"]["current"]`, so changing
the language or the entry mode never loses information.
"""
from __future__ import annotations

import base64
import html
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots
from rdkit import Chem
from rdkit.Chem import Crippen, Draw, Lipinski, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

from spectra import databases as db
from spectra.i18n import EXAMPLES, HELP, LANGS, UI, label
from spectra import nmr
from spectra.simulate import analyze

ICON = Path(__file__).parent / "launcher" / "favicon.png"
st.set_page_config(page_title="Spectrum Simulator – For Learning", page_icon=str(ICON) if ICON.exists() else "🧪", layout="wide")

SIM = "#2563eb"  # simulated
EXP = "#dc2626"  # experimental
DEFAULT_SMILES = "CC(=O)c1ccccc1"

st.markdown(
    """
    <style>
      .block-container {padding-top: 1.6rem;}
      .chip {display:inline-block;padding:2px 10px;margin:2px 4px 2px 0;border-radius:999px;
             background:rgba(37,99,235,.10);border:1px solid rgba(37,99,235,.35);font-size:.85rem}
      .legend-sim {color:#2563eb;font-weight:600} .legend-exp {color:#dc2626;font-weight:600}
      .note {font-size:.85rem;opacity:.8}
    </style>
    """,
    unsafe_allow_html=True,
)

# ---------------------------------------------------------------------------
# Persistent state
# ---------------------------------------------------------------------------
KEEP = st.session_state.setdefault("_keep", {
    "lang": "es", "mode": "example", "current": DEFAULT_SMILES,
    "name": "acetophenone", "smiles_text": DEFAULT_SMILES, "example": 0, "drawn": None,
    "want_exp": True, "overlay": "overlay", "labels": True, "section": "info", "nucleus": "1H", "mhz": 400, "integral": True, "struct_in_chart": True, "axes": {},
    "log_axis": True, "mirror": True, "ir_pick": 0, "not_found": None,
    "extract_name": "", "refs": [], "ref_error": None,
    "peaks": {"uid": "peaks", "on": True, "sim": False, "nist": True, "user": True, "ref": False, "n": 6, "prom": 5.0},
    "style": {"uid": "style", "title": "", "font": 14, "lw": 1.6, "grid": True, "legend": "bottom", "height": 480, "fmt": "png", "visband": True},
    "user": {"ir": [], "uv": [], "ms": [], "nmr": []}, "upload_n": 0, "uid": 0, "read_errors": [],
})


def wkey_for(name: str) -> str:
    # one widget per language: a label change then never resets or freezes a widget
    return f"w_{name}_{KEEP['lang']}"


def keep(name: str) -> str:
    """Return the widget key for a persisted setting, syncing it with the stored value."""
    wkey = wkey_for(name)
    if st.session_state.get(wkey, object()) != KEEP[name]:
        st.session_state[wkey] = KEEP[name]
    return wkey


def save(name: str):
    """Callback factory: copy a widget's value into the persistent store."""
    wkey = wkey_for(name)

    def _cb():
        value = st.session_state[wkey]
        if value is None:  # clicking the active option of a segmented control deselects it
            st.session_state[wkey] = KEEP[name]
        else:
            KEEP[name] = value
    return _cb


def T(key: str, **kw) -> str:
    text = UI[key][0 if KEEP["lang"] == "es" else 1]
    return text.format(**kw) if kw else text


L = lambda s: label(s, KEEP["lang"])


def item_widget(item, field, factory, *args, **kw):
    """Draw a widget bound to item[field]; its value persists even if the widget disappears."""
    wkey = f"u_{item['uid']}_{field}_{KEEP['lang']}"
    if st.session_state.get(wkey, object()) != item[field]:
        st.session_state[wkey] = item[field]

    def _cb():
        item[field] = st.session_state[wkey]
    return factory(*args, key=wkey, on_change=_cb, **kw)




# ---------------------------------------------------------------------------
# Cached data
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False, ttl=24 * 3600)
def cached_pubchem(query, kind):
    return db.pubchem_lookup(query, kind)


@st.cache_data(show_spinner=False, ttl=24 * 3600)
def cached_nist(cas):
    return {"IR": db.nist_ir_all(cas), "Mass": db.nist_spectrum(cas, "Mass"), "UVVis": db.nist_spectrum(cas, "UVVis")}


@st.cache_data(show_spinner=False)
def cached_analysis(smiles):
    return analyze(smiles)


def base_layout(fig, xtitle, ytitle, reverse_x=False, height=430):
    fig.update_layout(
        height=height, margin=dict(l=60, r=20, t=30, b=55), hovermode="x unified",
        legend=dict(orientation="h", y=1.08, x=0), template="plotly_white",
        xaxis=dict(title=xtitle, autorange="reversed" if reverse_x else True, showgrid=True, gridcolor="rgba(128,128,128,.15)"),
        yaxis=dict(title=ytitle, showgrid=True, gridcolor="rgba(128,128,128,.15)"),
        font=dict(family="Inter, Segoe UI, sans-serif", size=13),
    )
    return fig


# ---------------------------------------------------------------------------
# Callbacks that change the structure being studied
# ---------------------------------------------------------------------------
def on_name_search():
    KEEP["name"] = st.session_state[wkey_for("name")]
    q = KEEP["name"].strip()
    if not q:
        return
    kind = "inchikey" if len(q) == 27 and q.count("-") == 2 else "name"
    comp = cached_pubchem(q, kind)
    if comp and comp.smiles:
        KEEP["current"], KEEP["not_found"] = comp.smiles, None
    else:
        KEEP["not_found"] = q


def on_smiles_apply():
    KEEP["smiles_text"] = st.session_state[wkey_for("smiles_text")]
    if KEEP["smiles_text"].strip():
        KEEP["current"] = KEEP["smiles_text"].strip()


def on_example():
    KEEP["example"] = st.session_state[wkey_for("example")]
    KEEP["current"] = EXAMPLES[KEEP["example"]][2]


def on_add_ref():
    q = (st.session_state.get("ref_query") or "").strip()
    KEEP["ref_error"] = None
    if not q:
        return
    smiles, label = None, q
    looks_smiles = any(c in q for c in "=#()[]@/\\0123456789") or q.isupper()
    if looks_smiles and Chem.MolFromSmiles(q) is not None:
        smiles = q
        comp = cached_pubchem(q, "smiles")
        if comp:
            label = comp.name
    else:
        comp = cached_pubchem(q, "name")
        if comp and comp.smiles:
            smiles, label = comp.smiles, comp.name
        elif Chem.MolFromSmiles(q) is not None:
            smiles = q
    if smiles is None:
        KEEP["ref_error"] = q
        return
    KEEP["uid"] += 1
    KEEP["refs"].append({"uid": KEEP["uid"], "label": label, "smiles": smiles, "show": True})


def on_mode():
    KEEP["mode"] = st.session_state[wkey_for("mode")]
    KEEP["not_found"] = None
    # the text box starts from the structure currently on screen
    KEEP["smiles_text"] = KEEP["current"]


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
sb = st.sidebar
sb.title("🧪 Spectrum Simulator")
sb.segmented_control(T("language"), list(LANGS), format_func=LANGS.get, key=keep("lang"), on_change=save("lang"))
sb.caption(T("subtitle"))

MODES = {"name": "m_name", "draw": "m_draw", "smiles": "m_smiles", "example": "m_example", "extract": "m_extract"}
sb.radio(T("how"), list(MODES), format_func=lambda m: T(MODES[m]), key=keep("mode"), on_change=on_mode)
mode = KEEP["mode"]

if mode == "name":
    with sb.form("name_form", border=False):
        st.text_input(T("name_input"), key=keep("name"))
        st.form_submit_button(T("search"), on_click=on_name_search, use_container_width=True)
    if KEEP["not_found"]:
        sb.error(T("not_found", q=KEEP["not_found"]))
elif mode == "draw":
    try:
        from streamlit_ketcher import st_ketcher
        sb.info(T("draw_hint"))
        with st.expander(T("editor"), expanded=True):
            drawn = st_ketcher(KEEP["current"], height=430, key="ketcher")
        if drawn and drawn != KEEP["drawn"]:
            KEEP["drawn"] = drawn
            if drawn != KEEP["current"]:
                KEEP["current"] = drawn
                st.rerun()
    except ImportError:
        sb.error(T("no_ketcher"))
elif mode == "smiles":
    with sb.form("smiles_form", border=False):
        st.text_area(T("smiles_input"), key=keep("smiles_text"), height=110)
        st.form_submit_button(T("apply"), on_click=on_smiles_apply, use_container_width=True)
elif mode == "extract":
    sb.text_input(T("extract_name"), key=keep("extract_name"), on_change=save("extract_name"), placeholder=T("extract_default"))
    sb.caption(T("extract_hint"))
else:
    idx = 0 if KEEP["lang"] == "es" else 1
    sb.selectbox(T("compound"), range(len(EXAMPLES)), format_func=lambda i: EXAMPLES[i][idx], key=keep("example"), on_change=on_example)

with sb.expander(f"{T('refs')} ({len(KEEP['refs'])})", expanded=mode == "extract"):
    st.caption(T("refs_help"))
    with st.form("ref_form", border=False, clear_on_submit=True):
        st.text_input(T("ref_input"), key="ref_query")
        st.form_submit_button(T("add"), on_click=on_add_ref, use_container_width=True)
    if KEEP["ref_error"]:
        st.warning(T("ref_error", q=KEEP["ref_error"]))
    for r in list(KEEP["refs"]):
        rc = st.columns([3, 1])
        with rc[0]:
            item_widget(r, "show", st.checkbox, r["label"])
        with rc[1]:
            if st.button("✕", key=f"rmref_{r['uid']}"):
                KEEP["refs"].remove(r)
                st.rerun()

sb.divider()
want_exp = sb.toggle(T("want_exp"), key=keep("want_exp"), on_change=save("want_exp"))
sb.radio(T("compare"), ["overlay", "stacked"], format_func=lambda o: T(o), horizontal=True, key=keep("overlay"), on_change=save("overlay"))
overlay = KEEP["overlay"] == "overlay"
show_labels = sb.toggle(T("labels"), key=keep("labels"), on_change=save("labels"))

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
A, compound = None, None
exp = {"IR": [], "Mass": None, "UVVis": None}
if mode != "extract":
    if not KEEP["current"]:
        st.info(T("empty"))
        st.stop()
    try:
        A = cached_analysis(KEEP["current"])
    except ValueError:
        st.error(T("bad_structure"))
        st.stop()
    if A.inchikey:
        with st.spinner(T("finding")):
            compound = cached_pubchem(A.inchikey, "inchikey")
    if want_exp and compound and compound.cas:
        with st.spinner(T("downloading")):
            exp = cached_nist(compound.cas)
n_exp = len(exp["IR"]) + (exp["Mass"] is not None) + (exp["UVVis"] is not None)

REFS = []  # (label, Analysis) for the visible reference compounds
for r in KEEP["refs"]:
    if r["show"]:
        try:
            REFS.append((r["label"], cached_analysis(r["smiles"])))
        except ValueError:
            pass


def default_title():
    if mode == "extract":
        return KEEP["extract_name"].strip() or T("extract_default")
    return compound.name if compound else ""


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
if A is not None:
    c1, c2 = st.columns([1, 2.2], gap="large")
    with c1:
        st.image(Draw.MolToImage(A.mol, size=(420, 320)), use_container_width=True)
    with c2:
        st.subheader(compound.name if compound else T("unknown"))
        m1, m2, m3, m4 = st.columns(4)
        m1.markdown(f"**{T('formula')}**<br><span style='font-size:1.4rem'>{A.formula_html}</span>", unsafe_allow_html=True)
        m2.metric(T("molar_mass"), f"{A.mw:.2f}")
        m3.metric(T("exact_mass"), f"{A.exact_mass:.4f}")
        m4.metric(T("dbe"), f"{A.dbe:g}")
        st.code(A.smiles, language=None)
        if compound:
            links = [f"[PubChem CID {compound.cid}]({compound.url})"]
            if compound.cas:
                links.append(f"CAS {compound.cas} · [NIST WebBook]({db.nist_page(compound.cas)})")
            st.markdown(" · ".join(links))
            if want_exp and compound.cas:
                st.markdown(f"{T('found')}: **IR {len(exp['IR'])} · MS {int(exp['Mass'] is not None)} · UV-Vis {int(exp['UVVis'] is not None)}**")
            elif want_exp:
                st.markdown(f"<span class='note'>{T('no_cas')}</span>", unsafe_allow_html=True)
        else:
            st.markdown(f"<span class='note'>{T('new_structure')}</span>", unsafe_allow_html=True)
        legend = f"<span class='legend-sim'>━ {T('sim')}</span>"
        if n_exp:
            legend += f" &nbsp; <span class='legend-exp'>━ {T('exp')}</span>"
        st.markdown(legend, unsafe_allow_html=True)
else:
    st.subheader("🧫 " + default_title())
    st.caption(T("extract_hint"))
    if KEEP["refs"]:
        st.markdown(" ".join(f"<span class='chip'>{html.escape(r['label'])}</span>" for r in KEEP["refs"] if r["show"]), unsafe_allow_html=True)

SECTIONS = {"info": "tab_info", "ir": "tab_ir", "uv": "tab_uv", "ms": "tab_ms", "nmr": "tab_nmr", "help": "tab_help"}
st.segmented_control(T("section"), list(SECTIONS), format_func=lambda s: T(SECTIONS[s]), key=keep("section"),
                     on_change=save("section"), label_visibility="collapsed")
section = KEEP["section"]


# ---------------------------------------------------------------------------
# User instrument files (kept in KEEP["user"], so they survive language/mode changes)
# ---------------------------------------------------------------------------
USER_COLORS = ["#059669", "#d97706", "#7c3aed", "#db2777", "#0891b2", "#65a30d"]
YTYPES = {"ir": ["yt_T", "yt_T1", "yt_abs"], "uv": ["yt_abs", "yt_T", "yt_eps", "yt_logeps"], "ms": ["yt_int"], "nmr": ["yt_nmr"]}


def guess_ytype(sec, header, y):
    h = (header or "").lower()
    if sec in ("ms", "nmr"):
        return "yt_int" if sec == "ms" else "yt_nmr"
    if "log" in h:
        return "yt_logeps"
    if "abs" in h or h.strip() in ("a", "au", "a.u."):
        return "yt_abs"
    if "%" in h or "trans" in h or h.strip() in ("t", "%t"):
        return "yt_T1" if np.nanmax(y) <= 1.5 else "yt_T"
    if sec == "ir":
        return "yt_T1" if np.nanmax(y) <= 1.5 else "yt_T"
    return "yt_abs"


def user_panel(sec):
    items = KEEP["user"][sec]
    with st.expander(f"{T('my_data')} ({len(items)})", expanded=bool(items)):
        files = st.file_uploader(T("upload"), type=["csv", "txt", "tsv", "asc", "prn", "dat", "dx", "jdx", "jcamp"],
                                 accept_multiple_files=True, help=T("upload_help"), key=f"up_{sec}_{KEEP['upload_n']}")
        if files:
            KEEP["read_errors"] = []
            for f in files:
                try:
                    tab = db.read_user_file(f.getvalue(), f.name)
                except ValueError:
                    KEEP["read_errors"].append(f.name)
                    continue
                KEEP["uid"] += 1
                ncols = tab.data.shape[1]
                ycols = list(range(1, min(ncols, 5)))
                items.append({"uid": KEEP["uid"], "name": f.name, "cols": tab.columns, "data": tab.data.tolist(),
                              "x": 0, "y": ycols, "ytype": guess_ytype(sec, tab.columns[1] if ncols > 1 else "", tab.data[:, 1]),
                              "conc": 0.0, "path": 1.0, "show": True})
            KEEP["upload_n"] += 1  # new uploader key -> the widget is emptied, the data stays in KEEP
            st.rerun()
        for f in KEEP["read_errors"]:
            st.warning(T("read_error", f=f))
        for item in list(items):
            st.markdown(f"**{item['name']}** · {len(item['data'])} pts")
            cols = st.columns([1.1, 1.6, 1.4, .8, .8] if sec == "uv" else [1.1, 1.6, 1.4, .8])
            opts = list(range(len(item["cols"])))
            with cols[0]:
                item_widget(item, "x", st.selectbox, T("x_col"), opts, format_func=lambda i: item["cols"][i])
            with cols[1]:
                item_widget(item, "y", st.multiselect, T("y_cols"), opts, format_func=lambda i: item["cols"][i])
            with cols[2]:
                item_widget(item, "ytype", st.selectbox, T("y_type"), YTYPES[sec], format_func=T)
            if sec == "uv":
                with cols[3]:
                    item_widget(item, "conc", st.number_input, T("conc"), min_value=0.0, format="%.6f", step=1e-5)
            with cols[-1]:
                item_widget(item, "show", st.checkbox, T("show"))
                if st.button(T("remove"), key=f"rm_{item['uid']}"):
                    items.remove(item)
                    st.rerun()
        if sec == "uv" and items:
            st.caption(T("user_note_uv"))
    return [it for it in items if it["show"]]


def user_series(sec, sim_max=None):
    """Yield (label, x, y) for every visible user column, converted to the units of the section's plot."""
    for item in user_panel(sec):
        d = np.array(item["data"], dtype=float)
        for yc in item["y"]:
            if yc >= d.shape[1] or yc == item["x"]:
                continue
            x, y = d[:, item["x"]], d[:, yc]
            ok = ~(np.isnan(x) | np.isnan(y))
            x, y = x[ok], y[ok]
            if not len(x):
                continue
            t = item["ytype"]
            if sec == "ir":
                if np.nanmax(x) < 50:  # micrometres -> cm-1
                    x = 1e4 / x
                y = 100 * 10 ** (-y) if t == "yt_abs" else (100 * y if t == "yt_T1" else y)
            elif sec == "uv":
                if t == "yt_T":
                    y, t = -np.log10(np.clip(y, 1e-3, None) / 100), "yt_abs"
                if t == "yt_abs":
                    if item["conc"] > 0:
                        y = y / (item["conc"] * item["path"])  # Beer-Lambert: ε = A / (c·b)
                    elif sim_max is not None:
                        y = y / max(np.nanmax(y), 1e-9) * sim_max  # scale to the simulated maximum
                elif t == "yt_logeps":
                    y = 10 ** y
            elif sec == "nmr":
                y = y / max(np.nanmax(y), 1e-9)
            else:
                y = 100 * y / max(np.nanmax(y), 1e-9)
            yield f"{item['name']} · {item['cols'][yc]}", x, y


# ---------------------------------------------------------------------------
# Charts: series, peak marking, style
# ---------------------------------------------------------------------------
REF_COLORS = ["#9333ea", "#ea580c", "#0d9488", "#be123c", "#4f46e5", "#ca8a04"]


@dataclass
class Series:
    name: str
    x: np.ndarray
    y: np.ndarray
    kind: str          # sim, ref, nist, user
    color: str
    sticks: bool = False
    row: int = 1
    sign: int = 1


def find_peaks(x, y, n, prom_pct, minima=False, sticks=False):
    """Simple peak picking: local maxima (or minima) ranked by prominence. Returns [(x, y)]."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return []
    o = np.argsort(x)
    x, y = x[o], y[o]
    v = -y if minima else y
    span = float(v.max() - v.min()) or 1.0
    if sticks:
        idx = [i for i in np.argsort(-v) if v[i] - v.min() >= prom_pct / 100 * span][:n]
        return [(x[i], y[i]) for i in idx]
    cand = np.where((v[1:-1] > v[:-2]) & (v[1:-1] >= v[2:]))[0] + 1
    cand = cand[np.argsort(-v[cand])][:400]
    scored = []
    for i in cand:
        left = np.where(v[:i] > v[i])[0]
        right = np.where(v[i + 1:] > v[i])[0]
        lo = left[-1] if len(left) else 0
        hi = i + 1 + right[0] if len(right) else len(v) - 1
        prom = v[i] - max(v[lo:i + 1].min(), v[i:hi + 1].min())
        if prom >= prom_pct / 100 * span:
            scored.append((prom, i))
    scored.sort(reverse=True)
    return [(x[i], y[i]) for _, i in scored[:n]]


def style_layout(fig, xtitle, ytitle, reverse_x, stacked):
    stl = KEEP["style"]
    legend = {"bottom": dict(orientation="h", yref="container", y=0.01, yanchor="bottom", x=0.5, xanchor="center"),
              "top": dict(orientation="h", y=1.02, yanchor="bottom", x=0),
              "right": dict(orientation="v", y=1, x=1.02, xanchor="left")}.get(stl["legend"])
    title = (stl["title"] or default_title()).strip()
    grid = "rgba(128,128,128,.18)" if stl["grid"] else "rgba(0,0,0,0)"
    fig.update_layout(
        height=int(stl["height"]) + (180 if stacked else 0), template="simple_white",
        margin=dict(l=70, r=30 if stl["legend"] != "right" else 200, t=80 if title else 40,
                    b=60 + (24 + 22 * (len(fig.data) > 4)) * (stl["legend"] == "bottom")),
        title=dict(text=html.escape(title), x=0.5, xanchor="center", yref="container", y=1, yanchor="top", pad=dict(t=18),
                   font=dict(size=int(stl["font"]) + 6)) if title else None,
        showlegend=stl["legend"] != "none", legend=legend, hovermode="x unified",
        font=dict(family="Arial, Helvetica, sans-serif", size=int(stl["font"]), color="#111827"),
    )
    fig.update_xaxes(showgrid=stl["grid"], gridcolor=grid, ticks="outside", mirror=True, showline=True, linecolor="#111827")
    fig.update_yaxes(showgrid=stl["grid"], gridcolor=grid, ticks="outside", mirror=True, showline=True, linecolor="#111827",
                     title_text=ytitle)
    if stacked:
        fig.update_xaxes(title_text=xtitle, row=2, col=1)
    else:
        fig.update_xaxes(title_text=xtitle)
    if reverse_x:
        fig.update_xaxes(autorange="reversed")


def add_structure(fig, mol, signals, band=230):
    """Put the structure, with the signal letters/numbers on its atoms, above the plot (included in exports)."""
    m = Chem.Mol(mol)
    for sgn in signals:
        for i in sgn.atoms:
            m.GetAtomWithIdx(i).SetProp("atomNote", sgn.label)
    d2 = rdMolDraw2D.MolDraw2DCairo(720, 330)
    opts = d2.drawOptions()
    opts.setBackgroundColour((1, 1, 1, 0))
    opts.annotationFontScale = 0.85
    opts.bondLineWidth = 2
    d2.DrawMolecule(m)
    d2.FinishDrawing()
    uri = "data:image/png;base64," + base64.b64encode(d2.GetDrawingText()).decode()
    h0, t0, b0 = fig.layout.height, fig.layout.margin.t, fig.layout.margin.b
    fig.update_layout(height=h0 + band, margin_t=t0 + band)
    plot_h = h0 + band - (t0 + band) - b0
    fig.add_layout_image(source=uri, xref="paper", yref="paper", x=0.5, y=1.0 + 8 / plot_h, xanchor="center", yanchor="bottom",
                         sizex=1.0, sizey=(band - 16) / plot_h, sizing="contain", layer="above")


def build_fig(series, xtitle, ytitle, reverse_x=False, mirror=False):
    lw = float(KEEP["style"]["lw"])
    has_top = any(s.kind in ("sim", "ref") for s in series)
    has_bottom = any(s.kind in ("nist", "user") for s in series)
    stacked = (not overlay) and has_top and has_bottom and not mirror
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=.07,
                        subplot_titles=(T("sim"), T("exp_any"))) if stacked else go.Figure()
    for s in series:
        s.row = 2 if stacked and s.kind in ("nist", "user") else 1
        s.sign = -1 if mirror and s.kind in ("nist", "user") else 1
        y = s.sign * np.asarray(s.y, float)
        dash = "dash" if s.kind == "ref" else None
        kw = dict(row=s.row, col=1) if stacked else {}
        if s.sticks and len(s.x) <= 600:
            xs, ys = [], []
            for a, b in zip(s.x, y):
                xs += [a, a, None]
                ys += [0, b, None]
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name=s.name, hoverinfo="skip",
                                     line=dict(color=s.color, width=lw + .6, dash=dash)), **kw)
            fig.add_trace(go.Scatter(x=s.x, y=y, mode="markers", marker=dict(size=4, color=s.color), showlegend=False,
                                     hovertemplate=f"{html.escape(s.name)}<br>m/z %{{x}}<br>%{{y:.0f}}%<extra></extra>"), **kw)
        else:
            fig.add_trace(go.Scatter(x=s.x, y=y, mode="lines", name=s.name, line=dict(color=s.color, width=lw, dash=dash)), **kw)
    style_layout(fig, xtitle, ytitle, reverse_x, stacked)
    if mirror:
        fig.update_layout(hovermode="closest")
    fig._stacked = stacked
    return fig


def mark_peaks(fig, series, sec, fmt, minima=False, unit=""):
    """Annotate automatically picked peaks on the chosen series; returns rows for the peak table."""
    pk = KEEP["peaks"]
    rows = []
    if not pk["on"]:
        return rows
    for s in series:
        if not pk.get(s.kind, False):
            continue
        found = find_peaks(s.x, s.y, int(pk["n"]), float(pk["prom"]), minima=minima, sticks=s.sticks)
        if not found:
            continue
        top_y = max((yy for _, yy in found), key=lambda v: -v if minima else v)
        for xx, yy in sorted(found, key=lambda p: p[0]):
            is_max = yy == top_y and sec == "uv"
            text = ("λmax " if is_max else "") + fmt.format(xx) + (unit if sec == "uv" else "")
            kw = dict(row=s.row, col=1) if getattr(fig, "_stacked", False) else {}
            down = minima or s.sign < 0
            fig.add_annotation(x=xx, y=s.sign * yy, text=f"<b>{text}</b>" if is_max else text, showarrow=True, arrowhead=0,
                               arrowcolor=s.color, ax=0, ay=26 if down else -26,
                               font=dict(size=max(9, int(KEEP["style"]["font"]) - 2), color=s.color), **kw)
            rows.append((s.name, fmt.format(xx) + unit, f"{yy:.3g}"))
    return rows


def chart_controls(sec, rows):
    """Peak table and the expanders for peak marking, style and axes (all persisted)."""
    pk, stl = KEEP["peaks"], KEEP["style"]
    if rows:
        with st.expander(f"{T('picked')} ({len(rows)})"):
            df = pd.DataFrame(rows, columns=[T("series"), T("position"), T("intensity")])
            st.dataframe(df, hide_index=True, use_container_width=True)
            download(df, f"{sec}_peaks.csv")
    with st.expander(T("peaks")):
        c = st.columns([1.2, 1, 1, 1, 1])
        with c[0]:
            item_widget(pk, "on", st.toggle, T("mark_on"))
        for col, k in zip(c[1:], ("sim", "nist", "user", "ref")):
            with col:
                item_widget(pk, k, st.checkbox, T("k_" + k))
        c2 = st.columns(2)
        with c2[0]:
            item_widget(pk, "n", st.slider, T("n_peaks"), 1, 25)
        with c2[1]:
            item_widget(pk, "prom", st.slider, T("prom"), 0.5, 50.0, step=0.5)
    with st.expander(T("style")):
        item_widget(stl, "title", st.text_input, T("chart_title"), placeholder=default_title())
        c = st.columns(4)
        with c[0]:
            item_widget(stl, "font", st.slider, T("font"), 9, 24)
        with c[1]:
            item_widget(stl, "lw", st.slider, T("lw"), 0.5, 4.0, step=0.1)
        with c[2]:
            item_widget(stl, "height", st.slider, T("height"), 320, 900, step=20)
        with c[3]:
            item_widget(stl, "fmt", st.selectbox, T("fmt"), ["png", "svg", "jpeg"])
        c = st.columns(3)
        with c[0]:
            item_widget(stl, "legend", st.selectbox, T("legend_pos"), ["bottom", "top", "right", "none"], format_func=lambda v: T("lg_" + v))
        with c[1]:
            item_widget(stl, "grid", st.checkbox, T("grid"))
        with c[2]:
            if sec == "uv":
                item_widget(stl, "visband", st.checkbox, T("visband"))


PLOT_CONFIG = {"scrollZoom": True, "displaylogo": False, "modeBarButtonsToAdd": ["drawline", "eraseshape"]}


def show_chart(fig, sec, xdef, ydef=None, reverse=False, rows=()):
    """Draw a chart with zoom enabled and user-editable axis ranges (kept in KEEP["axes"])."""
    ax = KEEP["axes"].setdefault(sec, {"uid": "ax_" + sec, "auto": True, "x0": float(xdef[0]), "x1": float(xdef[1]),
                                       "y0": float(ydef[0]) if ydef else 0.0, "y1": float(ydef[1]) if ydef else 1.0})
    if ax["auto"]:
        xr, yr = list(xdef), list(ydef) if ydef else None
    else:
        xr, yr = sorted([ax["x0"], ax["x1"]]), sorted([ax["y0"], ax["y1"]])
    fig.update_xaxes(range=xr[::-1] if reverse else xr, autorange=False)
    if yr:
        fig.update_yaxes(range=yr, autorange=False)
    # uirevision keeps the user's zoom across reruns until the data or the ranges change
    fig.update_layout(dragmode="zoom", uirevision=f"{sec}|{KEEP['current']}|{mode}|{ax['auto']}|{xr}|{yr}|{KEEP['overlay']}")
    stl = KEEP["style"]
    name = re.sub(r"[^\w-]+", "_", (stl["title"] or default_title() or "spectrum"))[:40]
    st.plotly_chart(fig, use_container_width=True, config=PLOT_CONFIG | {"toImageButtonOptions": {
        "format": stl["fmt"], "scale": 3, "filename": f"{name}_{sec}"}})
    chart_controls(sec, list(rows))
    with st.expander(T("axes")):
        c = st.columns([1.2, 1, 1, 1, 1])
        with c[0]:
            item_widget(ax, "auto", st.checkbox, T("auto_range"))
        for col, key, nm in ((c[1], "x0", "X min"), (c[2], "x1", "X max"), (c[3], "y0", "Y min"), (c[4], "y1", "Y max")):
            with col:
                item_widget(ax, key, st.number_input, nm, disabled=ax["auto"], format="%.3f")
        st.caption(T("zoom_help"))


def download(df, name):
    st.download_button(T("download"), df.to_csv(index=False).encode(), file_name=name, mime="text/csv", key="dl_" + name)


def ref_series(get, kind="ref"):
    """Series for the reference compounds; get(analysis) -> (x, y)."""
    out = []
    for k, (lab, an) in enumerate(REFS):
        x, y = get(an)
        out.append(Series(lab, np.asarray(x), np.asarray(y), kind, REF_COLORS[k % 6]))
    return out


# ---------------------------------------------------------------------------
# IR
# ---------------------------------------------------------------------------
if section == "ir":
    series = []
    if A is not None:
        series.append(Series(T("ir_sim"), A.ir.x, A.ir.y, "sim", SIM))
    series += ref_series(lambda an: (an.ir.x, an.ir.y))
    if exp["IR"]:
        if KEEP["ir_pick"] >= len(exp["IR"]):
            KEEP["ir_pick"] = 0
        opts = [f"{i + 1}. {s.meta.get('STATE', T('state_nd'))} — {s.meta.get('ORIGIN', '')[:50]}" for i, s in enumerate(exp["IR"])]
        st.selectbox(T("exp_pick"), range(len(opts)), format_func=lambda i: opts[i], key=keep("ir_pick"), on_change=save("ir_pick"))
        ex, ey = db.ir_to_transmittance(exp["IR"][KEEP["ir_pick"]])
        series.append(Series(T("ir_exp"), ex, ey, "nist", EXP))
    for k, (nm, ux, uy) in enumerate(user_series("ir")):
        series.append(Series(nm, ux, uy, "user", USER_COLORS[k % 6]))
    fig = build_fig(series, T("wavenumber"), T("transmittance"), reverse_x=True)
    if show_labels and A is not None:
        kw = {"row": 1, "col": 1} if fig._stacked else {}
        for c, w, i, lab in A.ir.peaks[:14]:
            yv = float(np.interp(c, A.ir.x[::-1], A.ir.y[::-1]))
            fig.add_annotation(x=c, y=yv, text=f"{c:.0f}", showarrow=True, arrowhead=0, ay=28, ax=0, font=dict(size=10, color=SIM), **kw)
    rows = mark_peaks(fig, series, "ir", "{:.0f}", minima=True, unit=" cm⁻¹")
    show_chart(fig, "ir", (400, 4000), (0, 105), reverse=True, rows=rows)
    if A is not None:
        st.markdown(f"**{T('groups')}:** " + " ".join(f"<span class='chip'>{L(g)}</span>" for g in A.ir.groups), unsafe_allow_html=True)
        strength = lambda i: T("strong") if i > .7 else T("medium") if i > .4 else T("weak")
        df = pd.DataFrame([(round(c), L(lab), strength(i)) for c, w, i, lab in A.ir.peaks], columns=["ν (cm⁻¹)", T("assignment"), T("intensity")])
        st.dataframe(df, hide_index=True, use_container_width=True)
        download(pd.DataFrame({"wavenumber_cm-1": A.ir.x, "transmittance_%": A.ir.y}), "ir_simulated.csv")

# ---------------------------------------------------------------------------
# UV-Vis
# ---------------------------------------------------------------------------
elif section == "uv":
    raw = A is None and not REFS  # extract without references: plot the instrument units as they are
    log_scale = False
    if not raw:
        log_scale = st.toggle(T("log_axis"), key=keep("log_axis"), on_change=save("log_axis"))
    ytr = (lambda v: np.log10(np.maximum(v, 1))) if log_scale else (lambda v: v)
    series = []
    if A is not None:
        series.append(Series(T("uv_sim"), A.uv.x, ytr(A.uv.y), "sim", SIM))
    series += ref_series(lambda an: (an.uv.x, ytr(an.uv.y)))
    if exp["UVVis"] is not None:
        ux, uy, unit = db.uv_to_epsilon(exp["UVVis"])
        is_eps = unit == "ε"
        series.append(Series(f"{T('uv_exp')} ({unit if is_eps else T('absorbance')})", ux, ytr(uy) if is_eps else uy, "nist", EXP))
    maxima = ([float(A.uv.y.max())] if A is not None else []) + [float(an.uv.y.max()) for _, an in REFS]
    ref_max = max(maxima) if maxima else 1.0
    for k, (nm, ux, uy) in enumerate(user_series("uv", sim_max=None if raw else ref_max)):
        series.append(Series(nm, ux, ytr(np.maximum(uy, 1e-9)) if log_scale else uy, "user", USER_COLORS[k % 6]))
    ytitle = T("raw_abs") if raw else ("log ε" if log_scale else "ε (L·mol⁻¹·cm⁻¹)")
    fig = build_fig(series, T("wavelength"), ytitle)
    if KEEP["style"]["visband"]:
        for lo, hi, col in [(380, 450, "#7c3aed"), (450, 495, "#2563eb"), (495, 570, "#16a34a"), (570, 590, "#eab308"), (590, 620, "#f97316"), (620, 750, "#dc2626")]:
            fig.add_vrect(x0=lo, x1=hi, fillcolor=col, opacity=.07, line_width=0)
    if show_labels and A is not None:
        kw = {"row": 1, "col": 1} if fig._stacked else {}
        for l, e, lab in A.uv.bands:
            yv = float(ytr(np.array([max(np.interp(l, A.uv.x, A.uv.y), 1)]))[0])
            fig.add_annotation(x=l, y=yv, text=f"{l:.0f} nm", showarrow=True, ay=-26, ax=0, font=dict(size=10, color=SIM), **kw)
    rows = mark_peaks(fig, series, "uv", "{:.0f}", unit=" nm")
    xs = [float(np.nanmin(s.x)) for s in series if len(s.x)] + [float(np.nanmax(s.x)) for s in series if len(s.x)]
    xdef = (180, 750) if not raw or not xs else (min(xs), max(xs))
    show_chart(fig, "uv", xdef, rows=rows)
    if A is not None:
        uv = A.uv
        cc1, cc2 = st.columns([3, 1])
        with cc1:
            st.markdown(f"**{T('lmax')}:** {uv.lambda_max:.0f} nm")
            df = pd.DataFrame([(round(l), f"{e:,.0f}", f"{np.log10(max(e, 1)):.2f}", L(lab)) for l, e, lab in uv.bands], columns=["λ (nm)", "ε", "log ε", T("transition")])
            st.dataframe(df, hide_index=True, use_container_width=True)
        with cc2:
            st.markdown(f"**{T('color')}**")
            st.markdown(f"<div style='width:100%;height:90px;border-radius:12px;border:1px solid rgba(128,128,128,.4);background:{uv.color}'></div>", unsafe_allow_html=True)
            st.caption(T("color_note"))
        if exp["UVVis"] is not None:
            st.caption(f"{T('exp')}: " + exp["UVVis"].describe())
        download(pd.DataFrame({"wavelength_nm": uv.x, "epsilon": uv.y}), "uv_simulated.csv")

# ---------------------------------------------------------------------------
# MS
# ---------------------------------------------------------------------------
elif section == "ms":
    series = []
    if A is not None:
        series.append(Series(T("ms_sim"), np.array([p[0] for p in A.ms.peaks]), np.array([p[1] for p in A.ms.peaks]), "sim", SIM, sticks=True))
    for k, (lab, an) in enumerate(REFS):
        series.append(Series(lab, np.array([p[0] for p in an.ms.peaks]), np.array([p[1] for p in an.ms.peaks]), "ref", REF_COLORS[k % 6], sticks=True))
    mirror = False
    if (exp["Mass"] is not None or KEEP["user"]["ms"]) and (A is not None or REFS):
        mirror = st.toggle(T("mirror"), key=keep("mirror"), on_change=save("mirror"))
    if exp["Mass"] is not None:
        e = exp["Mass"]
        series.append(Series(T("ms_exp"), e.x, 100 * e.y / e.y.max(), "nist", EXP, sticks=True))
    for k, (nm, ux, uy) in enumerate(user_series("ms")):
        series.append(Series(nm, ux, uy, "user", USER_COLORS[k % 6], sticks=True))
    fig = build_fig(series, "m/z", T("rel_int"), mirror=mirror)
    fig.update_layout(hovermode="closest")
    if show_labels and A is not None:
        for m_, i_, lab in sorted(A.ms.peaks, key=lambda p: -p[1])[:8]:
            fig.add_annotation(x=m_, y=i_, text=str(m_), showarrow=False, yshift=10, font=dict(size=11, color=SIM))
    rows = mark_peaks(fig, series, "ms", "{:.0f}")
    xmax = max([float(np.nanmax(s.x)) for s in series if len(s.x)] + [50.0]) + 10
    show_chart(fig, "ms", (0, xmax), (-110 if mirror else 0, 115), rows=rows)
    if A is not None:
        ms = A.ms
        a, b = st.columns(2)
        with a:
            st.markdown(f"**{T('mol_ion')}:** m/z {ms.M} · " + T("n_rule", p=T("odd") if ms.M % 2 else T("even")))
            st.markdown(f"**{T('iso_pattern')}**")
            st.dataframe(pd.DataFrame([(m_, f"{r:.1f}") for m_, r in ms.molecular_ion], columns=["m/z", T("rel_pct")]), hide_index=True, use_container_width=True)
        with b:
            st.markdown(f"**{T('fragments')}**")
            top = sorted(ms.peaks, key=lambda p: -p[1])[:12]
            st.dataframe(pd.DataFrame([(m_, f"{i_:.0f}", str(ms.M - m_) if m_ < ms.M else "—", L(lab)) for m_, i_, lab in top],
                                      columns=["m/z", "%", T("loss"), T("interpretation")]), hide_index=True, use_container_width=True)
        download(pd.DataFrame([(m_, i_, L(lab)) for m_, i_, lab in ms.peaks], columns=["mz", "rel_intensity", "assignment"]), "ms_simulated.csv")

# ---------------------------------------------------------------------------
# NMR
# ---------------------------------------------------------------------------
elif section == "nmr":
    c_a, c_b, c_c = st.columns([1, 1.4, 1.4])
    with c_a:
        st.segmented_control(T("nucleus"), ["1H", "13C"], format_func=lambda n: {"1H": "¹H", "13C": "¹³C"}[n], key=keep("nucleus"), on_change=save("nucleus"))
    nucleus = KEEP["nucleus"]
    with c_b:
        if nucleus == "1H":
            st.select_slider(T("field"), [60, 90, 200, 300, 400, 500, 600, 800], key=keep("mhz"), on_change=save("mhz"))
    with c_c:
        if nucleus == "1H" and A is not None:
            st.toggle(T("show_integral"), key=keep("integral"), on_change=save("integral"))
        if A is not None:
            st.toggle(T("struct_in_chart"), key=keep("struct_in_chart"), on_change=save("struct_in_chart"))
    predict = (lambda m: nmr.predict_h(m)) if nucleus == "1H" else (lambda m: nmr.predict_c(m))
    spec = (lambda sg: nmr.spectrum_h(sg, KEEP["mhz"])) if nucleus == "1H" else nmr.spectrum_c
    series, signals = [], []
    if A is not None:
        signals = predict(A.mol)
        x, y = spec(signals)
        series.append(Series(f"{T('nmr_sim')} {'¹H' if nucleus == '1H' else '¹³C'}", x, y, "sim", SIM))
    series += ref_series(lambda an: spec(predict(an.mol)))
    for k, (nm, ux, uy) in enumerate(user_series("nmr")):
        series.append(Series(nm, ux, uy, "user", USER_COLORS[k % 6]))
    st.caption(T("nmr_exp_note"))
    fig = build_fig(series, T("shift"), T("rel_int").replace(" (%)", ""), reverse_x=True)
    kw = {"row": 1, "col": 1} if fig._stacked else {}
    if A is not None and nucleus == "1H" and KEEP["integral"] and signals:
        order = np.argsort(-x)
        cum = np.cumsum(y[order]) / max(np.sum(y), 1e-9)
        fig.add_trace(go.Scatter(x=x[order], y=0.25 + 0.7 * cum, name=T("integral_curve"), line=dict(color="#16a34a", width=1.2, dash="dot")), **kw)
    if show_labels and A is not None:
        # signals closer than this share one label (e.g. "4,5"), so the numbers do not overlap
        tol = 0.06 if nucleus == "1H" else 1.0
        groups = []
        for sgn in sorted(signals, key=lambda g: g.shift):
            if groups and sgn.shift - groups[-1][-1].shift < tol:
                groups[-1].append(sgn)
            else:
                groups.append([sgn])
        for grp in groups:
            xs_ = float(np.mean([g.shift for g in grp]))
            yv = max(float(np.interp(g.shift, x, y)) for g in grp)
            text = ",".join(g.label for g in sorted(grp, key=lambda g: g.label))
            fig.add_annotation(x=xs_, y=min(yv, 1.0), text=f"<b>{text}</b>", showarrow=False, yshift=12, font=dict(size=12, color=SIM), **kw)
    rows = mark_peaks(fig, series, "nmr", "{:.2f}" if nucleus == "1H" else "{:.1f}", unit=" ppm")
    if A is not None and KEEP["struct_in_chart"] and signals:
        add_structure(fig, A.mol, signals)
    show_chart(fig, "nmr_" + nucleus, (-0.5, 13) if nucleus == "1H" else (-5, 230), (-0.03, 1.15), reverse=True, rows=rows)

    if A is not None:
        with st.container():
            if nucleus == "1H":
                rows_t = [(s_.label, f"{s_.shift:.2f}", s_.n, s_.mult + (f" (J ≈ {s_.J:.0f} Hz)" if s_.mult not in ("s", "br s", "m") else ""), L(s_.desc)) for s_ in signals]
                df = pd.DataFrame(rows_t, columns=[T("signal"), T("shift"), T("integral"), T("mult"), T("type")])
            else:
                rows_t = [(s_.label, f"{s_.shift:.1f}", s_.n, s_.mult, L(s_.desc)) for s_ in signals]
                df = pd.DataFrame(rows_t, columns=[T("signal"), T("shift"), T("n_c"), "DEPT", T("type")])
            st.dataframe(df, hide_index=True, use_container_width=True)
            st.caption(T("solvent_note"))
        download(pd.DataFrame({"ppm": x[::10], "intensity": y[::10]}), f"nmr_{nucleus}_simulated.csv")

# ---------------------------------------------------------------------------
# Compound information and nomenclature
# ---------------------------------------------------------------------------
elif section == "info" and A is None:
    st.markdown(f"#### {default_title()}")
    st.info(T("no_info_extract"))
    if REFS:
        st.dataframe(pd.DataFrame([(lab, an.formula, f"{an.mw:.2f}", f"{an.exact_mass:.4f}", f"{an.dbe:g}", an.smiles) for lab, an in REFS],
                                  columns=[T("compound"), T("formula"), T("molar_mass"), T("exact_mass"), T("dbe"), "SMILES"]),
                     hide_index=True, use_container_width=True)
    else:
        st.caption(T("no_refs"))

elif section == "info":
    m = A.mol
    left, right = st.columns(2, gap="large")
    with left:
        st.markdown(f"#### {T('names')}")
        if compound:
            st.markdown(f"**{T('iupac')}:** {compound.iupac or '—'}")
            st.markdown(f"**{T('common')}:** {compound.name}")
            if compound.synonyms:
                st.markdown(f"**{T('synonyms')}:**")
                st.markdown(" ".join(f"<span class='chip'>{n}</span>" for n in compound.synonyms), unsafe_allow_html=True)
        else:
            st.info(T("no_name"))
        st.markdown(f"#### {T('identifiers')}")
        ids = [(T("formula"), A.formula), ("SMILES", A.smiles), ("InChI", A.inchi or "—"), ("InChIKey", A.inchikey or "—")]
        if compound:
            ids += [("PubChem CID", str(compound.cid)), ("CAS", compound.cas or "—")]
        st.dataframe(pd.DataFrame(ids, columns=[T("property"), T("value")]), hide_index=True, use_container_width=True)
        st.markdown(f"**{T('p_groups')}:** " + " ".join(f"<span class='chip'>{L(g)}</span>" for g in A.ir.groups if not g.startswith("C–H")), unsafe_allow_html=True)
    with right:
        st.markdown(f"#### {T('properties')}")
        hbd, hba, logp, mw = Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m), Crippen.MolLogP(m), A.mw
        props = [
            (T("molar_mass"), f"{A.mw:.3f}"), (T("exact_mass"), f"{A.exact_mass:.5f}"), (T("dbe"), f"{A.dbe:g}"),
            (T("p_heavy"), m.GetNumHeavyAtoms()), (T("p_atoms"), Chem.AddHs(m).GetNumAtoms()),
            (T("p_rings"), rdMolDescriptors.CalcNumRings(m)), (T("p_arom"), rdMolDescriptors.CalcNumAromaticRings(m)),
            (T("p_logp"), f"{logp:.2f}"), (T("p_mr"), f"{Crippen.MolMR(m):.2f}"), (T("p_tpsa"), f"{rdMolDescriptors.CalcTPSA(m):.1f}"),
            (T("p_hbd"), hbd), (T("p_hba"), hba), (T("p_rot"), rdMolDescriptors.CalcNumRotatableBonds(m)),
            (T("p_fsp3"), f"{rdMolDescriptors.CalcFractionCSP3(m):.2f}"),
            (T("p_stereo"), len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))),
            (T("p_charge"), Chem.GetFormalCharge(m)),
            (T("p_lipinski"), sum([mw > 500, logp > 5, hbd > 5, hba > 10])),
        ]
        st.dataframe(pd.DataFrame([(a, str(b)) for a, b in props], columns=[T("property"), T("value")]), hide_index=True,
                     use_container_width=True, height=36 * (len(props) + 1) + 4)
        st.markdown(f"**{T('p_composition')}**")
        from spectra.simulate import formula_counts
        cnt = formula_counts(m)
        masses = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "P": 30.974, "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.904, "Si": 28.085, "B": 10.81}
        tot = sum(masses.get(e, 0) * n for e, n in cnt.items())
        comp = [(e, n, f"{100 * masses.get(e, 0) * n / tot:.2f} %") for e, n in sorted(cnt.items(), key=lambda kv: (kv[0] != "C", kv[0] != "H", kv[0]))]
        st.dataframe(pd.DataFrame(comp, columns=["", "n", "%"]), hide_index=True, use_container_width=True)
        if compound and compound.props:
            st.markdown(f"**{T('pubchem_props')}**")
            st.dataframe(pd.DataFrame([(a, str(b)) for a, b in compound.props], columns=[T("property"), T("value")]), hide_index=True, use_container_width=True)

else:
    st.markdown(HELP[KEEP["lang"]])

st.caption(T("footer"))
