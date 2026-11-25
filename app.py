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

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Draw, Lipinski, rdMolDescriptors
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
    "want_exp": True, "overlay": "overlay", "labels": True, "section": "info", "nucleus": "1H", "mhz": 400, "integral": True, "axes": {},
    "log_axis": True, "mirror": True, "ir_pick": 0, "not_found": None,
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

MODES = {"name": "m_name", "draw": "m_draw", "smiles": "m_smiles", "example": "m_example"}
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
else:
    idx = 0 if KEEP["lang"] == "es" else 1
    sb.selectbox(T("compound"), range(len(EXAMPLES)), format_func=lambda i: EXAMPLES[i][idx], key=keep("example"), on_change=on_example)

sb.divider()
want_exp = sb.toggle(T("want_exp"), key=keep("want_exp"), on_change=save("want_exp"))
sb.radio(T("compare"), ["overlay", "stacked"], format_func=lambda o: T(o), horizontal=True, key=keep("overlay"), on_change=save("overlay"))
overlay = KEEP["overlay"] == "overlay"
show_labels = sb.toggle(T("labels"), key=keep("labels"), on_change=save("labels"))

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
if not KEEP["current"]:
    st.info(T("empty"))
    st.stop()
try:
    A = cached_analysis(KEEP["current"])
except ValueError:
    st.error(T("bad_structure"))
    st.stop()

compound = None
if A.inchikey:
    with st.spinner(T("finding")):
        compound = cached_pubchem(A.inchikey, "inchikey")

exp = {"IR": [], "Mass": None, "UVVis": None}
if want_exp and compound and compound.cas:
    with st.spinner(T("downloading")):
        exp = cached_nist(compound.cas)
n_exp = len(exp["IR"]) + (exp["Mass"] is not None) + (exp["UVVis"] is not None)

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
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


def item_widget(item, field, factory, *args, **kw):
    """Draw a widget bound to item[field]; its value persists even if the widget disappears."""
    wkey = f"u_{item['uid']}_{field}_{KEEP['lang']}"
    if st.session_state.get(wkey, object()) != item[field]:
        st.session_state[wkey] = item[field]

    def _cb():
        item[field] = st.session_state[wkey]
    return factory(*args, key=wkey, on_change=_cb, **kw)


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
                    y = y / (item["conc"] * item["path"]) if item["conc"] > 0 else y / max(np.nanmax(y), 1e-9) * sim_max
                elif t == "yt_logeps":
                    y = 10 ** y
            elif sec == "nmr":
                y = y / max(np.nanmax(y), 1e-9)
            else:
                y = 100 * y / max(np.nanmax(y), 1e-9)
            yield f"{item['name']} · {item['cols'][yc]}", x, y


PLOT_CONFIG = {"scrollZoom": True, "displaylogo": False, "modeBarButtonsToAdd": ["drawline", "eraseshape"],
               "toImageButtonOptions": {"format": "png", "scale": 3}}


def show_chart(fig, sec, xdef, ydef=None, reverse=False):
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
    fig.update_layout(dragmode="zoom", uirevision=f"{sec}|{KEEP['current']}|{ax['auto']}|{xr}|{yr}|{KEEP['overlay']}")
    st.plotly_chart(fig, use_container_width=True, config=PLOT_CONFIG | {"toImageButtonOptions": {"format": "png", "scale": 3, "filename": f"{sec}_spectrum"}})
    with st.expander(T("axes")):
        c = st.columns([1.2, 1, 1, 1, 1])
        with c[0]:
            item_widget(ax, "auto", st.checkbox, T("auto_range"))
        for col, key, name in ((c[1], "x0", "X min"), (c[2], "x1", "X max"), (c[3], "y0", "Y min"), (c[4], "y1", "Y max")):
            with col:
                item_widget(ax, key, st.number_input, name, disabled=ax["auto"], format="%.3f")
        st.caption(T("zoom_help"))


def stacked_or_overlay(traces_sim, traces_exp, xtitle, ytitle, reverse_x=False):
    if overlay or not traces_exp:
        return base_layout(go.Figure(traces_sim + traces_exp), xtitle, ytitle, reverse_x)
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=.06, subplot_titles=(T("sim"), T("exp_any")))
    for t in traces_sim:
        fig.add_trace(t, 1, 1)
    for t in traces_exp:
        fig.add_trace(t, 2, 1)
    base_layout(fig, None, ytitle, reverse_x, height=620)
    fig.update_xaxes(autorange="reversed" if reverse_x else True)
    fig.update_xaxes(title_text=xtitle, row=2, col=1)
    return fig


def download(df, name):
    st.download_button(T("download"), df.to_csv(index=False).encode(), file_name=name, mime="text/csv", key="dl_" + name)


# ---------------------------------------------------------------------------
# IR
# ---------------------------------------------------------------------------
if section == "ir":
    ir = A.ir
    sim = [go.Scatter(x=ir.x, y=ir.y, name=T("ir_sim"), line=dict(color=SIM, width=1.6))]
    exp_tr = []
    if exp["IR"]:
        if KEEP["ir_pick"] >= len(exp["IR"]):
            KEEP["ir_pick"] = 0

        opts = [f"{i + 1}. {s.meta.get('STATE', T('state_nd'))} — {s.meta.get('ORIGIN', '')[:50]}" for i, s in enumerate(exp["IR"])]
        st.selectbox(T("exp_pick"), range(len(opts)), format_func=lambda i: opts[i], key=keep("ir_pick"), on_change=save("ir_pick"))
        ex, ey = db.ir_to_transmittance(exp["IR"][KEEP["ir_pick"]])
        exp_tr = [go.Scatter(x=ex, y=ey, name=T("ir_exp"), line=dict(color=EXP, width=1.3))]
    for k, (nm, ux, uy) in enumerate(user_series("ir")):
        exp_tr.append(go.Scatter(x=ux, y=uy, name=nm, line=dict(color=USER_COLORS[k % 6], width=1.3)))
    fig = stacked_or_overlay(sim, exp_tr, T("wavenumber"), T("transmittance"), reverse_x=True)
    if show_labels:
        sub = {} if overlay or not exp_tr else {"row": 1, "col": 1}
        for c, w, i, lab in ir.peaks[:14]:
            yv = float(np.interp(c, ir.x[::-1], ir.y[::-1]))
            fig.add_annotation(x=c, y=yv, text=f"{c:.0f}", showarrow=True, arrowhead=0, ay=28, ax=0, font=dict(size=10, color=SIM), **sub)
    show_chart(fig, "ir", (400, 4000), (0, 105), reverse=True)
    st.markdown(f"**{T('groups')}:** " + " ".join(f"<span class='chip'>{L(g)}</span>" for g in ir.groups), unsafe_allow_html=True)
    strength = lambda i: T("strong") if i > .7 else T("medium") if i > .4 else T("weak")
    df = pd.DataFrame([(round(c), L(lab), strength(i)) for c, w, i, lab in ir.peaks], columns=["ν (cm⁻¹)", T("assignment"), T("intensity")])
    st.dataframe(df, hide_index=True, use_container_width=True)
    download(pd.DataFrame({"wavenumber_cm-1": ir.x, "transmittance_%": ir.y}), "ir_simulated.csv")

# ---------------------------------------------------------------------------
# UV-Vis
# ---------------------------------------------------------------------------
elif section == "uv":
    uv = A.uv
    log_scale = st.toggle(T("log_axis"), key=keep("log_axis"), on_change=save("log_axis"))
    ytr = (lambda v: np.log10(np.maximum(v, 1))) if log_scale else (lambda v: v)
    sim = [go.Scatter(x=uv.x, y=ytr(uv.y), name=T("uv_sim"), line=dict(color=SIM, width=2), fill="tozeroy", fillcolor="rgba(37,99,235,.08)")]
    exp_tr = []
    if exp["UVVis"] is not None:
        ux, uy, unit = db.uv_to_epsilon(exp["UVVis"])
        is_eps = unit == "ε"
        exp_tr = [go.Scatter(x=ux, y=ytr(uy) if is_eps else uy, name=f"{T('uv_exp')} ({unit if is_eps else T('absorbance')})", line=dict(color=EXP, width=1.6))]
    for k, (nm, ux, uy) in enumerate(user_series("uv", sim_max=float(uv.y.max()))):
        exp_tr.append(go.Scatter(x=ux, y=ytr(np.maximum(uy, 1e-9)) if log_scale else uy, name=nm, line=dict(color=USER_COLORS[k % 6], width=1.6)))
    fig = stacked_or_overlay(sim, exp_tr, T("wavelength"), "log ε" if log_scale else "ε (L·mol⁻¹·cm⁻¹)")
    for lo, hi, col in [(380, 450, "#7c3aed"), (450, 495, "#2563eb"), (495, 570, "#16a34a"), (570, 590, "#eab308"), (590, 620, "#f97316"), (620, 750, "#dc2626")]:
        fig.add_vrect(x0=lo, x1=hi, fillcolor=col, opacity=.07, line_width=0)
    if show_labels:
        for l, e, lab in uv.bands:
            yv = float(ytr(np.array([max(np.interp(l, uv.x, uv.y), 1)]))[0])
            fig.add_annotation(x=l, y=yv, text=f"{l:.0f} nm", showarrow=True, ay=-26, ax=0, font=dict(size=10, color=SIM))
    show_chart(fig, "uv", (180, 750))
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
    ms = A.ms
    mz = [p[0] for p in ms.peaks]
    it = [p[1] for p in ms.peaks]

    def sticks(x, y, color, name, neg=False):
        xs, ys = [], []
        for a, b in zip(x, y):
            xs += [a, a, None]
            ys += [0, -b if neg else b, None]
        return go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=color, width=2.2), name=name, hoverinfo="skip")

    traces = [sticks(mz, it, SIM, T("ms_sim")),
              go.Scatter(x=mz, y=it, mode="markers", marker=dict(size=6, color=SIM), showlegend=False,
                         text=[L(p[2]) for p in ms.peaks], hovertemplate="m/z %{x}<br>%{y:.0f}%<br>%{text}<extra></extra>")]
    mirror = False
    if exp["Mass"] is not None or KEEP["user"]["ms"]:
        mirror = st.toggle(T("mirror"), key=keep("mirror"), on_change=save("mirror"))
    if exp["Mass"] is not None:
        e = exp["Mass"]
        ey = 100 * e.y / e.y.max()
        traces.append(sticks(e.x, ey, EXP, T("ms_exp"), neg=mirror))
        traces.append(go.Scatter(x=e.x, y=-ey if mirror else ey, mode="markers", marker=dict(size=4, color=EXP), showlegend=False,
                                 hovertemplate="m/z %{x}<br>%{customdata:.0f}% (NIST)<extra></extra>", customdata=ey))
    for k, (nm, ux, uy) in enumerate(user_series("ms")):
        col = USER_COLORS[k % 6]
        if len(ux) <= 600:
            traces.append(sticks(ux, uy, col, nm, neg=mirror))
        else:
            traces.append(go.Scatter(x=ux, y=-uy if mirror else uy, name=nm, line=dict(color=col, width=1.2)))
    fig = base_layout(go.Figure(traces), "m/z", T("rel_int"))
    fig.update_layout(hovermode="closest")

    if show_labels:
        for m_, i_, lab in sorted(ms.peaks, key=lambda p: -p[1])[:8]:
            fig.add_annotation(x=m_, y=i_, text=str(m_), showarrow=False, yshift=10, font=dict(size=11, color=SIM))
    show_chart(fig, "ms", (0, ms.M + 10), (-110 if mirror else 0, 115))
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
        if nucleus == "1H":
            st.toggle(T("show_integral"), key=keep("integral"), on_change=save("integral"))
    signals = nmr.predict_h(A.mol) if nucleus == "1H" else nmr.predict_c(A.mol)
    x, y = nmr.spectrum_h(signals, KEEP["mhz"]) if nucleus == "1H" else nmr.spectrum_c(signals)
    sim = [go.Scatter(x=x, y=y, name=f"{T('nmr_sim')} {'¹H' if nucleus == '1H' else '¹³C'}", line=dict(color=SIM, width=1.3))]
    if nucleus == "1H" and KEEP["integral"] and signals:
        order = np.argsort(-x)
        cum = np.cumsum(y[order]) / max(np.sum(y), 1e-9)
        sim.append(go.Scatter(x=x[order], y=0.25 + 0.7 * cum, name=T("integral_curve"), line=dict(color="#16a34a", width=1.2, dash="dot")))
    exp_tr = [go.Scatter(x=ux, y=uy, name=nm, line=dict(color=USER_COLORS[k % 6], width=1.2)) for k, (nm, ux, uy) in enumerate(user_series("nmr"))]
    st.caption(T("nmr_exp_note"))
    fig = stacked_or_overlay(sim, exp_tr, T("shift"), T("rel_int").replace(" (%)", ""), reverse_x=True)
    if show_labels:
        sub = {} if overlay or not exp_tr else {"row": 1, "col": 1}
        for sgn in signals:
            yv = float(np.interp(sgn.shift, x, y))
            fig.add_annotation(x=sgn.shift, y=min(yv, 1.0), text=f"<b>{sgn.label}</b>", showarrow=False, yshift=12, font=dict(size=12, color=SIM), **sub)
    show_chart(fig, "nmr_" + nucleus, (-0.5, 13) if nucleus == "1H" else (-5, 230), (-0.03, 1.15), reverse=True)

    left, right = st.columns([1, 1.6], gap="large")
    with left:
        mol = Chem.Mol(A.mol)
        for sgn in signals:
            for i in sgn.atoms:
                mol.GetAtomWithIdx(i).SetProp("atomNote", sgn.label)
        d2 = rdMolDraw2D.MolDraw2DCairo(480, 380)
        d2.drawOptions().annotationFontScale = 0.9
        d2.DrawMolecule(mol)
        d2.FinishDrawing()
        st.image(d2.GetDrawingText(), use_container_width=True)
    with right:
        if nucleus == "1H":
            rows = [(s_.label, f"{s_.shift:.2f}", s_.n, s_.mult + (f" (J ≈ {s_.J:.0f} Hz)" if s_.mult not in ("s", "br s", "m") else ""), L(s_.desc)) for s_ in signals]
            df = pd.DataFrame(rows, columns=[T("signal"), T("shift"), T("integral"), T("mult"), T("type")])
        else:
            rows = [(s_.label, f"{s_.shift:.1f}", s_.n, s_.mult, L(s_.desc)) for s_ in signals]
            df = pd.DataFrame(rows, columns=[T("signal"), T("shift"), T("n_c"), "DEPT", T("type")])
        st.dataframe(df, hide_index=True, use_container_width=True)
        st.caption(T("solvent_note"))
    download(pd.DataFrame({"ppm": x[::10], "intensity": y[::10]}), f"nmr_{nucleus}_simulated.csv")

# ---------------------------------------------------------------------------
# Compound information and nomenclature
# ---------------------------------------------------------------------------
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
