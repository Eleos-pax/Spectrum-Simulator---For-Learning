"""Spectrum Simulator – For Learning
Simulador didáctico de espectros UV-Vis, IR y EI-MS con comparación frente a datos
experimentales (NIST Chemistry WebBook) e identificación del compuesto en PubChem.

Ejecutar:  streamlit run app.py
"""
from __future__ import annotations

import io

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

from spectra import databases as db
from spectra.simulate import analyze

st.set_page_config(page_title="Spectrum Simulator – For Learning", page_icon="🧪", layout="wide")

SIM = "#2563eb"   # simulado
EXP = "#dc2626"   # experimental
EXAMPLES = {
    "Acetofenona": "CC(=O)c1ccccc1",
    "Etanol": "CCO",
    "Ácido acético": "CC(=O)O",
    "Tolueno": "Cc1ccccc1",
    "2-Hexanona": "CCCCC(C)=O",
    "Benzaldehído": "O=Cc1ccccc1",
    "Acetato de etilo": "CCOC(C)=O",
    "Clorobenceno": "Clc1ccccc1",
    "Bromoetano": "CCBr",
    "Anilina": "Nc1ccccc1",
    "Nitrobenceno": "O=[N+]([O-])c1ccccc1",
    "Naftaleno": "c1ccc2ccccc2c1",
    "Acrilonitrilo": "C=CC#N",
    "β-Caroteno": "CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(C)CCCC2(C)C",
}

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
# Datos en caché
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
# Barra lateral: entrada de la estructura
# ---------------------------------------------------------------------------
st.sidebar.title("🧪 Spectrum Simulator")
st.sidebar.caption("UV-Vis · IR · Masas — para aprender")
mode = st.sidebar.radio("¿Cómo quieres indicar el compuesto?", ["Buscar por nombre", "Dibujar estructura", "Escribir SMILES / InChI", "Ejemplos de clase"], index=3)

smiles_in = None
compound = None
if mode == "Buscar por nombre":
    name = st.sidebar.text_input("Nombre (español o inglés), CAS o InChIKey", "acetophenone")
    kind = "inchikey" if len(name) == 27 and name.count("-") == 2 else "name"
    if name:
        with st.spinner("Consultando PubChem…"):
            compound = cached_pubchem(name.strip(), kind)
        if compound:
            smiles_in = compound.smiles
        else:
            st.sidebar.error("PubChem no encontró ese nombre. Prueba el nombre en inglés o dibuja la estructura.")
elif mode == "Dibujar estructura":
    try:
        from streamlit_ketcher import st_ketcher
        st.sidebar.info("Dibuja en el editor del panel principal y pulsa **Apply**.")
        with st.expander("✏️ Editor de estructuras (Ketcher)", expanded=True):
            smiles_in = st_ketcher(st.session_state.get("drawn", "CC(=O)c1ccccc1"), height=430)
            if smiles_in:
                st.session_state["drawn"] = smiles_in
    except ImportError:
        st.sidebar.error("Falta el paquete streamlit-ketcher (pip install streamlit-ketcher).")
elif mode == "Escribir SMILES / InChI":
    smiles_in = st.sidebar.text_area("SMILES, InChI o Molfile", "CC(=O)c1ccccc1", height=110)
else:
    ex = st.sidebar.selectbox("Compuesto", list(EXAMPLES))
    smiles_in = EXAMPLES[ex]

st.sidebar.divider()
want_exp = st.sidebar.toggle("Buscar espectros experimentales (NIST)", value=True)
overlay = st.sidebar.radio("Comparación", ["Superpuestos", "Uno sobre otro"], horizontal=True)
show_labels = st.sidebar.toggle("Etiquetar picos", value=True)

if not smiles_in:
    st.info("Indica un compuesto en la barra lateral para simular sus espectros.")
    st.stop()

try:
    A = cached_analysis(smiles_in)
except ValueError as e:
    st.error(str(e))
    st.stop()

# Identificación en PubChem (si el compuesto existe)
if compound is None and A.inchikey:
    with st.spinner("Buscando el compuesto en PubChem…"):
        compound = cached_pubchem(A.inchikey, "inchikey")

exp = {"IR": [], "Mass": None, "UVVis": None}
if want_exp and compound and compound.cas:
    with st.spinner("Descargando espectros experimentales del NIST WebBook…"):
        exp = cached_nist(compound.cas)

# ---------------------------------------------------------------------------
# Cabecera
# ---------------------------------------------------------------------------
c1, c2 = st.columns([1, 2.2], gap="large")
with c1:
    st.image(Draw.MolToImage(A.mol, size=(420, 320)), use_container_width=True)
with c2:
    title = compound.name if compound else "Compuesto dibujado (no encontrado en PubChem)"
    st.subheader(title)
    m1, m2, m3, m4 = st.columns(4)
    m1.markdown(f"**Fórmula**<br><span style='font-size:1.4rem'>{A.formula_html}</span>", unsafe_allow_html=True)
    m2.metric("Masa molar (g/mol)", f"{A.mw:.2f}")
    m3.metric("Masa exacta", f"{A.exact_mass:.4f}")
    m4.metric("Insaturaciones (DBE)", f"{A.dbe:g}")
    st.code(A.smiles, language=None)
    links = []
    if compound:
        links.append(f"[PubChem CID {compound.cid}]({compound.url})")
        if compound.cas:
            links.append(f"CAS {compound.cas} · [NIST WebBook]({db.nist_page(compound.cas)})")
    st.markdown(" · ".join(links) if links else "<span class='note'>Estructura nueva o ejercicio: solo se muestra el espectro simulado.</span>", unsafe_allow_html=True)
    n_exp = len(exp["IR"]) + (exp["Mass"] is not None) + (exp["UVVis"] is not None)
    if want_exp and compound and compound.cas:
        st.markdown(f"Espectros experimentales encontrados: **IR {len(exp['IR'])} · MS {int(exp['Mass'] is not None)} · UV-Vis {int(exp['UVVis'] is not None)}**")
    legend = "<span class='legend-sim'>━ Simulado</span>"
    if n_exp:
        legend += " &nbsp; <span class='legend-exp'>━ Experimental (NIST)</span>"
    st.markdown(legend, unsafe_allow_html=True)

tab_ir, tab_uv, tab_ms, tab_help = st.tabs(["📈 Infrarrojo (IR)", "🌈 UV-Visible", "⚛️ Espectrometría de masas", "📘 Cómo se simula"])


def stacked_or_overlay(traces_sim, traces_exp, xtitle, ytitle, reverse_x=False):
    if overlay == "Superpuestos" or not traces_exp:
        fig = go.Figure(traces_sim + traces_exp)
        return base_layout(fig, xtitle, ytitle, reverse_x)
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=.06, subplot_titles=("Simulado", "Experimental (NIST)"))
    for t in traces_sim:
        fig.add_trace(t, 1, 1)
    for t in traces_exp:
        fig.add_trace(t, 2, 1)
    base_layout(fig, None, ytitle, reverse_x, height=620)
    fig.update_xaxes(autorange="reversed" if reverse_x else True)
    fig.update_xaxes(title_text=xtitle, row=2, col=1)
    return fig


def download(df, name, label="Descargar datos (CSV)"):
    st.download_button(label, df.to_csv(index=False).encode(), file_name=name, mime="text/csv")


# ---------------------------------------------------------------------------
# IR
# ---------------------------------------------------------------------------
with tab_ir:
    ir = A.ir
    sim = [go.Scatter(x=ir.x, y=ir.y, name="IR simulado", line=dict(color=SIM, width=1.6))]
    exp_tr = []
    if exp["IR"]:
        opts = [f"{i + 1}. {s.meta.get('STATE', 'estado n.d.')} — {s.meta.get('ORIGIN', '')[:50]}" for i, s in enumerate(exp["IR"])]
        pick = st.selectbox("Espectro experimental", range(len(opts)), format_func=lambda i: opts[i])
        ex, ey = db.ir_to_transmittance(exp["IR"][pick])
        exp_tr = [go.Scatter(x=ex, y=ey, name="IR experimental", line=dict(color=EXP, width=1.3))]
    fig = stacked_or_overlay(sim, exp_tr, "Número de onda (cm⁻¹)", "Transmitancia (%)", reverse_x=True)
    if show_labels:
        for c, w, i, lab in ir.peaks[:14]:
            yv = float(np.interp(c, ir.x[::-1], ir.y[::-1]))
            fig.add_annotation(x=c, y=yv, text=f"{c:.0f}", showarrow=True, arrowhead=0, ay=28, ax=0, font=dict(size=10, color=SIM), row=1 if overlay != "Superpuestos" and exp_tr else None, col=1 if overlay != "Superpuestos" and exp_tr else None)
    fig.update_xaxes(range=[4000, 400])
    st.plotly_chart(fig, use_container_width=True)
    st.markdown("**Grupos funcionales detectados:** " + " ".join(f"<span class='chip'>{g}</span>" for g in ir.groups), unsafe_allow_html=True)
    df = pd.DataFrame([(round(c), lab, "fuerte" if i > .7 else "media" if i > .4 else "débil") for c, w, i, lab in ir.peaks], columns=["ν (cm⁻¹)", "Asignación", "Intensidad"])
    st.dataframe(df, hide_index=True, use_container_width=True)
    download(pd.DataFrame({"wavenumber_cm-1": ir.x, "transmittance_%": ir.y}), "ir_simulado.csv")

# ---------------------------------------------------------------------------
# UV-Vis
# ---------------------------------------------------------------------------
with tab_uv:
    uv = A.uv
    log_scale = st.toggle("Eje de absortividad en escala log ε", value=True)
    ytr = (lambda v: np.log10(np.maximum(v, 1))) if log_scale else (lambda v: v)
    sim = [go.Scatter(x=uv.x, y=ytr(uv.y), name="UV-Vis simulado", line=dict(color=SIM, width=2), fill="tozeroy", fillcolor="rgba(37,99,235,.08)")]
    exp_tr = []
    if exp["UVVis"] is not None:
        ux, uy, unit = db.uv_to_epsilon(exp["UVVis"])
        yy = ytr(uy) if unit == "ε" else uy
        exp_tr = [go.Scatter(x=ux, y=yy, name=f"UV-Vis experimental ({unit})", line=dict(color=EXP, width=1.6))]
    fig = stacked_or_overlay(sim, exp_tr, "Longitud de onda (nm)", "log ε" if log_scale else "ε (L·mol⁻¹·cm⁻¹)")
    # franja del espectro visible
    for lo, hi, col in [(380, 450, "#7c3aed"), (450, 495, "#2563eb"), (495, 570, "#16a34a"), (570, 590, "#eab308"), (590, 620, "#f97316"), (620, 750, "#dc2626")]:
        fig.add_vrect(x0=lo, x1=hi, fillcolor=col, opacity=.07, line_width=0)
    if show_labels:
        for l, e, lab in uv.bands:
            if l >= 180:
                fig.add_annotation(x=l, y=float(ytr(np.array([max(e, 1)]))[0]) if not log_scale else float(np.log10(max(np.interp(l, uv.x, uv.y), 1))),
                                   text=f"{l:.0f} nm", showarrow=True, ay=-26, ax=0, font=dict(size=10, color=SIM))
    fig.update_xaxes(range=[180, 750])
    st.plotly_chart(fig, use_container_width=True)
    cc1, cc2 = st.columns([3, 1])
    with cc1:
        st.markdown(f"**λmax simulada (≥ 200 nm):** {uv.lambda_max:.0f} nm")
        df = pd.DataFrame([(round(l), f"{e:,.0f}", f"{np.log10(max(e, 1)):.2f}", lab) for l, e, lab in uv.bands], columns=["λ (nm)", "ε", "log ε", "Transición"])
        st.dataframe(df, hide_index=True, use_container_width=True)
    with cc2:
        st.markdown("**Color aparente**")
        st.markdown(f"<div style='width:100%;height:90px;border-radius:12px;border:1px solid rgba(128,128,128,.4);background:{uv.color}'></div>", unsafe_allow_html=True)
        st.caption("Luz blanca menos lo que absorbe la molécula (estimación).")
    if exp["UVVis"] is not None:
        st.caption("Experimental: " + exp["UVVis"].describe())
    download(pd.DataFrame({"wavelength_nm": uv.x, "epsilon": uv.y}), "uv_simulado.csv")

# ---------------------------------------------------------------------------
# MS
# ---------------------------------------------------------------------------
with tab_ms:
    ms = A.ms
    mz = [p[0] for p in ms.peaks]
    it = [p[1] for p in ms.peaks]

    def sticks(x, y, color, name, neg=False):
        xs, ys = [], []
        for a, b in zip(x, y):
            xs += [a, a, None]
            ys += [0, -b if neg else b, None]
        return go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=color, width=2.2), name=name, hoverinfo="skip")

    hover_sim = go.Scatter(x=mz, y=it, mode="markers", marker=dict(size=6, color=SIM), name="m/z simulado",
                           text=[p[2] for p in ms.peaks], hovertemplate="m/z %{x}<br>%{y:.0f}%<br>%{text}<extra></extra>", showlegend=False)
    traces = [sticks(mz, it, SIM, "EI-MS simulado"), hover_sim]
    mirror = False
    if exp["Mass"] is not None:
        mirror = st.toggle("Vista espejo (experimental hacia abajo)", value=True)
        e = exp["Mass"]
        ey = 100 * e.y / e.y.max()
        traces.append(sticks(e.x, ey, EXP, "EI-MS experimental (NIST)", neg=mirror))
        traces.append(go.Scatter(x=e.x, y=-ey if mirror else ey, mode="markers", marker=dict(size=4, color=EXP), showlegend=False,
                                 hovertemplate="m/z %{x}<br>%{customdata:.0f}% (NIST)<extra></extra>", customdata=ey))
    fig = go.Figure(traces)
    base_layout(fig, "m/z", "Intensidad relativa (%)")
    fig.update_layout(hovermode="closest")
    fig.update_yaxes(range=[-110 if mirror else 0, 115])
    if show_labels:
        for m_, i_, lab in sorted(ms.peaks, key=lambda p: -p[1])[:8]:
            fig.add_annotation(x=m_, y=i_, text=str(m_), showarrow=False, yshift=10, font=dict(size=11, color=SIM))
    st.plotly_chart(fig, use_container_width=True)
    a, b = st.columns(2)
    with a:
        st.markdown(f"**Ion molecular M⁺·:** m/z {ms.M} · regla del nitrógeno: masa {'impar' if ms.M % 2 else 'par'}")
        st.markdown("**Patrón isotópico del ion molecular**")
        st.dataframe(pd.DataFrame([(m_, f"{r:.1f}") for m_, r in ms.molecular_ion], columns=["m/z", "% relativo"]), hide_index=True, use_container_width=True)
    with b:
        st.markdown("**Fragmentos principales (simulados)**")
        top = sorted(ms.peaks, key=lambda p: -p[1])[:12]
        st.dataframe(pd.DataFrame([(m_, f"{i_:.0f}", str(ms.M - m_) if m_ < ms.M else "—", lab) for m_, i_, lab in top], columns=["m/z", "%", "Pérdida", "Interpretación"]), hide_index=True, use_container_width=True)
    download(pd.DataFrame(ms.peaks, columns=["mz", "rel_intensity", "assignment"]), "ms_simulado.csv")

with tab_help:
    st.markdown(
        """
### Qué hace cada simulación
**Infrarrojo.** Busca grupos funcionales con patrones SMARTS (RDKit) y asigna bandas de una tabla de
frecuencias de grupo (C=O, O–H, N–H, C≡N, NO₂…). Las bandas C–H dependen de la hibridación (sp³ < 3000 cm⁻¹,
sp²/aromático > 3000 cm⁻¹). La conjugación baja la frecuencia del C=O unos 25 cm⁻¹. Cada banda es una
lorentziana; la suma se convierte a % de transmitancia. La región de huella (< 1500 cm⁻¹) es solo ilustrativa.

**UV-Visible.** Identifica los sistemas conjugados: número de enlaces π, anillos aromáticos fusionados y
auxocromos (–OH, –OR, –NH₂, halógenos, alquilos). Más conjugación desplaza λmax hacia el rojo
(efecto batocrómico), como en las reglas de Woodward-Fieser. Incluye transiciones n→π* débiles de C=O, NO₂ y N=N.
Si λmax entra en el visible, se estima el color que verías.

**Espectrometría de masas (EI, 70 eV).** Calcula el ion molecular con su patrón isotópico exacto
(M+1 por ¹³C, M+2 por ³⁷Cl/⁸¹Br/³⁴S) y fragmentos por reglas de estabilidad: ruptura α a carbonilo
(acilio), α a N/O (iminio, oxonio), ruptura bencílica (tropilio m/z 91), reordenamiento de McLafferty,
carbocationes terciarios > secundarios > primarios y pérdidas neutras (H₂O, CO, ·OCH₃, ·NO₂…).

### Datos experimentales
Si el compuesto existe, se identifica en **PubChem** por su InChIKey, se obtiene su número CAS y se
descargan sus espectros del **NIST Chemistry WebBook** (formato JCAMP-DX). Si dibujas una molécula que no está
en las bases de datos, se muestra solo el espectro simulado.

> Los espectros simulados son aproximaciones empíricas para aprender a interpretar. Para trabajo de
> investigación usa espectros medidos o cálculos de química cuántica (DFT/TD-DFT).
"""
    )

st.caption("Datos: PubChem (NCBI) y NIST Chemistry WebBook, SRD 69. Química computacional: RDKit. Gráficos: Plotly.")
