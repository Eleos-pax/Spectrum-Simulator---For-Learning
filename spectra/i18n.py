"""Bilingual (Spanish / English) user-interface text.

`UI` holds every interface string as {key: (es, en)}.
`label()` translates the band / fragment assignments produced by `simulate.py`
(written in Spanish) into English using an ordered phrase table.
"""
from __future__ import annotations

LANGS = {"es": "Español", "en": "English"}

UI = {
    "subtitle": ("UV-Vis · IR · Masas — para aprender", "UV-Vis · IR · Mass spec — for learning"),
    "language": ("Idioma", "Language"),
    "how": ("¿Cómo quieres indicar el compuesto?", "How do you want to enter the compound?"),
    "m_name": ("Buscar por nombre", "Search by name"),
    "m_draw": ("Dibujar estructura", "Draw structure"),
    "m_smiles": ("Escribir SMILES / InChI", "Type SMILES / InChI"),
    "m_extract": ("Extracto / muestra (solo mis datos)", "Extract / sample (my data only)"),
    "extract_name": ("Nombre del extracto o muestra (sale en la gráfica)", "Extract or sample name (shown on the chart)"),
    "extract_default": ("Extracto sin nombre", "Unnamed extract"),
    "extract_hint": ("Carga los archivos de tus equipos en cada pestaña (IR, UV-Vis, MS, RMN). Para la elucidación puedes superponer compuestos de referencia simulados.",
                     "Upload your instrument files in each tab (IR, UV-Vis, MS, NMR). For structure elucidation you can overlay simulated reference compounds."),
    "refs": ("🧪 Compuestos de referencia", "🧪 Reference compounds"),
    "refs_help": ("Se simulan y se dibujan con línea discontinua para compararlos con tus datos.", "They are simulated and drawn as dashed lines to compare with your data."),
    "ref_input": ("Nombre o SMILES", "Name or SMILES"),
    "add": ("Agregar", "Add"),
    "ref_error": ("No se encontró «{q}».", "“{q}” was not found."),
    "no_info_extract": ("Modo extracto: no hay una única estructura. Estos son los compuestos de referencia que agregaste:",
                        "Extract mode: there is no single structure. These are the reference compounds you added:"),
    "no_refs": ("Aún no hay compuestos de referencia. Agrégalos en la barra lateral.", "No reference compounds yet. Add them in the sidebar."),
    "raw_abs": ("Absorbancia (datos del equipo)", "Absorbance (instrument data)"),
    "peaks": ("📍 Marcar picos / λmax", "📍 Mark peaks / λmax"),
    "mark_on": ("Marcar picos automáticamente", "Mark peaks automatically"),
    "mark_in": ("Marcar en", "Mark on"),
    "k_sim": ("Simulado", "Simulated"), "k_nist": ("NIST", "NIST"), "k_user": ("Mis datos", "My data"), "k_ref": ("Referencias", "References"),
    "n_peaks": ("Picos por serie", "Peaks per series"),
    "prom": ("Prominencia mínima (%)", "Minimum prominence (%)"),
    "picked": ("Picos marcados", "Marked peaks"),
    "series": ("Serie", "Series"),
    "position": ("Posición", "Position"),
    "style": ("🎨 Estilo de la gráfica", "🎨 Chart style"),
    "chart_title": ("Título (vacío = nombre del compuesto o extracto)", "Title (empty = compound or extract name)"),
    "font": ("Tamaño de letra", "Font size"),
    "lw": ("Grosor de línea", "Line width"),
    "grid": ("Cuadrícula", "Grid"),
    "legend_pos": ("Leyenda", "Legend"),
    "lg_top": ("Arriba", "Top"), "lg_right": ("Derecha", "Right"), "lg_none": ("Oculta", "Hidden"),
    "height": ("Altura (px)", "Height (px)"),
    "fmt": ("Formato al guardar (icono de cámara)", "Format when saving (camera icon)"),
    "visband": ("Franja de colores del visible", "Visible colour band"),
    "m_example": ("Ejemplos de clase", "Class examples"),
    "name_input": ("Nombre (en inglés funciona mejor), CAS o InChIKey", "Name, CAS number or InChIKey"),
    "search": ("Buscar", "Search"),
    "querying": ("Consultando PubChem…", "Querying PubChem…"),
    "not_found": ("PubChem no encontró «{q}». Prueba el nombre en inglés o dibuja la estructura.",
                  "PubChem did not find “{q}”. Check the spelling or draw the structure."),
    "draw_hint": ("Dibuja en el editor y pulsa **Apply** para simular.", "Draw in the editor and press **Apply** to simulate."),
    "editor": ("✏️ Editor de estructuras (Ketcher)", "✏️ Structure editor (Ketcher)"),
    "no_ketcher": ("Falta el paquete streamlit-ketcher (pip install streamlit-ketcher).", "The streamlit-ketcher package is missing (pip install streamlit-ketcher)."),
    "smiles_input": ("SMILES, InChI o Molfile", "SMILES, InChI or Molfile"),
    "apply": ("Simular", "Simulate"),
    "compound": ("Compuesto", "Compound"),
    "want_exp": ("Buscar espectros experimentales (NIST)", "Fetch experimental spectra (NIST)"),
    "compare": ("Comparación", "Comparison"),
    "overlay": ("Superpuestos", "Overlaid"),
    "stacked": ("Uno sobre otro", "Stacked"),
    "labels": ("Etiquetar picos", "Label peaks"),
    "current": ("Estructura actual", "Current structure"),
    "empty": ("Indica un compuesto en la barra lateral para simular sus espectros.", "Enter a compound in the sidebar to simulate its spectra."),
    "bad_structure": ("No se pudo interpretar la estructura (SMILES, Molfile o InChI no válido).", "The structure could not be read (invalid SMILES, Molfile or InChI)."),
    "finding": ("Buscando el compuesto en PubChem…", "Looking the compound up in PubChem…"),
    "downloading": ("Descargando espectros experimentales del NIST WebBook…", "Downloading experimental spectra from the NIST WebBook…"),
    "unknown": ("Compuesto dibujado (no encontrado en PubChem)", "Drawn compound (not found in PubChem)"),
    "formula": ("Fórmula", "Formula"),
    "molar_mass": ("Masa molar (g/mol)", "Molar mass (g/mol)"),
    "exact_mass": ("Masa exacta", "Exact mass"),
    "dbe": ("Insaturaciones (DBE)", "Degrees of unsaturation"),
    "new_structure": ("Estructura nueva o ejercicio: solo se muestra el espectro simulado.", "New or exercise structure: only the simulated spectrum is shown."),
    "found": ("Espectros experimentales encontrados", "Experimental spectra found"),
    "no_cas": ("El compuesto existe en PubChem pero no tiene nº CAS: no se pueden buscar espectros en NIST.", "The compound is in PubChem but has no CAS number, so NIST spectra cannot be looked up."),
    "sim": ("Simulado", "Simulated"),
    "exp": ("Experimental (NIST)", "Experimental (NIST)"),
    "tab_ir": ("📈 Infrarrojo (IR)", "📈 Infrared (IR)"),
    "tab_uv": ("🌈 UV-Visible", "🌈 UV-Visible"),
    "tab_ms": ("⚛️ Espectrometría de masas", "⚛️ Mass spectrometry"),
    "tab_help": ("📘 Cómo se simula", "📘 How it is simulated"),
    "section": ("Espectro", "Spectrum"),
    "download": ("Descargar datos (CSV)", "Download data (CSV)"),
    # IR
    "ir_sim": ("IR simulado", "Simulated IR"),
    "ir_exp": ("IR experimental", "Experimental IR"),
    "exp_pick": ("Espectro experimental", "Experimental spectrum"),
    "state_nd": ("estado n.d.", "state n/a"),
    "wavenumber": ("Número de onda (cm⁻¹)", "Wavenumber (cm⁻¹)"),
    "transmittance": ("Transmitancia (%)", "Transmittance (%)"),
    "groups": ("Grupos funcionales detectados", "Functional groups detected"),
    "assignment": ("Asignación", "Assignment"),
    "intensity": ("Intensidad", "Intensity"),
    "strong": ("fuerte", "strong"), "medium": ("media", "medium"), "weak": ("débil", "weak"),
    # UV
    "uv_sim": ("UV-Vis simulado", "Simulated UV-Vis"),
    "uv_exp": ("UV-Vis experimental", "Experimental UV-Vis"),
    "log_axis": ("Eje de absortividad en escala log ε", "Absorptivity axis as log ε"),
    "wavelength": ("Longitud de onda (nm)", "Wavelength (nm)"),
    "lmax": ("λmax simulada (≥ 200 nm)", "Simulated λmax (≥ 200 nm)"),
    "transition": ("Transición", "Transition"),
    "color": ("Color aparente", "Apparent colour"),
    "color_note": ("Luz blanca menos lo que absorbe la molécula (estimación).", "White light minus what the molecule absorbs (estimate)."),
    "absorbance": ("absorbancia", "absorbance"),
    # MS
    "ms_sim": ("EI-MS simulado", "Simulated EI-MS"),
    "ms_exp": ("EI-MS experimental (NIST)", "Experimental EI-MS (NIST)"),
    "mirror": ("Vista espejo (experimental hacia abajo)", "Mirror view (experimental pointing down)"),
    "rel_int": ("Intensidad relativa (%)", "Relative intensity (%)"),
    "mol_ion": ("Ion molecular M⁺·", "Molecular ion M⁺·"),
    "n_rule": ("regla del nitrógeno: masa {p}", "nitrogen rule: {p} mass"),
    "odd": ("impar", "odd"), "even": ("par", "even"),
    "iso_pattern": ("Patrón isotópico del ion molecular", "Isotope pattern of the molecular ion"),
    "rel_pct": ("% relativo", "relative %"),
    "fragments": ("Fragmentos principales (simulados)", "Main fragments (simulated)"),
    "loss": ("Pérdida", "Loss"),
    "interpretation": ("Interpretación", "Interpretation"),
    # user files
    "my_data": ("📂 Mis datos del equipo", "📂 My instrument data"),
    "upload": ("Carga archivos exportados por el equipo (.csv, .txt, .dx/.jdx)", "Upload files exported by the instrument (.csv, .txt, .dx/.jdx)"),
    "upload_help": ("Se aceptan separadores , ; tabulador o espacios, coma o punto decimal y líneas de encabezado. Puedes cargar varios archivos y varias columnas Y.",
                    "Separators , ; tab or spaces, decimal comma or point and header lines are all accepted. You can load several files and several Y columns."),
    "read_error": ("No se pudo leer «{f}»: no se encontró una tabla de números.", "Could not read “{f}”: no numeric table was found."),
    "x_col": ("Columna X", "X column"),
    "y_cols": ("Columnas Y", "Y columns"),
    "y_type": ("Los valores Y son", "Y values are"),
    "remove": ("Quitar", "Remove"),
    "show": ("Mostrar", "Show"),
    "conc": ("c (mol/L) para calcular ε; 0 = normalizar", "c (mol/L) to compute ε; 0 = normalise"),
    "path": ("b (cm)", "b (cm)"),
    "yt_abs": ("Absorbancia", "Absorbance"),
    "yt_T": ("Transmitancia (%)", "Transmittance (%)"),
    "yt_T1": ("Transmitancia (0–1)", "Transmittance (0–1)"),
    "yt_eps": ("ε (L·mol⁻¹·cm⁻¹)", "ε (L·mol⁻¹·cm⁻¹)"),
    "yt_logeps": ("log ε", "log ε"),
    "yt_nmr": ("Intensidad (ppm en X)", "Intensity (ppm on X)"),
    "yt_int": ("Intensidad (cualquier escala)", "Intensity (any scale)"),
    "user_note_uv": ("Sin concentración, la absorbancia se escala al máximo del espectro simulado para comparar la forma.",
                     "Without a concentration, absorbance is scaled to the maximum of the simulated spectrum to compare band shapes."),
    "exp_any": ("Experimental", "Experimental"),
    # general information
    "tab_info": ("🧾 Compuesto", "🧾 Compound"),
    "names": ("Nomenclatura", "Nomenclature"),
    "iupac": ("Nombre IUPAC", "IUPAC name"),
    "common": ("Nombre común", "Common name"),
    "synonyms": ("Otros nombres y sinónimos", "Other names and synonyms"),
    "no_name": ("Estructura nueva: no figura en PubChem, así que no hay nombre IUPAC registrado. Usa la fórmula, el SMILES y el InChI como identificadores y propón el nombre con las reglas IUPAC.",
                "New structure: it is not in PubChem, so there is no registered IUPAC name. Use the formula, SMILES and InChI as identifiers and work the name out with the IUPAC rules."),
    "identifiers": ("Identificadores", "Identifiers"),
    "properties": ("Propiedades generales (calculadas con RDKit)", "General properties (computed with RDKit)"),
    "pubchem_props": ("Propiedades calculadas por PubChem", "Properties computed by PubChem"),
    "property": ("Propiedad", "Property"),
    "value": ("Valor", "Value"),
    "p_heavy": ("Átomos pesados", "Heavy atoms"), "p_atoms": ("Átomos totales (con H)", "Total atoms (with H)"),
    "p_rings": ("Anillos", "Rings"), "p_arom": ("Anillos aromáticos", "Aromatic rings"),
    "p_logp": ("logP (Crippen)", "logP (Crippen)"), "p_mr": ("Refractividad molar", "Molar refractivity"),
    "p_tpsa": ("Área polar (TPSA, Å²)", "Polar surface area (TPSA, Å²)"),
    "p_hbd": ("Donadores de puente H", "H-bond donors"), "p_hba": ("Aceptores de puente H", "H-bond acceptors"),
    "p_rot": ("Enlaces rotables", "Rotatable bonds"), "p_fsp3": ("Fracción sp³ de C", "Fraction of sp³ carbons"),
    "p_stereo": ("Estereocentros", "Stereocentres"), "p_charge": ("Carga formal", "Formal charge"),
    "p_lipinski": ("Regla de Lipinski (violaciones)", "Lipinski rule (violations)"),
    "p_groups": ("Grupos funcionales", "Functional groups"),
    "p_composition": ("Composición elemental (% en masa)", "Elemental composition (mass %)"),
    # NMR
    "tab_nmr": ("🧲 RMN", "🧲 NMR"),
    "nucleus": ("Núcleo", "Nucleus"),
    "field": ("Campo del espectrómetro (MHz, ¹H)", "Spectrometer field (MHz, ¹H)"),
    "solvent_note": ("Desplazamientos estimados en CDCl₃ por incrementos aditivos (±0.3 ppm en ¹H, ±5 ppm en ¹³C). Los H de OH/NH/COOH son intercambiables: su posición depende del disolvente y la concentración.",
                     "Shifts estimated in CDCl₃ with additive increments (±0.3 ppm for ¹H, ±5 ppm for ¹³C). OH/NH/COOH protons are exchangeable: their position depends on solvent and concentration."),
    "nmr_sim": ("RMN simulada", "Simulated NMR"),
    "shift": ("δ (ppm)", "δ (ppm)"),
    "integral": ("Integral (nH)", "Integral (nH)"),
    "mult": ("Multiplicidad", "Multiplicity"),
    "type": ("Tipo", "Type"),
    "n_c": ("nº C equivalentes", "Equivalent C"),
    "signal": ("Señal", "Signal"),
    "show_integral": ("Mostrar curva de integración", "Show integration curve"),
    "integral_curve": ("Integral", "Integral"),
    "nmr_exp_note": ("NIST no publica espectros de RMN; para comparar, carga el archivo exportado por tu equipo (ppm, intensidad).",
                     "NIST does not publish NMR spectra; to compare, upload the file exported by your instrument (ppm, intensity)."),
    # axes
    "axes": ("🔍 Ejes y zoom", "🔍 Axes and zoom"),
    "auto_range": ("Rango automático", "Automatic range"),
    "zoom_help": ("Rueda del ratón: acercar/alejar · arrastrar: zoom de caja · doble clic: restablecer · cámara: guardar imagen PNG.",
                  "Mouse wheel: zoom in/out · drag: box zoom · double click: reset · camera icon: save PNG image."),
    "footer": ("Datos: PubChem (NCBI) y NIST Chemistry WebBook, SRD 69. Química computacional: RDKit. Gráficos: Plotly.",
               "Data: PubChem (NCBI) and NIST Chemistry WebBook, SRD 69. Cheminformatics: RDKit. Charts: Plotly."),
}

HELP = {
    "es": """
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

**RMN ¹H y ¹³C.** Agrupa los átomos químicamente equivalentes (simetría) y estima los desplazamientos con
incrementos aditivos: tablas de Shoolery para ¹H alifático, benceno (7.27 ppm) con incrementos orto/meta/para para
¹H aromático, reglas de Grant-Paul para ¹³C alifático y 128.5 ppm con incrementos ipso/orto/meta/para para ¹³C
aromático. La multiplicidad sigue la regla n+1 (primer orden) y la integral es el número de H de cada grupo.

### Datos experimentales
Si el compuesto existe, se identifica en **PubChem** por su InChIKey, se obtiene su número CAS y se
descargan sus espectros del **NIST Chemistry WebBook** (formato JCAMP-DX). Si dibujas una molécula que no está
en las bases de datos, se muestra solo el espectro simulado.

> Los espectros simulados son aproximaciones empíricas para aprender a interpretar. Para trabajo de
> investigación usa espectros medidos o cálculos de química cuántica (DFT/TD-DFT).
""",
    "en": """
### What each simulation does
**Infrared.** Functional groups are found with SMARTS patterns (RDKit) and bands are taken from a
group-frequency table (C=O, O–H, N–H, C≡N, NO₂…). C–H stretches depend on hybridisation (sp³ < 3000 cm⁻¹,
sp²/aromatic > 3000 cm⁻¹). Conjugation lowers the C=O frequency by about 25 cm⁻¹. Each band is a
Lorentzian; the sum is converted to % transmittance. The fingerprint region (< 1500 cm⁻¹) is only illustrative.

**UV-Visible.** Conjugated systems are identified: number of π bonds, fused aromatic rings and
auxochromes (–OH, –OR, –NH₂, halogens, alkyl groups). More conjugation shifts λmax to longer wavelengths
(bathochromic shift), as in the Woodward-Fieser rules. Weak n→π* transitions of C=O, NO₂ and N=N are included.
If λmax falls in the visible range, the colour you would see is estimated.

**Mass spectrometry (EI, 70 eV).** The molecular ion is computed with its exact isotope pattern
(M+1 from ¹³C, M+2 from ³⁷Cl/⁸¹Br/³⁴S), and fragments follow stability rules: α-cleavage next to a carbonyl
(acylium), α to N/O (iminium, oxonium), benzylic cleavage (tropylium m/z 91), McLafferty rearrangement,
tertiary > secondary > primary carbocations, and neutral losses (H₂O, CO, ·OCH₃, ·NO₂…).

**¹H and ¹³C NMR.** Chemically equivalent atoms are grouped by symmetry and shifts are estimated with
additive increments: Shoolery-type tables for aliphatic ¹H, benzene (7.27 ppm) with ortho/meta/para increments for
aromatic ¹H, Grant-Paul rules for aliphatic ¹³C and 128.5 ppm with ipso/ortho/meta/para increments for aromatic
¹³C. Multiplicity follows the n+1 rule (first order) and the integral is the number of H in each group.

### Experimental data
If the compound exists, it is identified in **PubChem** by its InChIKey, its CAS number is retrieved and
its spectra are downloaded from the **NIST Chemistry WebBook** (JCAMP-DX format). If you draw a molecule that
is not in the databases, only the simulated spectrum is shown.

> Simulated spectra are empirical approximations meant for learning to interpret spectra. For research,
> use measured spectra or quantum-chemical calculations (DFT/TD-DFT).
""",
}

EXAMPLES = [  # (es, en, SMILES)
    ("Acetofenona", "Acetophenone", "CC(=O)c1ccccc1"),
    ("Etanol", "Ethanol", "CCO"),
    ("Ácido acético", "Acetic acid", "CC(=O)O"),
    ("Tolueno", "Toluene", "Cc1ccccc1"),
    ("2-Hexanona", "2-Hexanone", "CCCCC(C)=O"),
    ("Benzaldehído", "Benzaldehyde", "O=Cc1ccccc1"),
    ("Acetato de etilo", "Ethyl acetate", "CCOC(C)=O"),
    ("Clorobenceno", "Chlorobenzene", "Clc1ccccc1"),
    ("Bromoetano", "Bromoethane", "CCBr"),
    ("Anilina", "Aniline", "Nc1ccccc1"),
    ("Nitrobenceno", "Nitrobenzene", "O=[N+]([O-])c1ccccc1"),
    ("Naftaleno", "Naphthalene", "c1ccc2ccccc2c1"),
    ("Acrilonitrilo", "Acrylonitrile", "C=CC#N"),
    ("β-Caroteno", "β-Carotene", "CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(C)CCCC2(C)C"),
]

# Spanish -> English phrases for simulation labels (longest phrases first when applied)
_PHRASES = {
    "C aromático": "aromatic C",
    # IR groups
    "Alcohol / fenol": "Alcohol / phenol", "Ácido carboxílico": "Carboxylic acid", "Éster": "Ester",
    "Lactona": "Lactone", "Aldehído": "Aldehyde", "Cetona": "Ketone", "Amida primaria": "Primary amide",
    "Amida secundaria": "Secondary amide", "Amida terciaria": "Tertiary amide", "Haluro de acilo": "Acyl halide",
    "Anhídrido": "Anhydride", "Amina primaria": "Primary amine", "Amina secundaria": "Secondary amine",
    "Nitrilo": "Nitrile", "Alquino terminal": "Terminal alkyne", "Alquino interno": "Internal alkyne",
    "Alqueno": "Alkene", "Vinilo terminal": "Terminal vinyl", "Anillo aromático": "Aromatic ring",
    "Benceno monosustituido": "Monosubstituted benzene", "Nitro": "Nitro", "Éter": "Ether",
    "Fluoruro": "Fluoride", "Cloruro": "Chloride", "Bromuro": "Bromide", "Yoduro": "Iodide", "Tiol": "Thiol",
    "Sulfona / sulfonamida": "Sulfone / sulfonamide", "Sulfóxido": "Sulfoxide", "Isocianato": "Isocyanate",
    # IR assignments
    "(puente de H)": "(H-bonded)", "ácido (muy ancha)": "acid (very broad)", "ácido": "acid",
    "aldehído (Fermi)": "aldehyde (Fermi)", "lactona": "lactone", "éster": "ester", "cetona": "ketone",
    "aldehído": "aldehyde", "haluro de acilo": "acyl halide", "amida I": "amide I", "amida II": "amide II",
    "amida": "amide", "asim.": "asym.", "sim.": "sym.", "f.p.": "o.o.p.", "tijera": "scissoring",
    "paraguas": "umbrella", "balanceo": "rocking", "(débil)": "(weak)", "aromático": "aromatic",
    "sobretonos aromáticos": "aromatic overtones", "deformación del anillo": "ring deformation",
    "(monosust.)": "(monosubst.)", "(conjugado)": "(conjugated)",
    # UV
    "(banda E₂/K)": "(E₂/K band)", "banda B (estructura fina)": "B band (fine structure)", "banda p/B": "p/B band",
    "(banda E₁)": "(E₁ band)", "enlaces π conjugados": "conjugated π bonds", "aislado": "isolated",
    "(prohibida, débil)": "(forbidden, weak)", "transferencia de carga": "charge transfer",
    "(< 200 nm, no observable en UV convencional)": "(< 200 nm, not observable in routine UV)",
    # MS
    "ion molecular": "molecular ion", "ion acilio (ruptura α)": "acylium ion (α-cleavage)",
    "ion iminio (ruptura α a N)": "iminium ion (α-cleavage to N)", "ion oxonio (ruptura α a O)": "oxonium ion (α-cleavage to O)",
    "catión bencílico / tropilio": "benzyl / tropylium cation", "catión alílico": "allyl cation",
    "carbocatión primario": "primary carbocation", "carbocatión secundario": "secondary carbocation",
    "carbocatión terciario": "tertiary carbocation", "pérdida de": "loss of", "reordenamiento de McLafferty": "McLafferty rearrangement",
    "(deshidratación)": "(dehydration)", "(aldehído)": "(aldehyde)", "(fenol)": "(phenol)",
    "catión fenilo": "phenyl cation", "(fenilo − C₂H₂)": "(phenyl − C₂H₂)", "(tropilio − C₂H₂)": "(tropylium − C₂H₂)",
    "ion tropilio": "tropylium ion", "ion benzoílo": "benzoyl ion", "serie alquílica": "alkyl series",
    "(isótopos)": "(isotopes)",
}
_ORDER = sorted(_PHRASES, key=len, reverse=True)


def label(text: str, lang: str) -> str:
    """Translate a simulation label into the chosen language (labels are generated in Spanish)."""
    if lang == "es" or not text:
        return text
    for es in _ORDER:
        if es in text:
            text = text.replace(es, _PHRASES[es])
    return text
