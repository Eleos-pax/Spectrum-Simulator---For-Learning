"""Publication-quality export of the app's charts.

The interactive charts are Plotly figures; browsers can only save them as screen images (or SVG).
Journals usually ask for TIFF at 600-1200 dpi or vector PDF/EPS at an exact column width with
embedded fonts. This module redraws a Plotly figure with matplotlib so it can be saved that way:

* exact size in centimetres (journal column widths) and resolution in dpi;
* TIFF (LZW-compressed), PNG, PDF, SVG and EPS; TrueType fonts embedded (pdf/ps fonttype 42),
  editable text in SVG;
* lines, stick spectra, dashed reference spectra, peak labels, shaded bands, subplots (stacked view),
  reversed axes and the numbered NMR structure image are all carried over.
"""
from __future__ import annotations

import base64
import html
import io
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

FORMATS = {  # key: (extension, mime type)
    "tiff": ("tiff", "image/tiff"),
    "png": ("png", "image/png"),
    "pdf": ("pdf", "application/pdf"),
    "svg": ("svg", "image/svg+xml"),
    "eps": ("eps", "application/postscript"),
}
WIDTHS_CM = {"col1": 8.5, "col15": 11.4, "col2": 17.4}
DASH = {"dash": (0, (5, 3)), "dot": (0, (1, 2)), "dashdot": (0, (5, 2, 1, 2)), None: "-", "solid": "-"}
_TAG = re.compile(r"<[^>]+>")


_SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")
_SUB = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def _plain(text) -> str:
    """Strip HTML and turn Unicode super/subscripts (cm⁻¹, ¹³C, C₆H₅) into mathtext, which every font can draw."""
    t = html.unescape(_TAG.sub("", text or "")).strip().replace("$", r"\$")
    t = re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+", lambda m: "$^{" + m.group(0).translate(_SUP) + "}$", t)
    t = re.sub(r"[₀₁₂₃₄₅₆₇₈₉]+", lambda m: "$_{" + m.group(0).translate(_SUB) + "}$", t)
    return t


def _color(c, default="#111827"):
    if not c:
        return default
    m = re.match(r"rgba?\(([^)]+)\)", str(c))
    if m:
        v = [float(t) for t in m.group(1).split(",")]
        return (v[0] / 255, v[1] / 255, v[2] / 255, v[3] if len(v) > 3 else 1.0)
    return c


def _arr(v):
    return np.array([np.nan if t is None else t for t in (v if v is not None else [])], dtype=float)


def series_names(fig) -> list[str]:
    """Names of the curves shown in the legend (the user can choose which ones to export)."""
    out = []
    for t in fig.data:
        if t.showlegend is not False and t.name and t.name not in out:
            out.append(t.name)
    return out


def to_publication(fig, fmt="tiff", dpi=600, width_cm=17.4, height_cm=None, font_pt=8.0, line_pt=1.0,
                   include=None, title=True, legend="bottom", grid=False, font="Arial") -> tuple[bytes, str, str]:
    """Redraw the Plotly figure `fig` with matplotlib and return (file bytes, extension, mime)."""
    include = set(include) if include is not None else set(series_names(fig))
    groups = {t.legendgroup for t in fig.data if t.name in include and t.legendgroup}
    stacked = any((t.yaxis or "y") == "y2" for t in fig.data)
    image = next(iter(fig.layout.images), None) if fig.layout.images else None
    if height_cm is None:
        height_cm = width_cm * (0.6 if not stacked else 0.9) + (width_cm * 0.30 if image is not None else 0)

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [font, "Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": font_pt, "axes.titlesize": font_pt + 1, "axes.labelsize": font_pt, "xtick.labelsize": font_pt - 0.5,
        "ytick.labelsize": font_pt - 0.5, "legend.fontsize": font_pt - 0.5, "axes.linewidth": 0.6 * line_pt + 0.2,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.direction": "out", "ytick.direction": "out",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": dpi,
    })
    fig_m = plt.figure(figsize=(width_cm / 2.54, height_cm / 2.54))
    cm = 1 / height_cm  # one centimetre as a fraction of the figure height
    has_title = bool(title and fig.layout.title.text)
    top = 1 - (0.75 if has_title else 0.25) * cm * (font_pt / 8)
    if image is not None:
        img_frac = (width_cm * 0.30) / height_cm
        img_ax = fig_m.add_axes([0.08, top - img_frac, 0.84, img_frac])
        img_ax.axis("off")
        raw = base64.b64decode(image.source.split(",", 1)[1])
        img_ax.imshow(Image.open(io.BytesIO(raw)), interpolation="lanczos")
        top -= img_frac + 0.15 * cm
    bottom = ((1.0 if legend == "bottom" else 0.0) + 1.0) * cm * (font_pt / 8)
    rows = 2 if stacked else 1
    left = min(0.2, 1.5 / width_cm * (font_pt / 8))
    gs = fig_m.add_gridspec(rows, 1, left=left, right=0.97 if legend != "right" else 0.74, bottom=bottom, top=top, hspace=0.12)
    axes = [fig_m.add_subplot(gs[0])]
    if stacked:
        axes.append(fig_m.add_subplot(gs[1], sharex=axes[0]))
    pick = lambda t: axes[1] if stacked and (t.yaxis or "y") == "y2" else axes[0]

    handles = []
    for t in fig.data:
        if t.type != "scatter":
            continue
        own = t.name in include if t.showlegend is not False else (t.legendgroup in groups if t.legendgroup else True)
        if not own:
            continue
        ax = pick(t)
        x, y = _arr(t.x), _arr(t.y)
        mode = t.mode or "lines"
        lc = _color(t.line.color if t.line else None)
        if "text" in mode and t.text is not None:  # peak labels
            below = (t.textposition or "").startswith("bottom")
            tc = _color(t.textfont.color if t.textfont else None, lc)
            for xi, yi, s in zip(x, y, t.text):
                if np.isnan(xi) or np.isnan(yi):
                    continue
                ax.plot([xi, xi], [yi, yi], marker="|", color=tc, ms=3, mew=0.6)
                ax.annotate(_plain(s), (xi, yi), xytext=(0, -3 if below else 3), textcoords="offset points",
                            ha="center", va="top" if below else "bottom", fontsize=font_pt - 1.5, color=tc,
                            fontweight="bold" if "<b>" in str(s) else "normal", annotation_clip=False)
            continue
        if mode == "markers":  # hover helpers (tops of stick spectra) are not drawn
            continue
        lw = (t.line.width or 1.5) / 1.6 * line_pt if t.line else line_pt
        (h,) = ax.plot(x, y, color=lc, lw=lw, ls=DASH.get(t.line.dash if t.line else None, "-"), solid_capstyle="butt")
        if t.showlegend is not False and t.name:
            h.set_label(_plain(t.name))
            handles.append(h)

    for sh in fig.layout.shapes or []:  # shaded bands (visible spectrum strip)
        if sh.type == "rect" and sh.xref in ("x", None):
            for ax in axes:
                ax.axvspan(sh.x0, sh.x1, color=_color(sh.fillcolor), alpha=(sh.opacity or 0.1) * 0.9, lw=0, zorder=0)

    lay = fig.layout
    for i, ax in enumerate(axes):
        xa = lay.xaxis if i == 0 else (lay.xaxis2 or lay.xaxis)
        ya = lay.yaxis if i == 0 else (lay.yaxis2 or lay.yaxis)
        if xa.range:
            ax.set_xlim(xa.range[0], xa.range[1])
        elif xa.autorange == "reversed":
            ax.invert_xaxis()
        if ya.range:
            ax.set_ylim(ya.range[0], ya.range[1])
        ax.set_ylabel(_plain(ya.title.text) if ya.title else "")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(length=3)
        if grid:
            ax.grid(color="0.88", lw=0.4)
        if i < len(axes) - 1:
            ax.tick_params(labelbottom=False)
    xa_last = (lay.xaxis2 if stacked and lay.xaxis2 and lay.xaxis2.title.text else lay.xaxis)
    axes[-1].set_xlabel(_plain(xa_last.title.text) if xa_last.title else "")
    if stacked and lay.annotations:  # subplot titles
        for ax, ann in zip(axes, [a for a in lay.annotations if a.xref == "paper"][:2]):
            ax.set_title(_plain(ann.text), fontsize=font_pt, loc="left")
    if title and lay.title.text:
        fig_m.suptitle(_plain(lay.title.text), fontsize=font_pt + 2, fontweight="bold", y=1 - 0.1 * cm, va="top")
    if handles and legend != "none":
        if legend == "right":
            fig_m.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.76, top), frameon=False)
        else:
            fig_m.legend(handles=handles, loc="lower center", ncol=min(len(handles), 3), frameon=False,
                         bbox_to_anchor=(0.5, 0.05 * cm), handlelength=1.8, columnspacing=1.4)

    ext, mime = FORMATS[fmt]
    buf = io.BytesIO()
    kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if fmt == "tiff" else {}
    fig_m.savefig(buf, format=ext if fmt != "tiff" else "tiff", dpi=dpi, facecolor="white", **kw)
    plt.close(fig_m)
    return buf.getvalue(), ext, mime
