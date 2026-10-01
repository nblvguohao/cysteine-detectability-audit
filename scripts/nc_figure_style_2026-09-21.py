"""House style for the Nature Communications figures of manuscript R22 (2026-09-21).

Every rule here comes from the journal's own text (https://www.nature.com/ncomms/submit/how-to-submit, read 2026-09-21,
quoted in results/nc_figure_guidelines_quotes_2026-09-21.csv) or from looking at published NC figures:
  * "One- or two-column format figures are required."      -> drawn at the PRINT width, 89 mm or 183 mm, never rescaled,
                                                               so a 7 pt label is 7 pt on paper
  * "Figure lettering should be in a clear, sans-serif typeface (for example, Helvetica); ... the same typeface in
     approximately the same font size should be used for all figures"   -> Arial, 7 pt text, 8 pt bold panel letters
  * "All display items should be on a white background, and should avoid excessive boxing, unnecessary colour"
  * "The thinnest lines in the final figure should be no smaller than one point wide."   -> every line >= 1.0 pt
  * "prepare figures to fit the pdf page size (210 x 276 mm)"   -> height <= 240 mm leaves room for the legend
  * observed in published NC Articles (two read on 2026-09-21): Fig. 1 opens with a study-design schematic; panel
    letters are bold lowercase at the top left; panels carry short titles, and explanations live in the legend, not in
    the plot.
Colour: Okabe-Ito, which stays distinguishable under the common colour-vision deficiencies; one colour per artefact.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MM = 1 / 25.4
DOUBLE = 183 * MM
SINGLE = 89 * MM
MAX_H = 240 * MM
LW = 1.0
FS = 7
FS_PANEL = 8

INK, GREY, LIGHT = "#1a1a1a", "#8a8a8a", "#d4d4d4"
BLUE, ORANGE, GREEN, PINK, VERMIL, SKY, YELLOW = "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442"
ARTEFACT = {1: BLUE, 2: ORANGE, 3: GREEN, 4: PINK, 5: VERMIL}

RC = {
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS,
    "legend.fontsize": FS, "axes.linewidth": LW, "lines.linewidth": LW, "xtick.major.width": LW, "ytick.major.width": LW,
    "xtick.major.size": 3, "ytick.major.size": 3, "xtick.minor.width": LW, "ytick.minor.width": LW,
    "patch.linewidth": LW, "hatch.linewidth": LW, "errorbar.capsize": 2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titlelocation": "left", "axes.titlepad": 4,
    "legend.frameon": False, "legend.handlelength": 1.4, "legend.borderaxespad": 0.2,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
}


def apply():
    plt.rcParams.update(RC)


def panel_letter(fig, ax, letter, dx=-0.0, dy=0.0):
    """Bold lowercase letter at the top-left of the axes' bounding box (figure coordinates)."""
    fig.canvas.draw()
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bb.x0 + dx, bb.y1 + dy, letter, fontsize=FS_PANEL, fontweight="bold", va="bottom", ha="left")


def ci(ax, y, point, lo, hi, colour=INK, marker="o", filled=True, ms=3.6, horizontal=True, zorder=3):
    if horizontal:
        ax.plot([lo, hi], [y, y], color=colour, lw=LW, solid_capstyle="butt", zorder=zorder)
        ax.plot([point], [y], marker, color=colour, ms=ms, mfc=colour if filled else "white", mew=LW, zorder=zorder + 1)
    else:
        ax.plot([y, y], [lo, hi], color=colour, lw=LW, solid_capstyle="butt", zorder=zorder)
        ax.plot([y], [point], marker, color=colour, ms=ms, mfc=colour if filled else "white", mew=LW, zorder=zorder + 1)


def min_linewidth(fig):
    """Smallest line width actually present on the canvas, in points (gate: >= 1.0)."""
    widths = []
    for ax in fig.get_axes():
        for ln in ax.get_lines():
            if ln.get_visible() and ln.get_linestyle() not in ("None", "", " "):
                widths.append(ln.get_linewidth())
        for sp in ax.spines.values():
            if sp.get_visible():
                widths.append(sp.get_linewidth())
        for c in ax.collections:
            lw = c.get_linewidths()
            if len(lw) and c.get_visible():
                widths += [w for w in lw if w > 0]
        for p in ax.patches:
            if p.get_visible() and p.get_linewidth() > 0 and p.get_edgecolor()[3] > 0:
                widths.append(p.get_linewidth())
        for t in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
            if t.tick1line.get_visible():
                widths.append(t.tick1line.get_markeredgewidth())
    return min(widths) if widths else None


def font_sizes(fig):
    sizes = set()
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and (t.get_text() or "").strip():
            sizes.add(round(t.get_fontsize(), 2))
    return sorted(sizes)
