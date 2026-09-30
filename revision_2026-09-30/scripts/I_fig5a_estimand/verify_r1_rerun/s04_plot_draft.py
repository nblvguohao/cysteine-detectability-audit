"""Step 4 - DRAFT figure of the absolute AUCs behind Figure 5a (for the authors to restyle).

POST HOC revision analysis (2026-09-30), item I_fig5a_estimand. Not registered, not pre-specified.
Plots only stored values from t1_absolute_auc_by_arm.csv (written by s01_tables.py); computes nothing.

Encoding (matches the submitted Fig. 5a where it overlaps): colour and marker shape = negative set
(NEG_A ink circle, NEG_B blue square); fill = feature set (open = VIS10, filled = DIG25).
Identity is never carried by colour alone (shape + fill + legend).

Run:  python -B s04_plot_draft.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import common as C  # noqa: E402

MM = 1 / 25.4
LW, FS = 1.0, 7
INK, AXIS, MUTE, REF, BAND, BLUE = "#1f2326", "#4a4f55", "#6f757c", "#c3c8ce", "#f1f3f5", "#1f6aa5"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FS, "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": FS,
    "axes.linewidth": LW, "lines.linewidth": LW, "xtick.major.width": LW, "ytick.major.width": LW,
    "xtick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "axes.edgecolor": AXIS, "xtick.color": AXIS, "ytick.color": AXIS, "text.color": INK, "axes.labelcolor": INK,
    "pdf.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.titlesize": FS, "axes.titlelocation": "left", "axes.unicode_minus": True,
})
STYLE = {"NEG_A": dict(color=INK, marker="o"), "NEG_B": dict(color=BLUE, marker="s")}
XL = (0.25, 1.0)


def draw(ax, wide, caliber):
    y = 0.0
    ticks, labels = [], []
    for i, (did, lab) in enumerate(C.COHORTS):
        if i % 2 == 0:
            ax.add_patch(Rectangle((0, y - 0.5), 1, 1, transform=ax.get_yaxis_transform(),
                                   facecolor=BAND, edgecolor="none", zorder=0))
        for neg, dy in (("NEG_A", 0.2), ("NEG_B", -0.2)):
            r = wide[(wide.dataset_id == did) & (wide.negative_set == neg)].iloc[0]
            st = STYLE[neg]
            if r.usable != 1:
                ax.text(XL[0] + 0.01, y + dy, "no observed-unmodified class", color=MUTE, fontsize=FS,
                        va="center", ha="left", style="italic")
                continue
            for fs, off, filled in (("VIS10", 0.065, False), ("DIG25", -0.065, True)):
                p, lo, hi = r[f"{fs}_auc_{caliber}"], r[f"{fs}_auc_{caliber}_lo"], r[f"{fs}_auc_{caliber}_hi"]
                yy = y + dy + off
                ax.plot([max(lo, XL[0]), hi], [yy, yy], color=st["color"], lw=LW, solid_capstyle="butt", zorder=3)
                ax.plot([p], [yy], st["marker"], color=st["color"], ms=3.6, mew=LW,
                        mfc=st["color"] if filled else "white", zorder=4)
        ticks.append(y)
        labels.append(lab.replace(" (", "\n("))
        y -= 1.0
    ax.axvline(0.5, color=REF, lw=LW, ls=(0, (3, 2)), zorder=0.5)
    ax.set_xlim(*XL)
    ax.set_ylim(y + 0.5, 0.5)
    ax.set_xticks([0.3, 0.5, 0.7, 0.9])
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.spines["left"].set_visible(False)


def main():
    wide = pd.read_csv(C.OUT / "t1_absolute_auc_by_arm.csv")
    fig, axes = plt.subplots(1, 2, figsize=(183 * MM, 120 * MM), sharey=True,
                             gridspec_kw=dict(left=0.17, right=0.985, top=0.935, bottom=0.255, wspace=0.12))
    draw(axes[0], wide, "global")
    draw(axes[1], wide, "within")
    axes[1].tick_params(axis="y", labelleft=False)
    axes[0].set_xlabel("Out-of-fold AUC, all sites (global)")
    axes[1].set_xlabel("Out-of-fold AUC within protein or homology component")
    axes[0].set_title("Global", fontweight="bold")
    axes[1].set_title("Within protein", fontweight="bold")
    handles = [
        Line2D([], [], marker="o", color=INK, mfc="white", mew=LW, lw=LW, ms=3.6, label="VIS10, NEG_A"),
        Line2D([], [], marker="o", color=INK, mfc=INK, mew=LW, lw=LW, ms=3.6, label="DIG25, NEG_A"),
        Line2D([], [], marker="s", color=BLUE, mfc="white", mew=LW, lw=LW, ms=3.6, label="VIS10, NEG_B"),
        Line2D([], [], marker="s", color=BLUE, mfc=BLUE, mew=LW, lw=LW, ms=3.6, label="DIG25, NEG_B"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.17, 0.118), ncol=4, frameon=False)
    fig.text(0.17, 0.012, "VIS10: 10 peptide-visibility features; DIG25: all 25 theoretical-digest features "
             "(VIS10 plus 15 more).\n"
             "NEG_A: other cysteines of site-bearing proteins; NEG_B: cysteines the same "
             "experiment reported unmodified.\n"
             "Lines: 95% cluster-bootstrap intervals; dashed line: AUC 0.5.",
             ha="left", va="bottom", color=MUTE, fontsize=FS, linespacing=1.3)
    for ext in ("png", "pdf"):
        meta = {"CreationDate": None, "ModDate": None} if ext == "pdf" else {}
        fig.savefig(C.OUT / f"fig5a_absolute_auc_DRAFT.{ext}", dpi=300, metadata=meta)
    C.record_inputs("s04_plot_draft", [C.OUT / "t1_absolute_auc_by_arm.csv"])
    print("written", [str(C.OUT / f"fig5a_absolute_auc_DRAFT.{e}") for e in ("png", "pdf")])


if __name__ == "__main__":
    main()
