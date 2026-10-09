#!/usr/bin/env python3
"""Figure 5 for the revision (2026-09-30): panel a redrawn as ABSOLUTE out-of-fold AUCs, panels b and c unchanged.

WHAT CHANGED AND WHY
  As Supplemental Note 20 explains, the ratio plotted in the submitted panel a
  (AUC_VIS10 - 0.5)/(AUC_DIG25 - 0.5) is numerically right but cannot carry the reading "share of discrimination that
  is visibility rather than chemistry". The manuscript's Figure 5 legend now describes panel a as absolute AUCs:
  VIS10 (open symbols) and DIG25 (filled symbols), over all sites (left) and within proteins or homology components
  (right), NEG_A (circles) and NEG_B (squares), dashed line at AUC 0.5, NEG_B unavailable for one cohort (FAT-switch).
  Panel a draws exactly that. The ratio is no longer drawn; its Source Data rows stay, relabelled
  'ratio; moved to Supplemental Note 20'.

PANELS b AND c ARE CARRIED OVER, NOT REDRAWN
  The original drawing code (repo/scripts/plot_figs2to7_redesign_2026-09-23.py, fig5(), which calls
  plot_manuscript_f1_three_axes.collect_panel_b/c) cannot regenerate them from the local inputs: the local
  repo/results/phase2e_claim_retest.csv holds the panel-c specification superseded on 2026-09-24 (SNO-012
  0.4364/0.5402 and SNO-002 0.1555/0.1358, against 0.5636/0.4598 and 0.3403/0.3325 in the submitted panel), the local
  module's precision-tier wording is older than the submitted panel's, and it stacks panel b's two backgrounds in the
  other order. So b and c are taken from the submitted PDF (MCP/figures/Fig5_three_axes.pdf): a copy of its page has
  everything in the old panel-a region removed by PDF redaction (text and vector art, so the old ratio panel does not
  survive as hidden content), and what remains is placed 1:1, at its original coordinates, over the new panel a.
  The script then checks that every text span and vector path of b and c is unchanged and that both regions render
  pixel-identically to the submitted figure at 200 dpi.

DATA RULE
  No statistic is computed here. Every drawn point and interval end of panel a is one stored cell of
  inputs/repo_results/ptm_detectability_share.csv (sha256 pinned below; byte-identical to
  repo/results/ptm_detectability_share.csv), read as a string and converted to float for placement only. The only
  arithmetic is layout. Each drawn value is logged with its raw string and must appear verbatim in the written
  Source Data. Before drawing, every stored string is also checked against the I_fig5a_estimand tables
  (Source_Data_Fig5a_absolute_auc_ADDITION_PROPOSED.csv, t1_absolute_auc_long.csv); the 3-decimal values quoted in
  Supplemental Note 20 and the b/c Source Data rows against Supplemental Data 5 and 6 are checked and reported.

STYLE (identical to the submitted Fig. 5, whose b and c panels are kept)
  Arial 7 pt text, 8 pt bold panel letter, pdf.fonttype 42, 1.0 pt axes, ticks and interval lines, tick length
  2.5 pt, no top, right or left spines, pale alternating row bands, INK for NEG_A and BLUE for NEG_B as in the
  submitted panel a (scripts/nc_figure_style_v2_2026-09-23.py, copied here). Width 7.2 in (518.4 pt); height
  150 mm (425.2 pt, 5.91 in), as submitted. Panel a's axes box, letter and title sit where the submitted panel a had
  them, so its x axes, tick labels and axis labels share their lines with panel b.

OUTPUTS (this folder only): Fig5_three_axes.pdf, Fig5_three_axes.png (200 dpi, rendered from the PDF),
  Source_Data_Fig5_three_axes.csv. Nothing else is written; run with `python -B` so no bytecode is left.
USAGE: python -B fig5_rebuild.py [--old-pdf PATH] [--old-source-data PATH]
  The two options point at copies of the SUBMITTED files if MCP/ has been updated since; they are hash-checked.
Tested with Python 3.10.11, matplotlib 3.10.9, PyMuPDF 1.28.0, numpy 2.2.6.
"""
import argparse
import csv
import hashlib
import io
import pathlib
import sys
from decimal import ROUND_HALF_UP, Decimal

import fitz  # PyMuPDF
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent            # .../revision_2026-09-30/results/figures_rev/Fig5
REV = HERE.parents[2]                                      # .../revision_2026-09-30
WORK = REV.parents[1]                                      # .../_cys_repo_work
MCP = WORK.parent / "巯基化" / "MCP"                        # the submission package (read only)

SHARE = REV / "inputs" / "repo_results" / "ptm_detectability_share.csv"
SHARE_REPO = WORK / "repo" / "results" / "ptm_detectability_share.csv"
ADDITION = REV / "results" / "I_fig5a_estimand" / "Source_Data_Fig5a_absolute_auc_ADDITION_PROPOSED.csv"
T1_LONG = REV / "results" / "I_fig5a_estimand" / "t1_absolute_auc_long.csv"
NOTE20 = MCP / "supplemental" / "Supplemental_Note_20_fig5a_absolute_discrimination.md"
SD5 = MCP / "supplemental" / "Supplemental_Data_5_retest_round_d.csv"
SD6 = MCP / "supplemental" / "Supplemental_Data_6_retest_round_e.csv"
OLD_PDF_DEFAULT = MCP / "figures" / "Fig5_three_axes.pdf"
OLD_SD_DEFAULT = MCP / "source_data" / "Source_Data_Fig5_three_axes.csv"

PINNED = {  # inputs the figure depends on; any change refuses the run
    "share": "732fd67a04d0a5b9a261873cbe85196389ef98d79a1a25b1a54273dd80f9c249",
    "old_pdf": "0dabb1a0ea9e8dfcbe8134332f1e8e7f26708236201ad8cf33ced41ba615ea94",
    "old_sd": "624780e0dab500fae47e246b148c25a69ba6512515656e8456d21f99fe1d31cf",
    "addition": "af25af151888cbce5b561360b629bcf549128a67259fcab981a44cfa7f6e61fd",
    "t1_long": "7411fa53e01b3905298f30c107c41901794f6b0004b38d390e6038807520c909",
}
OUT_PDF = HERE / "Fig5_three_axes.pdf"
OUT_PNG = HERE / "Fig5_three_axes.png"
OUT_SD = HERE / "Source_Data_Fig5_three_axes.csv"
PNG_DPI = 200

# ---------------------------------------------------------------------------------------------- house style
# copied from repo/scripts/nc_figure_style_2026-09-21.py (v1 RC) and nc_figure_style_v2_2026-09-23.py (v2 RC,
# colours), plus the Arial mathtext setting of plot_figs_r34_2026-09-23.py
LW, FS, FS_PANEL = 1.0, 7, 8
INK, AXIS, MUTE, REF, BAND, BLUE = "#1f2326", "#4a4f55", "#6f757c", "#c3c8ce", "#f1f3f5", "#1f6aa5"
RC = {
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS,
    "legend.fontsize": FS, "axes.linewidth": LW, "lines.linewidth": LW, "xtick.major.width": LW,
    "ytick.major.width": LW, "xtick.minor.width": LW, "ytick.minor.width": LW, "patch.linewidth": LW,
    "hatch.linewidth": LW, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.borderaxespad": 0.2, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.edgecolor": AXIS, "xtick.color": AXIS, "ytick.color": AXIS, "text.color": INK, "axes.labelcolor": INK,
    "xtick.labelcolor": INK, "ytick.labelcolor": INK, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.major.pad": 2.0, "ytick.major.pad": 2.0, "axes.labelpad": 3.0, "legend.handlelength": 1.2,
    "legend.handletextpad": 0.45, "legend.columnspacing": 1.1, "legend.labelspacing": 0.35,
    "lines.markeredgewidth": LW, "lines.solid_capstyle": "butt", "axes.unicode_minus": True,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "mathtext.sf": "Arial",
}

# ---------------------------------------------------------------------------------------------- geometry (pt)
W_PT = 518.4                    # 7.2 in, as required
H_PT = 150 / 25.4 * 72          # 425.197 pt = 150 mm, the submitted height (5.91 in <= 6 in)
# the submitted panel a axes box (left, top, bottom; pt from the top-left corner), read from the submitted PDF's
# drawings; top and bottom are also panel b's axes box, which is asserted against that PDF below
AX_LEFT, AX_TOP, AX_BOTTOM = 90.78, 19.13, 212.60
AX_RIGHT, SUB_GAP = 304.0, 12.0  # right edge of the within-protein sub-panel; gap between the two sub-panels
LEGEND_TOP = 243.0               # top of the 2 x 2 key, under both sub-panels
# regions of the submitted page: everything of its panel a lies in REGION_A_OLD, nothing of b or c does
REGION_A_OLD = fitz.Rect(0, 0, 305, 268)
REGION_B = fitz.Rect(310, 0, W_PT, 240)
REGION_C = fitz.Rect(0, 270, W_PT, H_PT)

# ---------------------------------------------------------------------------------------------- panel a content
# row order and labels exactly as in the submitted panel a (labels from plot_manuscript_f1_three_axes.DATASET_LABEL,
# split before the species as the redesign script did)
COHORTS = [
    ("qtrp_S1_ph5", "QTRP S1 pH5 (human)"),
    ("qtrp_S2_ph5", "QTRP S2 pH5 (human)"),
    ("qpers_sid_tierB", "qPerS-SID tier B (human)"),
    ("cysboost2019_human_sno_hela", "Cys-BOOST HeLa (human)"),
    ("cysboost2019_human_sno_shsy5y", "Cys-BOOST SH-SY5Y (human)"),
    ("natcomm2023_ath_sno", "FAT-switch (Arabidopsis)"),
    ("abiotech2025_ath_sno", "PAT-switch (Arabidopsis)"),
    ("fps2020_ath_sulfenyl", "YAP1C reporter (Arabidopsis)"),
]
CONSTRUCTION = {"A": "NEG_A_not_observed", "B": "NEG_B_observed_unmodified"}
CALIBER = {"global": "auc_global", "within_protein": "auc_within_unit"}   # Source Data name -> stored column
FEATURE_SETS = ("VIS10", "DIG25", "CPL15")                                   # CPL15 goes to Source Data only
DRAWN_SETS = ("VIS10", "DIG25")
OFFSET = {("A", "VIS10"): 0.33, ("A", "DIG25"): 0.11, ("B", "VIS10"): -0.11, ("B", "DIG25"): -0.33}
COLOUR = {"A": INK, "B": BLUE}
MARKER = {"A": "o", "B": "s"}
MS = {"A": 3.3, "B": 3.0}          # a square reads larger than a circle of the same size
NOTE_X, NOTE_DY = 0.53, -0.22      # "NEG_B unavailable", centred on the two NEG_B sub-rows, right of AUC 0.5
XLIM, XTICKS = (0.25, 1.0), [0.3, 0.5, 0.7, 0.9]

SPEC_PLOT = "absolute out-of-fold AUC; plotted in panel a; post hoc revision 2026-09-30"
SPEC_CPL15 = "absolute out-of-fold AUC of CPL15 (the 15 columns DIG25 adds); not plotted; post hoc revision 2026-09-30"
SPEC_COUNT = "count; not plotted; post hoc revision 2026-09-30"
SPEC_NO_NEGB = ("count; no observed-unmodified class (stored usable = 0; negatives 0 < 30); "
                "no NEG_B symbols drawn; post hoc revision 2026-09-30")
SPEC_RATIO = "ratio; moved to Supplemental Note 20"
SHARE_TABLE = "analysis output: ptm_detectability_share.csv"
COUNT_FIELDS = [("n_positive", "n_positive"), ("n_negative", "n_negative"), ("n_units", "n_units"),
                ("n_units_with_positive", "n_units_with_positive"), ("n_units_both_classes", "n_informative_units")]


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def check(condition, message):
    if not condition:
        sys.exit("REFUSE: " + message)


# ---------------------------------------------------------------------------------------------- 1. stored values
def collect(share_rows):
    """One record per cohort; per construction either None (not usable) or the raw stored strings."""
    table = {(r["dataset_id"], r["negative_construction"], r["feature_set"]): r for r in share_rows}
    check(len(table) == len(share_rows), "duplicate (dataset, construction, feature set) in the stored table")
    records = []
    for dataset, label in COHORTS:
        rec = {"dataset_id": dataset, "label": label}
        for tag, construction in CONSTRUCTION.items():
            rows = {fs: table[(dataset, construction, fs)] for fs in FEATURE_SETS}
            usable = {rows[fs]["usable"] for fs in FEATURE_SETS}
            check(len(usable) == 1, "usable differs across feature sets for %s %s" % (dataset, construction))
            if usable == {"0"}:
                r0 = rows["VIS10"]
                rec[tag] = None
                rec[tag + "_excluded"] = {"n_negative": r0["n_negative"], "exclusion_reason": r0["exclusion_reason"]}
                continue
            counts = {}
            for out_field, col in COUNT_FIELDS:
                vals = {rows[fs][col] for fs in FEATURE_SETS}
                check(len(vals) == 1, "count %s differs across feature sets for %s" % (col, dataset))
                counts[out_field] = vals.pop()
            aucs = {}
            for fs in FEATURE_SETS:
                for cal, col in CALIBER.items():
                    raw = (rows[fs][col], rows[fs][col + "_lo"], rows[fs][col + "_hi"])
                    check(all(v.strip() for v in raw), "empty stored AUC for %s %s %s" % (dataset, tag, fs))
                    p, lo, hi = (float(v) for v in raw)
                    check(lo <= p <= hi, "point outside its interval for %s %s %s %s" % (dataset, tag, fs, cal))
                    aucs[(fs, cal)] = raw
            rec[tag] = {"counts": counts, "aucs": aucs}
        records.append(rec)
    check([r["dataset_id"] for r in records if r["B"] is None] == ["natcomm2023_ath_sno"],
          "expected exactly one cohort without NEG_B (FAT-switch)")
    return records


def source_rows(records):
    """Panel a Source Data rows in the order of the I_fig5a_estimand proposal (counts, then per feature set the
    global and within-protein AUC with their interval ends)."""
    out = []
    for rec in records:
        for tag in ("A", "B"):
            key = "%s|NEG_%s" % (rec["dataset_id"], tag)
            arm = rec[tag]
            if arm is None:
                out.append([key, "n_negative", rec[tag + "_excluded"]["n_negative"], SHARE_TABLE, SPEC_NO_NEGB])
                continue
            for out_field, _ in COUNT_FIELDS:
                out.append([key, out_field, arm["counts"][out_field], SHARE_TABLE, SPEC_COUNT])
            for fs in FEATURE_SETS:
                spec = SPEC_PLOT if fs in DRAWN_SETS else SPEC_CPL15
                for cal in CALIBER:
                    field = "auc_%s_%s" % (cal, fs)
                    p, lo, hi = arm["aucs"][(fs, cal)]
                    out += [[key, field, p, SHARE_TABLE, spec], [key, field + "_lo", lo, SHARE_TABLE, spec],
                            [key, field + "_hi", hi, SHARE_TABLE, spec]]
    return out


# ---------------------------------------------------------------------------------------------- 2. cross-checks
def check_against_item_tables(rows_a):
    """Hard checks: every stored string equals the I_fig5a_estimand tables."""
    mine = {(r[0], r[1]): r[2] for r in rows_a}
    addition = read_csv(ADDITION)
    check(len(addition) == 345, "the proposed addition should have 345 rows")
    for a in addition:
        check((a["row"], a["field"]) in mine, "proposed row missing: %s %s" % (a["row"], a["field"]))
        check(mine[(a["row"], a["field"])] == a["value"] and a["source_table"] == SHARE_TABLE,
              "value differs from the proposed addition: %s %s" % (a["row"], a["field"]))
    long_rows = [r for r in read_csv(T1_LONG) if r["in_fig5a"] == "1" and r["usable"] == "1"]
    check(len(long_rows) == 15 * 6, "t1_absolute_auc_long should hold 90 usable Fig. 5a rows")
    for r in long_rows:
        key = "%s|%s" % (r["dataset_id"], r["negative_set"])
        field = "auc_%s_%s" % (r["caliber"], r["feature_set"])
        for suffix, col in (("", "auc"), ("_lo", "auc_lo"), ("_hi", "auc_hi")):
            check(mine[(key, field + suffix)] == r[col], "t1_absolute_auc_long differs: %s %s" % (key, field))
    return {"proposed_addition_rows_matched": len(addition), "t1_long_rows_matched": len(long_rows)}


def three_decimals(raw):
    return str(Decimal(raw).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def check_note20(records):
    """Soft check, reported: the 3-decimal values of Table S14.1 in Supplemental Note 20 (the note the legend cites)."""
    if not NOTE20.exists():
        return {"checked": 0, "status": "Supplemental Note 20 not found"}
    lines = NOTE20.read_text(encoding="utf-8").splitlines()
    cols = {"VIS10 global": ("VIS10", "global"), "DIG25 global": ("DIG25", "global"),
            "VIS10 within protein": ("VIS10", "within_protein"), "DIG25 within protein": ("DIG25", "within_protein")}
    header = next((i for i, l in enumerate(lines) if l.startswith("| Cohort | Neg. set | sites / negatives")), None)
    if header is None:
        return {"checked": 0, "status": "Table S14.1 header not found"}
    names = [c.strip() for c in lines[header].strip("|").split("|")]
    by_label = {rec["label"].split(" (")[0]: rec for rec in records}
    checked, mismatches = 0, []
    for line in lines[header + 2:]:
        if not line.startswith("|"):
            break
        cells = dict(zip(names, [c.strip() for c in line.strip("|").split("|")]))
        rec = by_label.get(cells.get("Cohort"))
        tag = {"NEG_A": "A", "NEG_B": "B"}.get(cells.get("Neg. set"))
        if rec is None or tag is None or rec[tag] is None:
            continue
        for col, (fs, cal) in cols.items():
            p, lo, hi = rec[tag]["aucs"][(fs, cal)]
            expected = "%s [%s, %s]" % (three_decimals(p), three_decimals(lo), three_decimals(hi))
            checked += 1
            if cells.get(col) != expected:
                mismatches.append((cells.get("Cohort"), cells.get("Neg. set"), col, cells.get(col), expected))
    return {"checked": checked, "mismatches": mismatches, "sha256": sha256(NOTE20)}


def check_bc_against_supplemental(old_sd_rows):
    """Soft check, reported: panel b and c Source Data values against the supplemental tables they cite."""
    report = {}
    if SD5.exists():
        sd5 = {(r["claim_id"], r["specification"]): r for r in read_csv(SD5)}
        b = [r for r in old_sd_rows if r["panel"] == "b"]
        bad = [r for r in b if sd5.get(("PERS-008", r["row"]), {}).get(r["field"]) != r["value"]]
        report["b"] = {"rows": len(b), "mismatches": len(bad), "sha256_SD5": sha256(SD5)}
    if SD6.exists():
        sd6 = {(r["claim_id"], r["specification"]): r for r in read_csv(SD6)}
        c = [r for r in old_sd_rows if r["panel"] == "c"]
        bad = [r for r in c if sd6.get((r["row"], "primary"), {}).get(r["field"]) != r["value"]]
        report["c"] = {"rows": len(c), "mismatches": len(bad), "sha256_SD6": sha256(SD6)}
    return report


# ---------------------------------------------------------------------------------------------- 3. panel a
def band(ax, y0, y1):
    ax.add_patch(Rectangle((0, min(y0, y1)), 1, abs(y1 - y0), transform=ax.get_yaxis_transform(),
                           facecolor=BAND, edgecolor="none", lw=0, zorder=0))


def ci(ax, y, point, lo, hi, colour, marker, filled, ms):
    ax.plot([lo, hi], [y, y], color=colour, lw=LW, solid_capstyle="butt", zorder=3)
    ax.plot([point], [y], marker, color=colour, ms=ms, mfc=colour if filled else "white", mew=LW, zorder=4)


def draw_panel_a(records, letter_origin, title_origin):
    """Returns (PDF bytes of a page carrying panel a only, ledger of drawn values)."""
    plt.rcParams.update(RC)
    fig = plt.figure(figsize=(W_PT / 72, H_PT / 72))

    def fx(x):
        return x / W_PT

    def fy(y):          # y in pt from the top
        return 1 - y / H_PT

    w = (AX_RIGHT - AX_LEFT - SUB_GAP) / 2
    h = (AX_BOTTOM - AX_TOP) / H_PT
    ax_g = fig.add_axes([fx(AX_LEFT), fy(AX_BOTTOM), fx(w), h])
    ax_w = fig.add_axes([fx(AX_LEFT + w + SUB_GAP), fy(AX_BOTTOM), fx(w), h], sharey=ax_g)
    ledger = []
    for ax, cal in ((ax_g, "global"), (ax_w, "within_protein")):
        for i, rec in enumerate(records):
            y = -i
            if i % 2 == 0:
                band(ax, y + 0.5, y - 0.5)
            for tag in ("A", "B"):
                arm = rec[tag]
                if arm is None:
                    ax.text(NOTE_X, y + NOTE_DY, "NEG_%s unavailable" % tag, color=MUTE, ha="left", va="center",
                            fontsize=FS)
                    continue
                for fs in DRAWN_SETS:
                    raw = arm["aucs"][(fs, cal)]
                    p, lo, hi = (float(v) for v in raw)
                    ci(ax, y + OFFSET[(tag, fs)], p, lo, hi, COLOUR[tag], MARKER[tag], fs == "DIG25", MS[tag])
                    field = "auc_%s_%s" % (cal, fs)
                    for suffix, r, v in (("", raw[0], p), ("_lo", raw[1], lo), ("_hi", raw[2], hi)):
                        ledger.append(("%s|NEG_%s" % (rec["dataset_id"], tag), field + suffix, r, v))
        ax.axvline(0.5, color=REF, lw=LW, ls=(0, (3, 2)), zorder=0.5)
        ax.set_xlim(*XLIM)
        ax.set_xticks(XTICKS)
        ax.set_ylim(-len(records) + 0.5, 0.5)
        ax.tick_params(axis="y", length=0, pad=3)
        ax.spines["left"].set_visible(False)
    ax_g.set_yticks([-i for i in range(len(records))])
    ax_g.set_yticklabels([rec["label"].replace(" (", "\n(") for rec in records])
    ax_w.tick_params(axis="y", labelleft=False)
    ax_g.set_xlabel("Global AUC (all sites)")
    ax_w.set_xlabel("Within-protein AUC")
    # 2 x 2 key: columns = feature set (open VIS10, filled DIG25), rows = negative set (circle NEG_A, square NEG_B)
    handles = [Line2D([], [], marker=MARKER[t], color=COLOUR[t], mfc=COLOUR[t] if fs == "DIG25" else "white",
                      lw=LW, ms=MS[t], label="%s, NEG_%s" % (fs, t)) for fs in DRAWN_SETS for t in ("A", "B")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(fx((AX_LEFT + AX_RIGHT) / 2), fy(LEGEND_TOP)),
               ncol=2, columnspacing=2.0, frameon=False, borderaxespad=0)
    # letter and title at the submitted panel a's positions (read from the submitted PDF's text layer)
    fig.text(fx(letter_origin[0]), fy(letter_origin[1]), "a", fontsize=FS_PANEL, fontweight="bold", ha="left",
             va="baseline", color=INK)
    fig.text(fx(title_origin[0]), fy(title_origin[1]), "Negative-set construction", fontsize=FS, ha="left",
             va="baseline", color=INK)
    fonts = {round(t.get_fontsize(), 2) for t in fig.findobj(matplotlib.text.Text)
             if t.get_visible() and (t.get_text() or "").strip()}
    check(fonts <= {7.0, 8.0}, "panel a uses font sizes other than 7 and 8: %s" % sorted(fonts))
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    return buf.getvalue(), ledger


# ---------------------------------------------------------------------------------------------- 4. PDF helpers
def spans(page, clip=None):
    out = []
    for block in page.get_text("dict", clip=clip)["blocks"]:
        for line in block.get("lines", []):
            for s in line["spans"]:
                if s["text"].strip():
                    out.append((s["text"], s["font"], round(s["size"], 3), tuple(round(v, 2) for v in s["bbox"]),
                                tuple(round(v, 2) for v in s["origin"])))
    return sorted(out, key=lambda t: (t[3][1], t[3][0], t[0]))


def within(r, region):
    """r lies inside region. Written out because fitz treats a zero-height or zero-width rect (a line) as
    contained in every rect."""
    r = fitz.Rect(r)
    return region.x0 <= r.x0 and r.x1 <= region.x1 and region.y0 <= r.y0 and r.y1 <= region.y1


def touches(r, region):
    r = fitz.Rect(r)
    return r.x0 <= region.x1 and region.x0 <= r.x1 and r.y0 <= region.y1 and region.y0 <= r.y1


def drawings(page, clip):
    out = []
    for d in page.get_drawings():
        r = d["rect"]
        if r.width > 500 and r.height > 400:        # a full-page background
            continue
        if touches(r, clip):
            out.append((d["type"], tuple(round(v, 2) for v in r), tuple(round(v, 4) for v in (d.get("fill") or ())),
                        tuple(round(v, 4) for v in (d.get("color") or ())), d.get("width"), len(d["items"])))
    return sorted(out)


def pixels(page, dpi=PNG_DPI):
    pix = page.get_pixmap(dpi=dpi, alpha=False)
    return np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)


def px_slice(rect, dpi=PNG_DPI, margin=2):
    s = dpi / 72
    return (slice(int(rect.y0 * s) + margin, int(rect.y1 * s) - margin),
            slice(int(rect.x0 * s) + margin, int(rect.x1 * s) - margin))


def inside(bbox, rect):
    return within(bbox, rect)


def partition_old_page(old_page):
    """Every text span and path of the submitted page lies in exactly one of the three regions."""
    counts = {"a": [0, 0], "b": [0, 0], "c": [0, 0]}
    regions = {"a": REGION_A_OLD, "b": REGION_B, "c": REGION_C}
    for s in spans(old_page):
        hit = [k for k, r in regions.items() if inside(s[3], r)]
        check(len(hit) == 1, "submitted text span not in exactly one region: %r" % (s,))
        counts[hit[0]][0] += 1
    for d in old_page.get_drawings():
        r = d["rect"]
        if r.width > 500 and r.height > 400:
            continue
        hit = [k for k, reg in regions.items() if within(r, reg)]
        check(len(hit) == 1, "submitted path not in exactly one region: %r" % (r,))
        counts[hit[0]][1] += 1
    return counts


def clean_old_page(old_doc):
    """Copy of the submitted page with everything in REGION_A_OLD (and the full-page background) redacted away."""
    tmp = fitz.open()
    tmp.insert_pdf(old_doc)
    page = tmp[0]
    page.add_redact_annot(REGION_A_OLD, fill=False)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_TOUCHED,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    return fitz.open("pdf", tmp.tobytes(garbage=4, deflate=True))


# ---------------------------------------------------------------------------------------------- 5. main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--old-pdf", default=str(OLD_PDF_DEFAULT))
    ap.add_argument("--old-source-data", default=str(OLD_SD_DEFAULT))
    args = ap.parse_args()
    old_pdf, old_sd = pathlib.Path(args.old_pdf), pathlib.Path(args.old_source_data)

    for name, path in (("share", SHARE), ("old_pdf", old_pdf), ("old_sd", old_sd), ("addition", ADDITION),
                       ("t1_long", T1_LONG)):
        check(path.exists(), "missing input %s" % path)
        check(sha256(path) == PINNED[name], "input changed: %s (sha256 %s)" % (path, sha256(path)))
    if SHARE_REPO.exists():
        check(sha256(SHARE_REPO) == PINNED["share"], "repo copy of the stored table differs from the input copy")

    # 1. stored values and the rows they become
    records = collect(read_csv(SHARE))
    rows_a = source_rows(records)
    item_checks = check_against_item_tables(rows_a)
    note20 = check_note20(records)

    # 2. the submitted page: geometry assumptions, then a cleaned copy carrying b and c only
    old_doc = fitz.open(str(old_pdf))
    old_page = old_doc[0]
    check(abs(old_page.rect.height - H_PT) < 1e-3, "submitted page height is not 150 mm")
    partition = partition_old_page(old_page)
    old_spans = spans(old_page)
    letter = next(s for s in old_spans if s[0] == "a" and s[2] == 8.0 and inside(s[3], REGION_A_OLD))
    title = next(s for s in old_spans if s[0] == "Negative-set construction")
    b_letter = next(s for s in old_spans if s[0] == "b" and s[2] == 8.0)
    c_letter = next(s for s in old_spans if s[0] == "c" and s[2] == 8.0)
    b_axis = [d for d in old_page.get_drawings() if within(d["rect"], REGION_B) and d["type"] == "s"
              and abs(d["rect"].y0 - d["rect"].y1) < 1e-6 and d["rect"].width > 100]
    check(len(b_axis) == 1 and abs(b_axis[0]["rect"].y0 - AX_BOTTOM) < 0.01, "panel b x axis is not at AX_BOTTOM")
    b_band_top = min(d["rect"].y0 for d in old_page.get_drawings() if within(d["rect"], REGION_B))
    check(abs(b_band_top - AX_TOP) < 0.01, "panel b axes top is not at AX_TOP")
    cleaned = clean_old_page(old_doc)
    cpage = cleaned[0]
    check(not spans(cpage, REGION_A_OLD) and not drawings(cpage, REGION_A_OLD),
          "old panel a content survived the redaction")
    check(not any(d["rect"].width > 500 and d["rect"].height > 400 for d in cpage.get_drawings()),
          "the submitted page background survived the redaction")
    for region in (REGION_B, REGION_C):
        check(spans(cpage, region) == spans(old_page, region), "redaction changed text in b or c")
        check(drawings(cpage, region) == drawings(old_page, region), "redaction changed paths in b or c")
    old_px = pixels(old_page)
    clean_px = pixels(cpage)
    for region in (REGION_B, REGION_C):
        ys, xs = px_slice(region)
        check(np.array_equal(old_px[ys, xs], clean_px[ys, xs]), "redaction changed the rendering of b or c")

    # 3. panel a on its own page, then b and c placed over it 1:1
    a_pdf, ledger = draw_panel_a(records, letter[4], title[4])
    a_doc = fitz.open("pdf", a_pdf)
    a_page = a_doc[0]
    check(abs(a_page.rect.width - W_PT) < 1e-3 and abs(a_page.rect.height - H_PT) < 1e-3, "panel a page size")
    for region in (REGION_B, REGION_C):
        check(not spans(a_page, region) and not drawings(a_page, region), "panel a reaches into the b or c region")
    a_px = pixels(a_page)
    out = fitz.open("pdf", a_pdf)
    page = out[0]
    full = fitz.Rect(0, 0, W_PT, H_PT)
    page.show_pdf_page(full, cleaned, 0, clip=full)
    out.set_metadata({"title": "Figure 5. Three design axes determine whether a site-preference claim is testable",
                      "creator": "fig5_rebuild.py (matplotlib %s; panels b and c carried over with PyMuPDF %s)"
                                 % (matplotlib.__version__, fitz.VersionBind),
                      "producer": "PyMuPDF %s" % fitz.VersionBind, "creationDate": "", "modDate": ""})
    pdf_bytes = out.tobytes(garbage=4, deflate=True, no_new_id=True)   # no random trailer /ID: reproducible bytes

    # 4. checks on the finished page
    fin_doc = fitz.open("pdf", pdf_bytes)
    fin = fin_doc[0]
    check(abs(fin.rect.width - W_PT) < 1e-3 and fin.rect.height <= 6 * 72 + 1e-6, "final page size")
    for region, name in ((REGION_B, "b"), (REGION_C, "c")):
        check(spans(fin, region) == spans(old_page, region), "panel %s text differs from the submitted figure" % name)
        check(drawings(fin, region) == drawings(old_page, region), "panel %s paths differ from the submitted" % name)
    fin_spans = spans(fin)
    a_spans = spans(a_page)
    check(sorted(s for s in fin_spans if not (inside(s[3], REGION_B) or inside(s[3], REGION_C))) == sorted(a_spans),
          "text outside b and c is not exactly panel a's")
    small = [s for s in fin_spans if s[2] < 5.7]
    check(small == [s for s in old_spans if s[2] < 5.7] and all(s[0] == "2" and inside(s[3], REGION_B) for s in small),
          "text below 5.7 pt other than panel b's subscript in log2: %r" % small)
    fonts = {f[3].split("+")[-1] for f in fin.get_fonts(full=True)}
    check(fonts == {"ArialMT", "Arial-BoldMT"}, "fonts other than Arial: %s" % fonts)
    check(all(fitz.Rect(s[3]) in fin.rect for s in fin_spans), "text outside the page")
    new_letter = next(s for s in fin_spans if s[0] == "a" and s[2] == 8.0)
    check(abs(new_letter[4][1] - b_letter[4][1]) < 0.01 and abs(new_letter[4][0] - c_letter[4][0]) < 0.01,
          "letter a not aligned with b (baseline) and c (left edge)")
    label_box = fitz.Rect(0, 20, 89, 212)
    old_labels = [(s[0], s[3]) for s in old_spans if inside(s[3], label_box)]
    new_labels = [(s[0], s[3]) for s in fin_spans if inside(s[3], label_box)]
    check(len(old_labels) == len(new_labels) == 16 and [t for t, _ in old_labels] == [t for t, _ in new_labels],
          "cohort labels differ from the submitted panel a")
    label_shift = max(abs(a - b) for (_, ob), (_, nb) in zip(old_labels, new_labels) for a, b in zip(ob, nb))
    check(label_shift < 0.05, "cohort labels moved (%.3f pt)" % label_shift)
    fin_px = pixels(fin)
    for region in (REGION_B, REGION_C):
        ys, xs = px_slice(region)
        check(np.array_equal(fin_px[ys, xs], old_px[ys, xs]), "b or c does not render as submitted")
    mask = np.ones(fin_px.shape[:2], dtype=bool)
    for region in (REGION_B, REGION_C):
        ys, xs = px_slice(region, margin=-2)
        mask[max(ys.start, 0):ys.stop, max(xs.start, 0):xs.stop] = False
    check(np.array_equal(fin_px[mask], a_px[mask]), "outside b and c the page does not render as panel a alone")

    # 5. Source Data: new panel a rows, the old ratio rows relabelled, panels b and c byte-identical
    old_text = old_sd.read_bytes().decode("utf-8")
    check(old_text.endswith("\r\n") and "\n" not in old_text.replace("\r\n", ""), "submitted Source Data not CRLF")
    old_lines = old_text[:-2].split("\r\n")
    header = old_lines[0]
    check(header == "figure,panel,row,field,value,source_table,specification_label", "unexpected Source Data header")
    parsed = [next(csv.reader([ln])) for ln in old_lines[1:]]
    ratio = [p for p in parsed if p[1] == "a"]
    bc_lines = [ln for ln, p in zip(old_lines[1:], parsed) if p[1] in ("b", "c")]
    check(len(ratio) == 45 and len(bc_lines) == 28 and len(parsed) == 73, "unexpected submitted Source Data rows")
    check(all(p[6] == "" and p[3] in ("share_visibility", "lo", "hi") for p in ratio), "unexpected ratio rows")
    buf = io.StringIO()
    wr = csv.writer(buf, lineterminator="\r\n")
    for r in rows_a:
        wr.writerow(["Fig5", "a"] + r)
    for p in ratio:
        wr.writerow(p[:6] + [SPEC_RATIO])
    sd_text = header + "\r\n" + buf.getvalue() + "".join(ln + "\r\n" for ln in bc_lines)
    back = list(csv.DictReader(io.StringIO(sd_text, newline="")))
    index = {(r["panel"], r["row"], r["field"]): r for r in back}
    check(len(index) == len(back), "duplicate (panel, row, field) in the written Source Data")
    for key, field, raw, value in ledger:
        r = index.get(("a", key, field))
        check(r is not None and r["value"] == raw and float(r["value"]) == value and r["specification_label"]
              == SPEC_PLOT, "drawn value not in Source Data: %s %s" % (key, field))
    check(len(ledger) == 15 * 12 and len({(k, f) for k, f, _, _ in ledger}) == len(ledger), "expected 180 drawn values")
    check(sd_text.split("\r\n")[-1 - len(bc_lines):-1] == bc_lines, "b and c rows are not byte-identical")
    bc_report = check_bc_against_supplemental(read_csv(old_sd))

    # 6. write (this folder only)
    OUT_PDF.write_bytes(pdf_bytes)
    pix = fin.get_pixmap(dpi=PNG_DPI, alpha=False)
    pix.set_dpi(PNG_DPI, PNG_DPI)
    pix.save(str(OUT_PNG))
    OUT_SD.write_bytes(sd_text.encode("utf-8"))

    print("inputs: stored table sha256 %s...; submitted Fig. 5 %s...; submitted Source Data %s..."
          % (PINNED["share"][:12], PINNED["old_pdf"][:12], PINNED["old_sd"][:12]))
    print("submitted page partition (spans, paths): a %s, b %s, c %s" % (partition["a"], partition["b"], partition["c"]))
    print("panel a: %d drawn values, all logged and found verbatim in Source Data" % len(ledger))
    print("I_fig5a_estimand tables: %s" % item_checks)
    print("Supplemental Note 20 Table S14.1 (3 decimals): %d cells checked, %d mismatches %s"
          % (note20.get("checked", 0), len(note20.get("mismatches", [])), note20.get("mismatches", "")))
    print("panels b/c Source Data vs Supplemental Data 5/6: %s" % bc_report)
    print("b and c: text spans, paths and 200-dpi pixels identical to the submitted figure; cohort labels moved "
          "%.3f pt at most" % label_shift)
    print("page %.2f x %.2f pt (%.3f x %.3f in); fonts %s"
          % (fin.rect.width, fin.rect.height, fin.rect.width / 72, fin.rect.height / 72, sorted(fonts)))
    print("Source Data: %d rows (panel a absolute %d, ratio %d, b+c %d)" % (len(back), len(rows_a), len(ratio),
                                                                           len(bc_lines)))
    for p in (OUT_PDF, OUT_PNG, OUT_SD):
        print("wrote %s  sha256 %s" % (p.name, sha256(p)))


if __name__ == "__main__":
    main()
