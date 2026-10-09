#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main figure F1: the three axes the whole manuscript turns on, one panel each.

F1 is the first result the reader meets, so it has to carry the skeleton of the argument and
nothing else. The three axes (`reports/MANUSCRIPT_OUTLINE_MERGED_2026-09-17.md`, end of section 1)
are: (1) how the negative set is constructed, (2) how the background is chosen, (3) how precisely
the attribute was defined. Each gets one panel.

NOTHING IS RECOMPUTED HERE. Every point, every interval endpoint and every count drawn is one cell
of one stored CSV, read as a string and converted to float for placement only. The script creates
no statistic, no ratio, no difference, no mean. The only arithmetic is axis limits and marker
offsets. Two selections are made over stored cells and are declared as selections, not statistics:
the min and the max of a column (panel a's two range brackets) are themselves stored cells, named
with the dataset they come from in the audit.

PANEL SOURCES, FIXED BEFORE DRAWING
-----------------------------------
a) `results/ptm_detectability_share.csv` (the task brief named it
   `results/ptm_detectability_census_share.csv`; no such file exists - the stored table with those
   columns is `ptm_detectability_share.csv`, and the column that holds the share is
   `share_visibility` with `share_visibility_lo` / `share_visibility_hi`). One row per
   (dataset, negative construction), restricted to `feature_set == VIS10` and `usable == 1`, which
   is the only combination for which the share column is written at all.
   The quantity is the census's own registered definition, quoted from
   `scripts/run_ptm_detectability_census.py`:
       share_visibility = (AUC(VIS10) - 0.5) / (AUC(DIG25) - 0.5)
   on the primary (global) caliber - the fraction of the reachable above-chance discriminative
   power of the full 25-column in-silico digest description that the 10 pure-visibility columns
   already deliver. It is dimensionless, it is not clipped, and a value above 1 or below 0 means
   the visibility subset beat or fell below its own superset out of fold; both are drawn as stored.
b) `results/phase2d_claim_retest.csv`, the two calibers of `PERS-008` and nothing else:
   `primary` (background = every protein in the assembly) and `zf_background` (background = the
   authors' own eligible list). Three estimates per background - `baseline_log2_or`,
   `matched_log2_or`, `random_control_log2_or` - each with its stored interval.
c) `results/phase2e_claim_retest.csv`, `caliber == primary`, the five structural claims, with the
   attribute composition columns `attribute_rate_positive` / `attribute_rate_negative` and the
   verbatim `author_reported_composition` string. The grouping dimension - how precisely the
   original stated its attribute - is TRANSCRIBED from the table in section 1 of
   `reports/PHASE2E_STRUCTURAL_CLAIMS.md`, in that table's own row order, because that report is
   where the tiering was declared; the script does not rank anything itself. Each transcription is
   re-checked against a verbatim Chinese fragment of that report at run time (see `recheck`).
   The composition columns carry NO stored interval, so panel c draws points without intervals;
   computing one here would be inventing a statistic. The interval-bearing quantity for these five
   claims is the log2 odds ratio, and that is where they appear, in F2.

Two notes on the brief that produced this figure, so the difference is on the record:
  * it named `results/ptm_detectability_census_share.csv`; no such file exists and the stored table
    with those columns is `results/ptm_detectability_share.csv`, which is what is read.
  * it gave panel b's collapse as "3.4971 falls to 0.0727 [-0.9228, 0.9198]". Those numbers are
    correct and are drawn, but they are two BASELINES measured against two backgrounds, not a
    baseline and a control; the second background's own stored verdict is `undecidable`, because
    after matching (0.2428 [-0.9365, 1.4026]) the same-size random control (-0.0069
    [-1.2511, 1.3116]) also crosses zero at 474 proteins. The panel therefore draws all three
    estimates for BOTH backgrounds and prints each background's stored verdict word, rather than
    the two baselines alone.

WHY THE THREE PANELS MAY NOT SHARE AN AXIS
------------------------------------------
The three quantities are not commensurable and no renormalisation onto a common axis is honest:
  a  is a dimensionless RATIO OF TWO AUC EXCESSES on a per-dataset out-of-fold score set. Its zero
     means "the visibility columns carry none of the reachable signal"; its 1 means "they carry all
     of it". It has no sign convention tied to a claim and it is undefined when the denominator is
     small (the census writes NA below 0.02 excess).
  b  is a log2 ODDS RATIO of a 2x2 table (claimed attribute x author's positive label,
     Haldane-Anscombe 0.5), the estimand registered for every claim re-test round. Its zero means
     "no association"; it is unbounded and signed.
  c  is a PROPORTION OF SITES carrying the claimed attribute, separately inside the positive and
     the negative set. It lives in [0, 1] and its meaningful feature is the gap between the two
     markers and how that gap compares with the composition the authors reported.
  A shared axis would require dividing one by another or z-scoring them, which invents a statistic
  this figure is forbidden to invent, and would put "0" at three different meanings on one line.
  Each panel therefore states its own quantity in its own x-label.

WHY PANEL a's TWO NEGATIVE CONSTRUCTIONS MAY NOT BE POOLED INTO ONE RANGE
------------------------------------------------------------------------
NEG_A ("not observed") and NEG_B ("observed unmodified") are not two samples of one quantity; they
are two different questions asked of the same positives. Under NEG_A the negatives are every
cysteine of every positive-carrying protein minus the positives, so absence of a call cannot be
told from absence of detection, and the share measures how much of that readout is visibility.
Under NEG_B the negatives are only cysteines the same experiment reported unmodified, so detection
is held roughly fixed by construction and the share measures something else entirely. Quoting one
interval that spans both would be the exact error this manuscript is about. The panel therefore
(i) draws them with different markers, (ii) joins the pair belonging to one dataset with a thin
grey line so the per-dataset move is still readable, and (iii) prints TWO separate range brackets,
one per construction, each labelled with its own n of datasets. There is no pooled bracket.

One dataset has NO NEG_B at all - `natcomm2023_ath_sno`, with
`usable = 0` and `exclusion_reason = "negatives 0 < 30"` because the experiment reported zero
observed-unmodified cysteines. Its rows carry an explicit "NEG_B does not exist" note rather than
a blank, so an empty position is never read as a zero share.

DRAWING CONVENTIONS - inherited from scripts/plot_manuscript_f3_f4.py and
scripts/plot_manuscript_f2_claim_retest.py
------------------------------------------------------------------------
Font sizes, open spines, figure dpi, the GREY / INK / ACCENT palette roles, the
`interval_bars` helper (line with end caps plus a point marker), the verdict colour map and its
text, the italic in-axis group headers, the 600 dpi PNG plus vector PDF, and the audit JSON shape
including a `verification` key are all the same as those two scripts. Colour never carries
information the labels do not; no significance marker is drawn because no test is reported.

Interpreter: a separate plotting environment (python 3.9.6 / matplotlib 3.9.4), the same pair that drew F2,
F3 and F4; the project venv has no matplotlib.

Writes results/manuscript_f1_three_axes.png (600 dpi), .pdf (vector) and
results/manuscript_f1_audit.json. Overwrites no earlier product.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
REPORTS = os.path.join(ROOT, "reports")
SCRIPT = os.path.abspath(__file__)

SHARE_TABLE = "ptm_detectability_share.csv"
PHASE2D_TABLE = "phase2d_claim_retest.csv"
PHASE2E_TABLE = "phase2e_claim_retest.csv"
PHASE2E_REPORT = "PHASE2E_STRUCTURAL_CLAIMS.md"

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 200,
    "axes.unicode_minus": False,
})
GREY, INK, ACCENT = "#9a9a9a", "#222222", "#1f5fa9"
# same map as F2, so a verdict word means the same colour in both main figures
VERDICT_COLOUR = {
    "survives": ACCENT,
    "attenuated": "#6fa3d8",
    "null_broken_by_control": "#7b3f98",
    "vanishes": "#b26a00",
    "baseline_contradicts_claim": "#b0323a",
    "undecidable": GREY,
}
VERDICT_TEXT = {
    "survives": "survives",
    "attenuated": "attenuated",
    "null_broken_by_control": "null broken by control",
    "vanishes": "vanishes",
    "baseline_contradicts_claim": "baseline contradicts claim",
    "undecidable": "undecidable (power-limited)",
}

# Panel a: display names and row order, copied from scripts/plot_ptm_detectability_census_relabelled.py
# so the same dataset reads the same way in F1 and in the census figure.
DATASET_LABEL = {
    "qtrp_S1_ph5": "QTRP S1 pH5 (human)",
    "qtrp_S2_ph5": "QTRP S2 pH5 (human)",
    "qpers_sid_tierB": "qPerS-SID tier B (human)",
    "cysboost2019_human_sno_hela": "Cys-BOOST HeLa (human)",
    "cysboost2019_human_sno_shsy5y": "Cys-BOOST SH-SY5Y (human)",
    "natcomm2023_ath_sno": "FAT-switch (Arabidopsis)",
    "abiotech2025_ath_sno": "PAT-switch (Arabidopsis)",
    "fps2020_ath_sulfenyl": "YAP1C reporter (Arabidopsis)",
}
DATASET_ORDER = ["qtrp_S1_ph5", "qtrp_S2_ph5", "qpers_sid_tierB",
                 "cysboost2019_human_sno_hela", "cysboost2019_human_sno_shsy5y",
                 "natcomm2023_ath_sno", "abiotech2025_ath_sno", "fps2020_ath_sulfenyl"]
FAMILY_SHORT = {"persulfidation": "persulfidation",
                "s_nitrosylation": "S-nitrosylation",
                "sulfenylation": "sulfenylation"}
NEG_A, NEG_B = "NEG_A_not_observed", "NEG_B_observed_unmodified"
PANEL_A_XLIM = (-1.18, 1.14)

# Panel b: the two backgrounds of PERS-008, in the order the round table stores them.
PERS008_BACKGROUND = [
    ("primary", "background = every protein in the assembly",
     "as published: the whole assembled proteome"),
    ("zf_background", "background = the authors' own eligible list",
     "the list the paper itself says was eligible"),
]
ESTIMATE_ROWS = [
    ("baseline_log2_or", "baseline_ci_low", "baseline_ci_high",
     "baseline", "baseline, the paper's own caliber"),
    ("matched_log2_or", "matched_ci_low", "matched_ci_high",
     "after matching", "after detectability matching (registered primary caliber)"),
    ("random_control_log2_or", "random_control_ci_low", "random_control_ci_high",
     "random control", "same-size random control (power reference)"),
]

# Panel c: transcribed from the table in section 1 of reports/PHASE2E_STRUCTURAL_CLAIMS.md, in that
# table's own row order. `fragment` is a verbatim slice of that report and is checked at run time;
# `tier` is its English rendering for the figure. The script ranks nothing itself.
PRECISION_TIERS = [
    {"claim_id": "SFE-001", "fragment": "量与阈值都写明",
     "tier": "quantity and threshold both stated\n(relative accessible area > 25%)"},
    {"claim_id": "SNO-012", "fragment": "阈值写明但是",
     "tier": "threshold stated, but as absolute area\n(sulfur SASA \u2264 1.0 \u00c5\u00b2)"},
    {"claim_id": "SNO-002", "fragment": "阈值没写",
     "tier": "no threshold stated\n(Yang's 25% borrowed)"},
    {"claim_id": "SNO-001", "fragment": "类别写明但工具（DSSP）不可复现",
     "tier": "categories stated, assignment tool\nnot reproducible (DSSP vs P-SEA)"},
    {"claim_id": "SNO-009", "fragment": "窗口写明、统计量是逐位频率",
     "tier": "window stated, statistic is a\nper-position frequency"},
]
PANEL_C_DATA_MAX = 0.72          # end of the data region; text column starts to the right of it
PANEL_C_TEXT_X = 0.75
PANEL_C_XLIM = (-0.015, 1.52)


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(directory, name):
    with open(os.path.join(directory, name), encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def cell(row, key):
    """Return the raw string of a cell, or None when the column is absent or empty."""
    value = row.get(key)
    if value is None or value.strip() == "":
        return None
    return value.strip()


def interval_bars(ax, ys, points, lows, highs, colour=INK, marker="o", filled=True, ms=3.8):
    """The F3/F4 helper: a line with end caps and a point. Nothing is computed."""
    for y, point, low, high in zip(ys, points, lows, highs):
        ax.plot([low, high], [y, y], color=colour, lw=1.1, solid_capstyle="butt", zorder=2)
        ax.plot([low, low], [y - 0.11, y + 0.11], color=colour, lw=1.1, zorder=2)
        ax.plot([high, high], [y - 0.11, y + 0.11], color=colour, lw=1.1, zorder=2)
        if filled:
            ax.plot([point], [y], marker=marker, color=colour, ms=ms, zorder=4)
        else:
            ax.plot([point], [y], marker=marker, mfc="white", mec=colour, mew=1.0, ms=ms,
                    ls="none", zorder=4)


# --------------------------------------------------------------------------------------- panel a

def collect_panel_a():
    """Read every VIS10 row of the share table. No filtering beyond feature_set."""
    rows = [r for r in read(RESULTS, SHARE_TABLE) if r["feature_set"] == "VIS10"]
    by_key = {(r["dataset_id"], r["negative_construction"]): r for r in rows}
    assert len(by_key) == len(rows), "duplicate (dataset, construction) in the share table"
    datasets = sorted({r["dataset_id"] for r in rows})
    assert set(datasets) == set(DATASET_ORDER), (
        "dataset order list is out of sync with the table: %s" % sorted(
            set(datasets) ^ set(DATASET_ORDER)))
    records = []
    for dataset in DATASET_ORDER:
        entry = {"dataset_id": dataset, "label": DATASET_LABEL[dataset]}
        for tag, construction in (("A", NEG_A), ("B", NEG_B)):
            row = by_key.get((dataset, construction))
            if row is None:
                entry[tag] = None
                continue
            usable = row["usable"] == "1"
            share = cell(row, "share_visibility")
            if not usable or share is None:
                entry[tag] = {"defined": False,
                              "exclusion_reason": cell(row, "exclusion_reason"),
                              "n_positive": cell(row, "n_positive"),
                              "n_negative": cell(row, "n_negative")}
                continue
            entry[tag] = {
                "defined": True,
                "share__raw": share, "share": float(share),
                "lo__raw": cell(row, "share_visibility_lo"),
                "lo": float(cell(row, "share_visibility_lo")),
                "hi__raw": cell(row, "share_visibility_hi"),
                "hi": float(cell(row, "share_visibility_hi")),
                "n_positive": cell(row, "n_positive"),
                "n_negative": cell(row, "n_negative"),
                "family": row["chemistry_family"],
                "grouping_unit": row["grouping_unit"],
            }
        entry["family"] = entry["A"]["family"] if entry["A"] and entry["A"]["defined"] else None
        records.append(entry)
    return records


def range_bracket(records, tag):
    """Pick the smallest and the largest STORED share cell of one construction.

    This is a selection over stored cells, not a statistic: both endpoints are cells that exist in
    the table, and the dataset each one comes from is written into the audit. project notes rule 12
    requires any 'range' statement to be taken from the complete stored table, which is what the
    caller passes in.
    """
    defined = [(r["dataset_id"], r[tag]) for r in records if r[tag] and r[tag]["defined"]]
    lo_id, lo_rec = min(defined, key=lambda item: item[1]["share"])
    hi_id, hi_rec = max(defined, key=lambda item: item[1]["share"])
    return {"n_datasets": len(defined),
            "min_dataset": lo_id, "min__raw": lo_rec["share__raw"], "min": lo_rec["share"],
            "max_dataset": hi_id, "max__raw": hi_rec["share__raw"], "max": hi_rec["share"],
            "datasets": [item[0] for item in defined]}


def draw_panel_a(ax, records, drawn):
    ticks, labels, off_axis = [], [], []
    y = 0.0
    for record in records:
        a, b = record["A"], record["B"]
        if a and a["defined"] and b and b["defined"]:
            ax.plot([b["share"], a["share"]], [y, y], color="#c6c6c6", lw=0.9, zorder=1)
        if a and a["defined"]:
            interval_bars(ax, [y], [a["share"]], [a["lo"]], [a["hi"]], INK, "o", True, 4.0)
            for field, key in (("share_visibility", "share"), ("share_visibility_lo", "lo"),
                               ("share_visibility_hi", "hi")):
                drawn.append({"panel": "a", "table": "results/" + SHARE_TABLE,
                              "keys": {"dataset_id": record["dataset_id"],
                                       "negative_construction": NEG_A, "feature_set": "VIS10"},
                              "field": field, "value": a[key], "raw": a[key + "__raw"]})
        if b and b["defined"]:
            interval_bars(ax, [y], [b["share"]], [b["lo"]], [b["hi"]], ACCENT, "s", False, 4.2)
            for field, key in (("share_visibility", "share"), ("share_visibility_lo", "lo"),
                               ("share_visibility_hi", "hi")):
                drawn.append({"panel": "a", "table": "results/" + SHARE_TABLE,
                              "keys": {"dataset_id": record["dataset_id"],
                                       "negative_construction": NEG_B, "feature_set": "VIS10"},
                              "field": field, "value": b[key], "raw": b[key + "__raw"]})
            if b["lo"] < PANEL_A_XLIM[0]:
                # never let an interval stop at the axis edge without saying so
                ax.annotate("NEG_B interval runs to %s, beyond this axis" % b["lo__raw"],
                            xy=(PANEL_A_XLIM[0] + 0.02, y + 0.26), fontsize=6.0,
                            color="#5b6068", ha="left", va="bottom")
                off_axis.append({"dataset_id": record["dataset_id"], "field":
                                 "share_visibility_lo", "raw": b["lo__raw"]})
        elif b is not None:
            # the construction does not exist for this dataset: say so, never leave it blank
            ax.text(PANEL_A_XLIM[0] + 0.03, y, "NEG_B does not exist here (%s)"
                    % (b["exclusion_reason"] or "no reason stored"),
                    fontsize=6.0, color="#b26a00", ha="left", va="center", style="italic")
        ticks.append(y)
        labels.append(record["label"])
        y -= 1.0

    y -= 0.60
    brackets = {}
    for tag, construction, colour, name in (("A", NEG_A, INK, "NEG_A"),
                                            ("B", NEG_B, ACCENT, "NEG_B")):
        info = range_bracket(records, tag)
        brackets[name] = info
        ax.plot([info["min"], info["max"]], [y, y], color=colour, lw=3.2, alpha=0.30,
                solid_capstyle="butt", zorder=1)
        ax.plot([info["min"], info["min"]], [y - 0.17, y + 0.17], color=colour, lw=1.2)
        ax.plot([info["max"], info["max"]], [y - 0.17, y + 0.17], color=colour, lw=1.2)
        ticks.append(y)
        labels.append("%s range, %d datasets\n%s to %s"
                      % (name, info["n_datasets"], info["min__raw"], info["max__raw"]))
        y -= 1.15
    spine_bottom = y + 1.15 - 0.45

    ax.axvline(0.0, color=INK, lw=0.8, zorder=0)
    ax.axvline(1.0, color=GREY, ls=":", lw=0.8, zorder=0)
    ax.text(1.0, 0.70, "all of it", fontsize=6.0, color="#5b6068", ha="center", va="bottom")
    ax.text(0.0, 0.70, "none of it", fontsize=6.0, color="#5b6068", ha="center", va="bottom")
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_xlim(*PANEL_A_XLIM)
    ax.set_ylim(y - 2.05, 1.25)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_bounds(spine_bottom, 1.25)
    ax.set_xlabel("share of the reachable discriminative power carried by\n"
                  "visibility features alone:  (AUC$_{VIS10}$ - 0.5) / (AUC$_{DIG25}$ - 0.5)")
    ax.set_title("a  how the negative set is constructed", loc="left", pad=20)
    ax.text(0.0, 1.015, "the two constructions ask two different questions and are never pooled:\n"
                        "each carries its own range bracket, there is no combined one",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.4, color="#5b6068")
    handles = [
        plt.Line2D([], [], marker="o", color=INK, ls="none", ms=4.0,
                   label="NEG_A  negatives = cysteines never observed"),
        plt.Line2D([], [], marker="s", mfc="white", mec=ACCENT, mew=1.0, ls="none", ms=4.2,
                   label="NEG_B  negatives = seen unmodified in the same run"),
        plt.Line2D([], [], color=INK, lw=1.0,
                   label="95% interval, cluster bootstrap 5000, seed 20260915"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False, handletextpad=0.5,
              bbox_to_anchor=(-0.04, -0.012), fontsize=6.3, borderaxespad=0.0)
    return brackets, off_axis


# --------------------------------------------------------------------------------------- panel b

def collect_panel_b():
    rows = [r for r in read(RESULTS, PHASE2D_TABLE) if r["claim_id"] == "PERS-008"]
    by_caliber = {r["caliber"]: r for r in rows}
    assert set(by_caliber) == {c for c, _, _ in PERS008_BACKGROUND}, (
        "PERS-008 calibers on disk are %s" % sorted(by_caliber))
    out = []
    for caliber, title, subtitle in PERS008_BACKGROUND:
        row = by_caliber[caliber]
        entry = {"caliber": caliber, "title": title, "subtitle": subtitle,
                 "verdict": row["verdict"], "n_observations": cell(row, "n_observations"),
                 "n_positive": cell(row, "n_positive"), "n_attribute": cell(row, "n_attribute"),
                 "attribute_label": cell(row, "attribute_label"),
                 "baseline_reproduction": cell(row, "baseline_reproduction"), "estimates": []}
        for point_key, lo_key, hi_key, short, name in ESTIMATE_ROWS:
            entry["estimates"].append({
                "name": name, "short": short, "point_field": point_key,
                "point__raw": cell(row, point_key), "point": float(cell(row, point_key)),
                "lo_field": lo_key, "lo__raw": cell(row, lo_key), "lo": float(cell(row, lo_key)),
                "hi_field": hi_key, "hi__raw": cell(row, hi_key), "hi": float(cell(row, hi_key)),
            })
        out.append(entry)
    return out


def draw_panel_b(ax, blocks, drawn):
    ticks, labels = [], []
    y = 0.0
    for block in blocks:
        colour = VERDICT_COLOUR[block["verdict"]]
        ax.text(0.0, y, "%s\n%s  -  %s" % (block["title"], block["subtitle"],
                                           VERDICT_TEXT[block["verdict"]]),
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6.4,
                color=colour, style="italic", clip_on=False)
        y -= 1.15
        for index, estimate in enumerate(block["estimates"]):
            if index == 0:
                interval_bars(ax, [y], [estimate["point"]], [estimate["lo"]], [estimate["hi"]],
                              INK, "o", False, 4.2)
            elif index == 1:
                interval_bars(ax, [y], [estimate["point"]], [estimate["lo"]], [estimate["hi"]],
                              colour, "o", True, 4.6)
            else:
                interval_bars(ax, [y], [estimate["point"]], [estimate["lo"]], [estimate["hi"]],
                              "#5b6068", "x", True, 4.4)
            for field_key, value_key in (("point_field", "point"), ("lo_field", "lo"),
                                         ("hi_field", "hi")):
                drawn.append({"panel": "b", "table": "results/" + PHASE2D_TABLE,
                              "keys": {"claim_id": "PERS-008", "caliber": block["caliber"]},
                              "field": estimate[field_key], "value": estimate[value_key],
                              "raw": estimate[value_key + "__raw"]})
            ticks.append(y)
            labels.append(estimate["short"])
            y -= 1.0
        ax.text(0.0, y + 0.45, "%s proteins, %s positives, %s carry the attribute"
                % (block["n_observations"], block["n_positive"], block["n_attribute"]),
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6.0,
                color="#5b6068", clip_on=False)
        y -= 1.05
    spine_bottom = y + 1.05 - 0.45
    ax.axvline(0.0, color=INK, lw=0.8, zorder=1)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_ylim(y - 2.35, 1.0)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_bounds(spine_bottom, 1.0)
    ax.margins(x=0.10)
    ax.set_xlabel("effect size, log2 odds ratio\n(claimed attribute x author's positive label)")
    ax.set_title("b  how the background is chosen", loc="left", pad=20)
    ax.text(0.0, 1.015, "PERS-008 only: annotation to ubiquitin-dependent proteolysis,\n"
                        "one claim measured against two backgrounds",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.4, color="#5b6068")
    handles = [
        plt.Line2D([], [], marker="o", mfc="white", mec=INK, mew=1.0, ls="none", ms=4.2,
                   label="baseline, the paper's own caliber"),
        plt.Line2D([], [], marker="o", color=INK, ls="none", ms=4.6,
                   label="after detectability matching (primary caliber)"),
        plt.Line2D([], [], marker="x", color="#5b6068", ls="none", ms=4.4, mew=1.0,
                   label="same-size random control (power reference)"),
        plt.Line2D([], [], color=INK, lw=1.0,
                   label="95% interval, cluster bootstrap 5000, seed 20260915"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False, handletextpad=0.5,
              bbox_to_anchor=(-0.09, -0.012), fontsize=6.3, borderaxespad=0.0)


# --------------------------------------------------------------------------------------- panel c

def collect_panel_c():
    rows = [r for r in read(RESULTS, PHASE2E_TABLE) if r["caliber"] == "primary"]
    by_claim = {r["claim_id"]: r for r in rows}
    wanted = [tier["claim_id"] for tier in PRECISION_TIERS]
    assert set(by_claim) == set(wanted), (
        "phase2e primary rows on disk are %s" % sorted(by_claim))
    out = []
    for tier in PRECISION_TIERS:
        row = by_claim[tier["claim_id"]]
        out.append({
            "claim_id": tier["claim_id"], "tier": tier["tier"], "fragment": tier["fragment"],
            "verdict": row["verdict"], "unit": cell(row, "unit"),
            "n_observations": cell(row, "n_observations"), "n_positive": cell(row, "n_positive"),
            "attribute_label": cell(row, "attribute_label"),
            "author_reported_composition": cell(row, "author_reported_composition"),
            "pos__raw": cell(row, "attribute_rate_positive"),
            "pos": float(cell(row, "attribute_rate_positive")),
            "neg__raw": cell(row, "attribute_rate_negative"),
            "neg": float(cell(row, "attribute_rate_negative")),
        })
    return out


def draw_panel_c(ax, records, drawn):
    ticks, labels = [], []
    y = 0.0
    for record in records:
        colour = VERDICT_COLOUR[record["verdict"]]
        ax.text(0.012, y + 0.62, record["tier"], transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=6.3, color=colour, style="italic", clip_on=False)
        ax.plot([record["neg"], record["pos"]], [y, y], color=colour, lw=1.2, alpha=0.55, zorder=2)
        ax.plot([record["neg"]], [y], marker="o", mfc="white", mec=colour, mew=1.0, ms=4.6,
                ls="none", zorder=4)
        ax.plot([record["pos"]], [y], marker="o", color=colour, ms=4.6, zorder=4)
        for field, key in (("attribute_rate_positive", "pos"),
                           ("attribute_rate_negative", "neg")):
            drawn.append({"panel": "c", "table": "results/" + PHASE2E_TABLE,
                          "keys": {"claim_id": record["claim_id"], "caliber": "primary"},
                          "field": field, "value": record[key], "raw": record[key + "__raw"]})
        ax.text(PANEL_C_TEXT_X, y + 0.16, "authors reported  %s"
                % record["author_reported_composition"],
                fontsize=6.4, color=INK, ha="left", va="center")
        ax.text(PANEL_C_TEXT_X, y - 0.24, "verdict  %s   (%s sites, %s positive)"
                % (VERDICT_TEXT[record["verdict"]], record["n_observations"],
                   record["n_positive"]),
                fontsize=6.2, color=colour, ha="left", va="center")
        ticks.append(y)
        labels.append(record["claim_id"])
        y -= 1.38

    spine_bottom = y + 1.38 - 0.55
    ax.plot([PANEL_C_DATA_MAX, PANEL_C_DATA_MAX], [spine_bottom, 1.05], color="#dcdcdc", lw=0.8,
            zorder=0)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
    ax.set_xlim(*PANEL_C_XLIM)
    ax.set_ylim(spine_bottom - 0.20, 1.05)
    ax.tick_params(axis="y", length=0)
    ax.spines["bottom"].set_bounds(0.0, PANEL_C_DATA_MAX)
    ax.spines["left"].set_bounds(spine_bottom, 1.05)
    ax.set_xlabel("fraction of sites carrying the claimed attribute, as measured here")
    ax.set_title("c  how precisely the attribute was defined", loc="left", pad=20)
    ax.text(0.0, 1.015, "five structural claims, one feature table, one instrument, one control;\n"
                        "rows in the order of the tiering table in PHASE2E_STRUCTURAL_CLAIMS.md",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.4, color="#5b6068")
    handles = [
        plt.Line2D([], [], marker="o", color=INK, ls="none", ms=4.6,
                   label="positive sites"),
        plt.Line2D([], [], marker="o", mfc="white", mec=INK, mew=1.0, ls="none", ms=4.6,
                   label="negative sites"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, handletextpad=0.5,
              bbox_to_anchor=(1.0, 1.015), fontsize=6.4, ncol=2, columnspacing=1.4,
              borderaxespad=0.0)


# ------------------------------------------------------------------------------------ the figure

def figure_f1(panel_a, panel_b, panel_c):
    drawn = []
    fig = plt.figure(figsize=(7.6, 7.9))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.16, 1.0], height_ratios=[1.22, 1.0],
                            wspace=0.38, hspace=0.42)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    brackets, off_axis = draw_panel_a(ax_a, panel_a, drawn)
    draw_panel_b(ax_b, panel_b, drawn)
    draw_panel_c(ax_c, panel_c, drawn)

    fig.subplots_adjust(left=0.175, right=0.985, top=0.945, bottom=0.055)
    for ext, kwargs in (("png", {"dpi": 600}), ("pdf", {})):
        fig.savefig(os.path.join(RESULTS, "manuscript_f1_three_axes." + ext),
                    bbox_inches="tight", pad_inches=0.12, **kwargs)
    plt.close(fig)
    return drawn, brackets, off_axis


# --------------------------------------------------------------------------------------- recheck

def recheck(drawn, panel_a, panel_b, panel_c, brackets):
    """Find every drawn number back in its source table, cell by cell, by its key columns."""
    cache = {SHARE_TABLE: read(RESULTS, SHARE_TABLE),
             PHASE2D_TABLE: read(RESULTS, PHASE2D_TABLE),
             PHASE2E_TABLE: read(RESULTS, PHASE2E_TABLE)}
    misses, checked = [], 0
    for item in drawn:
        table = cache[os.path.basename(item["table"])]
        found = False
        for row in table:
            if any(row.get(key) != value for key, value in item["keys"].items()):
                continue
            raw = cell(row, item["field"])
            if raw is not None and raw == item["raw"] and float(raw) == item["value"]:
                found = True
                break
        checked += 1
        if not found:
            misses.append(item)

    # counts and ranges, recomputed here from the COMPLETE stored table, not from the drawn set
    share_rows = [r for r in read(RESULTS, SHARE_TABLE) if r["feature_set"] == "VIS10"]
    usable = {}
    for construction in (NEG_A, NEG_B):
        usable[construction] = sorted(
            r["dataset_id"] for r in share_rows
            if r["negative_construction"] == construction and r["usable"] == "1"
            and cell(r, "share_visibility") is not None)
    unusable_b = sorted(
        (r["dataset_id"], cell(r, "exclusion_reason"), cell(r, "n_negative"))
        for r in share_rows if r["negative_construction"] == NEG_B and r["usable"] != "1")
    extremes = {}
    for name, construction in (("NEG_A", NEG_A), ("NEG_B", NEG_B)):
        values = [(float(cell(r, "share_visibility")), r["dataset_id"],
                   cell(r, "share_visibility"))
                  for r in share_rows if r["negative_construction"] == construction
                  and r["usable"] == "1" and cell(r, "share_visibility") is not None]
        low, high = min(values), max(values)
        extremes[name] = {"min": low[0], "min_dataset": low[1], "min__raw": low[2],
                          "max": high[0], "max_dataset": high[1], "max__raw": high[2]}
    range_checks = []
    for name in ("NEG_A", "NEG_B"):
        drawn_bracket, truth = brackets[name], extremes[name]
        range_checks.append({
            "construction": name,
            "n_datasets_drawn": drawn_bracket["n_datasets"],
            "n_datasets_in_table": len(usable[NEG_A if name == "NEG_A" else NEG_B]),
            "min_drawn": drawn_bracket["min__raw"], "min_in_table": truth["min__raw"],
            "min_dataset_drawn": drawn_bracket["min_dataset"],
            "min_dataset_in_table": truth["min_dataset"],
            "max_drawn": drawn_bracket["max__raw"], "max_in_table": truth["max__raw"],
            "max_dataset_drawn": drawn_bracket["max_dataset"],
            "max_dataset_in_table": truth["max_dataset"],
            "agree": (drawn_bracket["min__raw"] == truth["min__raw"]
                      and drawn_bracket["max__raw"] == truth["max__raw"]
                      and drawn_bracket["min_dataset"] == truth["min_dataset"]
                      and drawn_bracket["max_dataset"] == truth["max_dataset"]
                      and drawn_bracket["n_datasets"]
                      == len(usable[NEG_A if name == "NEG_A" else NEG_B])),
        })

    # every verdict word printed on the figure must be the stored verdict of that stored row
    verdict_checks = []
    for row in read(RESULTS, PHASE2D_TABLE):
        if row["claim_id"] != "PERS-008":
            continue
        block = next(b for b in panel_b if b["caliber"] == row["caliber"])
        verdict_checks.append({"where": "b/PERS-008/" + row["caliber"],
                               "drawn": block["verdict"], "in_table": row["verdict"],
                               "agree": block["verdict"] == row["verdict"]})
    for row in read(RESULTS, PHASE2E_TABLE):
        if row["caliber"] != "primary":
            continue
        record = next(r for r in panel_c if r["claim_id"] == row["claim_id"])
        verdict_checks.append({"where": "c/" + row["claim_id"],
                               "drawn": record["verdict"], "in_table": row["verdict"],
                               "agree": record["verdict"] == row["verdict"]})

    # the precision tiering is transcribed from a report; confirm each source fragment is there
    with open(os.path.join(REPORTS, PHASE2E_REPORT), encoding="utf-8") as fh:
        report_text = fh.read()
    tier_checks = [{"claim_id": tier["claim_id"], "fragment": tier["fragment"],
                    "found_in_report": tier["fragment"] in report_text,
                    "report": "reports/" + PHASE2E_REPORT}
                   for tier in PRECISION_TIERS]

    # the author-reported composition strings are quoted verbatim, so quote-check them too
    quote_checks = []
    for row in read(RESULTS, PHASE2E_TABLE):
        if row["caliber"] != "primary":
            continue
        record = next(r for r in panel_c if r["claim_id"] == row["claim_id"])
        quote_checks.append({"claim_id": row["claim_id"],
                             "drawn": record["author_reported_composition"],
                             "in_table": cell(row, "author_reported_composition"),
                             "agree": record["author_reported_composition"]
                             == cell(row, "author_reported_composition")})

    # every COUNT printed as text on the figure, not just the plotted points (project notes 9.5:
    # a re-check has to cover each count that appears in the running text, not only the estimates)
    count_strings = []
    d_rows = {r["caliber"]: r for r in read(RESULTS, PHASE2D_TABLE) if r["claim_id"] == "PERS-008"}
    for block in panel_b:
        row = d_rows[block["caliber"]]
        for field in ("n_observations", "n_positive", "n_attribute"):
            count_strings.append({"where": "b/PERS-008/" + block["caliber"], "field": field,
                                  "printed": block[field], "in_table": cell(row, field),
                                  "agree": block[field] == cell(row, field)})
    e_rows = {r["claim_id"]: r for r in read(RESULTS, PHASE2E_TABLE) if r["caliber"] == "primary"}
    for record in panel_c:
        row = e_rows[record["claim_id"]]
        for field in ("n_observations", "n_positive"):
            count_strings.append({"where": "c/" + record["claim_id"], "field": field,
                                  "printed": record[field], "in_table": cell(row, field),
                                  "agree": record[field] == cell(row, field)})
    a_rows = {(r["dataset_id"], r["negative_construction"]): r
              for r in share_rows}
    for record in panel_a:
        if record["B"] and not record["B"]["defined"]:
            row = a_rows[(record["dataset_id"], NEG_B)]
            count_strings.append({"where": "a/" + record["dataset_id"] + "/NEG_B",
                                  "field": "exclusion_reason",
                                  "printed": record["B"]["exclusion_reason"],
                                  "in_table": cell(row, "exclusion_reason"),
                                  "agree": record["B"]["exclusion_reason"]
                                  == cell(row, "exclusion_reason")})

    panel_a_counts = {
        "datasets_in_table": len({r["dataset_id"] for r in share_rows}),
        "datasets_drawn": len(panel_a),
        "neg_a_usable": len(usable[NEG_A]),
        "neg_b_usable": len(usable[NEG_B]),
        "neg_b_not_usable": [{"dataset_id": d, "exclusion_reason": reason, "n_negative": n}
                             for d, reason, n in unusable_b],
        "neg_b_marked_nonexistent_on_figure": sorted(
            r["dataset_id"] for r in panel_a if r["B"] and not r["B"]["defined"]),
    }
    panel_a_counts["agree"] = (
        panel_a_counts["datasets_in_table"] == panel_a_counts["datasets_drawn"] == 9
        and panel_a_counts["neg_a_usable"] == 9 and panel_a_counts["neg_b_usable"] == 7
        and panel_a_counts["neg_b_marked_nonexistent_on_figure"]
        == sorted(d for d, _, _ in unusable_b))

    return {
        "values_drawn": checked,
        "values_not_found_in_sources": misses,
        "panel_a_counts": panel_a_counts,
        "panel_a_range_brackets": range_checks,
        "counts_printed_as_text": count_strings,
        "verdict_words": verdict_checks,
        "precision_tier_transcription": tier_checks,
        "author_composition_quotes": quote_checks,
        "all_pass": (not misses
                     and panel_a_counts["agree"]
                     and all(c["agree"] for c in range_checks)
                     and all(c["agree"] for c in count_strings)
                     and all(c["agree"] for c in verdict_checks)
                     and all(c["found_in_report"] for c in tier_checks)
                     and all(c["agree"] for c in quote_checks)),
    }


def main():
    started = time.time()
    panel_a = collect_panel_a()
    panel_b = collect_panel_b()
    panel_c = collect_panel_c()
    drawn, brackets, off_axis = figure_f1(panel_a, panel_b, panel_c)
    verification = recheck(drawn, panel_a, panel_b, panel_c, brackets)

    inputs = ["results/" + SHARE_TABLE, "results/" + PHASE2D_TABLE, "results/" + PHASE2E_TABLE,
              "reports/" + PHASE2E_REPORT]
    audit = {
        "script": "scripts/plot_manuscript_f1_three_axes.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "round": "main figure F1, the three axes, drawn from stored tables only",
        "computes_nothing": True,
        "arithmetic_performed": "none on data; min/max in panel a are SELECTIONS of stored cells "
                               "and both endpoints are named with the dataset they come from",
        "panels": {
            "a": {
                "axis": "how the negative set is constructed",
                "quantity": "share_visibility = (AUC(VIS10) - 0.5) / (AUC(DIG25) - 0.5), primary "
                            "(global) caliber, dimensionless, unclipped",
                "definition_source": "scripts/run_ptm_detectability_census.py docstring",
                "table": "results/" + SHARE_TABLE,
                "rows_used": "feature_set == VIS10 (the only rows the share column is written for)",
                "datasets": DATASET_ORDER,
                "constructions": {
                    "NEG_A_not_observed": "every cysteine of every positive-carrying protein, "
                                          "minus the positives",
                    "NEG_B_observed_unmodified": "only cysteines the same experiment reported "
                                                 "unmodified, minus the positives",
                },
                "pooling": "forbidden; two separate range brackets are drawn and labelled, there "
                           "is no pooled bracket",
                "range_brackets": brackets,
                "neg_b_absent": [{"dataset_id": r["dataset_id"],
                                  "exclusion_reason": r["B"]["exclusion_reason"],
                                  "n_negative": r["B"]["n_negative"],
                                  "shown_as": "NEG_B does not exist here"}
                                 for r in panel_a if r["B"] and not r["B"]["defined"]],
                "intervals_running_off_the_axis": off_axis,
                "chemistry_family_per_dataset": {r["dataset_id"]: FAMILY_SHORT.get(r["family"],
                                                                                   r["family"])
                                                 for r in panel_a},
                "grouping_unit_per_dataset": {r["dataset_id"]: r["A"]["grouping_unit"]
                                              for r in panel_a if r["A"] and r["A"]["defined"]},
            },
            "b": {
                "axis": "how the background is chosen",
                "quantity": "log2 odds ratio of the 2x2 table "
                            "(claimed attribute x author's positive label), Haldane-Anscombe 0.5",
                "table": "results/" + PHASE2D_TABLE,
                "claim": "PERS-008",
                "calibers": [{"caliber": b["caliber"], "verdict": b["verdict"],
                              "n_observations": b["n_observations"],
                              "n_positive": b["n_positive"], "n_attribute": b["n_attribute"],
                              "baseline_reproduction": b["baseline_reproduction"],
                              "attribute_label": b["attribute_label"],
                              "estimates": {e["name"]: {"point": e["point__raw"],
                                                        "ci": [e["lo__raw"], e["hi__raw"]]}
                                            for e in b["estimates"]}}
                             for b in panel_b],
                "why_the_random_control_is_drawn": "matching removes proteins and removing "
                                                   "proteins weakens effects on its own; without "
                                                   "the same-size random control a widened "
                                                   "interval cannot be told from a killed effect",
            },
            "c": {
                "axis": "how precisely the attribute was defined",
                "quantity": "fraction of sites carrying the claimed attribute, inside the "
                            "positive set and inside the negative set; a proportion in [0, 1]",
                "table": "results/" + PHASE2E_TABLE,
                "rows_used": "caliber == primary",
                "grouping_dimension": "definition precision, TRANSCRIBED verbatim from the table "
                                      "in section 1 of reports/" + PHASE2E_REPORT + ", in that "
                                      "table's own row order; this script ranks nothing itself",
                "claims": [{"claim_id": r["claim_id"], "tier_drawn": r["tier"].replace("\n", " "),
                            "tier_source_fragment": r["fragment"], "verdict": r["verdict"],
                            "attribute_label": r["attribute_label"],
                            "attribute_rate_positive": r["pos__raw"],
                            "attribute_rate_negative": r["neg__raw"],
                            "author_reported_composition": r["author_reported_composition"],
                            "n_observations": r["n_observations"],
                            "n_positive": r["n_positive"]}
                           for r in panel_c],
                "author_values_not_plotted": "the authors' reported compositions are bounds "
                                             "(>60% / <30%), approximations (~35%) and one "
                                             "difference (+10 percentage points); turning them "
                                             "into plotted numbers would create numbers, so they "
                                             "are quoted verbatim as text. SNO-012's quote gives "
                                             "the positive side only (~35% buried), so that row "
                                             "has no author counterpart for the negative marker.",
                "no_intervals_drawn": "phase2e stores no interval columns for "
                                      "attribute_rate_positive / attribute_rate_negative, so the "
                                      "two markers of a row are point values and the panel draws "
                                      "no interval. Computing one here would be a new statistic. "
                                      "The interval-bearing quantity for these five claims is the "
                                      "log2 odds ratio, which is where they appear in F2.",
            },
        },
        "axes_not_shared": "a is a ratio of two AUC excesses, b a signed unbounded log2 odds "
                           "ratio, c a proportion in [0, 1]; zero means three different things "
                           "and any renormalisation onto one axis would invent a statistic",
        "deliberately_not_drawn": {
            "DIG25 / CPL15 rows of the share table": "the share is defined only for VIS10",
            "auc_global / auc_within_unit": "F1 shows the share, not the AUCs it is built from",
            "phase2d / phase2e stratified_log2_or": "registered SECONDARY caliber; the primary "
                                                    "caliber is matching",
            "phase2e plddt70 caliber": "declared sensitivity caliber; panel c draws the primary "
                                       "one only and the sensitivity result stays in the round "
                                       "report",
            "the other 33 claims": "F2 is the figure that carries every claim; F1 carries one "
                                   "worked example per axis",
        },
        "verification": verification,
        "outputs": ["results/manuscript_f1_three_axes.png",
                    "results/manuscript_f1_three_axes.pdf"],
        "inputs": {name: sha256_of(os.path.join(ROOT, name)) for name in inputs},
        "versions": {"python": sys.version, "matplotlib": matplotlib.__version__},
        "elapsed_seconds": round(time.time() - started, 2),
    }
    out = os.path.join(RESULTS, "manuscript_f1_audit.json")
    with open(out + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    os.replace(out + ".tmp", out)
    print("wrote F1:", verification["values_drawn"], "values re-checked,",
          len(verification["values_not_found_in_sources"]), "not found;",
          "panel a counts agree:", verification["panel_a_counts"]["agree"], ";",
          "range brackets agree:", all(c["agree"] for c in verification["panel_a_range_brackets"]),
          "; all_pass:", verification["all_pass"])


if __name__ == "__main__":
    main()
