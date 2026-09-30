#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main figure F2: every re-tested claim, before and after the detectability control.

Replaces results/phase2_claim_retest_fig.* (A-tier, 8 claims only), which stays on disk
untouched. The row set of this figure is results/claim_verdict_tally_2026-09-17.csv, the
de-duplicated 35-claim tally (primary caliber only, latest round wins), and the effect sizes
come from the round tables that tally points at. NOTHING IS RECOMPUTED HERE: every point and
every interval endpoint drawn is one cell of one of those CSV files, read as a string and
converted to float. No statistic, no ratio, no difference is created by this script; the only
arithmetic is axis limits and marker offsets.

PANEL LAYOUT AND WHY IT IS SPLIT THAT WAY
-----------------------------------------
All four rounds registered the same estimand - "log2 odds ratio of the 2x2 table (claimed
attribute x author's positive label), Haldane-Anscombe 0.5" (phase2 audit, `predeclared.estimand`,
inherited verbatim by 2b/2d/2e). So the calibers ARE commensurable and a single axis is
defensible, with ONE declared exception: PERS-010, whose registered estimand is
`log2(observed shared / expected shared under a detection-independence model)`
(phase2 audit, `predeclared.pers010_estimand`). That is a different quantity that happens to be
on a log2 scale, so it is NOT drawn on the odds-ratio axis; it gets its own panel (c) with its
own axis label. This is the "facet per caliber" option of the two the brief allowed; the
alternative - renormalising everything to a retention ratio - was rejected because the stored
`retention_ratio` column is |controlled|/|baseline|, which is undefined in sign and explodes
when the baseline is near zero (values of 1.42, 2.16, 3.34 sit in the tables), so it would turn
the honest "undecidable" rows into spurious large numbers.

Panels a and b split the shared odds-ratio axis by UNIT OF OBSERVATION, not by caliber: site-level
claims (a) and protein-level claims (b), the `unit` column of the source tables. Same statistic,
same axis label, different axis range - the protein-level GO claims run to log2 OR ~ 4.9 while the
sequence-motif claims live inside +/-2, and one shared range would compress panel a to illegibility.
Both panels carry the same zero line and the same legend, so a reader can still compare across them.

Rows inside a panel are grouped by verdict, in the fixed order
survives / attenuated / null_broken_by_control / vanishes / baseline_contradicts_claim / undecidable,
with the verdict also carried by colour.

Panel d is the 7 claims with NO measured baseline (`has_measured_baseline == 0`):
6 out_of_instrument_scope and 1 blocked_material_unreachable. They are listed as text in grey,
never as points, and never inside the 28 that have a baseline - per the tally report, the
denominator for the five verdict categories is 28, not 35.

WHAT EACH ROW SHOWS (the three things the figure has to say)
-----------------------------------------------------------
1. open circle + thin dark interval = the baseline effect, the paper's own caliber
   (`baseline_log2_or`, `baseline_ci_low/high`);
   filled circle + coloured interval = after detectability matching, the registered PRIMARY
   caliber (`matched_log2_or`, `matched_ci_low/high`); an arrow joins the two, so the move the
   control produces is the visual unit.
2. grey cross + grey interval = the same-size random control
   (`random_control_log2_or`, `random_control_ci_low/high`), 20 draws at seeds 20260915..34.
   It is the only thing that separates "the control killed the effect" from "the sample got
   small"; without it a widened interval is unreadable.
3. colour + grouping = the verdict from the tally table.
   The secondary caliber (depth stratification, `stratified_*`) is deliberately NOT drawn: the
   registered rule makes matching primary and stratification secondary, and drawing four
   estimates per row makes the rows unreadable. The `*_strata.csv` files hold the per-stratum
   detail behind that secondary caliber and are likewise not used.

Row tags, all read from stored columns, none inferred:
  2 / 2B / 2D / 2E  the round that produced the drawn row. For the eight A-tier claims the tally
     credits the round as phase2b because the combined table supersedes the A-tier one, but the
     combined table's own `round` column still says phase2; the label follows that column, so an
     A-tier row reads "2".
  T  `is_transfer` is True - the baseline is a transfer test on another cohort, not the authors'
     own data (SFE-006, SFE-007, SFE-008, SFE-011, SNO-016).
  S  the tally records a sensitivity caliber whose verdict differs from the primary one
     (PERS-008 zf_background, SNO-012 plddt70). Both labels must be quoted together.
  R  the tally records a different verdict in an earlier round (PERS-008: 2B
     blocked_material_unreachable -> 2D survives, after the background changed).
  n  `claim_direction == null_no_preference` - the claim asserts NO preference, so "survives"
     and "broken by control" read in the opposite direction.

Interpretation discipline carried over from the round reports: a verdict on a row whose
`baseline_reproduction` is not the authors' own cohort means less than one on a row that is;
that column is written into the audit JSON for every row rather than encoded as a glyph.

Interpreter: the scratchpad plotenv (python 3.9.6 / matplotlib 3.9.4), the same pair that drew
F3 and F4; the project venv has no matplotlib. Font sizes, spine style, palette role
(INK / GREY / ACCENT) and the audit-JSON shape follow scripts/plot_manuscript_f3_f4.py.

Writes results/manuscript_f2_claim_retest.png (600 dpi), .pdf (vector) and
results/manuscript_f2_audit.json. Overwrites no earlier product.
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
SCRIPT = os.path.abspath(__file__)

TALLY = "claim_verdict_tally_2026-09-17.csv"
ROUND_TABLE = {
    "phase2": "phase2b_claim_retest_combined.csv",   # the combined table supersedes the A-tier one
    "phase2b": "phase2b_claim_retest_combined.csv",
    "phase2d": "phase2d_claim_retest.csv",
    "phase2e": "phase2e_claim_retest.csv",
}
A_TIER_TABLE = "phase2_claim_retest.csv"             # cross-checked, not plotted

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 200,
    "axes.unicode_minus": False,
})
GREY, INK, ACCENT = "#9a9a9a", "#222222", "#1f5fa9"
VERDICT_ORDER = ["survives", "attenuated", "null_broken_by_control", "vanishes",
                 "baseline_contradicts_claim", "undecidable"]
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
    "vanishes": "vanishes (transfer test)",
    "baseline_contradicts_claim": "baseline contradicts claim",
    "undecidable": "undecidable (power-limited)",
}
NO_BASELINE_TEXT = {
    "out_of_instrument_scope": "out of instrument scope",
    "blocked_material_unreachable": "material unreachable",
}
# PERS-010's registered estimand is not an odds ratio; see the docstring.
OWN_ESTIMAND = {"PERS-010": "log2(observed shared / expected shared)"}


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(name):
    with open(os.path.join(RESULTS, name), encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def cell(row, key):
    """Return the raw string of a cell, or None when the column is absent/empty."""
    value = row.get(key)
    if value is None or value.strip() == "":
        return None
    return value.strip()


def collect():
    """Join the tally (row set) to the round tables (effect sizes). Reads only; no arithmetic."""
    tally = read(TALLY)
    tables = {name: read(name) for name in sorted(set(ROUND_TABLE.values()))}
    rows, no_baseline, provenance = [], [], []
    for entry in tally:
        claim = entry["claim_id"]
        rnd = entry["round"]
        table_name = ROUND_TABLE[rnd]
        source = None
        for candidate in tables[table_name]:
            if candidate["claim_id"] != claim:
                continue
            if "caliber" in candidate and candidate["caliber"] != "primary":
                continue  # 2d/2e emit a sensitivity row as well; the tally counts primary only
            source = candidate
            break
        if source is None:
            raise SystemExit("no source row for %s in %s" % (claim, table_name))
        record = {
            "claim_id": claim,
            "round": rnd,
            "verdict": entry["verdict"],
            "has_measured_baseline": entry["has_measured_baseline"] == "1",
            "sensitivity_caliber": cell(entry, "sensitivity_caliber"),
            "sensitivity_verdict": cell(entry, "sensitivity_verdict"),
            "earlier_round": cell(entry, "earlier_round"),
            "earlier_verdict": cell(entry, "earlier_verdict"),
            "source_table": table_name,
            # the tally credits the A-tier claims to the combined table's round; the round that
            # actually produced the row is the source table's own `round` column where it exists
            "round_of_row": cell(source, "round") or rnd,
            "unit": cell(source, "unit"),
            "claim_direction": cell(source, "claim_direction"),
            "is_transfer": cell(source, "is_transfer") == "True",
            "transfer_cohort": cell(source, "transfer_cohort"),
            "baseline_reproduction": (cell(source, "baseline_reproduction")
                                      or cell(source, "baseline_reproduced")),
            "attribute_label": cell(source, "attribute_label"),
            "n_observations": cell(source, "n_observations"),
            "n_positive": cell(source, "n_positive"),
        }
        for key in ("baseline_log2_or", "baseline_ci_low", "baseline_ci_high",
                    "matched_log2_or", "matched_ci_low", "matched_ci_high",
                    "random_control_log2_or", "random_control_ci_low", "random_control_ci_high"):
            raw = cell(source, key)
            record[key] = None if raw is None else float(raw)
            record[key + "__raw"] = raw
        if record["has_measured_baseline"]:
            rows.append(record)
        else:
            no_baseline.append(record)
        provenance.append({"claim_id": claim, "round_in_tally": rnd,
                           "round_of_row": record["round_of_row"],
                           "table": "results/" + table_name, "verdict": entry["verdict"],
                           "has_measured_baseline": record["has_measured_baseline"]})
    return rows, no_baseline, provenance


def sort_key(record):
    return (VERDICT_ORDER.index(record["verdict"]), record["claim_id"])


def flags(record):
    marks = ""
    if record["is_transfer"]:
        marks += "T"
    if record["sensitivity_verdict"] and record["sensitivity_verdict"] != record["verdict"]:
        marks += "S"
    if record["earlier_verdict"] and record["earlier_verdict"] != record["verdict"]:
        marks += "R"
    if record["claim_direction"] == "null_no_preference":
        marks += "n"
    return marks


def label_of(record):
    mark = flags(record)
    tag = record["round_of_row"].replace("phase2", "2").upper()
    return "%s  %s%s" % (record["claim_id"], tag, ("  " + mark) if mark else "")


def draw_rows(ax, records, title, drawn):
    """One row per claim: baseline (open), matched (filled), random control (cross)."""
    ticks, labels = [], []
    y = 0.0
    for verdict in VERDICT_ORDER:
        block = [r for r in records if r["verdict"] == verdict]
        if not block:
            continue
        colour = VERDICT_COLOUR[verdict]
        ax.text(0.0, y, VERDICT_TEXT[verdict], transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=6.6, color=colour, style="italic",
                clip_on=False)
        y -= 1.0
        for record in sorted(block, key=sort_key):
            base, matched = record["baseline_log2_or"], record["matched_log2_or"]
            if base is not None and matched is not None:
                ax.annotate("", xy=(matched, y), xytext=(base, y),
                            arrowprops=dict(arrowstyle="-|>", color=colour, lw=0.8,
                                            alpha=0.45, shrinkA=3.0, shrinkB=3.0))
            if record["baseline_ci_low"] is not None:
                ax.plot([record["baseline_ci_low"], record["baseline_ci_high"]],
                        [y + 0.22] * 2, color=INK, lw=0.9, solid_capstyle="butt", zorder=2)
            if base is not None:
                ax.plot([base], [y + 0.22], marker="o", mfc="white", mec=INK, mew=0.9,
                        ms=4.2, zorder=3)
            if record["matched_ci_low"] is not None:
                ax.plot([record["matched_ci_low"], record["matched_ci_high"]],
                        [y] * 2, color=colour, lw=1.5, solid_capstyle="butt", zorder=2)
            if matched is not None:
                ax.plot([matched], [y], marker="o", color=colour, ms=4.6, zorder=4)
            if record["random_control_ci_low"] is not None:
                ax.plot([record["random_control_ci_low"], record["random_control_ci_high"]],
                        [y - 0.22] * 2, color=GREY, lw=0.8, solid_capstyle="butt", zorder=2)
            if record["random_control_log2_or"] is not None:
                ax.plot([record["random_control_log2_or"]], [y - 0.22], marker="x",
                        color="#5b6068", ms=4.0, mew=0.9, zorder=5)
            for key in ("baseline_log2_or", "baseline_ci_low", "baseline_ci_high",
                        "matched_log2_or", "matched_ci_low", "matched_ci_high",
                        "random_control_log2_or", "random_control_ci_low",
                        "random_control_ci_high"):
                if record[key] is not None:
                    drawn.append({"claim_id": record["claim_id"], "field": key,
                                  "value": record[key], "raw": record[key + "__raw"],
                                  "table": "results/" + record["source_table"]})
            ticks.append(y)
            labels.append(label_of(record))
            y -= 1.0
        y -= 0.25
    ax.axvline(0.0, color=INK, lw=0.8, zorder=1)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_ylim(y + 0.55, 0.75)
    ax.set_title(title, loc="left")
    ax.tick_params(axis="y", length=0)
    return ticks


def figure_f2(rows, no_baseline):
    drawn = []
    site = [r for r in rows if r["unit"] == "site"]
    protein = [r for r in rows if r["unit"] == "protein" and r["claim_id"] not in OWN_ESTIMAND]
    own = [r for r in rows if r["claim_id"] in OWN_ESTIMAND]

    n_site_rows = len(site) + len({r["verdict"] for r in site})
    n_prot_rows = len(protein) + len({r["verdict"] for r in protein})
    height = 0.255 * n_site_rows + 1.5
    fig = plt.figure(figsize=(7.6, height))
    grid = fig.add_gridspec(3, 2, width_ratios=[1.0, 0.92],
                            height_ratios=[n_prot_rows + 1, 3.0,
                                           max(n_site_rows - n_prot_rows - 4, 5)],
                            wspace=0.40, hspace=1.10)
    ax_site = fig.add_subplot(grid[:, 0])
    ax_prot = fig.add_subplot(grid[0, 1])
    ax_own = fig.add_subplot(grid[1, 1])
    ax_text = fig.add_subplot(grid[2, 1])
    ax_text.axis("off")

    draw_rows(ax_site, site, "a  site-level claims (%d)" % len(site), drawn)
    ax_site.set_xlabel("effect size, log2 odds ratio\n(claimed attribute x author's positive label)")

    draw_rows(ax_prot, protein, "b  protein-level claims (%d)" % len(protein), drawn)
    ax_prot.set_xlabel("effect size, log2 odds ratio")

    # (c) PERS-010: a different registered estimand, hence its own axis.
    draw_rows(ax_own, own, "c  PERS-010, separate estimand", drawn)
    ax_own.set_xlabel("log2(observed shared / expected shared)\nunder detection independence")

    # (d) the claims with no measured baseline; text only, never points.
    lines = []
    for record in sorted(no_baseline, key=lambda r: (r["verdict"], r["claim_id"])):
        lines.append("%s   %s" % (record["claim_id"], NO_BASELINE_TEXT[record["verdict"]]))
    ax_text.text(0.0, 1.0, "d  no measured baseline (%d claims, not part of the %d above)"
                 % (len(no_baseline), len(rows)), transform=ax_text.transAxes,
                 ha="left", va="top", fontsize=9, color=INK)
    ax_text.text(0.0, 0.86, "\n".join(lines), transform=ax_text.transAxes, ha="left",
                 va="top", fontsize=7, color="#5b6068", linespacing=1.5)

    handles = [plt.Line2D([], [], marker="o", color=VERDICT_COLOUR[v], ls="none", ms=4.6,
                          label=VERDICT_TEXT[v])
               for v in VERDICT_ORDER if any(r["verdict"] == v for r in rows)]
    handles += [
        plt.Line2D([], [], marker="o", mfc="white", mec=INK, mew=0.9, ls="none", ms=4.2,
                   label="baseline, the paper's own caliber"),
        plt.Line2D([], [], marker="o", color=INK, ls="none", ms=4.6,
                   label="after detectability matching (primary caliber)"),
        plt.Line2D([], [], marker="x", color="#5b6068", ls="none", ms=4.0, mew=0.9,
                   label="same-size random control (power reference)"),
        plt.Line2D([], [], color=INK, lw=1.0,
                   label="95% interval, cluster bootstrap 5000, seed 20260915"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.012),
               ncol=3, frameon=False, handletextpad=0.6, columnspacing=1.6)
    ax_text.text(0.0, 0.16, "row tags\n2 / 2B / 2D / 2E  the round that produced the row\n"
                 "T  transfer test, the baseline is another cohort\n"
                 "S  the sensitivity caliber gives a different verdict\n"
                 "R  the verdict changed between rounds\n"
                 "n  the claim asserts no preference",
                 transform=ax_text.transAxes, ha="left", va="top", fontsize=7,
                 color="#5b6068", linespacing=1.5)
    for axis in (ax_site, ax_prot, ax_own):
        axis.margins(x=0.06)
    fig.subplots_adjust(left=0.150, right=0.985, top=0.965, bottom=0.075)
    for ext, kwargs in (("png", {"dpi": 600}), ("pdf", {})):
        fig.savefig(os.path.join(RESULTS, "manuscript_f2_claim_retest." + ext),
                    bbox_inches="tight", pad_inches=0.12, **kwargs)
    plt.close(fig)
    return drawn


def recheck(drawn, rows, no_baseline):
    """Find every drawn number back in the source tables, cell by cell."""
    cache = {name: read(name) for name in sorted(set(ROUND_TABLE.values()))}
    misses, checked = [], 0
    for item in drawn:
        table = cache[os.path.basename(item["table"])]
        found = False
        for row in table:
            if row["claim_id"] != item["claim_id"]:
                continue
            if "caliber" in row and row["caliber"] != "primary":
                continue
            raw = cell(row, item["field"])
            if raw is not None and float(raw) == item["value"] and raw == item["raw"]:
                found = True
                break
        checked += 1
        if not found:
            misses.append(item)
    # counts printed on the figure must match the tally audit's own counts
    audit = json.load(open(os.path.join(RESULTS, "claim_verdict_tally_2026-09-17_audit.json"),
                           encoding="utf-8"))
    count_checks = {
        "n_with_measured_baseline_figure": len(rows),
        "n_with_measured_baseline_tally_audit": audit["n_with_measured_baseline"],
        "n_distinct_claims_figure": len(rows) + len(no_baseline),
        "n_distinct_claims_tally_audit": audit["n_distinct_claims"],
        "agree": (len(rows) == audit["n_with_measured_baseline"]
                  and len(rows) + len(no_baseline) == audit["n_distinct_claims"]),
    }
    # the A-tier table is superseded by the combined one; confirm they agree where both exist
    a_tier = {r["claim_id"]: r for r in read(A_TIER_TABLE)}
    combined = {r["claim_id"]: r for r in read("phase2b_claim_retest_combined.csv")}
    a_tier_diffs = []
    for claim, row in a_tier.items():
        other = combined.get(claim)
        if other is None:
            a_tier_diffs.append({"claim_id": claim, "issue": "absent from the combined table"})
            continue
        for key in ("verdict", "baseline_log2_or", "matched_log2_or", "random_control_log2_or"):
            left, right = cell(row, key), cell(other, key)
            if left != right:
                a_tier_diffs.append({"claim_id": claim, "field": key,
                                     "phase2_claim_retest.csv": left,
                                     "combined": right})
    return {"values_drawn": checked, "values_not_found_in_sources": misses,
            "counts": count_checks, "a_tier_vs_combined": a_tier_diffs}


def main():
    started = time.time()
    rows, no_baseline, provenance = collect()
    drawn = figure_f2(rows, no_baseline)
    verification = recheck(drawn, rows, no_baseline)
    inputs = [TALLY, A_TIER_TABLE] + sorted(set(ROUND_TABLE.values()))
    audit = {
        "script": "scripts/plot_manuscript_f2_claim_retest.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "round": "main figure F2, drawn from stored tables only",
        "computes_nothing": True,
        "row_set_source": "results/" + TALLY,
        "panels": {
            "a": {"content": "site-level claims, log2 odds ratio, grouped by verdict",
                  "unit": "site",
                  "claims": [r["claim_id"] for r in rows if r["unit"] == "site"]},
            "b": {"content": "protein-level claims, log2 odds ratio, grouped by verdict",
                  "unit": "protein",
                  "claims": [r["claim_id"] for r in rows
                             if r["unit"] == "protein" and r["claim_id"] not in OWN_ESTIMAND]},
            "c": {"content": "PERS-010 only; registered estimand is "
                             "log2(observed shared / expected shared), not an odds ratio, "
                             "so it may not share panel b's axis",
                  "claims": sorted(OWN_ESTIMAND)},
            "d": {"content": "claims with no measured baseline, listed as text",
                  "claims": [r["claim_id"] for r in no_baseline]},
        },
        "estimates_drawn_per_row": ["baseline (the paper's own caliber)",
                                    "detectability matching (registered primary caliber)",
                                    "same-size random control (power reference)"],
        "estimates_deliberately_not_drawn": {
            "stratified_log2_or": "registered SECONDARY caliber; four estimates per row is "
                                  "unreadable and the verdict rule keys on matching",
            "retention_ratio": "|controlled|/|baseline|, sign-blind and unstable near zero "
                               "(1.42/2.16/3.34 on disk); rejected as a common axis",
            "*_strata.csv": "per-stratum detail behind the secondary caliber; not used",
        },
        "claim_provenance": provenance,
        "baseline_reproduction_per_claim": {r["claim_id"]: r["baseline_reproduction"]
                                            for r in rows},
        "flags": {"T": "is_transfer", "S": "tally sensitivity verdict differs from primary",
                  "R": "tally earlier-round verdict differs", "n": "claim_direction is "
                       "null_no_preference"},
        "verification": verification,
        "outputs": ["results/manuscript_f2_claim_retest.png",
                    "results/manuscript_f2_claim_retest.pdf"],
        "inputs": {"results/" + name: sha256_of(os.path.join(RESULTS, name)) for name in inputs},
        "versions": {"python": sys.version, "matplotlib": matplotlib.__version__},
        "elapsed_seconds": round(time.time() - started, 2),
    }
    out = os.path.join(RESULTS, "manuscript_f2_audit.json")
    with open(out + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    os.replace(out + ".tmp", out)
    print("wrote F2:", len(rows), "claims with a baseline,", len(no_baseline), "without;",
          verification["values_drawn"], "values re-checked,",
          len(verification["values_not_found_in_sources"]), "not found;",
          "counts agree:", verification["counts"]["agree"], ";",
          "A-tier vs combined differences:", len(verification["a_tier_vs_combined"]))


if __name__ == "__main__":
    main()
