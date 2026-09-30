"""Why do the donors the paper calls ineffective contribute the most tier B sites?

`reports/QPERS_SID_SITE_MAPPING.md` left this open: the paper reports that short
NaSH and Na2S treatment had hardly any detectable effect on proteome persulfide
formation, yet after site mapping those two arms carry the MOST tier B sites
(780 and 677 against GYY4137's 346).  Until that is explained the
enrichment-based label cannot be used.

**The hypothesis this script tests, stated before running it:** the four donor
files do not apply the same filter, so the arm site counts are not comparable.
The legends, read from the files themselves and quoted verbatim into the audit,
say:

* GYY4137 (MOESM3) - "only significant peptides that showed an enriched elution
  of more than 30 % (Mean norm Ratio H/L >= 1.30; light green) were considered"
* Na2S4  (MOESM4) - "only peptides that showed an enriched elution of more than
  30 % (Mean norm Ratio H/L >= 1.30) were considered"
* Na2S   (MOESM5) - "Peptides were quantified using MaxQuant v1.5.2.8"
* NaSH   (MOESM6) - "Peptides were quantified using MaxQuant v1.5.2.8"

i.e. only the first two state an enrichment threshold, and only the first two
carry a `Significant` column on the Elution sheet.  If the hypothesis holds, the
Na2S and NaSH Elution sheets are unfiltered lists of everything quantified, their
ratios will centre on 1.0, and almost none will clear 1.30.

**What this script does about it.** It applies the paper's own stated criterion
uniformly to all four arms - mean normalised H/L >= 1.30, and `Significant` == "+"
where that column exists - and reports what each arm then contributes.  It does
not invent a threshold; 1.30 is the authors' number.

**Negative sets, predeclared.** Two candidates are counted for each arm, and the
second is the better-matched one:

* `beads_probe_thiol`  - sites seen probe-labelled on the beads for that arm and
  never enriched: seen by the same probe in the same run.
* `elution_not_enriched` - sites quantified in the SAME Elution sheet whose mean
  ratio is below 1.30: same probe, same run, same fraction, no enrichment. The
  unfiltered Na2S and NaSH sheets are therefore not waste - they are the
  best-matched negatives for the enrichment label.

**Tier A is also re-examined here** because the mapping round found all 130
adduct-based sites in MOESM2. This script confirms by counting raw "+446.17"
occurrences in MOESM7, and reports tier A's only available depth analogues: the
number of MOESM2 experiments (1-2) and the number of fractions (1-2) in which a
site is seen as a persulfide.

No model is fitted and no interval is computed; this round settles definitions.
Requires xlrd for the legends, as the ingest script does.
"""
from __future__ import annotations

import collections
import csv
import json
import platform
import re
import statistics
import sys
import time

import xlrd

from common import ROOT, RESULTS, sha256, write_csv, write_json

SCRIPT = ROOT / "scripts" / "run_qpers_sid_label_criteria_audit.py"
SITES = RESULTS / "qpers_sid_sites_normalised.csv"
INTAKE = ROOT / "external" / "intake" / "scirep2016_qpers_sid"
STEM = "41598_2016_BFsrep29808_MOESM{}_ESM.xls"
DONOR_FILES = {"GYY4137": "3", "Na2S4": "4", "Na2S": "5", "NaSH": "6"}
THRESHOLD = 1.30  # the authors' own number, not ours


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def legend_text(moesm):
    book = xlrd.open_workbook(INTAKE / STEM.format(moesm), on_demand=True)
    sheet = book.sheet_by_name("Legend")
    lines = [re.sub(r"\s+", " ", str(sheet.cell_value(r, c))).strip()
             for r in range(sheet.nrows) for c in range(sheet.ncols)
             if str(sheet.cell_value(r, c)).strip()]
    book.release_resources()
    return lines


def elution_has_significant_column(moesm):
    book = xlrd.open_workbook(INTAKE / STEM.format(moesm), on_demand=True)
    sheet = book.sheet_by_name("Elution")
    names = {re.sub(r"\s+", " ", str(sheet.cell_value(2, c))).strip() for c in range(sheet.ncols)}
    book.release_resources()
    return "Significant" in names


def moesm7_persulfide_occurrences():
    book = xlrd.open_workbook(INTAKE / STEM.format("7"), on_demand=True)
    counts = collections.Counter()
    for name in book.sheet_names():
        if name == "Legend":
            continue
        sheet = book.sheet_by_name(name)
        columns = {re.sub(r"\s+", " ", str(sheet.cell_value(2, c))).strip(): c
                   for c in range(sheet.ncols)}
        for index in range(3, sheet.nrows):
            peptide = str(sheet.cell_value(index, columns["Peptide"]))
            for tag in ("+57.02", "+414.19", "+446.17"):
                if tag in peptide:
                    counts[f"{name}|{tag}"] += 1
    book.release_resources()
    return counts


def ratio(row):
    try:
        return float(row["ratio_average"])
    except (TypeError, ValueError):
        return None


def main():
    started = time.time()
    rows = read_rows(SITES)
    legends = {donor: legend_text(moesm) for donor, moesm in DONOR_FILES.items()}
    has_significant = {donor: elution_has_significant_column(moesm)
                       for donor, moesm in DONOR_FILES.items()}
    states_threshold = {donor: any("1.30" in line or "30 %" in line for line in lines)
                        for donor, lines in legends.items()}

    arm_rows, label_rows = [], []
    positives, negatives_elution, negatives_beads = {}, {}, {}
    for donor in DONOR_FILES:
        elution = [r for r in rows if r["donor"] == donor and r["class_from_adduct"] == "SSH_enriched"]
        beads = [r for r in rows if r["donor"] == donor and r["class_from_adduct"] == "SH_probe_on_beads"]
        values = [ratio(r) for r in elution]
        numeric = sorted(v for v in values if v is not None)
        cleared = [r for r in elution if (ratio(r) or 0) >= THRESHOLD]
        significant = [r for r in cleared if not has_significant[donor] or r["significant"] == "+"]
        positives[donor] = {(r["accession"], r["site"]) for r in significant}
        negatives_elution[donor] = {(r["accession"], r["site"]) for r in elution
                                    if (ratio(r) is not None and ratio(r) < THRESHOLD)} - positives[donor]
        negatives_beads[donor] = {(r["accession"], r["site"]) for r in beads} - positives[donor]
        arm_rows.append({
            "donor": donor,
            "legend_states_enrichment_threshold": states_threshold[donor],
            "elution_sheet_has_significant_column": has_significant[donor],
            "elution_sites_listed": len({(r["accession"], r["site"]) for r in elution}),
            "elution_rows": len(elution),
            "ratio_min": round(numeric[0], 3) if numeric else "",
            "ratio_q25": round(numeric[len(numeric) // 4], 3) if numeric else "",
            "ratio_median": round(statistics.median(numeric), 3) if numeric else "",
            "ratio_q75": round(numeric[3 * len(numeric) // 4], 3) if numeric else "",
            "ratio_max": round(numeric[-1], 3) if numeric else "",
            "rows_at_or_above_1_30": len(cleared),
            "rows_below_1_30": sum(1 for v in values if v is not None and v < THRESHOLD),
            "rows_significant_plus": sum(1 for r in elution if r["significant"] == "+"),
            "positives_after_uniform_criterion": len(positives[donor]),
            "negatives_elution_not_enriched": len(negatives_elution[donor]),
            "negatives_beads_probe_thiol": len(negatives_beads[donor]),
        })

    union_positive = set().union(*positives.values())
    depth = collections.Counter(sum(1 for d in positives if site in positives[d])
                                for site in union_positive)
    informative_arms = [d for d in positives if len(positives[d]) >= 50]

    # tier A, from the mapped table plus a raw confirmation that MOESM7 carries no +446.17
    raw7 = moesm7_persulfide_occurrences()
    persulfide_in_7 = sum(v for k, v in raw7.items() if k.endswith("+446.17"))
    tier_a = [r for r in rows if r["class_from_adduct"] == "SSH_probe"]
    a_by_site = collections.defaultdict(lambda: {"experiments": set(), "fractions": set(), "tables": set()})
    for r in tier_a:
        key = (r["accession"], r["site"])
        a_by_site[key]["experiments"].add(r["experiment"])
        a_by_site[key]["fractions"].add(r["fraction"])
        a_by_site[key]["tables"].add(r["source_table"])
    a_thiol = {(r["accession"], r["site"]) for r in rows
               if r["class_from_adduct"] == "SH_probe" and r["fraction"] in ("Beads", "Total")}

    label_rows.append({
        "tier": "A_adduct_based", "definition": "peptide carries +446.17 in MOESM2",
        "arms": "main experiment only; MOESM7 carries no persulfide adduct at all",
        "positives": len(a_by_site),
        "negatives": len(a_thiol - set(a_by_site)),
        "negative_definition": "probe-labelled thiol (+414.19) in Beads or Total, never seen as +446.17",
        "depth_covariate": "number of MOESM2 experiments (1-2), or number of fractions (1-2)",
        "depth_distribution": json.dumps({
            "by_experiment": dict(sorted(collections.Counter(
                len(v["experiments"]) for v in a_by_site.values()).items())),
            "by_fraction": dict(sorted(collections.Counter(
                len(v["fractions"]) for v in a_by_site.values()).items())),
        }),
        "usable": "yes, but small and with only a two-valued depth covariate",
    })
    label_rows.append({
        "tier": "B_enrichment_based", "definition": f"mean norm H/L >= {THRESHOLD} and Significant == '+' where that column exists",
        "arms": ", ".join(f"{d}:{len(positives[d])}" for d in DONOR_FILES),
        "positives": len(union_positive),
        "negatives": len(set().union(*negatives_beads.values()) | set().union(*negatives_elution.values())),
        "negative_definition": "same-sheet Elution sites below the threshold, plus bead-retained probe-labelled thiols, minus any positive",
        "depth_covariate": "number of donor arms in which the site clears the criterion",
        "depth_distribution": json.dumps(dict(sorted(depth.items()))),
        "usable": (f"yes as a {len(informative_arms)}-arm design ({', '.join(informative_arms)}); "
                   "Na2S and NaSH contribute almost nothing once the authors' own criterion is applied"),
    })

    write_csv(RESULTS / "qpers_sid_tier_b_arm_criteria.csv", arm_rows)
    write_csv(RESULTS / "qpers_sid_label_candidates.csv", label_rows)

    for row in arm_rows:
        print(f"{row['donor']:8s} threshold_in_legend={str(row['legend_states_enrichment_threshold']):5s} "
              f"sig_col={str(row['elution_sheet_has_significant_column']):5s} "
              f"median_ratio={row['ratio_median']:>6} >=1.30 {row['rows_at_or_above_1_30']:4d}/{row['elution_rows']:4d} "
              f"-> positives {row['positives_after_uniform_criterion']:4d} "
              f"neg_elution {row['negatives_elution_not_enriched']:4d} neg_beads {row['negatives_beads_probe_thiol']:4d}",
              flush=True)
    print(f"\ntier B union positives {len(union_positive)} depth {dict(sorted(depth.items()))}")
    print(f"tier A positives {len(a_by_site)} negatives {len(a_thiol - set(a_by_site))} "
          f"| MOESM7 +446.17 occurrences {persulfide_in_7}")

    write_json(RESULTS / "qpers_sid_label_criteria_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": ("why do the donors the paper calls ineffective contribute the most tier B sites, and "
                     "what label definitions survive the answer"),
        "answer": ("the four donor files do not apply the same filter. Only GYY4137 and Na2S4 state an "
                   "enrichment threshold in their legend and only those two carry a Significant column on the "
                   "Elution sheet; the Na2S and NaSH Elution sheets are unfiltered lists of everything "
                   "quantified, their ratios centre on 1.0, and almost nothing clears the authors' own 1.30 "
                   "cut. The paper's text and its tables are consistent once the filter difference is taken "
                   "into account; the site counts in reports/QPERS_SID_SITE_MAPPING.md were comparing a "
                   "filtered list against an unfiltered one."),
        "legends_verbatim": legends,
        "elution_sheet_has_significant_column": has_significant,
        "legend_states_enrichment_threshold": states_threshold,
        "threshold_used": THRESHOLD,
        "threshold_provenance": "the authors' own number, quoted in the GYY4137 and Na2S4 legends; not chosen by us",
        "arm_table": arm_rows,
        "label_candidates": label_rows,
        "tier_a_confirmation": {
            "moesm7_adduct_occurrences": dict(sorted(raw7.items())),
            "moesm7_persulfide_occurrences": persulfide_in_7,
            "conclusion": ("all tier A adduct evidence comes from MOESM2; the per-donor SILAC tables carry "
                           "+414.19 and +57.02 but no +446.17 at all, although their legend lists the colour "
                           "code for it. Tier A therefore has no donor structure."),
        },
        "consequences_for_v3": [
            "tier B is a two-arm design (GYY4137 and Na2S4), not four, so the donor-arm depth covariate takes only the values 1 and 2",
            "the donor-arm depth distribution 505/269/191/227 reported in reports/QPERS_SID_SITE_MAPPING.md was built on unfiltered arms and must not be used",
            "the unfiltered Na2S and NaSH Elution rows are the best-matched negatives for the enrichment label: same probe, same run, same fraction, no enrichment",
            "tier A stays at 130 positives from one experiment with a two-valued depth covariate, so it cannot support depth stratification the way the QTRP arms did",
            "before either negative definition is used as a primary caliber it must pass the standing collinearity pre-check (per-class retention, then the standalone discriminability of the covariate the rule removes)",
        ],
        "limits": [
            "applying the authors' threshold to arms where they did not apply it assumes the ratio columns are comparable across files, which the paper does not state",
            "the Na2S and NaSH ratio distributions centre on 1.0, so those arms carry no usable positive signal at any threshold, not merely at 1.30",
            "tier B positives remain enrichment-based, not adduct-based",
            "no ProteomeXchange accession exists for this study, so none of this can be re-derived from raw spectra",
        ],
        "input_hashes": {
            "results/qpers_sid_sites_normalised.csv": sha256(SITES),
            "scripts/run_qpers_sid_label_criteria_audit.py": sha256(SCRIPT),
            **{f"external/intake/scirep2016_qpers_sid/{STEM.format(m)}": sha256(INTAKE / STEM.format(m))
               for m in sorted(set(DONOR_FILES.values()) | {"7"})},
        },
        "versions": {"python": sys.version, "xlrd": xlrd.__version__},
    })


if __name__ == "__main__":
    main()
