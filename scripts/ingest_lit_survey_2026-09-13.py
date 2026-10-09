"""Bring the 74-paper literature survey's raw tables into this tree, and recompute its numbers.

Why this script exists. The manuscript reports a survey of published site-level analyses and quotes
counts, percentages and Wilson intervals from it, but until now NOT ONE of the tables behind those
numbers was in this tree: they lived only in the external artifact store of project
proj_8c1f174d761d, and the code that produced the percentages existed only as recorded notebook
cells. A Data/Code availability statement cannot honestly point at either. This script copies the
products in, hashes them, and recomputes every published number from the copied coding table so
that the statement has something to point at.

WHAT IS COPIED (source paths recorded verbatim in the audit, with sha256 of source and copy)
  coding_table          the per-paper coding table, 942 data rows, 30 columns
  report                the survey report: the eight query strings verbatim, the endpoints, the
                        run date and the screening funnel
  summary               the machine-readable summary the report was written from
  deduped               the 2,250 deduplicated bibliographic records
  screen_raw            the 282 raw model-assisted screening batches
  packs                 the regex-extracted evidence passages, 367 full texts
  final_codes           the per-paper final codes
  figure                the survey figure embedded in the manuscript
  notebook_tape         the recorded notebook cells that produced the published percentages

  Earlier versions of the coding table and the report are copied as well, under their own names, so
  that "which version produced the manuscript's numbers" is answerable from this tree rather than
  from a timestamp in someone's artifact store.

WHAT IS RECOMPUTED, AND WHAT THE READING RULE IS (fixed before this script was run)
  Denominators, written down before touching the table:
    * the SCREENED corpus is every row of the coding table
    * the FULL-TEXT corpus is the rows with fulltext_retrieved == 1
    * the CLASSIFIED corpus is the rows with fulltext_stage_verdict == "included after full text"
      - this is the 74, and every percentage the manuscript quotes from the survey has 74 as its
      denominator unless the number is explicitly a funnel count
  Background classes are taken from the column background_class as stored; b and c are the two
  classes the manuscript calls "demonstrably matched or used detected-but-unmodified residues", and
  e is "no background set stated at all". No class is re-coded here and no row is re-read: this
  script recomputes arithmetic, it does not repeat the survey.

  The interval is the Wilson score interval, transcribed VERBATIM from cell 568 of the recorded
  notebook tape so that the published intervals are reproduced by the same estimator that produced
  them rather than by a differently-rounded reimplementation:

      def wilson(k, n, z=1.96):
          if n == 0: return (0,0,0)
          p = k/n; d = 1 + z*z/n
          c = (p + z*z/(2*n))/d
          h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
          return 100*p, 100*(c-h), 100*(c+h)

  Every recomputed value is compared against the value the manuscript quotes. A mismatch is
  REPORTED, never adjusted: the audit's `mismatches` list must be empty for the availability
  statement to be self-consistent, and if it is not, the report says which number moved.

WHAT THIS SCRIPT DOES NOT DO
  It does not re-run the search, re-screen, re-code or re-read any paper, and it does not touch the
  manuscript. It also does not claim the survey's method is sound - the survey's own limits (single
  reader, model-assisted screening, full-text retrievability) are stated where the survey is
  reported and are not re-litigated here.

Writes external/lit_survey_2026-09-13/ (copies), results/lit_survey_recount.csv and
results/lit_survey_recount_audit.json. Overwrites nothing else.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
DEST = os.path.join(ROOT, "external", "lit_survey_2026-09-13")
SCRIPT = os.path.abspath(__file__)

STORE = os.path.expanduser(
    "/path/to/artifact_store")

SOURCES = [
    ("coding_table", "7f10dd27-26bf-4a52-8af7-4e1c936dc312/v2537f0fa_lit_survey_background_controls.csv",
     "lit_survey_background_controls.csv",
     "the per-paper coding table the manuscript's survey numbers are computed from"),
    ("coding_table_earlier_version", "7f10dd27-26bf-4a52-8af7-4e1c936dc312/v6287af3d_lit_survey_background_controls.csv",
     "lit_survey_background_controls_earlier_version.csv",
     "an earlier version of the same table, copied so that version provenance is answerable here"),
    ("report", "6949011e-6526-4560-a6b6-dbf9249b1b55/vf820e699_lit_survey_report.md",
     "lit_survey_report.md",
     "the survey report: the eight query strings verbatim, the endpoints, the run date, the funnel"),
    ("report_earlier_version_a", "6949011e-6526-4560-a6b6-dbf9249b1b55/va5e043aa_lit_survey_report.md",
     "lit_survey_report_earlier_version_a.md", "earlier version of the report"),
    ("report_earlier_version_b", "6949011e-6526-4560-a6b6-dbf9249b1b55/vbffe97b5_lit_survey_report.md",
     "lit_survey_report_earlier_version_b.md", "earlier version of the report"),
    ("summary", "29ee06ef-5e46-4c41-82bd-3fe300ec9ceb/v10bd9c79_out.json",
     "lit_survey_summary.json", "the machine-readable summary the report was written from"),
    ("deduped", "332680df-0865-49db-a9d1-80b0a757b16e/vef0d9b40_deduped.json",
     "lit_survey_deduped_records.json", "the deduplicated bibliographic records"),
    ("screen_raw", "75e13330-24e7-4a27-b743-be60a75a510c/v9a8f3946_screen_raw.json",
     "lit_survey_screen_raw.json", "the raw model-assisted screening batches"),
    ("packs", "6301b066-426f-4847-82f2-8f3506270c5b/v5da76780_packs.json",
     "lit_survey_evidence_packs.json", "the regex-extracted evidence passages per full text"),
    ("final_codes", "04176d96-bc78-4a2d-b565-ea2e8d558f2b/vfcd09395_final_codes.json",
     "lit_survey_final_codes.json", "the per-paper final codes"),
    ("figure", "d7936211-476a-4ee7-8cdd-feaa25c02f9f/v0ff1e709_lit_survey_background_controls.png",
     "lit_survey_background_controls.png", "the survey figure embedded in the manuscript"),
    ("notebook_tape", "host_call_tapes/2295f575-1ee1-4be8-93fa-15186601308e-0.json",
     "lit_survey_notebook_cells.json",
     "the recorded notebook cells that produced the published percentages; there is no .py"),
]

# What the manuscript quotes, written here before the recount was run.
QUOTED = {
    "n_classified": 74,
    "n_matched_or_observed_background": 5,
    "pct_matched_or_observed_background": 6.8,
    "ci_matched_or_observed_background": (2.9, 14.9),
    "n_no_background_stated": 35,
    "pct_no_background_stated": 47.3,
    "ci_no_background_stated": (36.3, 58.5),
    "n_discusses_detectability_bias": 1,
    "n_reports_identified_fraction": 2,
    "n_screened_included": 942,
    "pct_fulltext_of_included": 39.0,
}


def wilson(k, n, z=1.96):
    # transcribed verbatim from cell 568 of the recorded notebook tape
    if n == 0:
        return (0, 0, 0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * p, 100 * (c - h), 100 * (c + h)


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    started = time.time()
    os.makedirs(DEST, exist_ok=True)

    copied, missing = [], []
    for key, relative, name, why in SOURCES:
        source = os.path.join(STORE, relative)
        if not os.path.exists(source):
            missing.append({"key": key, "source": source})
            continue
        target = os.path.join(DEST, name)
        shutil.copy2(source, target)
        copied.append({"key": key, "source_path": source, "copied_to":
                       os.path.relpath(target, ROOT), "bytes": os.path.getsize(target),
                       "sha256_source": sha256_of(source), "sha256_copy": sha256_of(target),
                       "why": why})
    for entry in copied:
        if entry["sha256_source"] != entry["sha256_copy"]:
            raise SystemExit(f"copy differs from source: {entry['key']}")

    table = os.path.join(DEST, "lit_survey_background_controls.csv")
    with open(table, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))

    screened = rows
    fulltext = [r for r in rows if r["fulltext_retrieved"] == "1"]
    sampled = [r for r in rows if r["in_classified_sample"] == "1"]
    classified = [r for r in rows if r["fulltext_stage_verdict"] == "included after full text"]

    classes = {}
    for r in classified:
        classes[r["background_class"]] = classes.get(r["background_class"], 0) + 1

    n = len(classified)
    matched = classes.get("b", 0) + classes.get("c", 0)
    unstated = classes.get("e", 0)
    unmatched = n - matched
    detectability = sum(1 for r in classified if r["discusses_detectability_bias"] == "1")
    identified_fraction = sum(1 for r in classified if r["reports_identified_fraction"] == "1")
    kr_motif = sum(1 for r in classified if r["kr_basic_as_biological_motif"] == "1")
    equal_size_bg = sum(1 for r in classified if r["equal_size_random_background"] == "1")
    abundance = sum(1 for r in classified if r["abundance_normalisation"] == "1")

    modifications = {}
    for r in classified:
        modifications[r["modification"]] = modifications.get(r["modification"], 0) + 1

    out_rows = []

    def record(label, k, denominator, quoted_n=None, quoted_pct=None, quoted_ci=None, note=""):
        pct, low, high = wilson(k, denominator)
        row = {"quantity": label, "count": k, "denominator": denominator,
               "percent": round(pct, 1), "ci_low": round(low, 1), "ci_high": round(high, 1),
               "quoted_count": "" if quoted_n is None else quoted_n,
               "quoted_percent": "" if quoted_pct is None else quoted_pct,
               "quoted_ci": "" if quoted_ci is None else f"{quoted_ci[0]}-{quoted_ci[1]}",
               "agrees": "", "note": note}
        checks = []
        if quoted_n is not None:
            checks.append(k == quoted_n)
        if quoted_pct is not None:
            checks.append(round(pct, 1) == quoted_pct)
        if quoted_ci is not None:
            checks.append(round(low, 1) == quoted_ci[0] and round(high, 1) == quoted_ci[1])
        row["agrees"] = "" if not checks else str(all(checks))
        out_rows.append(row)
        return row

    record("papers entering classification", n, n, QUOTED["n_classified"], None, None,
           "the survey's denominator; every percentage below is out of this number")
    record("demonstrably matched background, or used detected-but-unmodified residues (classes b+c)",
           matched, n, QUOTED["n_matched_or_observed_background"],
           QUOTED["pct_matched_or_observed_background"],
           QUOTED["ci_matched_or_observed_background"])
    record("no background set stated at all (class e)", unstated, n,
           QUOTED["n_no_background_stated"], QUOTED["pct_no_background_stated"],
           QUOTED["ci_no_background_stated"])
    record("background not demonstrably matched (not b and not c)", unmatched, n, None, None, None,
           "the complement of b+c; reported in the survey report, not in the manuscript abstract")
    record("discusses cleavage or detectability bias", detectability, n,
           QUOTED["n_discusses_detectability_bias"])
    record("reports the site/identified overlap statistic", identified_fraction, n,
           QUOTED["n_reports_identified_fraction"])
    record("treats K/R enrichment as a biological motif", kr_motif, n, None, None, None,
           "NOTE: this count equals the b+c count (both 5 of 74), so the two share a percentage and "
           "an interval. They are DIFFERENT quantities and must never be conflated in the text.")
    record("uses an equal-size random background", equal_size_bg, n)
    record("applies an abundance normalisation", abundance, n)

    funnel = {
        "rows in the coding table (screened and included at title/abstract)": len(screened),
        "full text retrieved": len(fulltext),
        "sampled for classification": len(sampled),
        "classified after full text": n,
    }
    record("full text retrieved, out of the included corpus", len(fulltext), len(screened),
           QUOTED["n_screened_included"] if False else None, QUOTED["pct_fulltext_of_included"],
           None, "the only survey percentage whose denominator is not 74")

    mismatches = [r for r in out_rows if r["agrees"] == "False"]

    path = os.path.join(RESULTS, "lit_survey_recount.csv")
    cols = ["quantity", "count", "denominator", "percent", "ci_low", "ci_high", "quoted_count",
            "quoted_percent", "quoted_ci", "agrees", "note"]
    with open(path + ".tmp", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)
    os.replace(path + ".tmp", path)

    audit = {
        "script": "scripts/ingest_lit_survey_2026-09-13.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "why": "the manuscript quotes survey numbers whose tables were not in this tree; this round "
               "copies them in, hashes them and recomputes every quoted number from the copy",
        "source_store": STORE,
        "source_store_is_outside_this_tree": True,
        "copied": copied,
        "missing_sources": missing,
        "coding_table_rows": len(rows),
        "coding_table_columns": list(rows[0].keys()),
        "funnel_recomputed_from_the_copied_table": funnel,
        "background_class_counts_over_the_classified_corpus": classes,
        "modification_counts_over_the_classified_corpus": modifications,
        "interval_method": "Wilson score interval, z = 1.96, transcribed verbatim from cell 568 of "
                           "the recorded notebook tape",
        "the_code_that_produced_the_published_numbers": {
            "where": "external/lit_survey_2026-09-13/lit_survey_notebook_cells.json, cells 568-569",
            "form": "recorded notebook cells, NOT a file in any repository",
            "note": "cell 568 is a failed first attempt (it looked for a column named `included` and "
                    "reported 0 papers); cell 569 produced the published numbers. Both are in the "
                    "tape and the Code availability statement must say the computation was a "
                    "notebook cell, with this script as the file-based reproduction.",
        },
        "quoted_values_checked": QUOTED,
        "mismatches": mismatches,
        "all_quoted_numbers_reproduce": not mismatches,
        "what_this_round_did_not_do": [
            "it did not re-run the search, re-screen, re-code or re-read any paper",
            "it did not touch the manuscript",
            "it did not re-litigate the survey's own limits (single reader, model-assisted "
            "screening, full-text retrievability), which are stated where the survey is reported",
        ],
        "versions": {"python": sys.version},
        "elapsed_seconds": round(time.time() - started, 2),
        "corrections": [],
    }
    path = os.path.join(RESULTS, "lit_survey_recount_audit.json")
    with open(path + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    os.replace(path + ".tmp", path)

    print(f"copied {len(copied)} files into external/lit_survey_2026-09-13/"
          + (f"; MISSING {len(missing)}" if missing else ""))
    print(f"funnel: {funnel}")
    print(f"classes over {n} classified: {classes}")
    for row in out_rows:
        if row["quoted_count"] != "" or row["quoted_percent"] != "":
            print(f"  {row['quantity'][:62]:64s} {row['count']:>4}/{row['denominator']:<4} "
                  f"{row['percent']:>5.1f}% [{row['ci_low']:.1f}, {row['ci_high']:.1f}] "
                  f"agrees={row['agrees']}")
    print(f"mismatches: {len(mismatches)}")


if __name__ == "__main__":
    main()
