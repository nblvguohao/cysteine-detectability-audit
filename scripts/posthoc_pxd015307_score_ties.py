"""POST-HOC, NOT PREDECLARED: when a spectrum admits both +32 chemistries, can its score tell them apart?

The registered round (`scripts/run_pxd015307_research_mass_accuracy.py`, registration
`protocols/pxd015307_research_preregistration_2026-09-17.json`) landed on its predeclared branch 4:
at 1% FDR the re-search of PXD015307 returns ZERO +32-class cysteine PSMs, in the proteome-scale
file and in all four targeted arms, so no rate may be quoted.

While reading that result a different pattern was visible in the output lines, and this script
measures it. **It is post-hoc exploration. It was not predeclared, it carries no verdict, and
nothing in it may be reported as a test.**

**What it measures.** Comet was asked for the top 5 lines per spectrum with both Sulfide
(+31.972071) and Dioxidation (+31.989829) available on cysteine. For every spectrum whose output
lines contain at least one +32-class assignment, this script groups the lines by their xcorr and
asks: among the lines sharing the BEST xcorr for that spectrum, how many distinct
sulfide/dioxidation combinations appear on the same peptide backbone, and how far apart are their
precursor mass errors?

A spectrum where several chemistries share one score is a spectrum whose FRAGMENTS do not decide
the chemistry; whatever decides it is then the precursor mass, the search space, or the tolerance.
That is the claim this project makes about site-level persulfidation data, here visible inside a
published deposit's own spectra.

**Columns.** `arm`, `scan`, `peptide`, `best_xcorr`, `n_lines_at_best_xcorr`,
`n_distinct_plus32_patterns_at_best_xcorr`, `chemistries_at_best_xcorr` (which of sulfide /
dioxidation appear), `ppm_min` / `ppm_max` / `ppm_span` over those lines, `e_value_best`.

**Limits.** These are low-confidence identifications - none of them passes the round's 1% FDR - and
this script does not filter by confidence, precisely because the point is about score ties rather
than about identification. The targeted arms are a bovine insulin model substrate. The deposit's
own search assigned the HUMAN insulin A-chain sequence (GIVEQCCTSICSLYQLENYCN) where our
database, which carries bovine insulin P01317, matches the bovine variant
(GIVEQCCASVCSLYQLENYCN); the two differ at three residues and this is recorded rather than
resolved. Ion-trap fragment binning (1.0005 Da) is used for those arms because their MS2 scans are
ion trap, which makes fragment-level discrimination of a 0.0178 Da difference impossible by
construction - that is the point, not a defect.

Writes results/pxd015307_posthoc_score_ties.csv and results/pxd015307_posthoc_score_ties_audit.json.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
SEARCH_DIR = os.path.join(ROOT, "external", "sulfhydrome_bta_tmt", "research")
SCRIPT = os.path.abspath(__file__)
SULFIDE = 31.972071
DIOXIDATION = 31.989829
TOL = 1e-4


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def plus32_pattern(row):
    """((position, kind), ...) for +32-class variable modifications on cysteine, or ()."""
    peptide = row.get("plain_peptide") or ""
    raw = (row.get("modifications") or "").strip()
    out = []
    if raw and raw != "-":
        for item in raw.split(","):
            parts = item.split("_")
            if len(parts) != 3 or parts[1] != "V":
                continue
            try:
                position, delta = int(parts[0]), float(parts[2])
            except ValueError:
                continue
            if not (1 <= position <= len(peptide)) or peptide[position - 1] != "C":
                continue
            if abs(delta - SULFIDE) < TOL:
                out.append((position, "sulfide"))
            elif abs(delta - DIOXIDATION) < TOL:
                out.append((position, "dioxidation"))
    return tuple(sorted(out))


def ppm_of(row):
    exp, calc = float(row["exp_neutral_mass"]), float(row["calc_neutral_mass"])
    return (exp - calc) / calc * 1e6


def main():
    started = time.time()
    rows_out, per_arm = [], {}
    for name in sorted(os.listdir(SEARCH_DIR)):
        if not name.endswith(".txt"):
            continue
        arm = name[:-4]
        path = os.path.join(SEARCH_DIR, name)
        with open(path, encoding="utf-8", newline="") as fh:
            fh.readline()
            rows = list(csv.DictReader(fh, delimiter="\t"))
        by_scan = {}
        for r in rows:
            by_scan.setdefault(r["scan"], []).append(r)

        tied = 0
        with_plus32 = 0
        for scan, lines in by_scan.items():
            if not any(plus32_pattern(r) for r in lines):
                continue
            with_plus32 += 1
            best = max(float(r["xcorr"]) for r in lines)
            at_best = [r for r in lines if float(r["xcorr"]) == best]
            backbones = {r["plain_peptide"] for r in at_best}
            patterns = {plus32_pattern(r) for r in at_best if plus32_pattern(r)}
            if not patterns:
                continue
            chemistries = sorted({kind for pattern in patterns for _, kind in pattern})
            ppms = [ppm_of(r) for r in at_best]
            if len(patterns) > 1:
                tied += 1
            rows_out.append({
                "arm": arm, "scan": scan,
                "peptide": sorted(backbones)[0] if len(backbones) == 1 else "|".join(sorted(backbones)),
                "n_backbones_at_best_xcorr": len(backbones),
                "best_xcorr": round(best, 4),
                "n_lines_at_best_xcorr": len(at_best),
                "n_distinct_plus32_patterns_at_best_xcorr": len(patterns),
                "chemistries_at_best_xcorr": "+".join(chemistries),
                "ppm_min": round(min(ppms), 3), "ppm_max": round(max(ppms), 3),
                "ppm_span": round(max(ppms) - min(ppms), 3),
                "e_value_best": min(r["e-value"] for r in at_best),
            })
        per_arm[arm] = {
            "n_spectra_with_any_plus32_line": with_plus32,
            "n_spectra_where_several_plus32_patterns_share_the_best_xcorr": tied,
            "input_sha256": sha256_of(path),
        }

    both = [r for r in rows_out if r["chemistries_at_best_xcorr"] == "dioxidation+sulfide"]
    audit = {
        "script": "scripts/posthoc_pxd015307_score_ties.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "status": "POST-HOC exploration, not predeclared, no verdict",
        "registered_round": "protocols/pxd015307_research_preregistration_2026-09-17.json",
        "per_arm": per_arm,
        "n_rows": len(rows_out),
        "n_spectra_with_both_chemistries_at_the_best_xcorr": len(both),
        "max_ppm_span_among_those": max((r["ppm_span"] for r in both), default=None),
        "versions": {"python": sys.version},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    if rows_out:
        cols = list(rows_out[0].keys())
        tmp = os.path.join(RESULTS, "pxd015307_posthoc_score_ties.csv") + ".tmp"
        with open(tmp, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=cols)
            writer.writeheader()
            for row in rows_out:
                writer.writerow(row)
        os.replace(tmp, os.path.join(RESULTS, "pxd015307_posthoc_score_ties.csv"))
    tmp = os.path.join(RESULTS, "pxd015307_posthoc_score_ties_audit.json") + ".tmp"
    json.dump(audit, open(tmp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    os.replace(tmp, os.path.join(RESULTS, "pxd015307_posthoc_score_ties_audit.json"))

    for arm, a in sorted(per_arm.items()):
        print(f"  {arm:26s} spectra_with_plus32={a['n_spectra_with_any_plus32_line']:5d}"
              f" tied={a['n_spectra_where_several_plus32_patterns_share_the_best_xcorr']:5d}")
    print("rows:", len(rows_out), "| both chemistries tied at best xcorr:", len(both),
          "| max ppm span among those:", audit["max_ppm_span_among_those"])


if __name__ == "__main__":
    main()
