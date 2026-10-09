"""Can a published persulfidome deposit's own data decide +S from +2O, once both are in the search space?

Registered before the search ran: `protocols/pxd015307_research_preregistration_2026-09-17.json`
(with four amendments, all recorded there). This script is the READING of that search. It fits
nothing and it changes no stored artefact.

**What was searched.** PXD015307 ("Quantitative profiling of protein sulfhydrome by the BTA-TMT"),
all five raw files, re-searched with Comet 2026.01 on amax against UniProt Bos taurus UP000009136
(51,611 entries, internal concatenated decoys), precursor +/-10 ppm, fragment bin 0.02, trypsin with
2 missed cleavages, `isotope_error = 0`, at most 3 variable modifications per peptide, and five
variable modifications: Sulfide +31.972071 (C), Dioxidation +31.989829 (C), iodoTMT6plex
+329.226595 (C), Carbamidomethyl +57.021464 (C) and Oxidation +15.9949 (M). **No static cysteine
term** - run 1 carried Comet's default `add_C_cysteine = 57.021464`, which forces every cysteine to
be carbamidomethylated and makes the question unanswerable by construction; that run is kept,
unmodified, under `discarded_run1/` on amax with its reason, and only run 2 is read here.

**Definitions fixed before this script was written.**

* `separation_da = 0.017758` - Dioxidation minus Sulfide. For a PSM it is expressed as
  `separation_ppm = separation_da / calc_neutral_mass * 1e6`.
* `ppm = (exp_neutral_mass - calc_neutral_mass) / calc_neutral_mass * 1e6`, uncorrected: no C13
  offsets were searched, so a 1 Da mis-assignment cannot be absorbed into this number.
* **Baseline** - the file's own instrument error, taken from PSMs that carry NO +32-class
  modification, summarised as median, SD and IQR. No `|ppm| < 10` gate is needed here because the search tolerance
  already is +/-10 ppm, and this script asserts that.
* **Resolvable** (the rule already registered in `scripts/audit_search_space_and_mass_accuracy.py`,
  unchanged): a +32-class PSM is RESOLVABLE if `|ppm - baseline_median| < separation_ppm / 2`.
* **An arm with no baseline PSMs cannot be read.** Added while implementing this script and BEFORE
  any result was seen (the first execution stopped on exactly this case): if an arm has zero
  1%-FDR PSMs without a +32-class modification, it has no instrument baseline, so its hits are
  listed with an empty deviation, its `fraction_resolvable` is null, and it is excluded from every
  rate. Recorded in the audit under `amendments`.
* 1% FDR: rank-1 PSMs only, ordered by e-value, decoys counted by the `decoy_prefix`; the threshold
  is the largest e-value whose cumulative `decoys / targets` stays at or below 0.01.

**Primary readout, predeclared.** The fraction of +32-class cysteine PSMs that are resolvable,
over the proteome-scale file (`Fig-1D`), reported with its count and with the assignment split
(sulfide versus dioxidation).

**Secondary readouts, predeclared.** (a) head-to-head: for a scan whose output lines contain BOTH
hypotheses for the SAME peptide backbone, the difference in xcorr and in `|ppm - baseline_median|`;
(b) the four targeted arms read the same way, which re-examines the deposit's own four Sulfide PSMs
with both candidates available; (c) per-protein counts of the two assignments.

**Outcome branches, predeclared** (from the registration):
  1. most +32-class PSMs NOT resolvable -> the sixth artefact gets a public, proteome-scale
     instance: chemical identity can be a product of the search space rather than the measurement;
  2. resolvable and sulfide-consistent -> the deposit's chemistry stands;
  3. resolvable and dioxidation-consistent -> directional evidence on published data;
  4. fewer than 20 +32-class cysteine PSMs at 1% FDR -> underpowered, and no fraction is quoted
     as a rate.

**Limits carried from the registration.** Different engine from the deposit (Comet versus Sequest
HT) with no engine-to-engine control; the deposit's persulfidome is defined by BTA-switch
enrichment chemistry, not by a +32 mass, so this round does not test their enrichment-based site
calls; the bovine proteome version is not matched to theirs; e-value target-decoy FDR only, no
Percolator; protein N-terminal acetylation was not searched. No published assignment is called
incorrect.

Writes results/pxd015307_research_psms.csv, results/pxd015307_research_summary.csv and
results/pxd015307_research_audit.json.
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
SEARCH_DIR = os.path.join(ROOT, "external", "sulfhydrome_bta_tmt", "research")
SCRIPT = os.path.abspath(__file__)

SEPARATION_DA = 0.017758
SULFIDE = 31.972071
DIOXIDATION = 31.989829
MASS_TOLERANCE = 1e-4          # how close a reported delta must be to a candidate to be that candidate
FDR = 0.01
DECOY_PREFIX = "DECOY_"
PROTEOME_FILE = "Fig-1D"
MIN_PSMS_FOR_A_RATE = 20


def sha256_of(path):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_comet_txt(path):
    with open(path, encoding="utf-8", newline="") as fh:
        header_line = fh.readline()
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return header_line.strip(), rows


def mod_list(row):
    """[(position, kind, delta)] parsed from Comet's modifications column."""
    raw = (row.get("modifications") or "").strip()
    if not raw or raw == "-":
        return []
    out = []
    for item in raw.split(","):
        parts = item.split("_")
        if len(parts) != 3:
            continue
        try:
            out.append((int(parts[0]), parts[1], float(parts[2])))
        except ValueError:
            continue
    return out


def plus32_kind(row):
    """'sulfide', 'dioxidation' or None, for a variable +32-class modification on a cysteine."""
    peptide = row.get("plain_peptide") or ""
    kinds = set()
    for position, kind, delta in mod_list(row):
        if kind != "V":
            continue
        if not (1 <= position <= len(peptide)) or peptide[position - 1] != "C":
            continue
        if abs(delta - SULFIDE) < MASS_TOLERANCE:
            kinds.add("sulfide")
        elif abs(delta - DIOXIDATION) < MASS_TOLERANCE:
            kinds.add("dioxidation")
    if not kinds:
        return None
    return "both" if len(kinds) > 1 else kinds.pop()


def ppm_of(row):
    exp = float(row["exp_neutral_mass"])
    calc = float(row["calc_neutral_mass"])
    return (exp - calc) / calc * 1e6


def is_decoy(row):
    proteins = (row.get("protein") or "").split(",")
    return all(p.startswith(DECOY_PREFIX) for p in proteins if p)


def fdr_threshold(rank1):
    ordered = sorted(rank1, key=lambda r: float(r["e-value"]))
    targets = decoys = 0
    threshold = None
    for row in ordered:
        if is_decoy(row):
            decoys += 1
        else:
            targets += 1
        if targets and decoys / targets <= FDR:
            threshold = float(row["e-value"])
    return threshold


def median(values):
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        return None
    mid = n // 2
    return ordered[mid] if n % 2 else (ordered[mid - 1] + ordered[mid]) / 2


def quantile(values, q):
    ordered = sorted(values)
    if not ordered:
        return None
    idx = (len(ordered) - 1) * q
    low, high = int(idx), min(int(idx) + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * (idx - low)


def sd(values):
    n = len(values)
    if n < 2:
        return None
    mean = sum(values) / n
    return (sum((v - mean) ** 2 for v in values) / (n - 1)) ** 0.5


def main():
    started = time.time()
    arms = sorted(f[:-4] for f in os.listdir(SEARCH_DIR) if f.endswith(".txt"))
    assert PROTEOME_FILE in arms, f"{PROTEOME_FILE}.txt not found in {SEARCH_DIR}"

    psm_rows, summary_rows, per_arm = [], [], {}
    for arm in arms:
        path = os.path.join(SEARCH_DIR, arm + ".txt")
        version_line, rows = read_comet_txt(path)
        rank1 = [r for r in rows if r.get("num") == "1"]
        threshold = fdr_threshold(rank1)
        if threshold is None:
            per_arm[arm] = {"n_rank1": len(rank1), "threshold": None,
                            "note": "no e-value threshold satisfies 1% FDR"}
            continue
        passing = [r for r in rank1 if float(r["e-value"]) <= threshold and not is_decoy(r)]
        for r in passing:
            assert abs(ppm_of(r)) <= 10.0001, f"{arm}: ppm outside the searched tolerance"

        baseline_ppm = [ppm_of(r) for r in passing if plus32_kind(r) is None]
        base_median = median(baseline_ppm)
        hits = [r for r in passing if plus32_kind(r) is not None]

        has_baseline = base_median is not None
        resolvable = 0
        for r in hits:
            calc = float(r["calc_neutral_mass"])
            separation_ppm = SEPARATION_DA / calc * 1e6
            observed = ppm_of(r)
            deviation = (observed - base_median) if has_baseline else None
            ok = (abs(deviation) < separation_ppm / 2) if has_baseline else None
            resolvable += int(bool(ok))
            psm_rows.append({
                "arm": arm, "scan": r["scan"], "charge": r["charge"],
                "plain_peptide": r["plain_peptide"], "modified_peptide": r["modified_peptide"],
                "assignment": plus32_kind(r),
                "protein": (r.get("protein") or "").split(",")[0],
                "calc_neutral_mass": round(calc, 4),
                "e_value": r["e-value"], "xcorr": r["xcorr"],
                "observed_ppm": round(observed, 4),
                "baseline_median_ppm": round(base_median, 4) if has_baseline else "",
                "deviation_from_baseline_ppm": round(deviation, 4) if has_baseline else "",
                "separation_ppm": round(separation_ppm, 4),
                "resolvable": int(ok) if has_baseline else "",
            })

        # head-to-head: the same scan and backbone carrying the other hypothesis on any output line
        by_scan = {}
        for r in rows:
            by_scan.setdefault(r["scan"], []).append(r)
        head_to_head = []
        for r in hits:
            mine = plus32_kind(r)
            for other in by_scan.get(r["scan"], []):
                theirs = plus32_kind(other)
                if theirs is None or theirs == mine or other["plain_peptide"] != r["plain_peptide"]:
                    continue
                if not has_baseline:
                    continue
                head_to_head.append({
                    "arm": arm, "scan": r["scan"], "peptide": r["plain_peptide"],
                    "kind_rank1": mine, "xcorr_rank1": float(r["xcorr"]),
                    "abs_dev_rank1": round(abs(ppm_of(r) - base_median), 4),
                    "kind_other": theirs, "xcorr_other": float(other["xcorr"]),
                    "abs_dev_other": round(abs(ppm_of(other) - base_median), 4),
                })

        counts = {"sulfide": 0, "dioxidation": 0, "both": 0}
        for r in hits:
            counts[plus32_kind(r)] += 1
        per_arm[arm] = {
            "comet_version_line": version_line,
            "n_rank1_psms": len(rank1),
            "e_value_threshold_at_1pct_fdr": threshold,
            "n_psms_at_1pct_fdr": len(passing),
            "n_baseline_psms_no_plus32": len(baseline_ppm),
            "baseline_median_ppm": round(base_median, 4) if has_baseline else None,
            "baseline_sd_ppm": round(sd(baseline_ppm), 4) if sd(baseline_ppm) is not None else None,
            "baseline_iqr_ppm": [round(quantile(baseline_ppm, 0.25), 4),
                                 round(quantile(baseline_ppm, 0.75), 4)] if has_baseline else None,
            "n_plus32_cys_psms": len(hits),
            "assignment_counts": counts,
            "n_resolvable": resolvable if has_baseline else None,
            "fraction_resolvable": (round(resolvable / len(hits), 4) if hits else None) if has_baseline else None,
            "readable": has_baseline,
            "head_to_head_pairs": head_to_head,
            "input_sha256": sha256_of(path),
        }
        summary_rows.append({
            "arm": arm,
            "n_psms_at_1pct_fdr": len(passing),
            "n_baseline_psms_no_plus32": len(baseline_ppm),
            "baseline_median_ppm": per_arm[arm]["baseline_median_ppm"],
            "baseline_sd_ppm": per_arm[arm]["baseline_sd_ppm"],
            "n_plus32_cys_psms": len(hits),
            "n_sulfide": counts["sulfide"], "n_dioxidation": counts["dioxidation"],
            "n_both_on_one_psm": counts["both"],
            "n_resolvable": resolvable if has_baseline else "",
            "fraction_resolvable": per_arm[arm]["fraction_resolvable"],
            "n_head_to_head_pairs": len(head_to_head),
        })

    proteome = per_arm[PROTEOME_FILE]
    n_hits = proteome["n_plus32_cys_psms"]
    if n_hits < MIN_PSMS_FOR_A_RATE:
        branch = "4_underpowered_no_rate_quoted"
    elif proteome["fraction_resolvable"] < 0.5:
        branch = "1_majority_not_resolvable_public_proteome_scale_instance"
    else:
        counts = proteome["assignment_counts"]
        if counts["dioxidation"] > counts["sulfide"]:
            branch = "3_resolvable_and_dioxidation_consistent"
        else:
            branch = "2_resolvable_and_sulfide_consistent"

    fields = list(psm_rows[0].keys()) if psm_rows else ["arm"]
    for path, rows_out, cols in (
            (os.path.join(RESULTS, "pxd015307_research_psms.csv"), psm_rows, fields),
            (os.path.join(RESULTS, "pxd015307_research_summary.csv"), summary_rows,
             list(summary_rows[0].keys()) if summary_rows else ["arm"])):
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=cols)
            writer.writeheader()
            for row in rows_out:
                assert set(row.keys()) == set(cols), row
                writer.writerow(row)
        os.replace(tmp, path)

    audit = {
        "script": "scripts/run_pxd015307_research_mass_accuracy.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "registration": "protocols/pxd015307_research_preregistration_2026-09-17.json",
        "accession": "PXD015307",
        "separation_da": SEPARATION_DA,
        "resolvable_rule": "abs(ppm - baseline_median) < separation_ppm / 2",
        "fdr": FDR, "decoy_prefix": DECOY_PREFIX,
        "search_ran_on": "amax:/home/amax/scratch/lgh_pxd015307 (run 2; run 1 kept under discarded_run1/ with its reason)",
        "amendments": [
            {"what": "an arm with no 1%-FDR PSM lacking a +32-class modification has no instrument baseline; its hits are listed with empty deviation, its fraction_resolvable is null, and it enters no rate",
             "when": "while implementing this script, before any result was seen - the first execution stopped on this case"}
        ],
        "per_arm": per_arm,
        "primary_file": PROTEOME_FILE,
        "branch": branch,
        "versions": {"python": sys.version},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    tmp = os.path.join(RESULTS, "pxd015307_research_audit.json") + ".tmp"
    json.dump(audit, open(tmp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    os.replace(tmp, os.path.join(RESULTS, "pxd015307_research_audit.json"))

    for arm in sorted(per_arm):
        a = per_arm[arm]
        if a.get("threshold") is None and "n_plus32_cys_psms" not in a:
            print(f"  {arm:26s} {a}")
            continue
        base = f"{a['baseline_median_ppm']:+.3f}" if a['baseline_median_ppm'] is not None else "  none"
        print(f"  {arm:26s} psms={a['n_psms_at_1pct_fdr']:6d} base={base}"
              f" sd={a['baseline_sd_ppm']} +32={a['n_plus32_cys_psms']:4d}"
              f" (S={a['assignment_counts']['sulfide']}, 2O={a['assignment_counts']['dioxidation']})"
              f" resolvable={a['n_resolvable']} h2h={len(a['head_to_head_pairs'])}")
    print("BRANCH:", branch)


if __name__ == "__main__":
    main()
