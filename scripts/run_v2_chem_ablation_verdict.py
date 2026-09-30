"""Read the three-arm basic-composition ablation with the already-registered instrument.

`scripts/run_v2_chem_basic_ablation.py` (sha256 48c38ab0486a6799bb7005456ef5b3e4aa05df47f320e82e7a03c0e096d35cf1)
refit the deployed chemistry stack on amax in three arms - `control` (1046 chem columns),
`drop_w710` (minus the four win7/win10 basic-composition columns) and `drop_all_basic` (minus all
twenty basic-composition columns) - and deliberately computed NO verdict statistic there, so that
the reading happens with the same local instrument that produced the published numbers.

**This script is that reading. It fits nothing and it changes no stored artefact.**

**The verdict statistic, fixed before this script was written** (quoted from the ablation's own
docstring): the top-100 DISTAL (6-12) cleavage-boundary log2 odds ratio under TRYPSIN, computed by
`scripts/run_self_audit_five_proteases.py`'s instrument on the out-of-fold scores the job returned.
Every statistical function here is imported from the registered modules rather than reimplemented:
`band_flags` / `sequences` (band definitions and proteome), `PROTEASES` / `boundaries_for`
(cleavage rules), and `make_weighted_table` / `log_odds_ratio` / `make_weighted_auc` /
`bootstrap_interval` with REPLICATES=5000 and SEED=20260915 (cluster bootstrap over homology
components). The comparison anchor is the stored trypsin distal top-100 value for the deployed
wet-lab blend, 3.7116 [2.3163, 5.5248] (`results/self_audit_five_proteases_topk.csv`).

**Reading rule, quoted from the ablation docstring, fixed before the run:**

  1. `control` reproduces the stored enrichment (its interval covers 3.7116) AND the ablated arms'
     enrichment intervals cross zero -> the ranking IS produced by the basic-composition features.
  2. `control` reproduces AND the ablated arms keep their enrichment -> the ranking survives
     removal; it comes from elsewhere in the feature set, and the declared-chemistry-feature
     explanation is refuted.
  3. `control` does NOT reproduce -> this machine is not comparable to the 2026-09-13 run; the
     ablation is read arm-to-arm only and NO comparison with any published number may be made.

**Declared as SECONDARY before running, not part of the verdict:**

* the same top-100 distal statistic under LysC - the only other protease whose distal band
  replicated in `reports/SELF_AUDIT_FIVE_PROTEASES.md`, so it is the one place where an ablation
  effect could be checked against an independent cleavage rule;
* the proximal (<=5) band under both proteases, which that round declared a CONTROL rather than a
  second test;
* the within-label-stratum AUC of score against the distal attribute (strata `all`,
  `label_positive`, `label_negative`), the clause the trypsin round's verdict rested on, because
  within a stratum the label is constant and an association cannot be inherited from it.

**Declared limits.**

* Intervals for different arms are NOT a paired difference test. No paired test between arms was
  predeclared and none is run; overlapping or non-overlapping intervals are reported as such and
  no difference is called significant on that basis.
* Task performance per arm is TRANSCRIBED from `results/v2_chem_ablation_audit.json` (the amax
  job's own `within_protein_auc_mean` for the blend), not recomputed here.
* The three arms were refit on amax (Python 3.8.10 / lightgbm 4.6.0) while the stored 2026-09-13
  scores came from the project venv on this machine; branch 3 exists precisely because that
  difference could break comparability.

Writes `results/v2_chem_ablation_verdict.csv` and `results/v2_chem_ablation_verdict_audit.json`.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from run_cross_protease_detectability_probe import PROTEASES
from run_phase2_claims_under_detectability_control import (
    REPLICATES, SEED, bootstrap_interval, log_odds_ratio, make_weighted_auc,
    make_weighted_table, sha256_of,
)
from run_self_audit_five_proteases import DISTAL, PROXIMAL, TOP_K, band_flags, read_csv, sequences

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
SCRIPT = os.path.abspath(__file__)
ARMS = ("control", "drop_w710", "drop_all_basic")
VERDICT_PROTEASE = "Trypsin"
SECONDARY_PROTEASE = "LysC"
STORED_ANCHOR = {"log2_or": 3.7116, "ci_low": 2.3163, "ci_high": 5.5248,
                 "source": "results/self_audit_five_proteases_topk.csv",
                 "model": "v2_chem_blend_wetlab", "protease": "Trypsin", "band": "distal_6_12"}


def arm_path(arm):
    return os.path.join(RESULTS, f"v2_chem_ablation_{arm}_oof.csv")


def main():
    started = time.time()
    seqs = sequences()
    arms = {arm: read_csv(arm_path(arm)) for arm in ARMS}

    reference = arms[ARMS[0]]
    keys = [(r["accession"], r["position"]) for r in reference]
    for arm, rows in arms.items():
        assert len(rows) == len(reference), f"{arm}: row count differs"
        assert [(r["accession"], r["position"]) for r in rows] == keys, f"{arm}: site order differs"

    y = np.asarray([int(r["label"]) for r in reference], dtype=int)
    groups = np.asarray([r["component"] for r in reference])
    scores = {arm: np.asarray([float(r["blend"]) for r in rows], dtype=float)
              for arm, rows in arms.items()}

    bands = {}
    for protease in (VERDICT_PROTEASE, SECONDARY_PROTEASE):
        rule = PROTEASES[protease]
        distal, proximal, peplen = [], [], []
        for r in reference:
            a, b, L = band_flags(seqs[r["accession"]], int(r["position"]), rule)
            distal.append(a); proximal.append(b); peplen.append(L)
        bands[protease] = {"distal_6_12": np.asarray(distal, dtype=int),
                           "proximal_le_5": np.asarray(proximal, dtype=int),
                           "peptide_length": np.asarray(peplen, dtype=float)}

    rows_out = []
    for protease, band_set in bands.items():
        role = "verdict" if protease == VERDICT_PROTEASE else "secondary"
        for band in ("distal_6_12", "proximal_le_5"):
            attr = band_set[band]
            band_role = role if band == "distal_6_12" else "control_band"
            for arm in ARMS:
                score = scores[arm]
                order = np.argsort(-score, kind="stable")
                in_top = np.zeros(len(score), dtype=int)
                in_top[order[:TOP_K]] = 1
                table = make_weighted_table(in_top, attr)
                lor = log_odds_ratio(*table(np.ones(len(score))))
                liv, _ = bootstrap_interval(lambda w: log_odds_ratio(*table(w)), groups)
                rows_out.append({
                    "statistic": "top_k_log2_odds_ratio", "role": band_role, "protease": protease,
                    "band": band, "arm": arm, "stratum": "all", "k": TOP_K, "n": len(score),
                    "rate_in_top_k": round(float(attr[order[:TOP_K]].mean()), 4),
                    "rate_rest": round(float(attr[order[TOP_K:]].mean()), 4),
                    "value": round(float(lor), 4),
                    "ci_low": round(liv[0], 4), "ci_high": round(liv[1], 4),
                    "mean_peptide_length_top_k": round(
                        float(band_set["peptide_length"][order[:TOP_K]].mean()), 2),
                    "mean_peptide_length_cohort": round(float(band_set["peptide_length"].mean()), 2),
                })
                if band != "distal_6_12":
                    continue
                for stratum, mask in (("all", np.ones(len(y), dtype=bool)),
                                      ("label_positive", y == 1), ("label_negative", y == 0)):
                    fn = make_weighted_auc(attr[mask], score[mask])
                    point = fn(np.ones(int(mask.sum())))
                    iv, _ = bootstrap_interval(fn, groups[mask])
                    rows_out.append({
                        "statistic": "within_stratum_auc_score_to_attribute", "role": "secondary",
                        "protease": protease, "band": band, "arm": arm, "stratum": stratum,
                        "k": "", "n": int(mask.sum()),
                        "rate_in_top_k": "", "rate_rest": round(float(attr[mask].mean()), 4),
                        "value": round(float(point), 4),
                        "ci_low": round(iv[0], 4), "ci_high": round(iv[1], 4),
                        "mean_peptide_length_top_k": "", "mean_peptide_length_cohort": "",
                    })

    fields = ["statistic", "role", "protease", "band", "arm", "stratum", "k", "n",
              "rate_in_top_k", "rate_rest", "value", "ci_low", "ci_high",
              "mean_peptide_length_top_k", "mean_peptide_length_cohort"]
    out_csv = os.path.join(RESULTS, "v2_chem_ablation_verdict.csv")
    tmp = out_csv + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows_out:
            assert set(row.keys()) == set(fields), row
            writer.writerow(row)
    os.replace(tmp, out_csv)

    verdict_rows = {r["arm"]: r for r in rows_out
                    if r["statistic"] == "top_k_log2_odds_ratio"
                    and r["protease"] == VERDICT_PROTEASE and r["band"] == "distal_6_12"}
    control = verdict_rows["control"]
    control_reproduces = bool(control["ci_low"] <= STORED_ANCHOR["log2_or"] <= control["ci_high"])
    ablated_cross_zero = {arm: bool(verdict_rows[arm]["ci_low"] <= 0.0 <= verdict_rows[arm]["ci_high"])
                          for arm in ARMS[1:]}
    if not control_reproduces:
        branch = "3_control_does_not_reproduce_read_arm_to_arm_only"
    elif all(ablated_cross_zero.values()):
        branch = "1_ranking_is_produced_by_the_basic_composition_features"
    elif not any(ablated_cross_zero.values()):
        branch = "2_ranking_survives_removal_declared_feature_explanation_refuted"
    else:
        branch = "mixed_not_covered_by_a_single_predeclared_branch"

    job_audit = json.load(open(os.path.join(RESULTS, "v2_chem_ablation_audit.json"), encoding="utf-8"))
    task = {}
    for entry in job_audit["arms"]:
        blend = [s for s in entry["summaries"] if s["model"] == "blend"][0]
        task[entry["arm"]] = {
            "n_columns": entry["n_columns"],
            "n_dropped": len(entry["dropped_columns"]),
            "within_protein_auc_mean": blend["within_protein_auc_mean"],
            "roc_auc": blend["roc_auc"], "top1": blend["top1"],
            "transcribed_from": "results/v2_chem_ablation_audit.json",
        }

    audit = {
        "script": "scripts/run_v2_chem_ablation_verdict.py",
        "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "round": "reading of the three-arm basic-composition ablation with the registered instrument",
        "fits_nothing": True,
        "verdict_statistic": "top-100 distal_6_12 log2 odds ratio under Trypsin, cluster bootstrap over homology components",
        "replicates": REPLICATES, "seed": SEED,
        "stored_anchor": STORED_ANCHOR,
        "control_interval_covers_stored_anchor": control_reproduces,
        "ablated_arm_interval_crosses_zero": ablated_cross_zero,
        "branch": branch,
        "verdict_rows": {arm: {k: verdict_rows[arm][k] for k in ("value", "ci_low", "ci_high",
                                                                 "rate_in_top_k", "rate_rest")}
                         for arm in ARMS},
        "task_performance_transcribed": task,
        "no_paired_test_between_arms": True,
        "inputs": {f"results/v2_chem_ablation_{arm}_oof.csv": sha256_of(arm_path(arm)) for arm in ARMS},
        "input_audit": {"results/v2_chem_ablation_audit.json":
                        sha256_of(os.path.join(RESULTS, "v2_chem_ablation_audit.json")),
                        "ablation_script_sha256_reported_by_job": job_audit["script_sha256"]},
        "instrument": {"scripts/run_self_audit_five_proteases.py":
                       sha256_of(os.path.join(ROOT, "scripts", "run_self_audit_five_proteases.py"))},
        "versions": {"python": sys.version, "numpy": np.__version__},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    out_json = os.path.join(RESULTS, "v2_chem_ablation_verdict_audit.json")
    tmp = out_json + ".tmp"
    json.dump(audit, open(tmp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    os.replace(tmp, out_json)

    print("control reproduces stored anchor:", control_reproduces, control["value"],
          [control["ci_low"], control["ci_high"]])
    for arm in ARMS:
        r = verdict_rows[arm]
        print(f"  {arm:15s} log2OR={r['value']:7.4f} [{r['ci_low']}, {r['ci_high']}] "
              f"top100={r['rate_in_top_k']} rest={r['rate_rest']} "
              f"withinProteinAUC={task[arm]['within_protein_auc_mean']:.4f}")
    print("BRANCH:", branch)
    print(f"done in {audit['elapsed_minutes']} min")


if __name__ == "__main__":
    main()
