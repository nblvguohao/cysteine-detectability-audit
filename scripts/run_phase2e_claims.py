"""Re-test the five structural claims under detectability control.

Same registered instrument as phase 2, 2b and 2d - `run_contingency_claim`, cluster bootstrap with
5000 replicates and a fixed seed, nearest-neighbour matching inside a 0.2 SD caliper, a size-matched
random control, and the verdict rule with its 0.5 retention threshold. The driver `run_one` is
IMPORTED from `run_phase2d_claims` rather than copied, and its per-claim reproduction-level table is
extended in place (`BASELINE_REPRODUCTION.update`) so one implementation governs both rounds.

**A fourth reproduction level is introduced and declared here**:
`author_cohort_attribute_substituted` - the cohort, negatives and threshold are the authors' own,
but the attribute is computed differently (AlphaFold + Shrake-Rupley instead of NetSurfP, homology
models, or experimental structures; P-SEA instead of DSSP). A verdict at this level says the claim
holds or fails *given a defensible reconstruction of the attribute*, not that the authors' own
measurement was reproduced.

**Each claim runs twice**: the primary caliber over every site with a structural value, and a
declared `plddt70` sensitivity restricted to residues with pLDDT >= 70. The restriction is a
sensitivity rather than the primary because it conditions on order and order correlates with the
accessibility being measured.

Writes results/phase2e_claim_retest.csv, _strata.csv and _audit.json.
"""
from __future__ import annotations

import csv
import json
import os
import platform
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import phase2e_claim_cohorts as cohorts_module
import run_phase2d_claims as p2d
from run_phase2_claims_under_detectability_control import (
    CALIPER_SD, REPLICATES, SEED, sha256_of,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
OUT_CSV = os.path.join(RESULTS, "phase2e_claim_retest.csv")
OUT_STRATA = os.path.join(RESULTS, "phase2e_claim_retest_strata.csv")
OUT_AUDIT = os.path.join(RESULTS, "phase2e_claim_retest_audit.json")
SCRIPT = os.path.abspath(__file__)

LEVELS = {}
for cid in ("SFE-001", "SNO-002", "SNO-012", "SNO-001", "SNO-009"):
    LEVELS[cid] = "author_cohort_attribute_substituted"
    LEVELS[f"{cid}_plddt70"] = "author_cohort_attribute_substituted"
p2d.BASELINE_REPRODUCTION.update(LEVELS)

EXTRA = ["n_sites_dropped_no_structure", "n_positives_lost"]


def author_composition_check(cohort, record):
    """The authors report a composition; this is the only check that the threshold is right."""
    y = np.asarray(cohort["y"]).astype(int)
    a = np.asarray(cohort["attribute"]).astype(int)
    return {"attribute_rate_positive": round(float(a[y == 1].mean()), 4) if (y == 1).any() else None,
            "attribute_rate_negative": round(float(a[y == 0].mean()), 4) if (y == 0).any() else None,
            "author_reported": cohort["author_statistic"].get("value")}


def main():
    started = time.time()
    only = sys.argv[1:] or None
    rows, strata_rows, records = [], [], {}
    plan = []
    for cid, builder in cohorts_module.COHORTS.items():
        plan.append((cid, builder, "primary", False))
        plan.append((cid, builder, "plddt70", True))
    for claim_id, builder, caliber, plddt_only in plan:
        if only and claim_id not in only:
            continue
        t0 = time.time()
        record, cohort = p2d.run_one(
            claim_id, (lambda b=builder, p=plddt_only: b(plddt_only=p)),
            caliber if caliber == "primary" else "plddt70")
        record["author_composition_check"] = author_composition_check(cohort, record)
        for k in EXTRA:
            record[k] = cohort.get(k)
        row = p2d.row_of(record)
        row.update({k: cohort.get(k) for k in EXTRA})
        row["attribute_rate_positive"] = record["author_composition_check"]["attribute_rate_positive"]
        row["attribute_rate_negative"] = record["author_composition_check"]["attribute_rate_negative"]
        row["author_reported_composition"] = record["author_composition_check"]["author_reported"]
        rows.append(row)
        records[f"{claim_id}|{caliber}"] = {k: v for k, v in record.items()
                                            if not k.startswith("_")}
        for s in (record.get("strata_detail") or []):
            strata_rows.append(dict(claim_id=claim_id, caliber=caliber,
                                    stratum_kind="propensity_quintile", **s))
        print(f"  {claim_id} [{caliber}] {record['verdict']} "
              f"baseline={record['baseline_log2_or']:.4f} matched={record.get('matched_log2_or')} "
              f"attr+={row['attribute_rate_positive']} attr-={row['attribute_rate_negative']} "
              f"({round(time.time()-t0, 1)}s)", flush=True)

    fieldnames = p2d.FIELDNAMES + EXTRA + [
        "attribute_rate_positive", "attribute_rate_negative", "author_reported_composition"]
    with open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    if strata_rows:
        with open(OUT_STRATA, "w", encoding="utf-8-sig", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(strata_rows[0]))
            w.writeheader()
            w.writerows(strata_rows)

    feature_audit = json.load(open(os.path.join(RESULTS, "structural_site_features_audit.json")))
    audit = {
        "script": "scripts/run_phase2e_claims.py", "script_sha256": sha256_of(SCRIPT),
        "cohort_module_sha256": sha256_of(os.path.join(ROOT, "scripts", "phase2e_claim_cohorts.py")),
        "driver_imported_from": "scripts/run_phase2d_claims.py",
        "driver_sha256": sha256_of(os.path.join(ROOT, "scripts", "run_phase2d_claims.py")),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.node(),
        "round": "phase 2e - the five structural claims with an explicit author threshold",
        "instrument": {"replicates": REPLICATES, "seed": SEED, "caliper_sd": CALIPER_SD,
                       "source": "run_phase2_claims_under_detectability_control (unchanged)"},
        "new_reproduction_level": {
            "author_cohort_attribute_substituted": (
                "cohort, negatives and threshold are the authors' own; the attribute is computed "
                "differently (AlphaFold + Shrake-Rupley rather than NetSurfP / homology models / "
                "experimental structures; P-SEA rather than DSSP)")},
        "plddt_caliber": ("plddt70 is a declared sensitivity, not the primary: restricting to "
                          "confident residues conditions on order, which correlates with the "
                          "accessibility being measured"),
        "structural_features": {
            "script": feature_audit["script"], "script_sha256": feature_audit["script_sha256"],
            "source": feature_audit["source"], "sasa_method": feature_audit["sasa_method"],
            "sse_method": feature_audit["sse_method"],
            "thresholds": feature_audit["thresholds"], "counts": feature_audit["counts"]},
        "verdicts": {r["claim_id"] + "|" + r["caliber"]: r["verdict"] for r in rows},
        "author_composition_checks": {
            r["claim_id"] + "|" + r["caliber"]: {
                "positive": r["attribute_rate_positive"], "negative": r["attribute_rate_negative"],
                "author_reported": r["author_reported_composition"]} for r in rows},
        "claims_still_not_attempted": feature_audit["claims_not_attempted"],
        "records": records,
        "versions": {"python": sys.version, "numpy": np.__version__},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    json.dump(audit, open(OUT_AUDIT, "w"), ensure_ascii=False, indent=2, default=str)
    print("rows", len(rows), f"| {audit['elapsed_minutes']} min")


if __name__ == "__main__":
    main()
