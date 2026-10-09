"""Re-test the claims unblocked by the 2026-09-16 material drop, under detectability control.

Same registered instrument as phase 2 and 2b (`run_contingency_claim`, cluster bootstrap with
5000 replicates and a fixed seed, nearest-neighbour matching inside a 0.2 SD caliper, a
size-matched random control, and the verdict rule with its 0.5 retention threshold). Nothing
about the estimand is re-tuned; this round only supplies new cohorts.

**Baseline-reproduction status is declared per claim BEFORE the run**, as in phase 2b: whether
the cohort is the authors' own data (`reproduced_on_author_cohort`), the authors' data under a
substituted negative rule (`author_positives_substituted_negatives`), or a background the authors
did not use (`author_background_substituted`). A verdict means much more in the first case, and
the report may not blur them.

**Two claims are recorded as not testable rather than forced into a table** - see NOT_TESTABLE.
**Twelve claims need structural attributes** and are listed in NEEDS_STRUCTURE with what each
requires; this round does not attempt them, and the coverage figure must say so.

**PERS-008 is run twice on purpose**: once against the authors' own whole-genome background
(DAVID's default) and once against the 472 total zinc-finger proteins. The second is the better
matched comparison, and the contrast between them is itself evidence about background choice -
which is this project's phase-1 result.

Writes results/phase2d_claim_retest.csv, _strata.csv and _audit.json.
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

import phase2d_claim_cohorts as cohorts_module
from run_phase2_claims_under_detectability_control import (
    CALIPER_SD, REPLICATES, SEED, _strata_detail, bootstrap_interval, log_odds_ratio,
    make_weighted_table, run_contingency_claim, sha256_of,
)
from run_phase2b_claims_backfill import (
    classify_null, crosses_zero, secondary_attribute_effect,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
OUT_CSV = os.path.join(RESULTS, "phase2d_claim_retest.csv")
OUT_STRATA = os.path.join(RESULTS, "phase2d_claim_retest_strata.csv")
OUT_AUDIT = os.path.join(RESULTS, "phase2d_claim_retest_audit.json")
SCRIPT = os.path.abspath(__file__)

BASELINE_REPRODUCTION = {
    "SFE-002": "author_positives_substituted_negatives",
    "SNO-004": "reproduced_on_author_cohort",
    "SNO-005": "reproduced_on_author_cohort",
    "SNO-006": "reproduced_on_author_cohort",
    "SNO-014": "author_positives_substituted_negatives",
    "SNO-008": "author_background_substituted",
    "SFI-001": "reproduced_on_author_cohort",
    "SFI-002": "reproduced_on_author_cohort",
    "PERS-008": "reproduced_on_author_cohort",
    "PERS-008_zf_background": "author_background_substituted",
}

NOT_TESTABLE = [
    {"claim_id": "SFE-003",
     "reason": "qualitative_cross_modification_comparison_no_author_statistic",
     "detail": ("作者把次磺酰化基序与他人报道的 SNO、亲电烷基化基序做定性对比，"
                "既没给基序间差异的统计量，也没有属于本文的阴性集——对照是别人论文的基序。"
                "强行装成 2×2 只能检验我们自己挑的对照，不是这条论断。")},
    {"claim_id": "PERS-007",
     "reason": "attribute_published_only_for_positives",
     "detail": ("锌指配位型（CCCC / C2H2 / CCCH-CCHC）在 Table S2 里只对 119 个持硫化锌指标注，"
                "472 个总锌指表没有型别列，所以 2×2 的阴性一侧在原始材料里不存在。"
                "自己按序列判型会变成检验我们的分类器而不是作者的论断。"
                "这一条本身就是本文论点的例证：属性只在被修饰的那一类上被公布。")},
]

NEEDS_STRUCTURE = [
    {"claim_id": "SFE-001", "attribute": "solvent_accessibility", "needs": "per-site SASA"},
    {"claim_id": "SFE-004", "attribute": "intrinsic_reactivity", "needs": "predicted thiol reactivity / pKa"},
    {"claim_id": "SNO-001", "attribute": "secondary_structure", "needs": "DSSP-class per site"},
    {"claim_id": "SNO-002", "attribute": "solvent_accessibility", "needs": "per-site SASA"},
    {"claim_id": "SNO-003", "attribute": "pKa_electrostatic", "needs": "predicted Cys pKa"},
    {"claim_id": "SNO-007", "attribute": "pKa_electrostatic", "needs": "predicted Cys pKa"},
    {"claim_id": "SNO-009", "attribute": "secondary_structure", "needs": "DSSP-class per site"},
    {"claim_id": "SNO-010", "attribute": "pKa_electrostatic", "needs": "predicted Cys pKa"},
    {"claim_id": "SNO-011", "attribute": "pKa_electrostatic", "needs": "predicted Cys pKa"},
    {"claim_id": "SNO-012", "attribute": "solvent_accessibility", "needs": "per-site SASA"},
    {"claim_id": "SNO-013", "attribute": "conservation", "needs": "orthologue alignment per site"},
    {"claim_id": "SNO-015", "attribute": "pKa_electrostatic", "needs": "predicted Cys pKa"},
]

FIELDNAMES = ["claim_id", "caliber", "verdict", "claim_direction", "baseline_reproduction",
              "unit", "n_observations", "n_positive", "n_attribute", "n_proteins",
              "covariate_set", "attribute_label",
              "baseline_log2_or", "baseline_ci_low", "baseline_ci_high", "baseline_fisher_p",
              "baseline_fold_enrichment", "baseline_direction_agrees",
              "propensity_auc", "propensity_auc_ci_low", "propensity_auc_ci_high",
              "stratified_log2_or", "stratified_ci_low", "stratified_ci_high",
              "matched_log2_or", "matched_ci_low", "matched_ci_high",
              "n_matched_pairs", "matched_positive_retention", "matched_propensity_auc",
              "random_control_log2_or", "random_control_ci_low", "random_control_ci_high",
              "retention_ratio", "precheck_matching_verdict", "precheck_permutation_verdict",
              "secondary_attribute_label", "secondary_baseline_log2_or",
              "secondary_matched_log2_or", "secondary_retention_ratio",
              "null_resolution_log2_or", "primary_reason"]


def _iv(record, key, which):
    v = record.get(key)
    if not v:
        return ""
    return v[0] if which == "low" else v[1]


def run_one(claim_id, builder, caliber="primary"):
    cohort = builder()
    record = run_contingency_claim(cohort)
    record["claim_direction"] = cohort.get("claim_direction", "positive_preference")
    record["n_proteins"] = cohort.get("n_proteins")

    from run_phase2_claims_under_detectability_control import (
        nearest_neighbour_match, propensity_score,
    )
    score, logit = propensity_score(cohort["covariates"], cohort["y"].astype(int))
    caliper = CALIPER_SD * float(np.std(logit))
    positive_index, negative_index = nearest_neighbour_match(
        logit, cohort["y"].astype(int), caliper)
    assert int(len(positive_index)) == int(record["n_matched_pairs"]), claim_id
    record["_matched_index"] = (np.concatenate([positive_index, negative_index])
                                if len(positive_index) else np.asarray([], dtype=int))
    if cohort.get("attribute_secondary") is not None:
        record["secondary_attribute"] = secondary_attribute_effect(
            cohort, record, np.asarray(cohort["attribute_secondary"]).astype(int))

    baseline = float(record["baseline_log2_or"])
    baseline_interval = record.get("baseline_interval")
    if crosses_zero(baseline_interval):
        record["baseline_direction_agrees"] = "no_baseline_effect"
    elif record["claim_direction"] == "null_no_preference":
        record["baseline_direction_agrees"] = "no"
    else:
        record["baseline_direction_agrees"] = "yes" if baseline > 0 else "no"

    if record["claim_direction"] == "null_no_preference":
        verdict, reason = classify_null(baseline_interval, record.get("matched_interval"))
        record["verdict_by_phase2_rule"] = record["verdict"]
        record["verdict"] = verdict
        record["primary_reason"] = reason
        iv = record.get("matched_interval") or [float("nan")] * 2
        record["null_resolution_log2_or"] = (
            float((iv[1] - iv[0]) / 2.0) if np.isfinite(iv[0]) and np.isfinite(iv[1]) else None)
    elif record["baseline_direction_agrees"] == "no":
        record["verdict_by_phase2_rule"] = record["verdict"]
        record["verdict"] = "baseline_contradicts_claim"
        record["primary_reason"] = (
            "补充 3：论断主张正向偏好，而本队列基线 log2 优势比 %.4f、区间 [%.4f, %.4f] 不跨零，"
            "方向相反——单列一类" % (baseline, baseline_interval[0], baseline_interval[1]))

    key = claim_id if caliber == "primary" else f"{claim_id}_{caliber}"
    record["baseline_reproduction"] = BASELINE_REPRODUCTION[key]
    record["caliber"] = caliber
    record["_cohort_notes"] = cohort["reproduction_notes"]
    record.pop("_matched_index", None)
    return record, cohort


def row_of(record):
    mt = record.get("matched_log2_or")
    b = record.get("baseline_log2_or")
    sec = record.get("secondary_attribute") or {}
    return {
        "claim_id": record["claim_id"], "caliber": record["caliber"],
        "verdict": record["verdict"], "claim_direction": record["claim_direction"],
        "baseline_reproduction": record["baseline_reproduction"],
        "unit": record["unit"], "n_observations": record["n_observations"],
        "n_positive": record["n_positive"], "n_attribute": record["n_attribute"],
        "n_proteins": record.get("n_proteins", ""),
        "covariate_set": record["covariate_set"], "attribute_label": record["attribute_label"],
        "baseline_log2_or": b,
        "baseline_ci_low": _iv(record, "baseline_interval", "low"),
        "baseline_ci_high": _iv(record, "baseline_interval", "high"),
        "baseline_fisher_p": record.get("baseline_fisher_p", ""),
        "baseline_fold_enrichment": record.get("baseline_fold_enrichment", ""),
        "baseline_direction_agrees": record.get("baseline_direction_agrees", ""),
        "propensity_auc": record.get("propensity_auc", ""),
        "propensity_auc_ci_low": _iv(record, "propensity_auc_interval", "low"),
        "propensity_auc_ci_high": _iv(record, "propensity_auc_interval", "high"),
        "stratified_log2_or": record.get("stratified_log2_or", ""),
        "stratified_ci_low": _iv(record, "stratified_interval", "low"),
        "stratified_ci_high": _iv(record, "stratified_interval", "high"),
        "matched_log2_or": mt,
        "matched_ci_low": _iv(record, "matched_interval", "low"),
        "matched_ci_high": _iv(record, "matched_interval", "high"),
        "n_matched_pairs": record.get("n_matched_pairs", ""),
        "matched_positive_retention": record.get("matched_positive_retention", ""),
        "matched_propensity_auc": record.get("matched_propensity_auc", ""),
        "random_control_log2_or": record.get("random_control_log2_or", ""),
        "random_control_ci_low": _iv(record, "random_control_interval", "low"),
        "random_control_ci_high": _iv(record, "random_control_interval", "high"),
        "retention_ratio": (abs(mt) / abs(b)) if (mt is not None and b) else "",
        "precheck_matching_verdict": (record.get("precheck_cleaning_collinearity") or {}).get("verdict", ""),
        "precheck_permutation_verdict": (record.get("precheck_permutation_degeneracy") or {}).get("verdict", ""),
        "secondary_attribute_label": sec.get("attribute_label", ""),
        "secondary_baseline_log2_or": sec.get("baseline_log2_or", ""),
        "secondary_matched_log2_or": sec.get("matched_log2_or", ""),
        "secondary_retention_ratio": sec.get("retention_ratio", ""),
        "null_resolution_log2_or": record.get("null_resolution_log2_or", ""),
        "primary_reason": record.get("primary_reason", ""),
    }


def main():
    started = time.time()
    only = sys.argv[1:] or None
    rows, strata_rows, records = [], [], {}
    plan = [(cid, fn, "primary") for cid, fn in cohorts_module.COHORTS.items()]
    plan.append(("PERS-008", None, "zf_background"))
    for claim_id, builder, caliber in plan:
        if only and claim_id not in only:
            continue
        if caliber == "zf_background":
            bg = cohorts_module.pers008_restricted_background()
            builder = lambda bg=bg: _pers008_with_background(bg)
        t0 = time.time()
        record, cohort = run_one(claim_id, builder, caliber)
        rows.append(row_of(record))
        records[f"{claim_id}|{caliber}"] = {
            k: v for k, v in record.items() if not k.startswith("_")}
        for s in (record.get("strata_detail") or []):
            strata_rows.append(dict(claim_id=claim_id, caliber=caliber,
                                    stratum_kind="propensity_quintile", **s))
        print(f"  {claim_id} [{caliber}] {record['verdict']} "
              f"baseline={record['baseline_log2_or']:.4f} "
              f"matched={record.get('matched_log2_or')} "
              f"({round(time.time()-t0, 1)}s)", flush=True)

    with open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(rows)
    if strata_rows:
        with open(OUT_STRATA, "w", encoding="utf-8-sig", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(strata_rows[0]))
            w.writeheader()
            w.writerows(strata_rows)

    audit = {
        "script": "scripts/run_phase2d_claims.py", "script_sha256": sha256_of(SCRIPT),
        "cohort_module_sha256": sha256_of(os.path.join(ROOT, "scripts", "phase2d_claim_cohorts.py")),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.node(),
        "round": "phase 2d - claims unblocked by the 2026-09-16 material drop",
        "instrument": {"replicates": REPLICATES, "seed": SEED, "caliper_sd": CALIPER_SD,
                       "source": "run_phase2_claims_under_detectability_control (unchanged)"},
        "claims_run": [r["claim_id"] + "|" + r["caliber"] for r in rows],
        "verdicts": {r["claim_id"] + "|" + r["caliber"]: r["verdict"] for r in rows},
        "baseline_reproduction": {r["claim_id"] + "|" + r["caliber"]: r["baseline_reproduction"]
                                  for r in rows},
        "not_testable": NOT_TESTABLE,
        "needs_structure": NEEDS_STRUCTURE,
        "cohort_notes": {k: v.get("_cohort_notes", []) if isinstance(v, dict) else []
                         for k, v in {f"{r['claim_id']}|{r['caliber']}": {} for r in rows}.items()},
        "records": records,
        "versions": {"python": sys.version, "numpy": np.__version__},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    json.dump(audit, open(OUT_AUDIT, "w"), ensure_ascii=False, indent=2, default=str)
    print("rows", len(rows), "| not testable", len(NOT_TESTABLE),
          "| needs structure", len(NEEDS_STRUCTURE),
          f"| {audit['elapsed_minutes']} min")


def _pers008_with_background(background):
    """PERS-008 against the 472 total zinc fingers instead of the whole proteome."""
    import re
    from phase2d_claim_cohorts import (
        GO_UBIQUITIN_CATABOLIC, _li_sheet, _protein_cohort, read_proteome_tsv,
    )
    ssh = {str(r.get("SSH ZFs") or "").strip() for r in _li_sheet("MEF SSH ZFs total")}
    ssh = {a for a in ssh if re.fullmatch(r"[A-Z0-9]{6,10}", a)}
    proteome = read_proteome_tsv("mmu")
    cohort = _protein_cohort(
        ssh, proteome, "PERS-008", [GO_UBIQUITIN_CATABOLIC],
        "注释到泛素依赖的蛋白分解（GO:0006511 及其后代）",
        {"kind": "pathway_enrichment", "value": None, "n_positive": len(ssh), "n_negative": None,
         "text": "GO 与 KEGG 均指向 ubiquitin-mediated proteolysis；阈值 0.05，正文未给具体 p 值"},
        ["次口径：背景改为 Table S2「MEF Total ZFs」的 472 个锌指蛋白，而不是全蛋白组",
         "这是匹配得多的对照——同为锌指蛋白，差别只在是否被检出持硫化",
         "作者用的是 DAVID 默认全基因组背景，那正是第一阶段判定会放大可检出性份额的构造"],
        "positive_preference", background=background)
    return cohort


if __name__ == "__main__":
    main()
