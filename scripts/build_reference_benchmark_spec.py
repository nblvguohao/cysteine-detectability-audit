"""Assemble a reference benchmark from the cohorts this tree can defend, and say what it cannot answer.

No model is fitted and no statistic is recomputed: every size is read from the artefact of the
round that DEFINED that cohort. Eligibility is declared before assembly as three separate
criteria, each reported on its own rather than folded into one boolean:

  C1 tier T1 or T2 under scripts/ptm_detectability_diagnostics.label_semantics_tier - a cohort
     whose negatives are "not observed" (T3) is excluded by construction, because a benchmark
     built on it scores detectability, which is the thing under audit.
  C2 the depth-stratified caliber must be identifiable: at least 30 grouping units contributing
     to the within-caliber statistic (threshold declared in protocols/v3_qtrp_analysis_plan.json).
  C3 folds group by homology component (results/v3_human_homology_components.csv).

Two defects found by the independent verification round (results/phase3_independent_verification.csv)
are fixed here and recorded rather than quietly patched:

  * the first version keyed the QTRP arms into the CENSUS round's datasets through a name map,
    which silently loaded a differently-constructed cohort - census `qtrp_S1_ph5` reports 277
    positives while the depth round's `S1_pH5.0|raw`, the cohort actually named in this table,
    has 272. Sizes now come from each cohort's own defining artefact and `size_source` names it.
  * the first version folded all three criteria into one column called `eligible` while
    implementing only C1, which overstated what had been checked. C2 was measured only for
    v3_pooled_ph5_raw (44 components); the other three now carry "not measured" rather than a pass.
  * `n_grouping_units` mixed homology components (v3) with proteins-carrying-a-positive (the
    others). Those are two different quantities and now have two columns.

What this benchmark is NOT: it is not a correct persulfidation benchmark. Every cohort comes from
ONE laboratory's chemistry per modification, the positive counts are in the hundreds, all are
human, and the effect sizes this project needs to resolve (0.02-0.04 AUC) are below what 314
positives can decide. It is a defensible starting point with its limits attached, not a standard.

Terminology: sulfenylation = 次磺酰化 (-SOH); sulfinylation = 亚磺酰化 (-SO2H).
"""
import csv, json, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
SCRIPT = Path(__file__).resolve()
sys.path.insert(0, str(ROOT / "scripts"))
from ptm_detectability_diagnostics import label_semantics_tier, sha256, read_csv

IDENTIFIABILITY_FLOOR = 30

v3 = json.loads((RESULTS / "v3_qtrp_stack_run_audit.json").read_text(encoding="utf-8"))
census = json.loads((RESULTS / "ptm_detectability_census_audit.json").read_text(encoding="utf-8"))
cen = census.get("datasets") or census.get("cohorts") or {}
comp = read_csv(RESULTS / "v3_human_homology_components.csv")
strata = read_csv(RESULTS / "qtrp_depth_strata.csv")

CANNOT_ANSWER = ("0.02-0.04 量级的 AUC 差值（阳性量级不足）；跨实验室泛化（每种修饰只有一个实验室的化学）；"
                 "植物体系（全部为人源）；以及任何需要化学倾向而非注释排序的绝对论断")

COHORTS = [
    {"cohort": "v3_pooled_ph5_raw", "modification": "persulfidation (持硫化)",
     "negative_class": "same_run_mutually_exclusive",
     "negative_definition": "同一 QTRP 运行内只检出游离巯基形式的位点（raw，保留高深度阴性）",
     "source_round": "reports/V3_QTRP_TRAINED_MODEL.md", "caliber": "detection-depth stratified",
     "grouping": "homology component (results/v3_human_homology_components.csv)",
     "prereg": "protocols/v3_qtrp_analysis_plan.json (registered, hash-pinned)",
     "size_source": "results/v3_qtrp_stack_run_audit.json"},
    {"cohort": "S1_pH5.0|raw", "modification": "persulfidation (持硫化)",
     "negative_class": "same_run_mutually_exclusive",
     "negative_definition": "主臂：同一运行只检出游离巯基（raw）",
     "source_round": "reports/QTRP_DEPTH_CONFOUND.md", "caliber": "detection-depth stratified",
     "grouping": "protein / homology component", "prereg": "n/a (pre-v3 round)",
     "size_source": "results/qtrp_depth_strata.csv"},
    {"cohort": "S2_pH5.0|raw", "modification": "persulfidation (持硫化)",
     "negative_class": "same_run_mutually_exclusive",
     "negative_definition": "重复臂：同一运行只检出游离巯基（raw）",
     "source_round": "reports/QTRP_DEPTH_CONFOUND.md", "caliber": "detection-depth stratified",
     "grouping": "protein / homology component", "prereg": "n/a (pre-v3 round)",
     "size_source": "results/qtrp_depth_strata.csv"},
    {"cohort": "qpers_sid_tierB", "modification": "persulfidation (持硫化)",
     "negative_class": "detected_background_proteome",
     "negative_definition": "同一供体臂内被检出但未达持硫化判据的位点（T2：位点层通道非互斥）",
     "source_round": "reports/QPERS_SID_LABEL_CRITERIA.md / reports/V3_QTRP_TRAINED_MODEL.md",
     "caliber": "donor-arm depth stratified", "grouping": "protein",
     "prereg": "results/v3_plan_registration_addendum_external.json",
     "size_source": "results/ptm_detectability_census_audit.json"},
]
CENSUS_KEY = {"S1_pH5.0|raw": "qtrp_S1_ph5", "S2_pH5.0|raw": "qtrp_S2_ph5"}


def from_strata(cohort):
    rows_ = [r for r in strata if r["cohort"] == cohort]
    return (sum(int(r["n_sites"]) for r in rows_), sum(int(r["n_positive"]) for r in rows_))


rows = []
for c in COHORTS:
    tier = label_semantics_tier(negative_class=c["negative_class"])
    e = dict(c)
    e["tier"] = tier["tier"]
    e["claims_licensed"] = tier["licence"]
    if c["cohort"] == "v3_pooled_ph5_raw":
        coh = v3["cohorts"]["v3_pooled_ph5_raw"]
        e["n_sites"], e["n_positive"] = coh.get("n_sites"), coh.get("n_positive")
        e["n_homology_components"] = coh.get("n_components")
        e["n_proteins_with_positive"] = coh.get("n_proteins")
        e["identifiable_units_in_within_caliber"] = 44
        e["identifiable_units_source"] = "results/v3_qtrp_stack_summary.csv components_contributing_to_within_caliber"
    elif c["cohort"] in CENSUS_KEY:
        e["n_sites"], e["n_positive"] = from_strata(c["cohort"])
        e["n_homology_components"] = ""
        e["n_proteins_with_positive"] = cen.get(CENSUS_KEY[c["cohort"]], {}).get("n_proteins_with_positive", "")
        e["identifiable_units_in_within_caliber"] = ""
        e["identifiable_units_source"] = "not measured for this cohort"
    else:
        src = cen.get("qpers_sid_tierB", {})
        e["n_sites"] = src.get("n_positive", 0) + src.get("n_neg_B_observed_unmodified", 0)
        e["n_positive"] = src.get("n_positive", "")
        e["n_homology_components"] = ""
        e["n_proteins_with_positive"] = src.get("n_proteins_with_positive", "")
        e["identifiable_units_in_within_caliber"] = ""
        e["identifiable_units_source"] = "not measured for this cohort"

    iu = e["identifiable_units_in_within_caliber"]
    c1 = e["tier"] in ("T1_direct_adduct", "T2_parallel_ambiguous")
    c2 = (iu >= IDENTIFIABILITY_FLOOR) if isinstance(iu, int) else None
    c3 = bool(e["n_homology_components"])
    e["criterion_C1_tier_T1_or_T2"] = c1
    e["criterion_C2_identifiable_units_ge_30"] = "not measured" if c2 is None else c2
    e["criterion_C3_component_grouped_folds"] = c3
    e["eligible_on_tier_only"] = c1
    e["eligible_on_all_three_criteria"] = bool(c1 and c2 is True and c3)
    e["cannot_answer"] = CANNOT_ANSWER
    rows.append(e)

out = RESULTS / "reference_benchmark_spec.csv"
with out.open("w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

audit = {
    "script": "scripts/build_reference_benchmark_spec.py", "script_sha256": sha256(SCRIPT),
    "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "new_computation": False,
    "declared_eligibility": {
        "C1_tier": "T1 or T2 only; T3 excluded by construction because a benchmark on 'not observed' negatives scores detectability",
        "C2_identifiability": f"at least {IDENTIFIABILITY_FLOOR} grouping units in the within-caliber statistic",
        "C3_folds": "homology component",
    },
    "identifiability_floor": {"value": IDENTIFIABILITY_FLOOR,
                              "declared_in": "protocols/v3_qtrp_analysis_plan.json",
                              "plan_sha256": sha256(ROOT / "protocols" / "v3_qtrp_analysis_plan.json")},
    "eligible_on_tier_only": [r["cohort"] for r in rows if r["eligible_on_tier_only"]],
    "eligible_on_all_three_criteria": [r["cohort"] for r in rows if r["eligible_on_all_three_criteria"]],
    "excluded_by_tier": [r["cohort"] for r in rows if not r["eligible_on_tier_only"]],
    "criteria_note": ("C1 is checked for every cohort; C2 was measured only for v3_pooled_ph5_raw (44 "
                      "components), so the other three carry 'not measured' rather than a pass, and "
                      "eligible_on_all_three_criteria is True for that cohort alone"),
    "tier_counts": {t: sum(1 for r in rows if r["tier"] == t) for t in {r["tier"] for r in rows}},
    "components_file_rows": len(comp),
    "corrections_from_independent_verification": {
        "source": "results/phase3_independent_verification.csv",
        "size_source_collision": ("the QTRP arms were keyed into the census round's datasets by a name map, "
                                  "loading a differently-constructed cohort: census qtrp_S1_ph5 has 277 "
                                  "positives, the depth round's S1_pH5.0|raw has 272. Sizes now come from "
                                  "each cohort's defining artefact and size_source names it."),
        "eligibility_overstated": ("one boolean called 'eligible' implemented only the tier test; the three "
                                   "declared criteria are now reported separately."),
        "mixed_grouping_quantity": ("n_grouping_units mixed homology components with "
                                    "proteins-carrying-a-positive; split into two columns."),
        "identifiability_floor_unpinned": "the floor and its source plan are now hash-pinned in this audit.",
    },
    "not_a_standard": (
        "This is not a correct persulfidation benchmark. One laboratory's chemistry per modification, "
        "positives in the hundreds, human only, and the 0.02-0.04 AUC effects this project needs are below "
        "what 314 positives can decide. Published as a defensible starting point with its limits attached."),
    "inputs": {p: sha256(RESULTS / p) for p in (
        "v3_qtrp_stack_run_audit.json", "ptm_detectability_census_audit.json",
        "v3_human_homology_components.csv", "qtrp_depth_strata.csv", "v3_qtrp_stack_summary.csv")},
    "outputs": {"results/reference_benchmark_spec.csv": sha256(out)},
    "versions": {"python": sys.version},
}
(RESULTS / "reference_benchmark_spec_audit.json").write_text(
    json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
for r in rows:
    print(f"  {r['cohort']:22s} {r['tier']:22s} C1={r['criterion_C1_tier_T1_or_T2']} "
          f"C2={r['criterion_C2_identifiable_units_ge_30']} C3={r['criterion_C3_component_grouped_folds']} "
          f"sites={r['n_sites']} pos={r['n_positive']} comps={r['n_homology_components'] or '-'} "
          f"prot={r['n_proteins_with_positive']} src={Path(r['size_source']).name}")
print("eligible on all three:", audit["eligible_on_all_three_criteria"])
