"""The registered leakage control is degenerate on this cohort. This adds the one that works.

`protocols/v3_qtrp_analysis_plan.json` requires a leakage control that permutes
labels WITHIN homology component and expects the depth-stratified AUC to return to
about 0.5.  That control was written by analogy with v2, whose cohort had 49,064
sites across 4,986 proteins - roughly ten sites per group, so a within-group
permutation genuinely scrambles the labels.

**On the human QTRP cohort it does not.** The components are tiny: of 1,232
components covering 2,049 candidate sites, 798 hold a single site and only 122 carry
both classes.  A within-component permutation therefore moves 198 labels out of
2,049 (9.7%) and only 99 of 324 positives; nine tenths of the label vector is left
exactly as it was.  The run in `logs/v3_qtrp_stack.log` shows the consequence: the
within-component permuted fit reaches an inner AUC around 0.67 rather than 0.5, which
is what a 90%-preserved label vector should give and says nothing about leakage.

This script therefore runs the control that answers the question the plan was asking -
does the pipeline, including its nested selection, manufacture signal from noise? -
by permuting labels GLOBALLY and running the identical pipeline on chem1046.  A global
permutation destroys all label structure, so a pipeline free of leakage must come back
to about 0.5 in both the plain and the depth-stratified statistic.

This is a declared ADDITION to the registered plan, not a replacement: the
within-component control still runs inside `scripts/run_v3_qtrp_stack.py` and its
number is reported next to this one, together with the degeneracy diagnostic below, so
a reader can see why the two differ.  Nothing else about the plan changes, and no
deployment criterion depends on this script.

Predeclared before the run: seed 20260915, the same component-grouped outer and inner
folds, the same configuration grid, the same reading rules.  The expectation is
recorded here in advance - plain and stratified AUC intervals containing 0.5 - so that
a failure would be visible as a failure.
"""
from __future__ import annotations

import collections
import json
import platform
import sys
import time

import numpy as np

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import v2_apply
from run_qtrp_persulfidation_test import SITES, load_proteome, probe_matrix, structural_matrix
from run_qtrp_trained_model import human_matrix
from run_qtrp_depth_conditioned import fast_auc, stratified_auc, verdict
from run_v3_qtrp_stack import (
    COMPONENTS, PLAN, SEED, pooled_labels, select_and_score, cluster_bootstrap,
)
from v2_stack import FEATURES

SCRIPT = ROOT / "scripts" / "run_v3_permutation_control.py"
COMPILED = RESULTS / "frbm_compiled_sites.csv"
PRIMARY = "v3_pooled_ph5_raw"


def main():
    started = time.time()
    def log(message):
        print(message, flush=True)

    component_of = {r["accession"]: r["component"] for r in read_csv(COMPONENTS)}
    sequences = load_proteome()
    compiled_ssh = {(r["accession"], int(r["site"]))
                    for r in read_csv(COMPILED) if r["compiled_ssh"] == "1"}
    cohorts, depth = pooled_labels(read_csv(SITES), compiled_ssh)
    positives, negatives = cohorts[PRIMARY]

    union_keys = sorted(positives | negatives)
    residue_ok = np.asarray([
        bool(sequences.get(a)) and 1 <= p <= len(sequences[a]) and sequences[a][p - 1] == "C"
        for a, p in union_keys])
    S_all, _, structure_stats = structural_matrix(union_keys, sequences)
    P_all, probe_ok = probe_matrix(union_keys, sequences)
    covered = np.isfinite(S_all).any(axis=1) & residue_ok & probe_ok
    keys = [k for k, keep in zip(union_keys, covered) if keep and k[0] in component_of]
    V_all, columns, _, _ = human_matrix(keys, sequences)

    y = np.asarray([1 if k in positives else 0 for k in keys])
    comps = np.asarray([component_of[a] for a, _ in keys])
    strata = np.asarray([depth[k] for k in keys])
    x = np.asarray(V_all[:, columns["chem"]], dtype=np.float64)

    sizes = collections.Counter(comps.tolist())
    mixed = [c for c in sizes if 0 < y[comps == c].sum() < (comps == c).sum()]
    rng = np.random.default_rng(SEED)
    y_within = y.copy()
    for component in np.unique(comps):
        mask = np.flatnonzero(comps == component)
        y_within[mask] = rng.permutation(y[mask])
    diagnostic = {
        "sites": int(len(y)), "positives": int(y.sum()), "components": len(sizes),
        "component_size_distribution": dict(sorted(collections.Counter(sizes.values()).items())),
        "components_carrying_both_classes": len(mixed),
        "sites_in_those_components": int(sum(int((comps == c).sum()) for c in mixed)),
        "labels_moved_by_within_component_permutation": int((y_within != y).sum()),
        "fraction_moved": round(float((y_within != y).mean()), 4),
        "positives_moved": int(((y == 1) & (y_within != y)).sum()),
        "reading": ("a within-component permutation leaves about nine tenths of this label vector "
                    "untouched, so it cannot return the pipeline to 0.5 and is uninformative about leakage "
                    "on this cohort"),
    }
    log(json.dumps(diagnostic, indent=2))

    rng_global = np.random.default_rng(SEED)
    y_global = rng_global.permutation(y)
    log(f"global permutation moved {int((y_global != y).sum())} of {len(y)} labels")
    scores, selection = select_and_score(x, y_global, comps, "chem1046_global_permuted", log)
    boot = cluster_bootstrap(y_global, comps, strata, {"chem1046_global_permuted": scores}, ())
    entry = boot["auc"]["chem1046_global_permuted"]

    rows = [{
        "control": "global_permutation",
        "design": "chem1046",
        "cohort": PRIMARY,
        "n_sites": int(len(y)), "n_positive": int(y_global.sum()),
        "labels_moved": int((y_global != y).sum()),
        "auc_plain": round(entry["plain"], 4),
        "plain_ci_low": round(entry["plain_interval"][0], 4),
        "plain_ci_high": round(entry["plain_interval"][1], 4),
        "auc_depth_stratified": round(entry["stratified"], 4),
        "stratified_ci_low": round(entry["stratified_interval"][0], 4),
        "stratified_ci_high": round(entry["stratified_interval"][1], 4),
        "verdict_stratified": verdict(entry["stratified_interval"]),
        "passes": bool(entry["stratified_interval"][0] <= 0.5 <= entry["stratified_interval"][1]),
    }]
    write_csv(RESULTS / "v3_permutation_control_summary.csv", rows)
    write_json(RESULTS / "v3_permutation_control_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "relationship_to_the_registered_plan": (
            "a declared ADDITION, not a replacement. The plan's within-component control still runs inside "
            "scripts/run_v3_qtrp_stack.py and its number is reported alongside this one. This script exists "
            "because that control is degenerate on this cohort, which the diagnostic below quantifies."),
        "expectation_recorded_before_the_run": (
            "a pipeline free of leakage must return plain and depth-stratified AUC intervals containing 0.5 "
            "after a global label permutation"),
        "within_component_degeneracy_diagnostic": diagnostic,
        "global_permutation_result": rows[0],
        "inner_selection": selection,
        "structure_coverage": structure_stats,
        "limits": [
            "a global permutation also destroys the component structure of the signal, so it tests pipeline leakage and not group-level confounding; the within-component control was meant to test the latter and cannot, on a cohort averaging 1.7 sites per component",
            "one permutation draw, not a permutation distribution; the interval here is the component bootstrap around that single draw",
        ],
        "input_hashes": {
            "results/qtrp_sites_normalised.csv": sha256(SITES),
            "results/v3_human_homology_components.csv": sha256(COMPONENTS),
            "features/v2_features.npz": sha256(FEATURES / "v2_features.npz"),
            "protocols/v3_qtrp_analysis_plan.json": sha256(PLAN),
            "scripts/run_v3_permutation_control.py": sha256(SCRIPT),
        },
        "versions": {"python": sys.version, "numpy": np.__version__,
                     "lightgbm": v2_apply.members.lgb.__version__},
    })
    log(f"done in {round((time.time() - started) / 60, 2)} min")


if __name__ == "__main__":
    main()
