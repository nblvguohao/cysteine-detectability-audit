"""Is the deployed chemistry model's cleavage-flank ranking produced by its basic-composition features?

`reports/SELF_AUDIT_CLEAVAGE_FLANK.md` found the wet-lab chemistry blend ranks cysteines by trypsin's
cleavage structure beyond what the label carries (within-label-stratum score-to-K/R AUC 0.6415 among
positives; top 100 98.0% K/R-flanked, log2 OR 3.4716). `reports/SELF_AUDIT_FIVE_PROTEASES.md` then
narrowed that to the lysine-cleaving proteases and recorded the limitation this script exists to close:
the chemistry feature set carries basic-residue composition at windows 1, 2, 3, 5, 7, 10, 15, 20, 30
and 50 residues, so the 6-12 association is fully reachable through a DECLARED chemistry feature, and
"the model reconstructed peptide-length logic" was not separable from "the model uses mid-range basic
composition, which for tryptic data IS a peptide-length proxy". Feature naming cannot separate them.
An ablation can.

Note what the chem feature set ALREADY excludes (`run_v2_stack.DETECTABILITY_TOKENS`): the whole `T:`
block, `nearest_K`, `nearest_R`, `pep_*` and `kr_within`. So the residual association cannot come from
an explicit K/R-distance or peptide feature - only from the composition features that remain, which is
why those are what this ablates.

**Three arms, predeclared. All refit from scratch with the same folds, the same member selection grids
and the same seeds; only the column set differs.**

  * `control`      - the chem set exactly as deployed (1046 columns). This arm is BOTH the ablation's
                     control AND the bridge to the stored 2026-09-13 out-of-fold scores, because this
                     round runs on a different machine and a different numpy/lightgbm build.
  * `drop_w710`    - chem minus the four win7/win10 basic-composition columns (frac and count).
  * `drop_all_basic` - chem minus all twenty basic-composition columns at every window. The maximal
                     version, to bound the effect: if `drop_w710` leaves the enrichment standing but
                     `drop_all_basic` removes it, the signal is spread across windows rather than
                     concentrated in the mid-range ones.

**Reading rule, fixed before the run.** The verdict statistic is the top-100 distal (6-12) K/R log2
odds ratio, computed AFTERWARDS by the already-registered local instrument
(`scripts/run_self_audit_five_proteases.py`) on the out-of-fold scores this job returns - not here, so
the audit stays byte-identical to the one that produced the published numbers.

  1. `control` reproduces the stored enrichment (its interval covers 3.7116) AND the ablated arms'
     enrichment intervals cross zero -> the ranking IS produced by the basic-composition features.
  2. `control` reproduces AND the ablated arms keep their enrichment -> the ranking survives removal;
     it comes from elsewhere in the feature set, and the declared-chemistry-feature explanation is
     refuted.
  3. `control` does NOT reproduce -> this machine is not comparable to the 2026-09-13 run; the ablation
     is read arm-to-arm only and no comparison with any published number may be made.

**And the question that decides whether this is actionable**: each arm's task performance is reported
(within-protein AUC, the project's headline metric). If the ablated arms lose the cleavage-flank
ranking at no cost in within-protein AUC, then dropping those columns is a defensible model
improvement rather than only a diagnosis. If it costs performance, the trade-off is the finding.

**Declared deviation from `run_v2_stack.main`.** The fold loop here calls exactly the same functions in
the same order with the same seeds (`fold_partitions`, `members.lgb_select`, `members.lgb_refit`,
`platt`/`apply_platt`, `members.fusion_select_and_refit`, `run_v2_stack.blend_weights`,
`common.summarize`). What is omitted is that script's reporting tail, which reads
`results/mlp_corrected_oof_by_seed.csv` and the v1 initial scores to build a cross-generation
comparison table; those are 26 MB of files unrelated to this ablation and are not uploaded. No
modelling step is skipped or altered.

Thread count is pinned by the caller (OMP_NUM_THREADS) and recorded, because LightGBM's histogram
construction is thread-count sensitive and all three arms must be mutually comparable.

Writes, per arm: results/v2_chem_ablation_<arm>_oof.csv and results/v2_chem_ablation_<arm>_folds.json;
plus one results/v2_chem_ablation_audit.json over all arms.
"""
from __future__ import annotations

import json
import os
import platform
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import v2_members as members
import run_v2_stack as base
from common import RESULTS, sha256, summarize
from pathlib import Path as _Path
from v2_stack import (
    apply_platt, fold_partitions, load_embeddings, load_v2_matrix, platt,
)

BASIC_RE = re.compile(r"^B:win(\d+)_(frac|count)_basic$")
ARMS = {
    "control": set(),
    "drop_w710": {"B:win7_frac_basic", "B:win7_count_basic",
                  "B:win10_frac_basic", "B:win10_count_basic"},
    "drop_all_basic": None,          # resolved against the real names below
}
SCRIPT = os.path.abspath(__file__)


def arm_columns(names, arm):
    chem = base.select_columns(names, "chem")
    if arm == "control":
        drop = set()
    elif arm == "drop_all_basic":
        drop = {str(names[i]) for i in chem if BASIC_RE.match(str(names[i]))}
    else:
        drop = ARMS[arm]
    keep = np.asarray([i for i in chem if str(names[i]) not in drop], dtype=int)
    return keep, sorted(drop)


def run_arm(arm, X_all, names, y, folds, proteins, components, positions, embeddings, centred,
            seeds, fold_list):
    started = time.time()
    columns, dropped = arm_columns(names, arm)
    X = np.ascontiguousarray(X_all[:, columns])
    print(json.dumps({"arm": arm, "n_columns": int(X.shape[1]),
                      "n_dropped": len(dropped), "dropped": dropped}), flush=True)

    member_oof = {name: np.full(len(y), np.nan) for name in base.MEMBERS}
    blend_oof = np.full(len(y), np.nan)
    records = []
    for outer_fold in fold_list:
        fold_started = time.time()
        subtrain, validation, outer_train, test, validation_fold = fold_partitions(
            folds, outer_fold, proteins, components)
        validation_scores, test_scores = {}, {}
        record = {"outer_fold": int(outer_fold), "inner_validation_fold": int(validation_fold),
                  "n_subtrain": int(len(subtrain)), "n_validation": int(len(validation)),
                  "n_outer_train": int(len(outer_train)), "n_test": int(len(test)),
                  "members": {}}
        for kind in ("rank", "bin"):
            name = f"lgb_{kind}"
            selection, valid_raw = members.lgb_select(
                kind, X, y, proteins, components, subtrain, validation)
            test_raw, _ = members.lgb_refit(kind, X, y, proteins, outer_train, test, selection)
            calibrator = platt(valid_raw, y[validation])
            validation_scores[name] = apply_platt(calibrator, valid_raw)
            test_scores[name] = apply_platt(calibrator, test_raw)
            record["members"][name] = {
                "selection": {k: v for k, v in selection.items() if k != "trace"},
                "trace_points": len(selection["trace"])}
            print(json.dumps({"arm": arm, "fold": int(outer_fold), "member": name,
                              "rounds": selection["rounds"],
                              "config_index": selection["config_index"],
                              "validation_objective": round(selection["validation_objective"], 5),
                              "seconds": round(time.time() - fold_started, 1)}), flush=True)

        fusion_records, fusion_validation, fusion_test = members.fusion_select_and_refit(
            X, embeddings, centred, y, proteins, components,
            subtrain, validation, outer_train, test, seeds)
        validation_scores["fusion_mlp"] = fusion_validation
        test_scores["fusion_mlp"] = fusion_test
        record["members"]["fusion_mlp"] = {
            "seeds": seeds,
            "selected_epochs": [r["epochs"] for r in fusion_records],
            "validation_objectives": [r["validation_objective"] for r in fusion_records],
            "device": fusion_records[0]["device"],
            "hyperparameters": fusion_records[0]["hyperparameters"]}

        weights, objective, evaluated = base.blend_weights(
            y, components, proteins, validation, validation_scores, step=0.1)
        blended = np.tensordot(weights, np.stack([test_scores[n] for n in base.MEMBERS]), axes=1)
        for name in base.MEMBERS:
            member_oof[name][test] = test_scores[name]
        blend_oof[test] = blended
        record["blend"] = {"members": list(base.MEMBERS),
                           "weights": [float(w) for w in weights],
                           "grid_step": 0.1, "grid_points_evaluated": evaluated,
                           "inner_validation_objective": objective}
        record["elapsed_seconds"] = time.time() - fold_started
        records.append(record)
        print(json.dumps({"arm": arm, "fold": int(outer_fold),
                          "weights": [round(float(w), 2) for w in weights],
                          "inner_objective": round(objective, 5),
                          "fold_seconds": round(record["elapsed_seconds"], 1)}), flush=True)

    mask = np.isin(folds, fold_list)
    if not np.isfinite(blend_oof[mask]).all():
        raise RuntimeError(f"arm {arm}: incomplete blended predictions")

    summaries = []
    for name in list(base.MEMBERS) + ["blend"]:
        scores = blend_oof if name == "blend" else member_oof[name]
        row = summarize(y[mask], scores[mask], proteins[mask], components[mask])
        row.update({"arm": arm, "model": name})
        summaries.append(row)

    out = os.path.join(RESULTS, f"v2_chem_ablation_{arm}_oof.csv")
    with open(out, "w", encoding="utf-8-sig", newline="") as fh:
        fh.write("accession,position,component,fold,label,blend,lgb_rank,lgb_bin,fusion_mlp\n")
        for i in range(len(y)):
            fh.write(f"{proteins[i]},{int(positions[i])},{components[i]},{int(folds[i])},"
                     f"{int(y[i])},{blend_oof[i]!r},{member_oof['lgb_rank'][i]!r},"
                     f"{member_oof['lgb_bin'][i]!r},{member_oof['fusion_mlp'][i]!r}\n")
    with open(os.path.join(RESULTS, f"v2_chem_ablation_{arm}_folds.json"), "w") as fh:
        json.dump({"arm": arm, "n_columns": int(X.shape[1]), "dropped_columns": dropped,
                   "folds": records, "summaries": summaries,
                   "elapsed_minutes": round((time.time() - started) / 60, 2)},
                  fh, ensure_ascii=False, indent=2)
    print(json.dumps({"arm": arm, "done_minutes": round((time.time() - started) / 60, 2),
                      "summaries": summaries}, ensure_ascii=False), flush=True)
    return {"arm": arm, "n_columns": int(X.shape[1]), "dropped_columns": dropped,
            "summaries": summaries, "elapsed_minutes": round((time.time() - started) / 60, 2)}


def main():
    started = time.time()
    arms = sys.argv[1:] or list(ARMS)
    fold_list = [1, 2, 3, 4, 5]
    X_all, names, y, folds, proteins, components, positions, meta = load_v2_matrix()
    embeddings, centred = load_embeddings(proteins, positions)
    print(json.dumps({"sites": int(len(y)), "positives": int(y.sum()),
                      "all_columns": int(X_all.shape[1]),
                      "chem_columns": int(len(base.select_columns(names, "chem"))),
                      "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
                      "arms": arms}), flush=True)
    results = [run_arm(arm, X_all, names, y, folds, proteins, components, positions,
                       embeddings, centred, base.SEEDS, fold_list) for arm in arms]
    audit = {
        "script": "scripts/run_v2_chem_basic_ablation.py", "script_sha256": sha256(_Path(SCRIPT)),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.node(),
        "round": ("ablation: does removing basic-residue composition remove the deployed chemistry "
                  "model's cleavage-flank ranking, and at what cost in task performance"),
        "arms": results,
        "seeds": base.SEEDS, "members": list(base.MEMBERS), "folds": fold_list,
        "detectability_tokens_already_excluded_by_chem": list(base.DETECTABILITY_TOKENS),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "verdict_computed_where": ("not here - the top-100 distal K/R log2 odds ratio is computed "
                                   "locally by scripts/run_self_audit_five_proteases.py on these "
                                   "out-of-fold scores, so the audit instrument stays unchanged"),
        "versions": {"python": sys.version, "numpy": np.__version__,
                     "lightgbm": members.lgb.__version__},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    try:
        import torch
        audit["versions"]["torch"] = torch.__version__
    except Exception as exc:
        audit["versions"]["torch"] = f"unavailable: {type(exc).__name__}"
    with open(os.path.join(RESULTS, "v2_chem_ablation_audit.json"), "w") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    print(json.dumps({"all_arms_done_minutes": audit["elapsed_minutes"]}), flush=True)


if __name__ == "__main__":
    main()
