"""One leakage-aware table holding every valid v1 and v2 OOF prediction.

The v1 model scores are read from their stored artefacts and are not refitted,
so the published v1 numbers reproduce exactly.  Stratification, weighted
bootstrap and metric definitions are imported from the v1 unified evaluation.
"""
from __future__ import annotations

import json
import platform
import sys
import time

import numpy as np

from common import (
    INPUTS,
    RESULTS,
    ROOT,
    load_frozen,
    load_initial_scores,
    read_csv,
    sha256,
    summarize,
    write_csv,
    write_json,
)
import run_unified_evaluation as v1
from v2_stack import FEATURES

V2_TAGS = ("v2_stack_chem", "v2_stack_full")


def load_v2(tag):
    archive = np.load(RESULTS / f"{tag}_oof.npz", allow_pickle=False)
    return {key: archive[key] for key in archive.files}


def main():
    started = time.time()
    X, y, folds, proteins, components, positions, meta = load_frozen()
    del X
    models = {
        "histgb_23_feature_current": load_initial_scores("hist_gradient_boosting_current"),
        "mlp_corrected_five_seed_mean": v1.load_corrected_ensemble(),
        "histgb_sequence_only": v1.load_score_column(INPUTS / "oof_sequence_only.csv"),
        "histgb_structure_only": v1.load_score_column(INPUTS / "oof_structure_only.csv"),
        "tabicl_v2_cuda": v1.load_complete_npz_scores(
            ROOT / "tabular" / "tabicl_cuda_oof_complete.npz",
            y=y, proteins=proteins, positions=positions,
        ),
        "protbert_frozen_linear_probe": v1.load_complete_npz_scores(
            ROOT / "plm" / "protbert_linear_probe_oof_complete.npz"
        ),
    }
    for tag in V2_TAGS:
        for member, scores in load_v2(tag).items():
            models[f"{tag}:{member}"] = scores
    for model, scores in models.items():
        if len(scores) != len(y) or not np.isfinite(scores).all():
            raise ValueError(f"Incomplete OOF predictions for {model}")

    metrics = []
    for model, scores in models.items():
        row = summarize(y, scores, proteins, components)
        row.update({
            "model": model,
            "label_scope": "method-dependent annotation versus unlabelled",
            "outer_split": "frozen five-fold homology components",
            "generation": "v2" if model.startswith("v2_") else "v1",
        })
        metrics.append(row)
    metrics.sort(key=lambda row: -row["within_protein_auc_mean"])
    write_csv(RESULTS / "v2_unified_model_metrics.csv", metrics)

    strata = v1.build_strata(meta, proteins)
    stratified = v1.stratified_auc_ap(y, models, strata, folds)
    write_csv(RESULTS / "v2_unified_stratified_auc_ap.csv", stratified)

    comparisons = {}
    pairs = [
        ("v2_stack_full:blend", "mlp_corrected_five_seed_mean", 20260913),
        ("v2_stack_chem:blend", "mlp_corrected_five_seed_mean", 20260914),
        ("v2_stack_full:blend", "v2_stack_chem:blend", 20260915),
        ("v2_stack_chem:blend", "histgb_23_feature_current", 20260916),
        ("v2_stack_full:blend", "tabicl_v2_cuda", 20260917),
        ("v2_stack_chem:blend", "protbert_frozen_linear_probe", 20260918),
        ("v2_stack_full:blend", "v2_stack_full:lgb_bin", 20260919),
        ("v2_stack_chem:blend", "v2_stack_chem:lgb_bin", 20260920),
        ("v2_stack_full:blend", "v2_stack_full:fusion_mlp", 20260921),
    ]
    for candidate, reference, seed in pairs:
        label = f"{candidate}__vs__{reference}"
        print("bootstrap:", label, flush=True)
        comparisons[label] = v1.paired_component_bootstrap_all(
            y, models[reference], models[candidate], proteins, components, seed=seed
        )
    write_json(RESULTS / "v2_unified_paired_component_intervals.json", comparisons)

    write_json(RESULTS / "v2_unified_evaluation_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "n_models": len(models),
        "models": sorted(models),
        "v1_scores_refitted": False,
        "metric_definitions": "imported from scripts/common.py and scripts/run_unified_evaluation.py",
        "label_semantics": "method-dependent annotation versus unlabelled Cys; not chemical positive versus confirmed negative",
        "hardware": {"platform": platform.platform()},
        "versions": {"python": sys.version, "numpy": np.__version__},
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "frozen_benchmark_data.npz",
                FEATURES / "v2_features.npz",
                RESULTS / "v2_stack_full_oof.npz",
                RESULTS / "v2_stack_chem_oof.npz",
                RESULTS / "mlp_corrected_oof_by_seed.csv",
            ]
        },
    })
    print(json.dumps([
        {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()
         if k in ("model", "roc_auc", "average_precision", "within_protein_auc_mean",
                  "top1", "top2", "mean_reciprocal_rank")}
        for row in metrics
    ], indent=1))


if __name__ == "__main__":
    main()
