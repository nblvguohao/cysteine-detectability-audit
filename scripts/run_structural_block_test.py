"""Does the new structural chemistry block help on the annotation benchmark?

53 three-dimensional descriptors are now available for 39,240 of the 49,064
cohort cysteines (`features/structural_chemistry_features.npz`).  This program
tests them on the frozen annotation benchmark with the project's own discipline:
frozen five-fold homology components, boosting rounds selected on the cyclic
inner fold only, one fixed LightGBM configuration for every design so that
differences come from the feature set.

The expectation from everything measured so far is a small effect, because the
annotation label is explained by peptide visibility rather than chemistry.  The
test is run anyway so the claim is measured rather than assumed, and because the
structure-only design is the cleanest available answer to "is there any
chemistry signal in this label at all".
"""
from __future__ import annotations

import json
import platform
import sys
import time

import numpy as np

from common import (
    INPUTS, RESULTS, ROOT, load_initial_scores, paired_component_bootstrap,
    sha256, summarize, write_csv, write_json,
)
from v2_stack import FEATURES, fold_partitions, load_v2_matrix
import v2_members as members

STRUCT = FEATURES / "structural_chemistry_features.npz"
FIXED_CONFIG = members.LGB_BIN_GRID[0]
DETECTABILITY_TOKENS = ("T:", "nearest_K", "nearest_R", "pep_", "kr_within")


def main():
    started = time.time()
    X, names, y, folds, proteins, components, positions, meta = load_v2_matrix()
    archive = np.load(STRUCT, allow_pickle=False)
    assert np.array_equal(archive["proteins"].astype(str), proteins)
    assert np.array_equal(archive["positions"].astype(int), positions)
    S = archive["X"]
    S_names = archive["names"].astype(str)
    covered = archive["covered"]
    print(json.dumps({"structural_features": int(S.shape[1]),
                      "sites_with_structure": int(covered.sum())}), flush=True)

    chem_columns = np.asarray([
        i for i, name in enumerate(names)
        if not any(token in name for token in DETECTABILITY_TOKENS)
    ])
    designs = {
        "v2_full": np.ascontiguousarray(X),
        "v2_full_plus_structure53": np.hstack([X, S]),
        "v2_chem": np.ascontiguousarray(X[:, chem_columns]),
        "v2_chem_plus_structure53": np.hstack([X[:, chem_columns], S]),
        "structure53_only": np.ascontiguousarray(S),
    }

    rows = []
    predictions = {}
    for label, matrix in designs.items():
        scores = np.full(len(y), np.nan)
        selections = []
        for outer_fold in range(1, 6):
            subtrain, validation, outer_train, test, validation_fold = fold_partitions(
                folds, outer_fold, proteins, components
            )
            selection, _ = members.lgb_select(
                "bin", matrix, y, proteins, components, subtrain, validation,
                grid=(FIXED_CONFIG,),
            )
            fold_scores, _ = members.lgb_refit(
                "bin", matrix, y, proteins, outer_train, test, selection
            )
            scores[test] = fold_scores
            selections.append(selection["rounds"])
        if not np.isfinite(scores).all():
            raise RuntimeError(f"Incomplete scores for {label}")
        predictions[label] = scores
        summary = summarize(y, scores, proteins, components)
        summary.update({
            "design": label,
            "n_features": int(matrix.shape[1]),
            "rounds_per_fold": ",".join(map(str, selections)),
        })
        rows.append(summary)
        print(f"{label:28s} n={matrix.shape[1]:5d} AUC {summary['roc_auc']:.4f} "
              f"within {summary['within_protein_auc_mean']:.4f} "
              f"Top-1 {summary['top1']:.4f}  {round(time.time()-started)}s", flush=True)

        # structure-covered subset, where the block is actually informative
        if label in ("v2_full", "v2_full_plus_structure53", "structure53_only"):
            sub = summarize(y[covered], scores[covered], proteins[covered], components[covered])
            sub.update({"design": f"{label}|structure_covered_subset",
                        "n_features": int(matrix.shape[1]), "rounds_per_fold": ""})
            rows.append(sub)

    intervals = {}
    for candidate, reference in (
        ("v2_full_plus_structure53", "v2_full"),
        ("v2_chem_plus_structure53", "v2_chem"),
    ):
        intervals[f"{candidate}__minus__{reference}"] = paired_component_bootstrap(
            y, predictions[reference], predictions[candidate], proteins, components
        )
    intervals["structure53_only__minus__current_histgb"] = paired_component_bootstrap(
        y, load_initial_scores("hist_gradient_boosting_current"),
        predictions["structure53_only"], proteins, components
    )

    write_csv(RESULTS / "structural_block_test.csv", rows)
    np.savez_compressed(RESULTS / "structural_block_test_oof.npz", **predictions)
    write_json(RESULTS / "structural_block_test_intervals.json", intervals)
    write_json(RESULTS / "structural_block_test_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": "does a three-dimensional chemistry block help on the annotation benchmark",
        "structural_features": int(S.shape[1]),
        "structural_feature_names": list(S_names),
        "sites_with_structure": int(covered.sum()),
        "structure_coverage": round(float(covered.mean()), 4),
        "uncovered_rows_are_nan": "LightGBM handles missing values natively; no imputation is applied",
        "splits": "frozen five-fold homology components; rounds selected on the cyclic inner fold",
        "fixed_lightgbm_config": dict(FIXED_CONFIG),
        "designs": {k: int(v.shape[1]) for k, v in designs.items()},
        "limits": [
            "AlphaFold single conformer, no dynamics",
            "no explicit pKa; electrostatics are geometric proxies",
            "20% of sites have no model, and those rows stay missing rather than imputed",
        ],
        "versions": {"python": sys.version, "numpy": np.__version__,
                     "lightgbm": members.lgb.__version__},
        "input_hashes": {
            "features/v2_features.npz": sha256(FEATURES / "v2_features.npz"),
            "features/structural_chemistry_features.npz": sha256(STRUCT),
        },
    })
    print(json.dumps(intervals, indent=2)[:900])


if __name__ == "__main__":
    main()
