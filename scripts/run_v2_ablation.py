"""Predeclared feature-block ablation for the v2 stack, on all frozen folds.

The purpose is attribution, not model selection: it measures how much of the
annotation-ranking signal comes from the in-silico tryptic-digest block, from
the new sequence-context block and from the within-protein contextual
transforms.  One fixed LightGBM configuration is used for every variant so that
differences come from the feature set only; boosting rounds are still selected
on the inner validation fold and never on the outer test fold.
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from common import (
    INPUTS,
    RESULTS,
    ROOT,
    load_initial_scores,
    paired_component_bootstrap,
    sha256,
    summarize,
    write_csv,
    write_json,
)
from v2_stack import FEATURES, fold_partitions, load_v2_matrix
import v2_members as members

DETECTABILITY_TOKENS = ("T:", "nearest_K", "nearest_R", "pep_", "kr_within")

VARIANTS = {
    "A_frozen23": lambda name: name.startswith("A:"),
    "T_digest_only": lambda name: "T:" in name and not name.startswith("C:"),
    "A_plus_T": lambda name: name.startswith("A:") or ("T:" in name and not name.startswith("C:")),
    "A_plus_B_chem": lambda name: (
        not name.startswith("C:")
        and not any(token in name for token in DETECTABILITY_TOKENS)
    ),
    "A_plus_B_all": lambda name: not name.startswith("C:"),
    "chem_no_detectability": lambda name: not any(
        token in name for token in DETECTABILITY_TOKENS
    ),
    "full": lambda name: True,
}

FIXED_CONFIG = members.LGB_BIN_GRID[0]


def main():
    started = time.time()
    X, names, y, folds, proteins, components, positions, meta = load_v2_matrix()
    reference_hgb = load_initial_scores("hist_gradient_boosting_current")

    rows = []
    predictions = {}
    records = []
    for variant, predicate in VARIANTS.items():
        columns = np.asarray([i for i, name in enumerate(names) if predicate(name)])
        Xv = np.ascontiguousarray(X[:, columns])
        for kind in ("bin",):
            scores = np.full(len(y), np.nan)
            fold_notes = []
            for outer_fold in range(1, 6):
                subtrain, validation, outer_train, test, validation_fold = fold_partitions(
                    folds, outer_fold, proteins, components
                )
                selection, _ = members.lgb_select(
                    kind, Xv, y, proteins, components, subtrain, validation,
                    grid=(FIXED_CONFIG,),
                )
                fold_scores, _ = members.lgb_refit(
                    kind, Xv, y, proteins, outer_train, test, selection
                )
                scores[test] = fold_scores
                fold_notes.append({
                    "outer_fold": outer_fold,
                    "inner_validation_fold": validation_fold,
                    "rounds": selection["rounds"],
                    "validation_objective": selection["validation_objective"],
                })
            if not np.isfinite(scores).all():
                raise RuntimeError(f"Incomplete scores for {variant}/{kind}")
            label = f"{variant}|lgb_{kind}"
            predictions[label] = scores
            summary = summarize(y, scores, proteins, components)
            summary.update({
                "variant": variant,
                "objective": f"lgb_{kind}",
                "n_features": int(len(columns)),
                "digest_block_present": bool(
                    any("T:" in names[i] for i in columns)
                ),
                "within_protein_context_present": bool(
                    any(names[i].startswith("C:") for i in columns)
                ),
            })
            rows.append(summary)
            records.append({"variant": variant, "objective": f"lgb_{kind}", "folds": fold_notes})
            print(json.dumps({
                "variant": variant, "objective": f"lgb_{kind}",
                "n_features": int(len(columns)),
                "roc_auc": round(summary["roc_auc"], 4),
                "average_precision": round(summary["average_precision"], 4),
                "top1": round(summary["top1"], 4),
                "elapsed": round(time.time() - started, 1),
            }), flush=True)

    write_csv(RESULTS / "v2_ablation_summary.csv", rows)
    write_json(RESULTS / "v2_ablation_fold_records.json", records)
    np.savez_compressed(RESULTS / "v2_ablation_oof.npz", **predictions)

    intervals = {}
    pairs = [
        ("A_plus_T|lgb_bin", "A_frozen23|lgb_bin", "digest_block_added_to_frozen_features"),
        ("A_plus_B_chem|lgb_bin", "A_frozen23|lgb_bin", "chemistry_sequence_block_added_to_frozen_features"),
        ("full|lgb_bin", "chem_no_detectability|lgb_bin", "detectability_block_added_to_chemistry_model"),
        ("chem_no_detectability|lgb_bin", "A_plus_B_chem|lgb_bin", "within_protein_context_added_to_chemistry_model"),
    ]
    for candidate, reference, label in pairs:
        intervals[label] = paired_component_bootstrap(
            y, predictions[reference], predictions[candidate], proteins, components
        )
    intervals["chemistry_model_vs_current_hist_gradient_boosting"] = paired_component_bootstrap(
        y, reference_hgb, predictions["chem_no_detectability|lgb_bin"], proteins, components
    )
    write_json(RESULTS / "v2_ablation_paired_intervals.json", intervals)
    write_json(RESULTS / "v2_ablation_run_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "purpose": "attribution of the v2 gain to feature blocks; not model selection",
        "objective_used": "LightGBM binary only, to keep the ablation cost bounded; the stack run reports the ranking objective",
        "fixed_lightgbm_config": dict(FIXED_CONFIG),
        "variants": {name: sorted({
            n.split(":")[0] for i, n in enumerate(names) if predicate(n)
        }) for name, predicate in VARIANTS.items()},
        "detectability_tokens": list(DETECTABILITY_TOKENS),
        "split_discipline": "frozen outer folds; rounds selected on the cyclic inner fold only",
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "lightgbm": members.lgb.__version__,
        },
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "frozen_benchmark_data.npz",
                FEATURES / "v2_features.npz",
            ]
        },
    })
    print(json.dumps(rows, indent=2, default=float)[:400])


if __name__ == "__main__":
    main()
