"""Run the v2 nested stack on the homology folds and audit the result."""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time

import numpy as np
import sklearn

from common import (
    INPUTS,
    RESULTS,
    ROOT,
    load_initial_scores,
    paired_component_bootstrap,
    sha256,
    summarize,
    validation_objective,
    write_csv,
    write_json,
)
from v2_stack import (
    FEATURES,
    PLM,
    apply_platt,
    fold_partitions,
    load_embeddings,
    load_v2_matrix,
    platt,
    simplex_grid,
)
import v2_members as members

SEEDS = [11, 29, 47, 71, 97]
MEMBERS = ("lgb_rank", "lgb_bin", "fusion_mlp")

# Predeclared feature sets.  ``full`` is the complete v2 matrix.  ``chem``
# removes the in-silico tryptic-digest block and the two trypsin-proximity
# sequence features, so that the part of the annotation signal attributable to
# mass-spectrometry peptide detectability is dropped rather than exploited.
# Residual detectability correlation cannot be fully removed by any such filter.
DETECTABILITY_TOKENS = ("T:", "nearest_K", "nearest_R", "pep_", "kr_within")


def select_columns(names, feature_set):
    if feature_set == "full":
        return np.arange(len(names))
    if feature_set == "chem":
        return np.asarray([
            i for i, name in enumerate(names)
            if not any(token in name for token in DETECTABILITY_TOKENS)
        ])
    raise SystemExit(f"Unknown feature set {feature_set}")


def blend_weights(y, components, proteins, validation, member_scores, step=0.1):
    """Pick non-negative simplex weights on the inner fold only."""
    stacked = np.stack([member_scores[name] for name in MEMBERS])
    best_weights = np.full(len(MEMBERS), 1.0 / len(MEMBERS))
    best_objective = -np.inf
    evaluated = 0
    for weights in simplex_grid(len(MEMBERS), step=step):
        blended = np.tensordot(weights, stacked, axes=1)
        objective = validation_objective(
            y[validation], blended, proteins[validation], components[validation]
        )
        evaluated += 1
        if objective > best_objective:
            best_objective, best_weights = objective, weights
    return best_weights, float(best_objective), evaluated


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", type=int, nargs="*", default=[1, 2, 3, 4, 5])
    parser.add_argument("--seeds", type=int, nargs="*", default=SEEDS)
    parser.add_argument("--tag", default="v2_stack")
    parser.add_argument("--feature-set", choices=("full", "chem"), default="full")
    parser.add_argument("--blend-step", type=float, default=0.1)
    args = parser.parse_args()

    started = time.time()
    X, names, y, folds, proteins, components, positions, meta = load_v2_matrix()
    columns = select_columns(names, args.feature_set)
    kept_names = [names[i] for i in columns]
    X = np.ascontiguousarray(X[:, columns])
    embeddings, centred = load_embeddings(proteins, positions)

    member_oof = {name: np.full(len(y), np.nan) for name in MEMBERS}
    blend_oof = np.full(len(y), np.nan)
    records = []

    for outer_fold in args.folds:
        fold_started = time.time()
        subtrain, validation, outer_train, test, validation_fold = fold_partitions(
            folds, outer_fold, proteins, components
        )
        validation_scores = {}
        test_scores = {}
        fold_record = {
            "outer_fold": int(outer_fold),
            "inner_validation_fold": int(validation_fold),
            "n_subtrain": int(len(subtrain)),
            "n_validation": int(len(validation)),
            "n_outer_train": int(len(outer_train)),
            "n_test": int(len(test)),
            "members": {},
        }

        for kind in ("rank", "bin"):
            name = f"lgb_{kind}"
            selection, valid_raw = members.lgb_select(
                kind, X, y, proteins, components, subtrain, validation
            )
            test_raw, _ = members.lgb_refit(
                kind, X, y, proteins, outer_train, test, selection
            )
            calibrator = platt(valid_raw, y[validation])
            validation_scores[name] = apply_platt(calibrator, valid_raw)
            test_scores[name] = apply_platt(calibrator, test_raw)
            fold_record["members"][name] = {
                "selection": {k: v for k, v in selection.items() if k != "trace"},
                "trace_points": len(selection["trace"]),
                "platt_coefficient": calibrator["coefficient"],
                "platt_intercept": calibrator["intercept"],
                "platt_centre": calibrator["centre"],
                "platt_spread": calibrator["spread"],
            }
            print(json.dumps({
                "fold": outer_fold, "member": name,
                "rounds": selection["rounds"],
                "config_index": selection["config_index"],
                "validation_objective": round(selection["validation_objective"], 5),
                "seconds": round(time.time() - fold_started, 1),
            }), flush=True)

        fusion_records, fusion_validation, fusion_test = members.fusion_select_and_refit(
            X, embeddings, centred, y, proteins, components,
            subtrain, validation, outer_train, test, args.seeds,
        )
        validation_scores["fusion_mlp"] = fusion_validation
        test_scores["fusion_mlp"] = fusion_test
        fold_record["members"]["fusion_mlp"] = {
            "seeds": args.seeds,
            "selected_epochs": [record["epochs"] for record in fusion_records],
            "validation_objectives": [record["validation_objective"] for record in fusion_records],
            "device": fusion_records[0]["device"],
            "hyperparameters": fusion_records[0]["hyperparameters"],
        }
        print(json.dumps({
            "fold": outer_fold, "member": "fusion_mlp",
            "epochs": [record["epochs"] for record in fusion_records],
            "seconds": round(time.time() - fold_started, 1),
        }), flush=True)

        weights, objective, evaluated = blend_weights(
            y, components, proteins, validation, validation_scores, step=args.blend_step
        )
        blended_test = np.tensordot(
            weights, np.stack([test_scores[name] for name in MEMBERS]), axes=1
        )
        for name in MEMBERS:
            member_oof[name][test] = test_scores[name]
        blend_oof[test] = blended_test
        fold_record["blend"] = {
            "members": list(MEMBERS),
            "weights": [float(w) for w in weights],
            "grid_step": args.blend_step,
            "grid_points_evaluated": evaluated,
            "inner_validation_objective": objective,
            "weights_estimated_on": "inner validation fold stage-1 out-of-sample predictions",
        }
        fold_record["elapsed_seconds"] = time.time() - fold_started
        records.append(fold_record)
        print(json.dumps({
            "fold": outer_fold, "weights": [round(float(w), 2) for w in weights],
            "inner_objective": round(objective, 5),
            "fold_seconds": round(fold_record["elapsed_seconds"], 1),
        }), flush=True)

    complete = all(np.isfinite(blend_oof[np.isin(folds, args.folds)]).all() for _ in [0])
    if not complete:
        raise RuntimeError("Incomplete blended predictions")

    evaluated_mask = np.isin(folds, args.folds)
    summaries = []
    for name in MEMBERS:
        row = summarize(
            y[evaluated_mask], member_oof[name][evaluated_mask],
            proteins[evaluated_mask], components[evaluated_mask],
        )
        row.update({"model": f"{args.tag}:{name}", "folds": ",".join(map(str, args.folds))})
        summaries.append(row)
    blend_summary = summarize(
        y[evaluated_mask], blend_oof[evaluated_mask],
        proteins[evaluated_mask], components[evaluated_mask],
    )
    blend_summary.update({"model": f"{args.tag}:blend", "folds": ",".join(map(str, args.folds))})
    summaries.append(blend_summary)

    reference_hgb = load_initial_scores("hist_gradient_boosting_current")
    rows = []
    import csv as _csv
    with (RESULTS / "mlp_corrected_oof_by_seed.csv").open(encoding="utf-8-sig", newline="") as handle:
        for row in _csv.DictReader(handle):
            rows.append((int(row["seed"]), row["accession"], int(row["position"]), float(row["score"])))
    by_seed = {}
    for seed, accession, position, score in rows:
        by_seed.setdefault(seed, []).append(score)
    reference_mlp = np.mean(np.stack([by_seed[seed] for seed in sorted(by_seed)]), axis=0)

    comparisons = {}
    if evaluated_mask.all():
        for label, reference in (
            ("v2_blend_vs_corrected_mlp", reference_mlp),
            ("v2_blend_vs_current_hist_gradient_boosting", reference_hgb),
        ):
            comparisons[label] = paired_component_bootstrap(
                y, reference, blend_oof, proteins, components
            )
        for name in MEMBERS:
            comparisons[f"v2_{name}_vs_corrected_mlp"] = paired_component_bootstrap(
                y, reference_mlp, member_oof[name], proteins, components
            )

    prediction_rows = []
    for i in np.flatnonzero(evaluated_mask):
        record = {
            "accession": proteins[i],
            "position": int(positions[i]),
            "component": components[i],
            "fold": int(folds[i]),
            "label": int(y[i]),
            "label_scope": "method_dependent_annotation_vs_unlabelled",
            "blend": float(blend_oof[i]),
        }
        for name in MEMBERS:
            record[name] = float(member_oof[name][i])
        prediction_rows.append(record)

    write_csv(RESULTS / f"{args.tag}_oof.csv", prediction_rows)
    write_csv(RESULTS / f"{args.tag}_summary.csv", summaries)
    write_json(RESULTS / f"{args.tag}_fold_records.json", records)
    if comparisons:
        write_json(RESULTS / f"{args.tag}_paired_intervals.json", comparisons)
    np.savez_compressed(
        RESULTS / f"{args.tag}_oof.npz",
        blend=blend_oof,
        **{name: member_oof[name] for name in MEMBERS},
    )
    write_json(RESULTS / f"{args.tag}_run_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "folds_evaluated": list(args.folds),
        "members": list(MEMBERS),
        "feature_set": args.feature_set,
        "feature_set_definition": {
            "full": "every v2 column",
            "chem": "v2 columns whose name contains none of " + ", ".join(DETECTABILITY_TOKENS),
        }[args.feature_set],
        "n_features_used": int(X.shape[1]),
        "n_features_dropped": int(len(names) - X.shape[1]),
        "feature_blocks": {
            "frozen_23": True,
            "new_sequence_context": True,
            "within_protein_rank_and_centring": True,
            "protbert_frozen_embedding": ["lgb: not used", "fusion_mlp: raw and protein-centred"],
            "n_tabular_features": int(X.shape[1]),
            "blocks_present": sorted({name.split(":")[0] for name in kept_names}),
            "digest_features_present": any("T:" in name for name in kept_names),
        },
        "split_discipline": {
            "outer": "frozen five-fold homology components",
            "inner": "cyclic remaining frozen component fold (k % 5 + 1)",
            "test_fold_used_for_preprocessing": False,
            "test_fold_used_for_hyperparameters": False,
            "test_fold_used_for_epochs": False,
            "test_fold_used_for_calibration": False,
            "test_fold_used_for_blend_weights": False,
            "calibration_and_blend_source": "stage-1 out-of-sample predictions on the inner validation fold",
        },
        "seeds": args.seeds,
        "lgb_rank_grid": [dict(config) for config in members.LGB_RANK_GRID],
        "lgb_bin_grid": [dict(config) for config in members.LGB_BIN_GRID],
        "lgb_max_rounds": members.MAX_ROUNDS,
        "lgb_round_block": members.BLOCK,
        "lgb_patience_blocks": members.PATIENCE_BLOCKS,
        "validation_objective": "unweighted mean(within-protein AUC, Top-1, Top-2, MRR)",
        "label_semantics": "method-dependent annotation versus unlabelled Cys",
        "hardware": {"platform": platform.platform()},
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "sklearn": sklearn.__version__,
            "lightgbm": members.lgb.__version__,
            "torch": __import__("torch").__version__,
        },
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "training_matrix.npz",
                INPUTS / "training_site_folds.csv",
                FEATURES / "v2_features.npz",
                PLM / "protbert_cys_embeddings.npz",
            ]
        },
        "summaries": summaries,
    })
    print(json.dumps(summaries, indent=2, default=float))


if __name__ == "__main__":
    main()
