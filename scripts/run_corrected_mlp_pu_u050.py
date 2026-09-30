"""Five-seed corrected MLP with a fixed 0.50 weight for unlabelled sites."""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import scipy
import sklearn
from threadpoolctl import threadpool_limits

from common import (
    INPUTS,
    RESULTS,
    ROOT,
    load_frozen,
    load_initial_scores,
    sha256,
    summarize,
    validation_objective,
    write_csv,
    write_json,
)
from run_corrected_mlp import (
    MAX_EPOCHS,
    MIN_DELTA,
    PATIENCE,
    SEEDS,
    fit_preprocessor,
    make_mlp,
    transform,
)
from run_unified_evaluation import load_corrected_ensemble, paired_component_bootstrap_all


UNLABELLED_WEIGHT = 0.50


def weights(y):
    return np.where(y == 1, 1.0, UNLABELLED_WEIGHT)


def select_epochs(X, y, proteins, components, folds, outer_fold, seed):
    validation_fold = outer_fold % 5 + 1
    subtrain = (folds != outer_fold) & (folds != validation_fold)
    validation = folds == validation_fold
    outer_test = folds == outer_fold
    for left, right in ((subtrain, validation), (subtrain, outer_test), (validation, outer_test)):
        if set(components[left]) & set(components[right]):
            raise ValueError("Component leakage")
        if set(proteins[left]) & set(proteins[right]):
            raise ValueError("Protein leakage")
    imputer, scaler, X_train = fit_preprocessor(X[subtrain])
    X_valid = transform(imputer, scaler, X[validation])
    model = make_mlp(seed)
    best_epoch, best_objective, stale = 1, -np.inf, 0
    trace = []
    for epoch in range(1, MAX_EPOCHS + 1):
        model.partial_fit(
            X_train,
            y[subtrain],
            classes=np.asarray([0, 1]),
            sample_weight=weights(y[subtrain]),
        )
        valid_scores = model.predict_proba(X_valid)[:, 1]
        objective = validation_objective(
            y[validation], valid_scores, proteins[validation], components[validation]
        )
        trace.append({"epoch": epoch, "validation_objective": objective})
        if objective > best_objective + MIN_DELTA:
            best_epoch, best_objective, stale = epoch, objective, 0
        else:
            stale += 1
        if stale >= PATIENCE:
            break
    return best_epoch, best_objective, validation_fold, trace


def refit_predict(X, y, folds, outer_fold, seed, epochs):
    train = folds != outer_fold
    test = folds == outer_fold
    imputer, scaler, X_train = fit_preprocessor(X[train])
    X_test = transform(imputer, scaler, X[test])
    model = make_mlp(seed)
    for _ in range(epochs):
        model.partial_fit(
            X_train,
            y[train],
            classes=np.asarray([0, 1]),
            sample_weight=weights(y[train]),
        )
    return test, model.predict_proba(X_test)[:, 1]


def main():
    started = time.time()
    X, y, folds, proteins, components, positions, _ = load_frozen()
    per_seed = []
    rows = []
    summaries = []
    selections = []
    for seed in SEEDS:
        scores = np.full(len(y), np.nan)
        for outer_fold in range(1, 6):
            fold_started = time.time()
            epoch, objective, validation_fold, trace = select_epochs(
                X, y, proteins, components, folds, outer_fold, seed
            )
            test, predicted = refit_predict(X, y, folds, outer_fold, seed, epoch)
            scores[test] = predicted
            selections.append({
                "seed": seed,
                "outer_fold": outer_fold,
                "inner_validation_fold": validation_fold,
                "selected_epochs": epoch,
                "best_inner_validation_objective": objective,
                "trace": trace,
                "elapsed_seconds": time.time() - fold_started,
            })
            print(
                f"PU0.50 seed={seed} fold={outer_fold} epochs={epoch} complete",
                flush=True,
            )
        per_seed.append(scores)
        summary = summarize(y, scores, proteins, components)
        summary.update({"model": "mlp_corrected_pu_u050", "seed": seed})
        summaries.append(summary)
        for i, score in enumerate(scores):
            rows.append({
                "model": "mlp_corrected_pu_u050",
                "seed": seed,
                "accession": proteins[i],
                "position": int(positions[i]),
                "component": components[i],
                "fold": int(folds[i]),
                "annotation_label": int(y[i]),
                "score": float(score),
                "unlabelled_weight": UNLABELLED_WEIGHT,
            })
    ensemble = np.mean(np.stack(per_seed), axis=0)
    ensemble_summary = summarize(y, ensemble, proteins, components)
    ensemble_summary.update({
        "model": "mlp_corrected_pu_u050_five_seed_mean",
        "seed": "five_seed_mean",
    })
    summaries.append(ensemble_summary)
    write_csv(RESULTS / "mlp_corrected_pu_u050_oof_by_seed.csv", rows)
    write_csv(RESULTS / "mlp_corrected_pu_u050_summary.csv", summaries)
    write_json(RESULTS / "mlp_corrected_pu_u050_epoch_selection.json", selections)
    comparisons = {
        "pu_u050_vs_corrected_unweighted_mlp": paired_component_bootstrap_all(
            y,
            load_corrected_ensemble(),
            ensemble,
            proteins,
            components,
            reps=2000,
            seed=20260930,
        ),
        "pu_u050_vs_current_hgb": paired_component_bootstrap_all(
            y,
            load_initial_scores("hist_gradient_boosting_current"),
            ensemble,
            proteins,
            components,
            reps=2000,
            seed=20260931,
        ),
    }
    write_json(RESULTS / "mlp_corrected_pu_u050_paired_intervals.json", comparisons)
    write_json(RESULTS / "mlp_corrected_pu_u050_run_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "unlabelled_weight": UNLABELLED_WEIGHT,
        "positive_weight": 1.0,
        "class_prior_estimated": False,
        "calibrated_probability_claimed": False,
        "selection_rule": "weight fixed after HGB sensitivity signal; no MLP outer-test tuning; inner component fold selects epoch only",
        "outer_split": "unchanged frozen five-fold homology components",
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sklearn": sklearn.__version__,
        },
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "frozen_benchmark_data.npz",
                INPUTS / "primary_site_folds.csv",
                RESULTS / "mlp_corrected_oof_by_seed.csv",
            ]
        },
    })


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
