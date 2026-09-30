"""Five-seed MLP with homology-isolated inner validation and fixed-epoch refit.

For each outer fold, one of the four remaining frozen homology folds is used only
to choose the epoch count. The model is then reinitialized and refit for that
many epochs on every component in the outer-training partition. The outer test
fold is never used for preprocessing, epoch selection, or threshold selection.
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import sys
import time
from pathlib import Path

import numpy as np
import scipy
import sklearn
from sklearn.impute import SimpleImputer
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from common import (
    INPUTS,
    RESULTS,
    ROOT,
    load_frozen,
    load_initial_five_seed_mean,
    load_initial_scores,
    paired_component_bootstrap,
    sha256,
    summarize,
    validation_objective,
    write_csv,
    write_json,
)


SEEDS = [11, 29, 47, 71, 97]
MAX_EPOCHS = 200
PATIENCE = 15
MIN_DELTA = 1e-4
HIDDEN = (64, 32)


def make_mlp(seed):
    return MLPClassifier(
        hidden_layer_sizes=HIDDEN,
        activation="relu",
        solver="adam",
        alpha=1e-3,
        batch_size=256,
        learning_rate_init=1e-3,
        max_iter=1,
        shuffle=True,
        early_stopping=False,
        random_state=seed,
        tol=0.0,
    )


def fit_preprocessor(X):
    imputer = SimpleImputer(strategy="median", add_indicator=True)
    scaler = StandardScaler()
    return imputer, scaler, scaler.fit_transform(imputer.fit_transform(X))


def transform(imputer, scaler, X):
    return scaler.transform(imputer.transform(X))


def select_epochs(
    X,
    y,
    proteins,
    components,
    folds,
    outer_fold,
    seed,
):
    validation_fold = outer_fold % 5 + 1
    subtrain = (folds != outer_fold) & (folds != validation_fold)
    validation = folds == validation_fold
    outer_test = folds == outer_fold
    for left, right in ((subtrain, validation), (subtrain, outer_test), (validation, outer_test)):
        assert not set(components[left]).intersection(components[right])
        assert not set(proteins[left]).intersection(proteins[right])

    imputer, scaler, X_train = fit_preprocessor(X[subtrain])
    X_valid = transform(imputer, scaler, X[validation])
    model = make_mlp(seed)
    best_epoch = 1
    best_objective = -np.inf
    stale = 0
    trace = []
    for epoch in range(1, MAX_EPOCHS + 1):
        model.partial_fit(X_train, y[subtrain], classes=np.asarray([0, 1]))
        valid_scores = model.predict_proba(X_valid)[:, 1]
        objective = validation_objective(
            y[validation], valid_scores, proteins[validation], components[validation]
        )
        trace.append({"epoch": epoch, "validation_objective": objective})
        if objective > best_objective + MIN_DELTA:
            best_objective = objective
            best_epoch = epoch
            stale = 0
        else:
            stale += 1
        if stale >= PATIENCE:
            break
    return best_epoch, best_objective, validation_fold, trace


def refit_and_predict(X, y, folds, outer_fold, seed, epochs):
    outer_train = folds != outer_fold
    outer_test = folds == outer_fold
    imputer, scaler, X_train = fit_preprocessor(X[outer_train])
    X_test = transform(imputer, scaler, X[outer_test])
    model = make_mlp(seed)
    for _ in range(epochs):
        model.partial_fit(X_train, y[outer_train], classes=np.asarray([0, 1]))
    return outer_test, model.predict_proba(X_test)[:, 1]


def main():
    started = time.time()
    RESULTS.mkdir(parents=True, exist_ok=True)
    X, y, folds, proteins, components, positions, _ = load_frozen()
    reference = load_initial_scores("hist_gradient_boosting_current")
    flawed_mean, flawed_seeds = load_initial_five_seed_mean()
    seed_predictions = []
    seed_summaries = []
    selection_records = []

    for seed in SEEDS:
        scores = np.full(len(y), np.nan)
        for outer_fold in range(1, 6):
            fold_started = time.time()
            epochs, objective, validation_fold, trace = select_epochs(
                X, y, proteins, components, folds, outer_fold, seed
            )
            test_mask, test_scores = refit_and_predict(
                X, y, folds, outer_fold, seed, epochs
            )
            scores[test_mask] = test_scores
            selection_records.append({
                "seed": seed,
                "outer_fold": outer_fold,
                "inner_validation_fold": validation_fold,
                "selected_epochs": epochs,
                "best_inner_validation_objective": objective,
                "trace": trace,
                "elapsed_seconds": time.time() - fold_started,
                "n_subtrain_components": int(len(np.unique(components[(folds != outer_fold) & (folds != validation_fold)]))),
                "n_validation_components": int(len(np.unique(components[folds == validation_fold]))),
                "n_outer_test_components": int(len(np.unique(components[folds == outer_fold]))),
            })
            print(
                f"seed={seed} outer_fold={outer_fold} validation_fold={validation_fold} "
                f"epochs={epochs} complete",
                flush=True,
            )
        if not np.isfinite(scores).all():
            raise RuntimeError(f"Incomplete scores for seed {seed}")
        seed_predictions.append(scores)
        seed_summary = summarize(y, scores, proteins, components)
        seed_summary["model"] = "mlp_corrected"
        seed_summary["seed"] = seed
        seed_summaries.append(seed_summary)

    ensemble = np.mean(np.stack(seed_predictions), axis=0)
    ensemble_summary = summarize(y, ensemble, proteins, components)
    ensemble_summary.update({"model": "mlp_corrected_five_seed_mean", "seed": "five_seed_mean"})
    seed_summaries.append(ensemble_summary)

    rows = []
    for seed, scores in zip(SEEDS, seed_predictions):
        for i, score in enumerate(scores):
            rows.append({
                "model": "mlp_corrected",
                "seed": seed,
                "accession": proteins[i],
                "position": int(positions[i]),
                "component": components[i],
                "fold": int(folds[i]),
                "label": int(y[i]),
                "score": float(score),
                "label_scope": "method_dependent_annotation_vs_unlabelled",
            })
    write_csv(RESULTS / "mlp_corrected_oof_by_seed.csv", rows)
    write_csv(RESULTS / "mlp_corrected_summary.csv", seed_summaries)
    write_json(RESULTS / "mlp_corrected_epoch_selection.json", selection_records)

    comparisons = {
        "corrected_mlp_vs_current_hist_gradient_boosting": paired_component_bootstrap(
            y, reference, ensemble, proteins, components
        ),
        "corrected_mlp_vs_initial_leaky_early_stopping_mlp": paired_component_bootstrap(
            y, flawed_mean, ensemble, proteins, components
        ),
    }
    write_json(RESULTS / "mlp_corrected_paired_intervals.json", comparisons)
    write_json(
        RESULTS / "mlp_corrected_run_audit.json",
        {
            "completed": True,
            "elapsed_seconds": time.time() - started,
            "algorithm": "sklearn MLPClassifier partial_fit; inner component holdout selects epoch; reinitialize and refit on all outer-train components",
            "outer_split": "unchanged frozen five-fold homology components",
            "inner_split": "cyclic one of the remaining frozen component folds",
            "validation_objective": "unweighted mean(within-protein AUC, Top-1, Top-2, MRR)",
            "test_fold_used_for_selection": False,
            "threshold_selected": False,
            "seeds": SEEDS,
            "initial_flawed_mlp_seeds": flawed_seeds,
            "max_epochs": MAX_EPOCHS,
            "patience": PATIENCE,
            "min_delta": MIN_DELTA,
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
                    INPUTS / "initial_oof_predictions.csv",
                    INPUTS / "initial_mlp_repeated_oof.csv",
                ]
            },
            "ensemble_summary": ensemble_summary,
        },
    )


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
