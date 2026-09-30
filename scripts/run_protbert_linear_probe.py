"""Leakage-aware five-seed linear probe on frozen ProtBERT Cys embeddings."""
from __future__ import annotations

import json
import platform
import sys
import time

import numpy as np
import torch

from common import (
    INPUTS,
    RESULTS,
    ROOT,
    load_frozen,
    load_initial_scores,
    paired_component_bootstrap,
    sha256,
    summarize,
    validation_objective,
    write_csv,
    write_json,
)


PLM = ROOT / "plm"
SEEDS = [11, 29, 47, 71, 97]
MAX_EPOCHS = 80
PATIENCE_EVALUATIONS = 8
EVALUATE_EVERY = 2
MIN_DELTA = 1e-4
LEARNING_RATE = 0.03
WEIGHT_DECAY = 1e-3


def standardize(train_x, other_x):
    mean = train_x.mean(axis=0, dtype=np.float64).astype(np.float32)
    scale = train_x.std(axis=0, dtype=np.float64).astype(np.float32)
    scale[scale < 1e-6] = 1.0
    return (train_x - mean) / scale, (other_x - mean) / scale


def make_model(seed, width):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    model = torch.nn.Linear(width, 1).to("cuda:0")
    torch.nn.init.normal_(model.weight, mean=0.0, std=0.01)
    torch.nn.init.zeros_(model.bias)
    return model


def tensors(train_x, train_y, eval_x):
    return (
        torch.from_numpy(np.ascontiguousarray(train_x)).to("cuda:0"),
        torch.from_numpy(np.asarray(train_y, dtype=np.float32)).to("cuda:0"),
        torch.from_numpy(np.ascontiguousarray(eval_x)).to("cuda:0"),
    )


def train_epoch(model, optimizer, train_x, train_y):
    optimizer.zero_grad(set_to_none=True)
    logits = model(train_x).squeeze(1)
    loss = torch.nn.functional.binary_cross_entropy_with_logits(logits, train_y)
    loss.backward()
    optimizer.step()
    return float(loss.detach().cpu())


def predict(model, x):
    with torch.inference_mode():
        return torch.sigmoid(model(x).squeeze(1)).cpu().numpy()


def select_epochs(X, y, proteins, components, folds, outer_fold, seed):
    validation_fold = outer_fold % 5 + 1
    subtrain = (folds != outer_fold) & (folds != validation_fold)
    validation = folds == validation_fold
    outer_test = folds == outer_fold
    for left, right in ((subtrain, validation), (subtrain, outer_test), (validation, outer_test)):
        assert not set(components[left]).intersection(components[right])
        assert not set(proteins[left]).intersection(proteins[right])
    train_x, valid_x = standardize(X[subtrain], X[validation])
    train_x, train_y, valid_x = tensors(train_x, y[subtrain], valid_x)
    model = make_model(seed, X.shape[1])
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )
    best_epoch, best_objective, stale = 1, -np.inf, 0
    trace = []
    for epoch in range(1, MAX_EPOCHS + 1):
        loss = train_epoch(model, optimizer, train_x, train_y)
        if epoch % EVALUATE_EVERY != 0:
            continue
        scores = predict(model, valid_x)
        objective = validation_objective(
            y[validation], scores, proteins[validation], components[validation]
        )
        trace.append({"epoch": epoch, "loss": loss, "validation_objective": objective})
        if objective > best_objective + MIN_DELTA:
            best_epoch, best_objective, stale = epoch, objective, 0
        else:
            stale += 1
        if stale >= PATIENCE_EVALUATIONS:
            break
    del train_x, train_y, valid_x, model, optimizer
    torch.cuda.empty_cache()
    return best_epoch, best_objective, validation_fold, trace


def refit_and_predict(X, y, folds, outer_fold, seed, epochs):
    train = folds != outer_fold
    test = folds == outer_fold
    train_x, test_x = standardize(X[train], X[test])
    train_x, train_y, test_x = tensors(train_x, y[train], test_x)
    model = make_model(seed, X.shape[1])
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )
    for _ in range(epochs):
        train_epoch(model, optimizer, train_x, train_y)
    scores = predict(model, test_x)
    del train_x, train_y, test_x, model, optimizer
    torch.cuda.empty_cache()
    return test, scores


def main():
    started = time.time()
    torch.cuda.reset_peak_memory_stats()
    X_table, y, folds, proteins, components, positions, _ = load_frozen()
    archive = np.load(PLM / "protbert_cys_embeddings.npz")
    X = archive["embeddings"].astype(np.float32)
    assert X.shape == (len(y), 1024)
    assert np.array_equal(archive["proteins"].astype(str), proteins)
    assert np.array_equal(archive["positions"].astype(int), positions)
    reference_hgb = load_initial_scores("hist_gradient_boosting_current")
    predictions, selection = [], []
    for seed in SEEDS:
        scores = np.full(len(y), np.nan)
        for outer_fold in range(1, 6):
            fold_started = time.time()
            epochs, objective, validation_fold, trace = select_epochs(
                X, y, proteins, components, folds, outer_fold, seed
            )
            mask, fold_scores = refit_and_predict(
                X, y, folds, outer_fold, seed, epochs
            )
            scores[mask] = fold_scores
            selection.append({
                "seed": seed,
                "outer_fold": outer_fold,
                "inner_validation_fold": validation_fold,
                "selected_epochs": epochs,
                "best_inner_validation_objective": objective,
                "trace": trace,
                "elapsed_seconds": time.time() - fold_started,
            })
            print(json.dumps({"seed": seed, "outer_fold": outer_fold, "epochs": epochs}), flush=True)
        if not np.isfinite(scores).all():
            raise RuntimeError(f"Incomplete scores for seed {seed}")
        predictions.append(scores)

    ensemble = np.mean(np.stack(predictions), axis=0)
    summaries = []
    for seed, scores in zip(SEEDS, predictions):
        row = summarize(y, scores, proteins, components)
        row.update({"model": "protbert_frozen_linear_probe", "seed": seed})
        summaries.append(row)
    ensemble_summary = summarize(y, ensemble, proteins, components)
    ensemble_summary.update({"model": "protbert_frozen_linear_probe_five_seed_mean", "seed": "five_seed_mean"})
    summaries.append(ensemble_summary)
    rows = []
    for seed, scores in zip(SEEDS, predictions):
        for i, score in enumerate(scores):
            rows.append({
                "model": "protbert_frozen_linear_probe",
                "seed": seed,
                "accession": proteins[i],
                "position": int(positions[i]),
                "component": components[i],
                "fold": int(folds[i]),
                "label": int(y[i]),
                "score": float(score),
                "label_scope": "method_dependent_annotation_vs_unlabelled",
            })
    write_csv(RESULTS / "protbert_linear_probe_oof_by_seed.csv", rows)
    write_csv(RESULTS / "protbert_linear_probe_summary.csv", summaries)
    write_json(RESULTS / "protbert_linear_probe_epoch_selection.json", selection)
    write_json(RESULTS / "protbert_linear_probe_paired_intervals.json", {
        "protbert_probe_vs_current_hist_gradient_boosting": paired_component_bootstrap(
            y, reference_hgb, ensemble, proteins, components, reps=5000, seed=20260913
        )
    })
    # Site order is asserted against the embedding archive above and recorded in
    # the run audit. Keep this compact derivative to the score vector only.
    np.savez(PLM / "protbert_linear_probe_oof_complete.npz", scores=ensemble)
    write_json(RESULTS / "protbert_linear_probe_run_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "representation": "frozen ProtBERT last-layer Cys token; 15 residues per flank",
        "probe": "single linear layer trained with full-batch unweighted BCE and AdamW",
        "outer_split": "unchanged frozen five-fold homology components",
        "inner_split": "cyclic one of the remaining frozen component folds",
        "validation_objective": "unweighted mean(within-protein AUC, Top-1, Top-2, MRR)",
        "test_fold_used_for_selection": False,
        "pretraining_exposure_controlled": False,
        "pretraining_caveat": "downstream homology isolation does not prove that related sequences were absent from ProtBERT pretraining",
        "seeds": SEEDS,
        "hyperparameters": {
            "max_epochs": MAX_EPOCHS,
            "patience_evaluations": PATIENCE_EVALUATIONS,
            "evaluate_every": EVALUATE_EVERY,
            "learning_rate": LEARNING_RATE,
            "weight_decay": WEIGHT_DECAY,
            "min_delta": MIN_DELTA,
        },
        "hardware": {
            "platform": platform.platform(),
            "gpu": torch.cuda.get_device_name(0),
            "capability": list(torch.cuda.get_device_capability(0)),
            "peak_gpu_memory_bytes": int(torch.cuda.max_memory_allocated()),
        },
        "versions": {"python": sys.version, "numpy": np.__version__, "torch": torch.__version__},
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "frozen_benchmark_data.npz",
                INPUTS / "primary_site_folds.csv",
                INPUTS / "initial_oof_predictions.csv",
                PLM / "protbert_cys_embeddings.npz",
            ]
        },
        "ensemble_summary": ensemble_summary,
    })


if __name__ == "__main__":
    main()
