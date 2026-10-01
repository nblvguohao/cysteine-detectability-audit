"""Shared, leakage-aware evaluation utilities for the Cys annotation ranker."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "inputs"
RESULTS = ROOT / "results"


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_frozen():
    archive = np.load(INPUTS / "frozen_benchmark_data.npz", allow_pickle=False)
    meta = read_csv(INPUTS / "primary_site_folds.csv")
    X = archive["X"]
    y = archive["y"].astype(int)
    folds = archive["folds"].astype(int)
    proteins = archive["proteins"].astype(str)
    components = archive["components"].astype(str)
    positions = np.asarray([int(row["position"]) for row in meta], dtype=int)
    assert X.shape == (49064, 23)
    assert len(meta) == len(y)
    assert np.array_equal(proteins, np.asarray([row["accession"] for row in meta]))
    assert np.array_equal(folds, np.asarray([int(row["fold"]) for row in meta]))
    assert np.array_equal(components, np.asarray([row["component"] for row in meta]))
    assert np.array_equal(y, np.asarray([int(row["label"]) for row in meta]))
    return X, y, folds, proteins, components, positions, meta


def protein_rows(y, scores, proteins, components):
    rows = []
    for accession in np.unique(proteins):
        idx = np.flatnonzero(proteins == accession)
        yy, ss = y[idx], scores[idx]
        if not 0 < yy.sum() < len(yy):
            continue
        tied = ss == ss.max()
        order = np.argsort(-ss, kind="stable")
        positive_ranks = np.flatnonzero(yy[order] == 1) + 1
        rows.append({
            "accession": str(accession),
            "component": str(components[idx[0]]),
            "auc": float(roc_auc_score(yy, ss)),
            "top1": float(yy[tied].mean()),
            "top2": float(yy[order[:2]].max()),
            "mrr": float(1.0 / positive_ranks.min()),
            "random_top1": float(yy.mean()),
            "random_top2": float(min(1.0, 2.0 * yy.sum() / len(yy))),
        })
    return rows


def summarize(y, scores, proteins, components):
    if not np.isfinite(scores).all():
        raise ValueError("Scores contain missing or non-finite values")
    rows = protein_rows(y, scores, proteins, components)
    top1 = float(np.mean([row["top1"] for row in rows]))
    random_top1 = float(np.mean([row["random_top1"] for row in rows]))
    return {
        "roc_auc": float(roc_auc_score(y, scores)),
        "average_precision": float(average_precision_score(y, scores)),
        "within_protein_auc_mean": float(np.mean([row["auc"] for row in rows])),
        "top1": top1,
        "top2": float(np.mean([row["top2"] for row in rows])),
        "mean_reciprocal_rank": float(np.mean([row["mrr"] for row in rows])),
        "random_top1": random_top1,
        "random_top2": float(np.mean([row["random_top2"] for row in rows])),
        "top1_enrichment": float(top1 / random_top1),
        "n_sites": int(len(y)),
        "n_positive_annotations": int(y.sum()),
        "n_proteins": int(len(rows)),
        "n_components": int(len(np.unique(components))),
    }


def validation_objective(y, scores, proteins, components):
    """Predeclared ranking objective used only inside the outer training data."""
    summary = summarize(y, scores, proteins, components)
    return float(np.mean([
        summary["within_protein_auc_mean"],
        summary["top1"],
        summary["top2"],
        summary["mean_reciprocal_rank"],
    ]))


def paired_component_bootstrap(
    y,
    reference,
    candidate,
    proteins,
    components,
    reps=5000,
    seed=20260912,
):
    ref = protein_rows(y, reference, proteins, components)
    cand = protein_rows(y, candidate, proteins, components)
    assert [r["accession"] for r in ref] == [r["accession"] for r in cand]
    assert [r["component"] for r in ref] == [r["component"] for r in cand]
    unique_components = np.unique([r["component"] for r in ref])
    component_index = {component: i for i, component in enumerate(unique_components)}
    protein_component = np.asarray([component_index[r["component"]] for r in ref])
    metrics = ("top1", "top2", "auc", "mrr")
    delta = {
        metric: np.asarray([c[metric] - r[metric] for r, c in zip(ref, cand)])
        for metric in metrics
    }
    rng = np.random.default_rng(seed)
    draws = {metric: np.empty(reps, dtype=float) for metric in metrics}
    for replicate in range(reps):
        component_weights = np.bincount(
            rng.integers(0, len(unique_components), len(unique_components)),
            minlength=len(unique_components),
        )
        protein_weights = component_weights[protein_component]
        for metric in metrics:
            draws[metric][replicate] = np.average(delta[metric], weights=protein_weights)
    return {
        "unit": "homology_component",
        "replicates": int(reps),
        "seed": int(seed),
        "estimand": "paired difference conditional on fixed outer-fold predictions",
        "does_not_cover": [
            "future experiment label uncertainty",
            "alternative homology graphs",
            "model refitting within each bootstrap replicate",
        ],
        "differences": {
            metric: {
                "candidate_minus_reference": float(delta[metric].mean()),
                "percentile_95_interval": [
                    float(value) for value in np.quantile(draws[metric], [0.025, 0.975])
                ],
            }
            for metric in metrics
        },
    }


def load_initial_scores(model_name):
    rows = read_csv(INPUTS / "initial_oof_predictions.csv")
    selected = [row for row in rows if row["model"] == model_name]
    if len(selected) != 49064:
        raise ValueError(f"Expected 49,064 predictions for {model_name}, got {len(selected)}")
    return np.asarray([float(row["score"]) for row in selected])


def load_initial_five_seed_mean():
    rows = read_csv(INPUTS / "initial_mlp_repeated_oof.csv")
    seeds = sorted({row["seed"] for row in rows}, key=int)
    blocks = []
    for seed in seeds:
        block = [row for row in rows if row["seed"] == seed]
        if len(block) != 49064:
            raise ValueError(f"Expected 49,064 initial MLP predictions for seed {seed}")
        blocks.append(np.asarray([float(row["score"]) for row in block]))
    return np.mean(np.stack(blocks), axis=0), seeds
