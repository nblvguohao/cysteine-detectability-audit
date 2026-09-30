"""Unified leakage-aware evaluation of every currently valid OOF prediction."""
from __future__ import annotations

import csv
import json
import sys
import time
from collections import defaultdict

import numpy as np
import scipy
import sklearn
from sklearn.metrics import average_precision_score, roc_auc_score

from common import (
    INPUTS,
    RESULTS,
    ROOT,
    load_frozen,
    load_initial_scores,
    protein_rows,
    read_csv,
    sha256,
    summarize,
    write_csv,
    write_json,
)


def load_score_column(path):
    rows = read_csv(path)
    return np.asarray([float(row["score"]) for row in rows])


def load_corrected_ensemble():
    rows = read_csv(RESULTS / "mlp_corrected_oof_by_seed.csv")
    by_seed = defaultdict(list)
    keys = defaultdict(list)
    for row in rows:
        seed = int(row["seed"])
        by_seed[seed].append(float(row["score"]))
        keys[seed].append((row["accession"], int(row["position"])))
    seeds = sorted(by_seed)
    if seeds != [11, 29, 47, 71, 97]:
        raise ValueError(f"Unexpected seeds: {seeds}")
    if any(keys[seed] != keys[seeds[0]] for seed in seeds[1:]):
        raise ValueError("Corrected MLP site order differs across seeds")
    return np.mean(np.stack([by_seed[seed] for seed in seeds]), axis=0)


def load_complete_npz_scores(path, y=None, proteins=None, positions=None):
    archive = np.load(path)
    scores = archive["scores"]
    if y is not None and "y" in archive:
        assert np.array_equal(archive["y"].astype(int), y)
    if proteins is not None and "proteins" in archive:
        assert np.array_equal(archive["proteins"].astype(str), proteins)
    if positions is not None and "positions" in archive:
        assert np.array_equal(archive["positions"].astype(int), positions)
    return scores


def parse_bool(value):
    return str(value).strip().lower() == "true"


def build_strata(meta, proteins):
    structure = np.asarray([parse_bool(row["structure_available"]) for row in meta])
    disulfide = np.asarray([parse_bool(row["predicted_disulfide"]) for row in meta])
    pka_masked = np.asarray([parse_bool(row["pka_masked"]) for row in meta])
    counts = defaultdict(int)
    positive_counts = defaultdict(int)
    for row in meta:
        counts[row["accession"]] += 1
        positive_counts[row["accession"]] += int(row["label"])
    site_counts = np.asarray([counts[protein] for protein in proteins])
    observed_counts = np.asarray([positive_counts[protein] for protein in proteins])
    strata = {
        "all": np.ones(len(meta), dtype=bool),
        "structure_available": structure,
        "structure_missing": ~structure,
        "predicted_disulfide": disulfide,
        "not_predicted_disulfide": ~disulfide,
        "pka_masked": pka_masked,
        "pka_observed": structure & ~pka_masked,
        "protein_2_to_5_candidate_cys": (site_counts >= 2) & (site_counts <= 5),
        "protein_6_to_10_candidate_cys": (site_counts >= 6) & (site_counts <= 10),
        "protein_over_10_candidate_cys": site_counts > 10,
        "protein_single_observed_annotation": observed_counts == 1,
        "protein_multiple_observed_annotations": observed_counts > 1,
    }
    return strata


def stratified_auc_ap(y, models, strata, folds):
    rows = []
    expanded = dict(strata)
    for fold in range(1, 6):
        expanded[f"outer_fold_{fold}"] = folds == fold
    for model, scores in models.items():
        for stratum, mask in expanded.items():
            yy, ss = y[mask], scores[mask]
            if len(np.unique(yy)) < 2:
                auc = ap = None
            else:
                auc = float(roc_auc_score(yy, ss))
                ap = float(average_precision_score(yy, ss))
            rows.append({
                "model": model,
                "stratum": stratum,
                "n_sites": int(mask.sum()),
                "n_annotation_positive": int(yy.sum()),
                "annotation_rate": float(yy.mean()),
                "roc_auc_annotation_vs_unlabelled": auc,
                "average_precision_annotation_vs_unlabelled": ap,
                "interpretation": "method-dependent annotation discrimination; not chemical specificity",
            })
    return rows


def score_groups(scores):
    """Return fixed tie groups so weighted AUC/AP avoids repeated sorting."""
    _, inverse = np.unique(scores, return_inverse=True)
    return inverse, int(inverse.max()) + 1


def weighted_auc_ap(y, weights, inverse, n_groups):
    positive = np.bincount(inverse, weights=weights * y, minlength=n_groups)
    negative = np.bincount(inverse, weights=weights * (1 - y), minlength=n_groups)
    total_positive = positive.sum()
    total_negative = negative.sum()
    negative_before = np.cumsum(negative) - negative
    auc = np.sum(positive * (negative_before + 0.5 * negative)) / (
        total_positive * total_negative
    )
    positive_desc = positive[::-1]
    negative_desc = negative[::-1]
    cumulative_positive = np.cumsum(positive_desc)
    cumulative_total = np.cumsum(positive_desc + negative_desc)
    ap = np.sum(
        (positive_desc / total_positive)
        * np.divide(
            cumulative_positive,
            cumulative_total,
            out=np.zeros_like(cumulative_positive),
            where=cumulative_total > 0,
        )
    )
    return float(auc), float(ap)


def paired_component_bootstrap_all(
    y, reference, candidate, proteins, components, reps=2000, seed=20260912
):
    ref_protein = protein_rows(y, reference, proteins, components)
    cand_protein = protein_rows(y, candidate, proteins, components)
    if [row["accession"] for row in ref_protein] != [row["accession"] for row in cand_protein]:
        raise ValueError("Protein order mismatch")
    unique_components, site_component = np.unique(components, return_inverse=True)
    component_lookup = {component: i for i, component in enumerate(unique_components)}
    protein_component = np.asarray([
        component_lookup[row["component"]] for row in ref_protein
    ])
    delta_protein = {
        "within_protein_auc": np.asarray([
            candidate_row["auc"] - reference_row["auc"]
            for reference_row, candidate_row in zip(ref_protein, cand_protein)
        ]),
        "top1": np.asarray([
            candidate_row["top1"] - reference_row["top1"]
            for reference_row, candidate_row in zip(ref_protein, cand_protein)
        ]),
        "top2": np.asarray([
            candidate_row["top2"] - reference_row["top2"]
            for reference_row, candidate_row in zip(ref_protein, cand_protein)
        ]),
        "mrr": np.asarray([
            candidate_row["mrr"] - reference_row["mrr"]
            for reference_row, candidate_row in zip(ref_protein, cand_protein)
        ]),
    }
    draws = {
        metric: np.empty(reps) for metric in [
            "roc_auc", "average_precision", *delta_protein
        ]
    }
    reference_groups = score_groups(reference)
    candidate_groups = score_groups(candidate)
    rng = np.random.default_rng(seed)
    for replicate in range(reps):
        component_weights = np.bincount(
            rng.integers(0, len(unique_components), len(unique_components)),
            minlength=len(unique_components),
        )
        site_weights = component_weights[site_component]
        protein_weights = component_weights[protein_component]
        reference_auc, reference_ap = weighted_auc_ap(
            y, site_weights, *reference_groups
        )
        candidate_auc, candidate_ap = weighted_auc_ap(
            y, site_weights, *candidate_groups
        )
        draws["roc_auc"][replicate] = candidate_auc - reference_auc
        draws["average_precision"][replicate] = candidate_ap - reference_ap
        for metric, delta in delta_protein.items():
            draws[metric][replicate] = np.average(delta, weights=protein_weights)
    point = {
        "roc_auc": float(roc_auc_score(y, candidate) - roc_auc_score(y, reference)),
        "average_precision": float(
            average_precision_score(y, candidate) - average_precision_score(y, reference)
        ),
        **{metric: float(delta.mean()) for metric, delta in delta_protein.items()},
    }
    return {
        "bootstrap_unit": "homology_component",
        "replicates": reps,
        "seed": seed,
        "estimand": "paired difference conditional on frozen folds and fixed OOF predictions",
        "does_not_cover": [
            "model refitting in each replicate",
            "alternative homology graphs",
            "future experimental label uncertainty",
        ],
        "differences": {
            metric: {
                "candidate_minus_reference": point[metric],
                "percentile_95_interval": [
                    float(value) for value in np.quantile(draws[metric], [0.025, 0.975])
                ],
            }
            for metric in draws
        },
    }


def partial_tabicl_audit(y, folds, proteins, components, reference):
    archive = np.load(INPUTS / "tabicl_oof_partial.npz")
    scores = archive["scores"]
    completed = [int(value) for value in archive["folds_completed"]]
    mask = np.isfinite(scores)
    expected = np.isin(folds, completed)
    if not np.array_equal(mask, expected):
        raise ValueError("TabICL finite predictions do not match recorded completed folds")
    return {
        "status": "incomplete_not_eligible_for_overall_model_selection",
        "completed_outer_folds": completed,
        "n_scored_sites": int(mask.sum()),
        "tabicl_partial_metrics": summarize(
            y[mask], scores[mask], proteins[mask], components[mask]
        ),
        "same_sites_hgb_metrics": summarize(
            y[mask], reference[mask], proteins[mask], components[mask]
        ),
        "warning": "Fold-1-only results are descriptive and cannot be compared with complete five-fold OOF results.",
    }


def main():
    started = time.time()
    X, y, folds, proteins, components, positions, meta = load_frozen()
    del X
    models = {
        "histgb_23_feature_current": load_initial_scores("hist_gradient_boosting_current"),
        "mlp_corrected_five_seed_mean": load_corrected_ensemble(),
        "histgb_sequence_only": load_score_column(INPUTS / "oof_sequence_only.csv"),
        "histgb_structure_only": load_score_column(INPUTS / "oof_structure_only.csv"),
        "tabicl_v2_cuda": load_complete_npz_scores(
            ROOT / "tabular" / "tabicl_cuda_oof_complete.npz",
            y=y,
            proteins=proteins,
            positions=positions,
        ),
        "protbert_frozen_linear_probe": load_complete_npz_scores(
            ROOT / "plm" / "protbert_linear_probe_oof_complete.npz"
        ),
    }
    if any(len(scores) != len(y) or not np.isfinite(scores).all() for scores in models.values()):
        raise ValueError("Incomplete valid-model OOF predictions")
    metrics = []
    for model, scores in models.items():
        row = summarize(y, scores, proteins, components)
        row.update({
            "model": model,
            "label_scope": "method-dependent annotation versus unlabelled",
            "outer_split": "frozen five-fold homology components",
        })
        metrics.append(row)
    strata = build_strata(meta, proteins)
    stratified = stratified_auc_ap(y, models, strata, folds)
    comparisons = {}
    print("bootstrap: corrected MLP vs HGB", flush=True)
    comparisons["mlp_corrected_vs_histgb_current"] = paired_component_bootstrap_all(
            y,
            models["histgb_23_feature_current"],
            models["mlp_corrected_five_seed_mean"],
            proteins,
            components,
        )
    print("bootstrap: 23-feature HGB vs sequence-only", flush=True)
    comparisons["histgb_23_feature_vs_sequence_only"] = paired_component_bootstrap_all(
            y,
            models["histgb_sequence_only"],
            models["histgb_23_feature_current"],
            proteins,
            components,
            seed=20260913,
        )
    print("bootstrap: 23-feature HGB vs structure-only", flush=True)
    comparisons["histgb_23_feature_vs_structure_only"] = paired_component_bootstrap_all(
            y,
            models["histgb_structure_only"],
            models["histgb_23_feature_current"],
            proteins,
            components,
            seed=20260914,
        )
    print("bootstrap: TabICL vs HGB", flush=True)
    comparisons["tabicl_vs_histgb_current"] = paired_component_bootstrap_all(
        y,
        models["histgb_23_feature_current"],
        models["tabicl_v2_cuda"],
        proteins,
        components,
        seed=20260915,
    )
    print("bootstrap: TabICL vs corrected MLP", flush=True)
    comparisons["tabicl_vs_corrected_mlp"] = paired_component_bootstrap_all(
        y,
        models["mlp_corrected_five_seed_mean"],
        models["tabicl_v2_cuda"],
        proteins,
        components,
        seed=20260916,
    )
    print("bootstrap: ProtBERT probe vs HGB", flush=True)
    comparisons["protbert_probe_vs_histgb_current"] = paired_component_bootstrap_all(
        y,
        models["histgb_23_feature_current"],
        models["protbert_frozen_linear_probe"],
        proteins,
        components,
        seed=20260917,
    )
    print("bootstrap: ProtBERT probe vs corrected MLP", flush=True)
    comparisons["protbert_probe_vs_corrected_mlp"] = paired_component_bootstrap_all(
        y,
        models["mlp_corrected_five_seed_mean"],
        models["protbert_frozen_linear_probe"],
        proteins,
        components,
        seed=20260918,
    )
    status = {
        "complete_valid_oof_models": list(models),
        "invalidated": {
            "initial_mlp_random_site_early_stopping": "protein/homology leakage in inner validation; superseded by corrected MLP",
        },
        "completed_gpu_branches": {
            "tabicl": {
                "status": "complete_five_fold_oof",
                "execution": "home RTX 5070 Ti",
                "audit": "tabular/tabicl_cuda_run_audit.json",
            },
            "frozen_plm": {
                "status": "complete_five_seed_five_fold_oof",
                "execution": "home RTX 5070 Ti",
                "audit": "results/protbert_linear_probe_run_audit.json",
            },
        },
        "incomplete": {
            "tabpfn": {
                "status": "not_run",
                "blocker": "licensed/auditable local model weights absent",
                "substitution_used": False,
            },
        },
        "historical_partial_result": {
            "tabicl_fold_1_only": partial_tabicl_audit(
                y, folds, proteins, components, models["histgb_23_feature_current"]
            )
        },
        "task_separation": {
            "site_localization": "current cohort contains known-annotation proteins and evaluates within-protein ranking",
            "protein_occurrence_detection": "not evaluated; requires proteins with defensible protein-level positives and negatives/unlabelled design",
        },
    }
    write_csv(RESULTS / "unified_model_metrics.csv", metrics)
    write_csv(RESULTS / "unified_stratified_auc_ap.csv", stratified)
    write_json(RESULTS / "unified_paired_component_intervals.json", comparisons)
    write_json(RESULTS / "model_completion_status.json", status)
    write_json(RESULTS / "unified_evaluation_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sklearn": sklearn.__version__,
        },
        "n_sites": len(y),
        "n_annotation_positive": int(y.sum()),
        "n_proteins": len(np.unique(proteins)),
        "n_homology_components": len(np.unique(components)),
        "label_scope": "method-dependent annotation versus unlabelled; not a chemical positive/negative dataset",
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "frozen_benchmark_data.npz",
                INPUTS / "primary_site_folds.csv",
                INPUTS / "initial_oof_predictions.csv",
                INPUTS / "oof_sequence_only.csv",
                INPUTS / "oof_structure_only.csv",
                INPUTS / "tabicl_oof_partial.npz",
                RESULTS / "mlp_corrected_oof_by_seed.csv",
                ROOT / "tabular" / "tabicl_cuda_oof_complete.npz",
                ROOT / "tabular" / "tabicl_cuda_run_audit.json",
                ROOT / "plm" / "protbert_linear_probe_oof_complete.npz",
                RESULTS / "protbert_linear_probe_run_audit.json",
            ]
        },
    })


if __name__ == "__main__":
    main()
