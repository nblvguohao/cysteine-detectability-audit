"""Predeclared PU-cost and evidence-specificity sensitivity branches.

This is not an identifiable chemical-probability PU estimator: the positive
class prior and selection propensity are unknown.  It tests whether within-
protein annotation ranking is robust when unlabelled examples are downweighted
and ambiguous protein-group mappings receive lower relative positive weight.
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import scipy
import sklearn
from sklearn.ensemble import HistGradientBoostingClassifier

from common import INPUTS, RESULTS, ROOT, load_frozen, read_csv, sha256, summarize, write_csv, write_json
from run_unified_evaluation import paired_component_bootstrap_all


def evidence_specificity(meta, y):
    evidence_rows = read_csv(INPUTS / "site_evidence_index.csv")
    index = {
        (row["protein_accession"], int(row["position"])): row
        for row in evidence_rows
    }
    specificity = np.ones(len(y), dtype=float)
    joined = 0
    ambiguous = 0
    max_probability = []
    repeated_rows = 0
    for i, row in enumerate(meta):
        if not y[i]:
            continue
        evidence = index.get((row["accession"], int(row["position"])))
        if evidence is None:
            raise ValueError(f"Positive site absent from evidence index: {row}")
        joined += 1
        protein_groups = {value for value in evidence["protein_groups"].split(";") if value}
        specificity[i] = 1.0 / len(protein_groups)
        ambiguous += len(protein_groups) > 1
        max_probability.append(float(evidence["max_site_probability"]))
        repeated_rows += len(evidence["source_excel_rows"].split(";")) > 1
    return specificity, {
        "n_positive_annotations": int(y.sum()),
        "n_joined": joined,
        "n_ambiguous_protein_group_mapping": ambiguous,
        "n_multiple_source_rows": repeated_rows,
        "max_site_probability_equal_one": int(np.sum(np.asarray(max_probability) == 1.0)),
        "evidence_rule": "relative positive weight = 1 / number of mapped protein-group accessions",
        "normalization": "positive weights normalized to mean one inside each outer-training fold",
        "not_used": [
            "source row multiplicity because rows are not established independent biological replicates",
            "max localization probability because nearly all joined values equal 1 and it is effectively degenerate",
        ],
    }


def parameters():
    return {
        "max_iter": 300,
        "learning_rate": 0.06,
        "max_leaf_nodes": 31,
        "max_depth": None,
        "min_samples_leaf": 20,
        "l2_regularization": 0.0,
        "max_bins": 255,
        "early_stopping": False,
        "random_state": 0,
        "loss": "log_loss",
        "categorical_features": None,
    }


def training_weights(y, evidence, train, unlabeled_weight, use_evidence):
    weights = np.where(y[train] == 1, 1.0, float(unlabeled_weight))
    if use_evidence:
        positive = y[train] == 1
        raw = evidence[train][positive]
        weights[positive] = raw / raw.mean()
    return weights


def main():
    started = time.time()
    X, y, folds, proteins, components, positions, meta = load_frozen()
    evidence, evidence_audit = evidence_specificity(meta, y)
    variants = [
        {"model": "hgb_unweighted_reproduction", "unlabelled_weight": 1.0, "evidence": False},
        {"model": "hgb_pu_sensitivity_u025", "unlabelled_weight": 0.25, "evidence": False},
        {"model": "hgb_pu_sensitivity_u050", "unlabelled_weight": 0.50, "evidence": False},
        {"model": "hgb_pu_sensitivity_u075", "unlabelled_weight": 0.75, "evidence": False},
        {"model": "hgb_evidence_specificity", "unlabelled_weight": 1.0, "evidence": True},
        {"model": "hgb_pu_u050_evidence_specificity", "unlabelled_weight": 0.50, "evidence": True},
    ]
    predictions = {}
    summaries = []
    prediction_rows = []
    for variant in variants:
        scores = np.full(len(y), np.nan)
        for fold in range(1, 6):
            train = folds != fold
            test = folds == fold
            if set(components[train]) & set(components[test]):
                raise ValueError("Homology component crossed folds")
            model = HistGradientBoostingClassifier(**parameters())
            weights = training_weights(
                y, evidence, train, variant["unlabelled_weight"], variant["evidence"]
            )
            model.fit(X[train], y[train], sample_weight=weights)
            scores[test] = model.predict_proba(X[test])[:, 1]
            print(variant["model"], "fold", fold, "complete", flush=True)
        predictions[variant["model"]] = scores
        summary = summarize(y, scores, proteins, components)
        summary.update({
            **variant,
            "branch_scope": "annotation-ranking sensitivity; not calibrated chemical PU probability",
        })
        summaries.append(summary)
        for i, score in enumerate(scores):
            prediction_rows.append({
                "model": variant["model"],
                "accession": proteins[i],
                "position": int(positions[i]),
                "component": components[i],
                "fold": int(folds[i]),
                "annotation_label": int(y[i]),
                "score": float(score),
                "unlabelled_weight": variant["unlabelled_weight"],
                "uses_evidence_specificity": variant["evidence"],
            })
    reference = predictions["hgb_unweighted_reproduction"]
    comparisons = {}
    for offset, variant in enumerate(variants[1:], 1):
        name = variant["model"]
        print("bootstrap", name, flush=True)
        comparisons[name + "_vs_unweighted"] = paired_component_bootstrap_all(
            y,
            reference,
            predictions[name],
            proteins,
            components,
            reps=1000,
            seed=20260920 + offset,
        )
    write_csv(RESULTS / "pu_evidence_oof.csv", prediction_rows)
    write_csv(RESULTS / "pu_evidence_summary.csv", summaries)
    write_json(RESULTS / "pu_evidence_paired_intervals.json", comparisons)
    write_json(RESULTS / "pu_evidence_run_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "design": "predeclared cost-sensitivity grid and evidence-specificity ablation; no test-fold model selection",
        "outer_split": "unchanged frozen five-fold homology components",
        "class_prior_estimated": False,
        "selection_propensity_identified": False,
        "probability_calibration_claimed": False,
        "chemical_negative_class_claimed": False,
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sklearn": sklearn.__version__,
        },
        "evidence_audit": evidence_audit,
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "frozen_benchmark_data.npz",
                INPUTS / "primary_site_folds.csv",
                INPUTS / "site_evidence_index.csv",
            ]
        },
    })


if __name__ == "__main__":
    main()
