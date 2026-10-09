"""Retrospective external observed-site coverage audit with the corrected MLP.

PXD044043 and PXD072089 appeared in prior project work, so this is explicitly
retrospective and is not a prospective blind test.  Cysteines absent from the
reported enriched-site lists remain unlabelled.  Consequently this program does
not calculate site-level ROC-AUC/AP or within-protein AUC on these datasets.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path

import joblib
import numpy as np
import sklearn
from threadpoolctl import threadpool_limits

from common import INPUTS, RESULTS, ROOT, load_training_matrix, sha256, write_csv, write_json
from run_corrected_mlp import SEEDS, fit_preprocessor, make_mlp, transform


EXTERNAL = ROOT / "external"
DEFAULT_SOURCE_PROJECT = Path(
    "/path/to/source_project"
)
SOURCE_PROJECT = Path(os.environ.get("PERSULFIDATION_SOURCE_PROJECT", DEFAULT_SOURCE_PROJECT))
SOURCE_PIPELINE = SOURCE_PROJECT / "outputs/20260911_模型与Q1c推进/model_pipeline.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    sys.path.insert(0, str(path.parent))
    spec.loader.exec_module(module)
    return module


def read_external_inputs():
    human_path = INPUTS / "PXD044043_mapped_human.json"
    rice_path = INPUTS / "PXD072089_mapped_rice.json"
    human = json.loads(human_path.read_text())["hsa_persulfidation"]["proteins"]
    rice = json.loads(rice_path.read_text())["osa_persulfidation"]["proteins"]
    rows = []
    for dataset, records in (("human_PXD044043", human), ("rice_PXD072089", rice)):
        for record in records:
            sequence = record["seq"]
            observed = sorted(set(map(int, record["modified"])))
            assert observed and all(sequence[position - 1] == "C" for position in observed)
            rows.append({
                "dataset": dataset,
                "accession": record["acc"],
                "sequence": sequence,
                "observed": observed,
            })
    return rows


def read_overlap_screen():
    path = INPUTS / "retrospective_overlap_screen.csv"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        (row["dataset"], row["accession"]): row["strict_novel_protein"].lower() == "true"
        for row in rows
    }


class FrozenPreprocessedMLP:
    def __init__(self, imputer, scaler, model):
        self.imputer = imputer
        self.scaler = scaler
        self.model = model

    def predict_proba(self, X):
        return self.model.predict_proba(transform(self.imputer, self.scaler, X))


def final_epoch_by_seed():
    records = json.loads((RESULTS / "mlp_corrected_epoch_selection.json").read_text())
    selected = defaultdict(list)
    for record in records:
        selected[int(record["seed"])].append(int(record["selected_epochs"]))
    epochs = {}
    for seed in SEEDS:
        values = selected[seed]
        if len(values) != 5:
            raise ValueError(f"Expected five outer-fold selections for seed {seed}")
        epochs[seed] = int(statistics.median(values))
    return epochs


def fit_corrected_full_models(X, y, epochs):
    models = []
    for seed in SEEDS:
        imputer, scaler, Xt = fit_preprocessor(X)
        model = make_mlp(seed)
        for _ in range(epochs[seed]):
            model.partial_fit(Xt, y, classes=np.asarray([0, 1]))
        models.append(FrozenPreprocessedMLP(imputer, scaler, model))
        print(f"full-fit corrected MLP seed={seed} epochs={epochs[seed]}", flush=True)
    return models


def centered_window(sequence, position, radius=7):
    return sequence[max(0, position - 1 - radius): min(len(sequence), position + radius)]


def score_sites(proteins, encoder, hist_model, mlp_models):
    rows = []
    for protein in proteins:
        positions = [i + 1 for i, aa in enumerate(protein["sequence"]) if aa == "C"]
        observed = set(protein["observed"])
        X = np.asarray([
            encoder.vector(protein["sequence"], position, None, "all_basic_removed")
            for position in positions
        ])
        hist = hist_model.predict_proba(X)[:, 1]
        per_seed = np.stack([model.predict_proba(X)[:, 1] for model in mlp_models])
        corrected = per_seed.mean(axis=0)
        hist_order = np.lexsort((np.asarray(positions), -hist))
        mlp_order = np.lexsort((np.asarray(positions), -corrected))
        hist_rank = {int(index): rank + 1 for rank, index in enumerate(hist_order)}
        mlp_rank = {int(index): rank + 1 for rank, index in enumerate(mlp_order)}
        for index, position in enumerate(positions):
            rows.append({
                "dataset": protein["dataset"],
                "accession": protein["accession"],
                "position": int(position),
                "observed_in_retrospective_dataset": int(position in observed),
                "unobserved_status": "not_applicable" if position in observed else "unlabelled_not_negative",
                "histgb_score": float(hist[index]),
                "histgb_rank": hist_rank[index],
                "mlp_corrected_score": float(corrected[index]),
                "mlp_corrected_seed_sd": float(per_seed[:, index].std(ddof=1)),
                "mlp_corrected_rank": mlp_rank[index],
                "context_15aa_or_terminal_shorter": centered_window(protein["sequence"], position),
                "structure_features": "missing_for_external_sequence",
                "score_scope": "full-fit annotation-ranking score; not a chemical probability",
            })
    return rows


def random_mrr(n: int, m: int) -> float:
    denominator = math.comb(n, m)
    return sum(
        (1 / rank) * math.comb(n - rank, m - 1) / denominator
        for rank in range(1, n - m + 2)
    )


def per_protein_metrics(site_rows, score_key, model_name):
    grouped = defaultdict(list)
    for row in site_rows:
        grouped[(row["dataset"], row["accession"])].append(row)
    result = []
    for (dataset, accession), rows in sorted(grouped.items()):
        ordered = sorted(rows, key=lambda row: (-row[score_key], row["position"]))
        observed = np.asarray([row["observed_in_retrospective_dataset"] for row in ordered])
        n, m = len(observed), int(observed.sum())
        top_score = ordered[0][score_key]
        tied_top = [row for row in ordered if row[score_key] == top_score]
        top1 = float(np.mean([row["observed_in_retrospective_dataset"] for row in tied_top]))
        k = min(2, n)
        result.append({
            "dataset": dataset,
            "accession": accession,
            "model": model_name,
            "n_cys": n,
            "n_observed_sites": m,
            "top1_observed_coverage": top1,
            "top2_observed_coverage": float(observed[:k].max()),
            "reciprocal_rank_first_observed": float(1 / (np.flatnonzero(observed)[0] + 1)),
            "random_top1": float(m / n),
            "random_top2": float(1 - (math.comb(n - m, k) / math.comb(n, k) if n - m >= k else 0)),
            "random_mrr": float(random_mrr(n, m)),
            "best_observed_rank": int(np.flatnonzero(observed)[0] + 1),
            "interpretation": "coverage of a previously observed site; not future experimental success",
        })
    return result


def bootstrap_mean(values, seed, reps=5000):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.empty(reps)
    for i in range(reps):
        draws[i] = values[rng.integers(0, len(values), len(values))].mean()
    return [float(value) for value in np.quantile(draws, [0.025, 0.975])]


def summarize(protein_rows, site_rows, strict_lookup):
    site_lookup = defaultdict(list)
    for row in site_rows:
        site_lookup[(row["dataset"], row["accession"])].append(row)
    summaries = []
    comparisons = []
    for dataset in sorted({row["dataset"] for row in site_rows}):
        for subset in ("all_retrospective", "strict_novel_screen"):
            allowed = {
                (row["dataset"], row["accession"])
                for row in site_rows
                if row["dataset"] == dataset
                and (subset == "all_retrospective" or strict_lookup[(row["dataset"], row["accession"])])
            }
            n_cys = sum(len(site_lookup[key]) for key in allowed)
            n_observed = sum(
                sum(row["observed_in_retrospective_dataset"] for row in site_lookup[key])
                for key in allowed
            )
            by_model = {}
            for model in ("histgb", "mlp_corrected"):
                rows = [
                    row for row in protein_rows
                    if row["model"] == model and (row["dataset"], row["accession"]) in allowed
                ]
                by_model[model] = {row["accession"]: row for row in rows}
                metric_map = {
                    "top1": "top1_observed_coverage",
                    "top2": "top2_observed_coverage",
                    "mrr": "reciprocal_rank_first_observed",
                }
                summary = {
                    "dataset": dataset,
                    "subset": subset,
                    "model": model,
                    "n_proteins": len(rows),
                    "n_cys": n_cys,
                    "n_observed_sites": n_observed,
                    "unobserved_cys_treatment": "unlabelled; excluded from negative-class metrics",
                }
                for offset, (label, key) in enumerate(metric_map.items()):
                    random_key = "random_" + label
                    value = float(np.mean([row[key] for row in rows]))
                    random = float(np.mean([row[random_key] for row in rows]))
                    gain = [row[key] - row[random_key] for row in rows]
                    summary[label] = value
                    summary[random_key] = random
                    summary[label + "_gain_over_random"] = float(np.mean(gain))
                    summary[label + "_gain_95_protein_bootstrap"] = bootstrap_mean(
                        gain, 20260912 + offset
                    )
                summaries.append(summary)
            accessions = sorted(set(by_model["histgb"]) & set(by_model["mlp_corrected"]))
            comparison = {
                "dataset": dataset,
                "subset": subset,
                "candidate": "mlp_corrected",
                "reference": "histgb",
                "bootstrap_unit": "protein",
                "replicates": 5000,
                "paired_differences": {},
            }
            for offset, (label, key) in enumerate({
                "top1": "top1_observed_coverage",
                "top2": "top2_observed_coverage",
                "mrr": "reciprocal_rank_first_observed",
            }.items()):
                delta = [
                    by_model["mlp_corrected"][accession][key]
                    - by_model["histgb"][accession][key]
                    for accession in accessions
                ]
                comparison["paired_differences"][label] = {
                    "candidate_minus_reference": float(np.mean(delta)),
                    "percentile_95_interval": bootstrap_mean(delta, 20260922 + offset),
                }
            comparisons.append(comparison)
    return summaries, comparisons


def main():
    started = time.time()
    EXTERNAL.mkdir(parents=True, exist_ok=True)
    encoder = load_module("source_model_pipeline", SOURCE_PIPELINE)
    X, y, _, _, _, _, _ = load_training_matrix()
    epochs = final_epoch_by_seed()
    models = fit_corrected_full_models(X, y, epochs)
    payload = joblib.load(INPUTS / "annotation_ranker_histgb.joblib")
    if payload["variant"] != "all_basic_removed":
        raise ValueError("Unexpected HGB feature variant")
    proteins = read_external_inputs()
    site_rows = score_sites(proteins, encoder, payload["model"], models)
    protein_rows = (
        per_protein_metrics(site_rows, "histgb_score", "histgb")
        + per_protein_metrics(site_rows, "mlp_corrected_score", "mlp_corrected")
    )
    strict_lookup = read_overlap_screen()
    if set(strict_lookup) != {(row["dataset"], row["accession"]) for row in proteins}:
        raise ValueError("Overlap-screen keys do not match retrospective proteins")
    summaries, comparisons = summarize(protein_rows, site_rows, strict_lookup)
    write_csv(EXTERNAL / "retrospective_site_scores.csv", site_rows)
    write_csv(EXTERNAL / "retrospective_protein_metrics.csv", protein_rows)
    write_csv(EXTERNAL / "retrospective_summary.csv", summaries)
    write_json(EXTERNAL / "retrospective_paired_intervals.json", comparisons)
    audit = {
        "status": "complete_retrospective_external_observed_site_coverage_audit",
        "prospective_blind_test": False,
        "datasets": ["human_PXD044043", "rice_PXD072089"],
        "label_semantics": {
            "observed": "reported enriched/mapped Cys site in the archived dataset",
            "unobserved": "unlabelled, not a confirmed negative",
        },
        "forbidden_interpretations": [
            "Top-2 is a future experimental success probability",
            "unobserved cysteines are chemical negatives",
            "scores are calibrated probabilities of persulfidation",
            "these datasets are untouched prospective blind tests",
        ],
        "negative_class_metrics_computed": False,
        "site_roc_auc_computed": False,
        "site_average_precision_computed": False,
        "within_protein_auc_computed": False,
        "model": {
            "mlp": "five corrected full-fit models; median outer-fold-selected epoch per seed",
            "epochs_by_seed": epochs,
            "histgb": "previous frozen all_basic_removed full-fit ranker",
            "structure_features": "missing for external sequences",
        },
        "strict_screen_reused": "versioned prior sequence/homology/window screen; its exact input is hashed",
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "sklearn": sklearn.__version__,
            "joblib": joblib.__version__,
        },
        "elapsed_seconds": time.time() - started,
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "PXD044043_mapped_human.json",
                INPUTS / "PXD072089_mapped_rice.json",
                INPUTS / "retrospective_overlap_screen.csv",
                INPUTS / "annotation_ranker_histgb.joblib",
                INPUTS / "training_matrix.npz",
                RESULTS / "mlp_corrected_epoch_selection.json",
            ]
        },
        "source_encoder": {
            "path": str(SOURCE_PIPELINE),
            "sha256": sha256(SOURCE_PIPELINE),
        },
    }
    write_json(EXTERNAL / "retrospective_run_audit.json", audit)
    output_paths = sorted(EXTERNAL.glob("retrospective_*"))
    write_json(EXTERNAL / "retrospective_output_hashes.json", {
        str(path.relative_to(ROOT)): sha256(path)
        for path in output_paths
        if path.name != "retrospective_output_hashes.json"
    })
    print(json.dumps(summaries, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
