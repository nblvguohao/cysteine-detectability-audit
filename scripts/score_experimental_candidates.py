"""Score supplied proteins with corrected MLP and PU sensitivity ensembles."""
from __future__ import annotations

import csv
import importlib.util
import json
import os
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from common import INPUTS, RESULTS, ROOT, load_frozen, read_csv, sha256, write_csv, write_json
from run_corrected_mlp import SEEDS, fit_preprocessor, make_mlp, transform


DEFAULT_SOURCE_PROJECT = Path(
    "/Users/lyuguohao/.codex/.chatgpt-projects/"
    "g-p-6aa27e760cf88191b6a228fe3c180b34"
)
SOURCE_PROJECT = Path(os.environ.get("PERSULFIDATION_SOURCE_PROJECT", DEFAULT_SOURCE_PROJECT))
SOURCE_PIPELINE = SOURCE_PROJECT / "outputs/20260911_模型与Q1c推进/model_pipeline.py"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    sys.path.insert(0, str(path.parent))
    spec.loader.exec_module(module)
    return module


def read_fasta(path):
    records, name, parts = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(parts)
            name, parts = line[1:].strip(), []
        else:
            parts.append(line.strip())
    if name is not None:
        records[name] = "".join(parts)
    return records


def median_epochs(path):
    records = json.loads(path.read_text())
    values = defaultdict(list)
    for record in records:
        values[int(record["seed"])].append(int(record["selected_epochs"]))
    return {seed: int(statistics.median(values[seed])) for seed in SEEDS}


class WrappedMLP:
    def __init__(self, imputer, scaler, model):
        self.imputer, self.scaler, self.model = imputer, scaler, model

    def predict_proba(self, X):
        return self.model.predict_proba(transform(self.imputer, self.scaler, X))


def fit_ensemble(X, y, epochs, unlabelled_weight):
    imputer, scaler, Xt = fit_preprocessor(X)
    sample_weight = np.where(y == 1, 1.0, unlabelled_weight)
    models = []
    for seed in SEEDS:
        model = make_mlp(seed)
        for _ in range(epochs[seed]):
            model.partial_fit(
                Xt,
                y,
                classes=np.asarray([0, 1]),
                sample_weight=sample_weight,
            )
        models.append(WrappedMLP(imputer, scaler, model))
    return models


def bool_value(value):
    return str(value).strip().lower() == "true"


def ranks(scores, positions):
    order = np.lexsort((np.asarray(positions), -scores))
    return order, {int(index): rank + 1 for rank, index in enumerate(order)}


def stability(base_rank, pu_rank, hgb_rank, base_sd, pu_sd):
    span = max(base_rank, pu_rank, hgb_rank) - min(base_rank, pu_rank, hgb_rank)
    if max(base_rank, pu_rank) <= 2 and span <= 2 and max(base_sd, pu_sd) <= 0.08:
        return "high"
    if max(base_rank, pu_rank) <= 4 and span <= 4:
        return "moderate"
    return "low"


def main():
    encoder = load_module("candidate_source_pipeline", SOURCE_PIPELINE)
    X, y, _, _, _, _, _ = load_frozen()
    base_epochs = median_epochs(RESULTS / "mlp_corrected_epoch_selection.json")
    pu_epochs = median_epochs(RESULTS / "mlp_corrected_pu_u050_epoch_selection.json")
    base_models = fit_ensemble(X, y, base_epochs, 1.0)
    pu_models = fit_ensemble(X, y, pu_epochs, 0.5)
    source_sequences, structures, _, _ = encoder.read_inputs()
    reverse = defaultdict(list)
    for accession, sequence in source_sequences.items():
        reverse[sequence].append(accession)
    targets = read_fasta(INPUTS / "candidate_proteins.fasta")
    hist_rows = read_csv(INPUTS / "candidate_histgb_scores.csv")
    hist = {(row["protein"], int(row["position"])): row for row in hist_rows}
    all_rows = []
    summaries = []
    for protein, sequence in targets.items():
        if "X" in sequence:
            summaries.append({
                "protein": protein,
                "status": "withheld_unknown_residues",
                "n_cys": sequence.count("C"),
                "model_top1_position": "",
                "model_top1_peptide": "",
                "model_top1_stability": "",
                "model_top2_position": "",
                "model_top2_peptide": "",
                "model_top2_stability": "",
                "assay_priority1_position": "",
                "assay_priority1_peptide": "",
                "assay_priority1_mutant": "",
                "assay_priority2_position": "",
                "assay_priority2_peptide": "",
                "assay_priority2_mutant": "",
                "assay_rule": "sequence contains unknown residues; no silent imputation",
            })
            continue
        matches = sorted(reverse.get(sequence, []))
        usable = [
            accession for accession in matches
            if structures.get(accession, {}).get("status") == "ok"
        ]
        feature_accession = usable[0] if usable else (matches[0] if matches else "")
        by_position = {
            int(row["pos"]): row for row in structures[feature_accession]["cys"]
        } if usable else {}
        positions = [i + 1 for i, aa in enumerate(sequence) if aa == "C"]
        target_X = np.asarray([
            encoder.vector(sequence, position, by_position.get(position), "all_basic_removed")
            for position in positions
        ])
        base_seed = np.stack([model.predict_proba(target_X)[:, 1] for model in base_models])
        pu_seed = np.stack([model.predict_proba(target_X)[:, 1] for model in pu_models])
        base_mean, pu_mean = base_seed.mean(axis=0), pu_seed.mean(axis=0)
        base_order, base_rank = ranks(base_mean, positions)
        pu_order, pu_rank = ranks(pu_mean, positions)
        hgb_scores = np.asarray([float(hist[(protein, position)]["full_fit_annotation_score"]) for position in positions])
        _, hgb_rank = ranks(hgb_scores, positions)
        protein_rows = []
        for index, position in enumerate(positions):
            old = hist[(protein, position)]
            base_sd = float(base_seed[:, index].std(ddof=1))
            pu_sd = float(pu_seed[:, index].std(ddof=1))
            row = {
                "protein": protein,
                "position": position,
                "corrected_mlp_score": float(base_mean[index]),
                "corrected_mlp_seed_sd": base_sd,
                "corrected_mlp_rank": base_rank[index],
                "pu_u050_mlp_score": float(pu_mean[index]),
                "pu_u050_mlp_seed_sd": pu_sd,
                "pu_u050_mlp_rank": pu_rank[index],
                "histgb_score": float(hgb_scores[index]),
                "histgb_rank": hgb_rank[index],
                "rank_span_three_models": max(base_rank[index], pu_rank[index], hgb_rank[index]) - min(base_rank[index], pu_rank[index], hgb_rank[index]),
                "rank_stability": stability(base_rank[index], pu_rank[index], hgb_rank[index], base_sd, pu_sd),
                "structure_present": old["structure_present"],
                "relative_sasa": old["rel_sasa"],
                "plddt": old["plddt"],
                "predicted_disulfide": old["predicted_disulfide"],
                "in_theoretical_digest_window": old["in_theoretical_digest_window"],
                "peptide": old["peptide"],
                "peptide_cys_count": old["peptide_cys_count"],
                "context": old["context"],
                "primary_mutant": f"C{position}S",
                "orthogonal_mutant": f"C{position}A",
                "mutation_warning": "interpret with expression, folding and activity controls",
                "scope": "full-fit experimental prioritization; not independent validation or chemical probability",
            }
            all_rows.append(row)
            protein_rows.append(row)
        model_top = [protein_rows[index] for index in base_order[:2]]
        assay_candidates = sorted(
            [
                row for row in protein_rows
                if bool_value(row["in_theoretical_digest_window"])
                and int(row["peptide_cys_count"]) == 1
            ],
            key=lambda row: (row["corrected_mlp_rank"], row["position"]),
        )
        assay_top = assay_candidates[:2]
        summaries.append({
            "protein": protein,
            "status": "scored",
            "n_cys": len(positions),
            "model_top1_position": model_top[0]["position"],
            "model_top1_peptide": model_top[0]["peptide"],
            "model_top1_stability": model_top[0]["rank_stability"],
            "model_top2_position": model_top[1]["position"] if len(model_top) > 1 else "",
            "model_top2_peptide": model_top[1]["peptide"] if len(model_top) > 1 else "",
            "model_top2_stability": model_top[1]["rank_stability"] if len(model_top) > 1 else "",
            "assay_priority1_position": assay_top[0]["position"] if assay_top else "",
            "assay_priority1_peptide": assay_top[0]["peptide"] if assay_top else "",
            "assay_priority1_mutant": assay_top[0]["primary_mutant"] if assay_top else "",
            "assay_priority2_position": assay_top[1]["position"] if len(assay_top) > 1 else "",
            "assay_priority2_peptide": assay_top[1]["peptide"] if len(assay_top) > 1 else "",
            "assay_priority2_mutant": assay_top[1]["primary_mutant"] if len(assay_top) > 1 else "",
            "assay_rule": "corrected MLP rank restricted to theoretical digest peptides with one Cys; model ranks remain unchanged",
        })
    write_csv(RESULTS / "experimental_candidate_all_cys.csv", all_rows)
    write_csv(RESULTS / "experimental_candidate_top2_and_assay.csv", summaries)
    write_json(RESULTS / "experimental_candidate_run_audit.json", {
        "completed": True,
        "primary_model": "corrected unweighted MLP five-seed full-fit ensemble",
        "sensitivity_model": "corrected PU-cost MLP, unlabelled weight 0.50, five-seed full-fit ensemble",
        "reference_model": "frozen 23-feature HGB full-fit",
        "base_epochs_by_seed": base_epochs,
        "pu_epochs_by_seed": pu_epochs,
        "unknown_residue_policy": "withhold entire protein",
        "candidate_scores_used_for_validation": False,
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "candidate_proteins.fasta",
                INPUTS / "candidate_histgb_scores.csv",
                INPUTS / "frozen_benchmark_data.npz",
                RESULTS / "mlp_corrected_epoch_selection.json",
                RESULTS / "mlp_corrected_pu_u050_epoch_selection.json",
            ]
        },
        "source_encoder": {"path": str(SOURCE_PIPELINE), "sha256": sha256(SOURCE_PIPELINE)},
    })


if __name__ == "__main__":
    main()
