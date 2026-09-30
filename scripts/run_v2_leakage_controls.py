"""Negative controls for the v2 feature set.

Control 1, within-protein label permutation.  Annotation labels are permuted
among the cysteines of each protein, which preserves every protein's annotation
count and every feature value.  If any feature carried label information about
which cysteine inside a protein is annotated, the within-protein metrics would
stay above chance after permutation.  They must collapse to the random
baseline.

Control 2, protein-level label permutation.  Whole label vectors are exchanged
between proteins of equal cysteine count, breaking the link between a protein's
own features and its labels.  This is deliberately a weak control and is
reported as such: the swap preserves which cysteine *index* inside a protein
carries an annotation, so features describing a cysteine's rank or position
among the cysteines of its protein remain informative and the global ROC-AUC
stays above chance.  Only control 1 is a clean negative control, and it targets
the within-protein metrics that the project selects models on.

Both controls use the frozen outer folds and select boosting rounds on the
inner fold only, exactly as the real runs do.
"""
from __future__ import annotations

import json
import sys
import time
from collections import defaultdict

import numpy as np

from common import RESULTS, ROOT, sha256, summarize, write_csv, write_json
from v2_stack import FEATURES, fold_partitions, load_v2_matrix
import v2_members as members

SEEDS = (20260913, 20260914, 20260915)
OUTER_FOLD = 1
CONFIG = members.LGB_BIN_GRID[0]
MAX_ROUNDS = 300


def permute_within_protein(y, proteins, seed):
    rng = np.random.default_rng(seed)
    permuted = y.copy()
    order = np.argsort(proteins, kind="stable")
    starts = np.flatnonzero(np.r_[True, proteins[order][1:] != proteins[order][:-1]])
    for group in np.split(order, starts[1:]):
        permuted[group] = rng.permutation(y[group])
    return permuted


def permute_between_proteins(y, proteins, seed):
    """Swap whole label vectors between proteins with the same cysteine count."""
    rng = np.random.default_rng(seed)
    permuted = y.copy()
    order = np.argsort(proteins, kind="stable")
    starts = np.flatnonzero(np.r_[True, proteins[order][1:] != proteins[order][:-1]])
    groups = np.split(order, starts[1:])
    by_size = defaultdict(list)
    for group in groups:
        by_size[len(group)].append(group)
    for size, members_of_size in by_size.items():
        if len(members_of_size) < 2:
            continue
        shuffled = rng.permutation(len(members_of_size))
        for target, source in enumerate(shuffled):
            permuted[members_of_size[target]] = y[members_of_size[source]]
    return permuted


def main():
    started = time.time()
    members.MAX_ROUNDS = MAX_ROUNDS
    members.PATIENCE_BLOCKS = 4
    X, names, y, folds, proteins, components, positions, meta = load_v2_matrix()
    subtrain, validation, outer_train, test, validation_fold = fold_partitions(
        folds, OUTER_FOLD, proteins, components
    )
    rows = []

    def evaluate(label, labels, seed):
        selection, _ = members.lgb_select(
            "bin", X, labels, proteins, components, subtrain, validation,
            grid=(CONFIG,),
        )
        scores, _ = members.lgb_refit(
            "bin", X, labels, proteins, outer_train, test, selection
        )
        summary = summarize(labels[test], scores, proteins[test], components[test])
        summary.update({
            "control": label,
            "permutation_seed": seed,
            "outer_fold": OUTER_FOLD,
            "rounds": selection["rounds"],
            "top1_minus_random": summary["top1"] - summary["random_top1"],
            "top2_minus_random": summary["top2"] - summary["random_top2"],
            "within_protein_auc_minus_half": summary["within_protein_auc_mean"] - 0.5,
            "roc_auc_minus_half": summary["roc_auc"] - 0.5,
        })
        rows.append(summary)
        print(json.dumps({
            "control": label, "seed": seed,
            "roc_auc": round(summary["roc_auc"], 4),
            "within_protein_auc": round(summary["within_protein_auc_mean"], 4),
            "top1": round(summary["top1"], 4),
            "random_top1": round(summary["random_top1"], 4),
            "seconds": round(time.time() - started, 1),
        }), flush=True)

    evaluate("observed_labels", y, 0)
    for seed in SEEDS:
        evaluate("within_protein_permutation", permute_within_protein(y, proteins, seed), seed)
    for seed in SEEDS:
        evaluate("between_protein_permutation", permute_between_proteins(y, proteins, seed), seed)

    within = [row for row in rows if row["control"] == "within_protein_permutation"]
    between = [row for row in rows if row["control"] == "between_protein_permutation"]
    observed = rows[0]
    verdict = {
        "observed_within_protein_auc": observed["within_protein_auc_mean"],
        "permuted_within_protein_auc_mean": float(
            np.mean([row["within_protein_auc_mean"] for row in within])
        ),
        "permuted_within_protein_auc_max": float(
            max(row["within_protein_auc_mean"] for row in within)
        ),
        "permuted_top1_minus_random_max": float(
            max(row["top1_minus_random"] for row in within)
        ),
        "between_protein_permuted_roc_auc_mean": float(
            np.mean([row["roc_auc"] for row in between])
        ),
        "within_protein_label_information_detected": bool(
            max(row["within_protein_auc_mean"] for row in within) > 0.55
        ),
        "primary_control": "within_protein_permutation",
        "interpretation": "a permuted within-protein AUC near 0.5 means no feature encodes which cysteine of a protein is annotated; it does not address ProtBERT pretraining exposure, which the v2 tabular models do not use",
        "between_protein_control_caveat": "the between-protein swap keeps the annotated cysteine index pattern, which position and rank features legitimately capture, so its ROC-AUC staying above chance is expected and is not evidence of label leakage",
    }
    write_csv(RESULTS / "v2_leakage_control_summary.csv", rows)
    write_json(RESULTS / "v2_leakage_control_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "outer_fold": OUTER_FOLD,
        "inner_validation_fold": int(validation_fold),
        "model": "LightGBM binary, fixed configuration",
        "max_rounds": MAX_ROUNDS,
        "permutation_seeds": list(SEEDS),
        "verdict": verdict,
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "lightgbm": members.lgb.__version__,
        },
        "input_hashes": {
            "features/v2_features.npz": sha256(FEATURES / "v2_features.npz")
        },
    })
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
