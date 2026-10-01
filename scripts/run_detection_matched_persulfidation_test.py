"""The decisive test: does chemistry signal survive when detectability is matched?

Every negative set used in this project so far has been "cysteine without an
annotation", which conflates "not modified" with "never seen".  PXD076162
(Schaffer et al., Cell Chem Biol 2026, PMID 42673962) removes that conflation.
Its competitive iodoacetamide-alkyne ABPP experiment quantifies 27,897 human
cysteine sites, each of which was therefore demonstrably visible to mass
spectrometry, in an independent laboratory, cell line and chemistry.

Intersecting those sites with the human persulfidome the project already uses
(PXD044043) gives a negative class that is *detected but not annotated*:

* matched arm   - only cysteines that the ABPP experiment quantified
* unmatched arm - every cysteine of the same proteins, the conventional design

Four models are compared, all fitted on tomato only and applied zero-shot, so
nothing here is trained on human data:

* v2 chem      - the detectability-ablated v2 model
* v2 full      - the full v2 model
* probe digest - the 25 in-silico digest features
* probe visib. - the 10 features that carry only peptide visibility

The contrast between the two arms is the result: if the visibility models lose
their discrimination once detection is matched while the chemistry model keeps
it, that is chemistry signal beyond being seen.  If the visibility models keep
discriminating, the confound survives detection matching.
"""
from __future__ import annotations

import collections
import json
import platform
import re
import sys
import time

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

from common import INPUTS, RESULTS, ROOT, load_frozen, sha256, summarize, write_csv, write_json
import v2_apply
import v2_features as features_module
from run_cross_protease_detectability_probe import (
    FEATURE_NAMES, HGB_KWARGS, PROTEASES, VISIBILITY_ONLY, boundaries_for, features,
)
from v2_stack import FEATURES

ABPP = ROOT / "external/lysis_oxidant_abpp/abpp_detected_cys_sites.json"
MODELS = ROOT / "models_v2"
RULE = "Trypsin"


def load_persulfidome():
    payload = json.loads((INPUTS / "PXD044043_mapped_human.json").read_text())
    proteins = payload["hsa_persulfidation"]["proteins"]
    records = {}
    for entry in proteins:
        records[entry["acc"]] = {
            "sequence": entry["seq"].upper(),
            "observed": {int(p) for p in entry["modified"]},
            "cysteines": {int(p) for p in entry["cys"]},
        }
    return records


def tomato_probe_matrix():
    _, y, _, proteins, _, positions, _ = load_frozen()
    sequences = {
        key: features_module._normalise(value)
        for key, value in features_module._sequences().items()
    }
    rule = PROTEASES[RULE]
    cache = {}
    rows = []
    for accession, position in zip(proteins, positions):
        sequence = sequences[accession]
        if accession not in cache:
            cache[accession] = boundaries_for(sequence, rule)
        rows.append(features(sequence, int(position) - 1, cache[accession], rule))
    return np.asarray(rows, dtype=np.float64), y


def probe_matrix(records, keys):
    rule = PROTEASES[RULE]
    cache = {}
    rows = []
    for accession, position in keys:
        sequence = records[accession]["sequence"]
        if accession not in cache:
            cache[accession] = boundaries_for(sequence, rule)
        rows.append(features(sequence, position - 1, cache[accession], rule))
    return np.asarray(rows, dtype=np.float64)


def subset_columns(subset):
    return [
        i for i, name in enumerate(FEATURE_NAMES)
        if subset == "all" or name in VISIBILITY_ONLY
    ]


def main():
    started = time.time()
    records = load_persulfidome()
    abpp = {(accession, int(position)) for accession, position in json.loads(ABPP.read_text())}

    # the candidate universe: every cysteine of the persulfidome proteins
    keys = []
    detected_flag = []
    labels = []
    mismatched = 0
    for accession, record in sorted(records.items()):
        sequence = record["sequence"]
        for position in sorted(record["cysteines"]):
            if not (1 <= position <= len(sequence)) or sequence[position - 1] != "C":
                mismatched += 1
                continue
            keys.append((accession, position))
            detected_flag.append((accession, position) in abpp)
            labels.append(1 if position in record["observed"] else 0)
    keys = list(keys)
    detected = np.asarray(detected_flag)
    y = np.asarray(labels, dtype=int)
    groups = np.asarray([accession for accession, _ in keys])
    print(json.dumps({
        "universe_sites": int(len(y)),
        "annotated": int(y.sum()),
        "abpp_detected": int(detected.sum()),
        "abpp_detected_and_annotated": int((detected & (y == 1)).sum()),
        "position_residue_mismatches_dropped": mismatched,
    }), flush=True)

    # ---- zero-shot scores from the exported v2 models
    reference_names = np.load(FEATURES / "v2_features.npz", allow_pickle=False)[
        "names"
    ].astype(str)
    encoder = v2_apply.load_source_encoder()
    matrix, matrix_proteins, matrix_positions = v2_apply.build_matrix(
        [
            {"accession": accession, "sequence": record["sequence"]}
            for accession, record in sorted(records.items())
        ],
        encoder,
        reference_names,
    )
    row_of = {
        (protein, int(position)): i
        for i, (protein, position) in enumerate(zip(matrix_proteins, matrix_positions))
    }
    index = np.asarray([row_of[key] for key in keys])
    scores = {}
    model_audit = {}
    for feature_set in ("chem", "full"):
        booster = v2_apply.members.lgb.Booster(
            model_file=str(MODELS / f"v2_lgb_bin_{feature_set}.txt")
        )
        columns = np.asarray(json.loads(
            (MODELS / f"v2_feature_columns_{feature_set}.json").read_text()
        )["column_indices_into_v2_features_npz"])
        scores[f"v2_{feature_set}"] = booster.predict(
            np.ascontiguousarray(matrix[np.ix_(index, columns)])
        )
        model_audit[f"v2_{feature_set}"] = {
            "source": f"models_v2/v2_lgb_bin_{feature_set}.txt",
            "sha256": sha256(MODELS / f"v2_lgb_bin_{feature_set}.txt"),
            "training": "tomato cohort only, zero-shot here",
        }

    # ---- tomato-fitted probe models, so all four models are zero-shot
    tomato_x, tomato_y = tomato_probe_matrix()
    human_x = probe_matrix(records, keys)
    for label, subset in (("probe_digest25", "all"), ("probe_visibility10", "visibility_only")):
        columns = subset_columns(subset)
        model = HistGradientBoostingClassifier(**HGB_KWARGS)
        model.fit(tomato_x[:, columns], tomato_y)
        scores[label] = model.predict_proba(human_x[:, columns])[:, 1]
        model_audit[label] = {
            "source": "fitted here on the frozen tomato cohort",
            "n_features": len(columns),
            "classifier": HGB_KWARGS,
            "training": "tomato cohort only, zero-shot here",
        }
        print(json.dumps({"fitted": label, "n_features": len(columns),
                          "seconds": round(time.time() - started, 1)}), flush=True)

    # ---- arms
    def within_protein_auc(mask):
        order = np.argsort(groups[mask], kind="stable")
        selected_groups = groups[mask][order]
        values = {name: [] for name in scores}
        n_proteins = 0
        start = 0
        while start < len(order):
            stop = start
            while stop < len(order) and selected_groups[stop] == selected_groups[start]:
                stop += 1
            block = np.flatnonzero(mask)[order[start:stop]]
            labels_block = y[block]
            if 0 < labels_block.sum() < len(labels_block):
                n_proteins += 1
                for name, vector in scores.items():
                    values[name].append(roc_auc_score(labels_block, vector[block]))
            start = stop
        return {name: float(np.mean(v)) if v else None for name, v in values.items()}, n_proteins

    rows = []
    arms = {
        "unmatched_all_cysteines": np.ones(len(y), dtype=bool),
        "detection_matched_abpp": detected,
    }
    per_arm = {}
    for arm, mask in arms.items():
        within, n_proteins = within_protein_auc(mask)
        for name, vector in scores.items():
            summary = summarize(y[mask], vector[mask], groups[mask], groups[mask])
            row = {
                "arm": arm,
                "model": name,
                "training_data": "tomato only; zero-shot on human",
                "n_sites": int(mask.sum()),
                "n_annotated": int(y[mask].sum()),
                "annotation_rate": round(float(y[mask].mean()), 4),
                "n_proteins": int(len(set(groups[mask]))),
                "roc_auc": round(float(roc_auc_score(y[mask], vector[mask])), 4),
                "average_precision": round(
                    float(average_precision_score(y[mask], vector[mask])), 4
                ),
                "within_protein_auc": round(within[name], 4) if within[name] else "",
                "within_protein_n_proteins": n_proteins,
                "top1": round(summary["top1"], 4),
                "random_top1": round(summary["random_top1"], 4),
                "mean_reciprocal_rank": round(summary["mean_reciprocal_rank"], 4),
            }
            rows.append(row)
            per_arm.setdefault(arm, {})[name] = row
            print(f"{arm:26s} {name:20s} AUC {row['roc_auc']:.4f} "
                  f"within {row['within_protein_auc']} Top-1 {row['top1']:.4f}", flush=True)

    # paired protein bootstrap for chemistry minus visibility, matched arm only
    def paired_protein_bootstrap(mask, left, right, replicates=5000, seed=20260914):
        proteins_in_arm = np.asarray(sorted(set(groups[mask])))
        by_protein = {}
        for protein in proteins_in_arm:
            block = np.flatnonzero(mask & (groups == protein))
            labels_block = y[block]
            if not (0 < labels_block.sum() < len(labels_block)):
                continue
            by_protein[protein] = (
                roc_auc_score(labels_block, scores[left][block])
                - roc_auc_score(labels_block, scores[right][block])
            )
        values = np.asarray(list(by_protein.values()))
        rng = np.random.default_rng(seed)
        draws = np.asarray([
            values[rng.integers(0, len(values), len(values))].mean()
            for _ in range(replicates)
        ])
        return {
            "n_proteins": int(len(values)),
            "mean_difference": float(values.mean()),
            "percentile_95_interval": [
                float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))
            ],
        }

    comparisons = {}
    for left, right in (
        ("v2_chem", "probe_visibility10"),
        ("v2_chem", "probe_digest25"),
        ("v2_full", "probe_visibility10"),
    ):
        comparisons[f"{left}__minus__{right}__within_protein_detection_matched"] = (
            paired_protein_bootstrap(detected, left, right)
        )

    write_csv(RESULTS / "detection_matched_persulfidation_test.csv", rows)
    write_json(RESULTS / "detection_matched_persulfidation_test_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "question": "does the chemistry signal survive when the negative class is restricted to cysteines proven visible to mass spectrometry",
        "detected_universe": {
            "source": "PXD076162 Table.xlsx, competitive IAA-alkyne ABPP, human PANC-1",
            "paper": "Schaffer et al., Cell Chem Biol 2026, PMID 42673962, doi 10.1016/j.chembiol.2026.08.005",
            "unique_sites_in_that_dataset": len(abpp),
        },
        "label_source": "inputs/PXD044043_mapped_human.json, observed persulfidation sites",
        "arms": {
            "unmatched_all_cysteines": "every cysteine of the persulfidome proteins, the conventional design",
            "detection_matched_abpp": "only cysteines that the ABPP experiment quantified",
        },
        "models": model_audit,
        "results": per_arm,
        "paired_protein_bootstrap": comparisons,
        "limits": [
            "the ABPP experiment used PANC-1 cells and iodoacetamide-alkyne chemistry, so it matches generic mass-spectrometry visibility and not visibility inside the persulfidome experiment",
            "IAA-ABPP detection is itself biased toward reactive and abundant cysteines; that bias applies to both classes",
            "unannotated cysteines remain unannotated, not chemically negative",
            "all four models are zero-shot from tomato, so absolute values understate what a human-trained model would reach",
        ],
        "input_hashes": {
            "external/lysis_oxidant_abpp/abpp_detected_cys_sites.json": sha256(ABPP),
            "inputs/PXD044043_mapped_human.json": sha256(INPUTS / "PXD044043_mapped_human.json"),
            "features/v2_features.npz": sha256(FEATURES / "v2_features.npz"),
        },
        "hardware": {"platform": platform.platform()},
        "versions": {
            "python": sys.version, "numpy": np.__version__,
            "lightgbm": v2_apply.members.lgb.__version__,
        },
    })
    print(json.dumps(comparisons, indent=2))


if __name__ == "__main__":
    main()
