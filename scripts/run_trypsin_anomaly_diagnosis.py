"""Why does visibility-only beat the full digest block on PhosphoLingo Trypsin?

`results/detectability_audit_instrument.csv` records one row where the ten
visibility features score *above* the 25-feature block: Trypsin with unfiltered
negatives, 0.7925 against 0.7629, a visibility share of 1.1126. The confidence
intervals barely overlap, so noise alone is an unlikely explanation, and this is
the only such row among the six unfiltered protease sets.

Three candidate explanations are tested on Zuallaert's own train/test split:

1. a particular feature group inside the 25 hurts, found by adding each group to
   the visibility block and by removing each group from the full block;
2. the harmful features shift between train and test, measured as a standardised
   mean difference per feature;
3. ordinary overfitting, visible as a large gap between in-sample training AUC
   and held-out AUC for the full block but not for the visibility block.

The paired difference between the two score vectors is bootstrapped on the test
set so the reversal itself carries an interval.
"""
from __future__ import annotations

import json
import platform
import sys
import time

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

from common import RESULTS, ROOT, sha256, write_csv, write_json
from run_cross_protease_detectability_probe import (
    FEATURE_NAMES, HGB_KWARGS, PROTEASES, VISIBILITY_ONLY, build, read_fasta,
)

DATA = ROOT / "external/phospholingo"
DATASET = "Trypsin"
VARIANT = "NPF"
REFERENCE = "Chymotrypsin"  # the unfiltered set with the lowest visibility share

GROUPS = {
    "visibility": sorted(VISIBILITY_ONLY),
    "peptide_composition_and_position": [
        "pep_cys_count", "pep_offset_in_peptide", "pep_relative_offset",
        "pep_is_protein_n_terminal", "pep_is_protein_c_terminal",
    ],
    "missed_cleavage_lengths": [
        "pep_mc1_left_len", "pep_mc1_right_len", "pep_mc1_min_len",
        "pep_mc2_left_len", "pep_mc2_right_len", "pep_mc2_min_len",
    ],
    "cleavage_context": [
        "pep_n_peptides_in_protein", "pep_cleavage_sites_within_20",
        "pep_dist_to_peptide_n_term", "pep_dist_to_peptide_c_term",
    ],
}


def check_group_cover():
    covered = sorted({name for names in GROUPS.values() for name in names})
    if covered != sorted(FEATURE_NAMES):
        missing = sorted(set(FEATURE_NAMES) - set(covered))
        extra = sorted(set(covered) - set(FEATURE_NAMES))
        raise SystemExit(f"Group cover mismatch. missing={missing} extra={extra}")


def columns_for(names):
    index = {name: i for i, name in enumerate(FEATURE_NAMES)}
    return [index[name] for name in names]


def load(dataset, variant):
    folder = "single_protease" if dataset in PROTEASES else "."
    base = DATA / folder / dataset / variant if folder != "." else DATA / dataset / variant
    train = read_fasta(base / "train.fasta")
    test = read_fasta(base / "test.fasta")
    return train, test, base


def fit_and_score(x_train, y_train, x_test, y_test):
    model = HistGradientBoostingClassifier(**HGB_KWARGS)
    model.fit(x_train, y_train)
    test_scores = model.predict_proba(x_test)[:, 1]
    train_scores = model.predict_proba(x_train)[:, 1]
    return {
        "test_roc_auc": float(roc_auc_score(y_test, test_scores)),
        "test_average_precision": float(average_precision_score(y_test, test_scores)),
        "train_roc_auc": float(roc_auc_score(y_train, train_scores)),
    }, test_scores


def paired_bootstrap(y, a, b, replicates=2000, seed=20260914):
    """Interval for AUC(a) - AUC(b) on the same test sites."""
    rng = np.random.default_rng(seed)
    draws = np.empty(replicates)
    n = len(y)
    for i in range(replicates):
        index = rng.integers(0, n, n)
        yy = y[index]
        if yy.min() == yy.max():
            draws[i] = np.nan
            continue
        draws[i] = roc_auc_score(yy, a[index]) - roc_auc_score(yy, b[index])
    draws = draws[np.isfinite(draws)]
    return {
        "difference": float(roc_auc_score(y, a) - roc_auc_score(y, b)),
        "percentile_95_interval": [
            float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))
        ],
        "replicates": int(len(draws)),
    }


def standardised_shift(x_train, x_test):
    out = []
    for j, name in enumerate(FEATURE_NAMES):
        train_column = x_train[:, j]
        test_column = x_test[:, j]
        pooled = np.sqrt(0.5 * (train_column.var() + test_column.var()))
        out.append({
            "feature": name,
            "train_mean": float(train_column.mean()),
            "test_mean": float(test_column.mean()),
            "standardised_mean_difference": (
                float((test_column.mean() - train_column.mean()) / pooled)
                if pooled > 1e-12 else 0.0
            ),
        })
    return sorted(out, key=lambda row: -abs(row["standardised_mean_difference"]))


def main():
    started = time.time()
    check_group_cover()
    rows = []
    audit = {}
    score_cache = {}

    for dataset in (DATASET, REFERENCE):
        train_records, test_records, base = load(dataset, VARIANT)
        x_train_full, y_train = build(train_records, dataset, "all")
        x_test_full, y_test = build(test_records, dataset, "all")

        combinations = [("visibility_only", GROUPS["visibility"])]
        for name, names in GROUPS.items():
            if name == "visibility":
                continue
            combinations.append((f"visibility_plus_{name}", GROUPS["visibility"] + names))
        combinations.append(("all_25", list(FEATURE_NAMES)))
        for name in GROUPS:
            if name == "visibility":
                continue
            combinations.append((
                f"all_minus_{name}",
                [feature for feature in FEATURE_NAMES if feature not in GROUPS[name]],
            ))

        for label, names in combinations:
            keep = columns_for(names)
            result, test_scores = fit_and_score(
                x_train_full[:, keep], y_train, x_test_full[:, keep], y_test
            )
            score_cache[(dataset, label)] = test_scores
            rows.append({
                "dataset": dataset,
                "variant": VARIANT,
                "feature_set": label,
                "n_features": len(keep),
                "n_train_sites": int(len(y_train)),
                "n_test_sites": int(len(y_test)),
                "test_positive_rate": round(float(y_test.mean()), 4),
                "test_roc_auc": round(result["test_roc_auc"], 4),
                "test_average_precision": round(result["test_average_precision"], 4),
                "train_roc_auc": round(result["train_roc_auc"], 4),
                "train_minus_test_auc": round(
                    result["train_roc_auc"] - result["test_roc_auc"], 4
                ),
            })
            print(json.dumps({
                "dataset": dataset, "feature_set": label, "n": len(keep),
                "test_auc": round(result["test_roc_auc"], 4),
                "train_auc": round(result["train_roc_auc"], 4),
                "seconds": round(time.time() - started, 1),
            }), flush=True)

        audit[dataset] = {
            "split_source": str(base.relative_to(ROOT)),
            "hashes": {
                str((base / name).relative_to(ROOT)): sha256(base / name)
                for name in ("train.fasta", "test.fasta")
            },
            "reversal_interval": paired_bootstrap(
                y_test,
                score_cache[(dataset, "visibility_only")],
                score_cache[(dataset, "all_25")],
            ),
            "top_shifted_features": standardised_shift(x_train_full, x_test_full)[:8],
        }

    write_csv(RESULTS / "trypsin_anomaly_diagnosis.csv", rows)

    trypsin = {row["feature_set"]: row for row in rows if row["dataset"] == DATASET}
    baseline = trypsin["visibility_only"]["test_roc_auc"]
    harmful = sorted(
        (
            (name, trypsin[f"visibility_plus_{name}"]["test_roc_auc"] - baseline)
            for name in GROUPS if name != "visibility"
        ),
        key=lambda item: item[1],
    )
    verdict = {
        "reversal_confirmed": trypsin["visibility_only"]["test_roc_auc"]
        > trypsin["all_25"]["test_roc_auc"],
        "reversal_interval_excludes_zero": (
            audit[DATASET]["reversal_interval"]["percentile_95_interval"][0] > 0
        ),
        "group_effect_on_visibility_block": {name: round(delta, 4) for name, delta in harmful},
        "most_harmful_group": harmful[0][0],
        "train_minus_test_gap_visibility": trypsin["visibility_only"]["train_minus_test_auc"],
        "train_minus_test_gap_all_25": trypsin["all_25"]["train_minus_test_auc"],
        "reference_dataset_reverses_too": (
            [row for row in rows if row["dataset"] == REFERENCE
             and row["feature_set"] == "visibility_only"][0]["test_roc_auc"]
            > [row for row in rows if row["dataset"] == REFERENCE
               and row["feature_set"] == "all_25"][0]["test_roc_auc"]
        ),
    }
    write_json(RESULTS / "trypsin_anomaly_diagnosis_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "question": "why the ten visibility features beat the 25-feature digest block on PhosphoLingo Trypsin with unfiltered negatives",
        "datasets": [DATASET, REFERENCE],
        "variant": VARIANT,
        "groups": GROUPS,
        "classifier": HGB_KWARGS,
        "splits": "Zuallaert's own train.fasta and test.fasta, unchanged",
        "verdict": verdict,
        "per_dataset": audit,
        "limits": [
            "one classifier family and one hyperparameter setting",
            "standardised mean difference describes marginal shift only, not joint shift",
            "in-sample training AUC is a coarse overfitting proxy on an imbalanced set",
        ],
        "hardware": {"platform": platform.platform()},
        "versions": {"python": sys.version, "numpy": np.__version__},
    })
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
