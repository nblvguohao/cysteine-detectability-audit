"""Aggregate LightGBM split-gain importance into feature families.

The purpose is to make visible which kinds of evidence the deployed model uses,
in particular how much of the total gain comes from in-silico tryptic
detectability.  Importance is descriptive: it is computed on the training fit
and is not a causal or per-site explanation.
"""
from __future__ import annotations

import json
import sys
import time
from collections import defaultdict

import numpy as np

from common import RESULTS, ROOT, sha256, write_csv, write_json
import v2_apply
from v2_stack import FEATURES

FAMILY_RULES = (
    ("frozen_23_engineered", lambda name: name.startswith("A:")),
    ("digest_detectability", lambda name: "T:" in name),
    ("trypsin_proximity", lambda name: "nearest_K" in name or "nearest_R" in name),
    ("sequence_window_composition", lambda name: "win" in name and "frac_" in name),
    ("sequence_window_scales", lambda name: "win" in name and "mean_" in name),
    ("sequence_window_counts", lambda name: "win" in name and ("count_" in name or "net_charge" in name)),
    ("neighbour_position_profile", lambda name: "pos_" in name and "_is_" not in name and "index" not in name and "relative" not in name and "terminus" not in name and "term" not in name),
    ("neighbour_position_class", lambda name: "_is_" in name),
    ("cysteine_topology", lambda name: "cys_" in name),
    ("other_residue_distance", lambda name: "nearest_" in name),
    ("position_in_protein", lambda name: any(
        token in name for token in ("pos_index", "pos_from", "pos_min_terminus", "pos_relative")
    )),
    ("protein_level", lambda name: "protein_" in name),
)


def family(name):
    core = name
    for prefix in ("C:rank_", "C:centred_"):
        if core.startswith(prefix):
            core = core[len(prefix):]
            break
    for label, predicate in FAMILY_RULES:
        if predicate(core):
            return label
    return "unassigned"


def transform(name):
    if name.startswith("C:rank_"):
        return "within_protein_rank"
    if name.startswith("C:centred_"):
        return "within_protein_centred"
    return "raw_value"


def main():
    started = time.time()
    rows = []
    audit = {}
    for feature_set in ("full", "chem"):
        booster, names, columns, deployment = v2_apply.fit_deployed_model(feature_set)
        used = [names[i] for i in columns]
        gains = booster.feature_importance(importance_type="gain")
        splits = booster.feature_importance(importance_type="split")
        total_gain = float(gains.sum())
        by_family = defaultdict(float)
        by_transform = defaultdict(float)
        by_family_split = defaultdict(int)
        for name, gain, split in zip(used, gains, splits):
            by_family[family(name)] += float(gain)
            by_transform[transform(name)] += float(gain)
            by_family_split[family(name)] += int(split)
        for label, gain in sorted(by_family.items(), key=lambda item: -item[1]):
            rows.append({
                "feature_set": feature_set,
                "grouping": "family",
                "group": label,
                "gain": gain,
                "gain_share": gain / total_gain,
                "splits": by_family_split[label],
                "n_features": sum(1 for name in used if family(name) == label),
            })
        for label, gain in sorted(by_transform.items(), key=lambda item: -item[1]):
            rows.append({
                "feature_set": feature_set,
                "grouping": "transform",
                "group": label,
                "gain": gain,
                "gain_share": gain / total_gain,
                "splits": sum(
                    int(split) for name, split in zip(used, splits)
                    if transform(name) == label
                ),
                "n_features": sum(1 for name in used if transform(name) == label),
            })
        order = np.argsort(-gains)[:30]
        audit[feature_set] = {
            "rounds": deployment["rounds"],
            "config": deployment["config"],
            "n_features": int(len(used)),
            "total_gain": total_gain,
            "top_30_features": [
                {"feature": used[i], "gain_share": float(gains[i] / total_gain)}
                for i in order
            ],
        }
        print(json.dumps({
            "feature_set": feature_set,
            "detectability_gain_share": round(
                (by_family["digest_detectability"] + by_family["trypsin_proximity"]) / total_gain, 4
            ),
            "seconds": round(time.time() - started, 1),
        }), flush=True)

    write_csv(RESULTS / "v2_feature_importance_families.csv", rows)
    write_json(RESULTS / "v2_feature_importance_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "importance_type": "LightGBM split gain on the full-fit deployed model",
        "interpretation_limit": "descriptive attribution of split gain; not causal and not a per-site explanation",
        "families": audit,
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "lightgbm": v2_apply.members.lgb.__version__,
        },
        "input_hashes": {
            str(FEATURES.relative_to(ROOT) / "v2_features.npz"): sha256(
                FEATURES / "v2_features.npz"
            )
        },
    })
    print(json.dumps(rows[:14], indent=1, default=float))


if __name__ == "__main__":
    main()
