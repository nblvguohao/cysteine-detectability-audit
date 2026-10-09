"""The decisive test for the structural block: does it help on the chemical label?

The same 53 three-dimensional descriptors that added nothing on the annotation
benchmark (`results/structural_block_test.csv`) are computed here for the human
ABPP cysteines and tested against the reactivity label, trained directly on it
with protein-grouped folds.

This is where three-dimensional chemistry should matter if it matters anywhere:
oxidant accessibility, burial, the electrostatic environment of the thiol and
its hydrogen-bond partners are the standard determinants of cysteine
nucleophilicity, and none of them is visible to a sequence window.

Baselines on the identical 18,458 sites and folds
(`results/abpp_reactivity_learnability.csv`, `results/reactivity_model_comparison.csv`):

    probe visibility 10            Spearman 0.0733
    v2 chemistry features 1,046    0.3246
    v2 all features 1,107          0.3426
    no-oxidant intensity alone     0.3033
    v2 all + intensity             0.4675
    ESM-2 sequence only            0.4394
    ESM-2 + intensity              0.5020
"""
from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import build_structural_chemistry_features as builder

PDB_DIR = ROOT / "external/alphafold_human/pdb"
SITES = RESULTS / "abpp_reactivity_sites.csv"
N_SPLITS = 5
REGRESSOR = dict(
    max_iter=400, learning_rate=0.06, max_leaf_nodes=31, min_samples_leaf=20,
    l2_regularization=0.0, early_stopping=False, random_state=0,
)


def build_structural(keys):
    by_protein = {}
    for index, (accession, position) in enumerate(keys):
        by_protein.setdefault(accession, []).append((index, position))
    names = None
    matrix = None
    parsed = 0
    missing_model = 0
    missing_site = 0
    for accession, entries in sorted(by_protein.items()):
        path = PDB_DIR / f"AF-{accession}-F1-model_v6.pdb"
        if not path.exists():
            missing_model += 1
            continue
        structure = builder.read_pdb(path)
        if structure is None:
            missing_model += 1
            continue
        parsed += 1
        for index, position in entries:
            values = builder.site_features(structure, position)
            if values is None:
                missing_site += 1
                continue
            if names is None:
                names = sorted(values)
                matrix = np.full((len(keys), len(names)), np.nan, dtype=np.float32)
            matrix[index] = [values[name] for name in names]
        if parsed % 400 == 0:
            print(json.dumps({"proteins_parsed": parsed}), flush=True)
    return matrix, names, {
        "proteins_parsed": parsed,
        "proteins_without_model": missing_model,
        "sites_without_sulfur_in_model": missing_site,
    }


def grouped_predictions(x, y, groups):
    predictions = np.full(len(y), np.nan)
    for train, test in GroupKFold(n_splits=N_SPLITS).split(x, y, groups):
        model = HistGradientBoostingRegressor(**REGRESSOR)
        model.fit(x[train], y[train])
        predictions[test] = model.predict(x[test])
    return predictions


def main():
    started = time.time()
    rows = read_csv(SITES)
    keys = [(r["accession"], int(r["position"])) for r in rows]
    engagement = np.asarray([float(r["engagement_pervan"]) for r in rows])
    intensity = np.asarray([float(r["log2_no_ox"]) for r in rows])
    groups = np.asarray([r["accession"] for r in rows])
    print(json.dumps({"sites": len(keys), "proteins": int(len(set(groups)))}), flush=True)

    S, S_names, counts = build_structural(keys)
    if S is None:
        raise SystemExit("No human structures parsed")
    covered = np.isfinite(S).any(axis=1)
    print(json.dumps({**counts, "sites_with_structure": int(covered.sum()),
                      "coverage": round(float(covered.mean()), 4),
                      "minutes": round((time.time() - started) / 60, 1)}), flush=True)

    designs = {
        "structure53_only": S,
        "structure53_plus_intensity": np.hstack([S, intensity.reshape(-1, 1)]),
    }
    # the v2 tabular block, so the comparison is like-for-like on these sites
    v2_path = RESULTS / "abpp_reactivity_learnability.csv"
    baseline = {r["design"]: r for r in read_csv(v2_path)}

    low, high = np.quantile(engagement, [0.2, 0.8])
    tail = (engagement <= low) | (engagement >= high)
    tail_label = (engagement[tail] <= low).astype(int)

    summary = []
    for label, x in designs.items():
        predictions = grouped_predictions(x, engagement, groups)
        finite = np.isfinite(predictions)
        rho = stats.spearmanr(predictions[finite], engagement[finite])
        within = []
        order = np.argsort(groups, kind="stable")
        starts = np.flatnonzero(np.r_[True, groups[order][1:] != groups[order][:-1]])
        for block in np.split(order, starts[1:]):
            block = block[np.isfinite(predictions[block])]
            if len(block) >= 5 and np.std(engagement[block]) > 0 and np.std(predictions[block]) > 0:
                value = stats.spearmanr(predictions[block], engagement[block]).statistic
                # a protein whose predictions are constant yields nan; drop it rather
                # than letting it poison the median
                if np.isfinite(value):
                    within.append(float(value))
        row = {
            "design": label,
            "n_features": int(x.shape[1]),
            "n_sites": int(finite.sum()),
            "spearman_out_of_fold": round(float(rho.statistic), 4),
            "spearman_p": f"{rho.pvalue:.3g}",
            "within_protein_spearman_median": round(float(np.median(within)), 4) if within else None,
            "within_protein_n_proteins": int(len(within)),
            "quintile_auc": round(float(roc_auc_score(
                tail_label[np.isfinite(predictions[tail])],
                -predictions[tail][np.isfinite(predictions[tail])])), 4),
        }
        summary.append(row)
        print(f"{label:30s} rho {row['spearman_out_of_fold']:+.4f} "
              f"within {row['within_protein_spearman_median']:+.4f} "
              f"quintile {row['quintile_auc']:.4f}", flush=True)

    for design in ("v2_all_features", "v2_chem_features", "v2_all_plus_intensity",
                   "no_oxidant_intensity_only", "probe_visibility10"):
        if design in baseline:
            row = baseline[design]
            summary.append({
                "design": f"{design} (baseline)",
                "n_features": row["n_features"],
                "n_sites": row["n_sites"],
                "spearman_out_of_fold": row["spearman_out_of_fold"],
                "spearman_p": row["spearman_p"],
                "within_protein_spearman_median": row["within_protein_spearman_median"],
                "within_protein_n_proteins": row["within_protein_n_proteins"],
                "quintile_auc": row["quintile_auc"],
            })

    write_csv(RESULTS / "structural_reactivity_test.csv", summary)
    write_json(RESULTS / "structural_reactivity_test_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": "does three-dimensional chemistry predict the ABPP reactivity label",
        "structural_features": int(S.shape[1]),
        "structural_feature_names": list(S_names),
        "counts": counts,
        "sites_with_structure": int(covered.sum()),
        "structure_coverage": round(float(covered.mean()), 4),
        "splits": f"GroupKFold by protein, {N_SPLITS} folds",
        "regressor": REGRESSOR,
        "results": summary,
        "limits": [
            "AlphaFold single conformer; engagement is measured in a lysate, where conformations differ",
            "no explicit pKa; the electrostatic descriptors are geometric proxies",
            "rows without a model stay missing and LightGBM-style native handling is not available to HistGradientBoosting for the Spearman on uncovered rows, so uncovered rows are dropped from the correlation",
            "engagement is reactivity toward one lysis-derived oxidant, not persulfidation",
        ],
        "versions": {"python": sys.version, "numpy": np.__version__},
        "input_hashes": {
            "results/abpp_reactivity_sites.csv": sha256(SITES),
        },
    })
    print(json.dumps(summary, indent=2, default=str)[:900])


if __name__ == "__main__":
    main()
