"""Finish the QTRP round: uncertainty intervals, and the annotation models on this label.

`reports/QTRP_PERSULFIDATION_TEST.md` reported point estimates only, on 272 and
80 positives.  At that size the two headline claims - struct53 discriminates,
pure visibility does not - cannot be read without intervals, and the replication
arm's 0.3617 / 0.4295 needs a test of whether it is really below chance rather
than noise around it.  This script adds the intervals and the one missing arm of
the comparison.

Two things are computed, both predeclared here before the run:

1. **Protein cluster bootstrap** over every (arm, negative definition, design):
   the AUC with a 95 percent percentile interval, plus the predeclared paired
   differences below.  The unit is the protein, matching every other external
   analysis in this project.  Out-of-fold scores are held fixed, so the estimand
   is "paired difference conditional on the fixed protein-grouped fold
   predictions", the same convention as `common.paired_component_bootstrap`.

2. **The annotation-trained models on the QTRP label.**  Every model the project
   actually deploys was fitted on the tomato capture-workflow annotation, and it
   has never been evaluated against a persulfidation label with a defensible
   negative set.  `results/detection_matched_persulfidation_test.csv` ran the
   same four models against an *annotation* label with detection-matched
   negatives; this arm runs them against a *persulfidation* label with same-run
   negatives, so the two tables can be read side by side.  All four are
   zero-shot: fitted on tomato, applied to human, no QTRP label involved.

Predeclared paired differences (in-domain designs, fitted here with
protein-grouped folds):

* struct53 - probe_visibility10          the headline contrast
* struct53 - probe_digest25              structure against detectability
* struct53_plus_digest25 - struct53      what the digest block adds
* probe_digest25 - probe_visibility10    whether cleavage context adds anything

Predeclared paired differences involving the zero-shot annotation-trained models:

* v2_chem_zeroshot - v2_full_zeroshot             the project's two rankers
* v2_chem_zeroshot - probe_visibility10_zeroshot
* struct53 - v2_chem_zeroshot             in-domain chemistry against the
                                          deployed chemistry ranker; not a fair
                                          comparison, since one is fitted on
                                          this label and one is not, and
                                          reported only to size the gap

Predeclared reading rule: an AUC interval whose lower bound stays above 0.5 is
called discrimination; one that contains 0.5 is called undetermined at this
sample size; one whose upper bound stays below 0.5 is called inverted.  No
threshold is chosen after seeing the numbers.
"""
from __future__ import annotations

import json
import platform
import sys
import time

import numpy as np
from scipy.stats import rankdata
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import v2_apply
from run_qtrp_persulfidation_test import (
    ARMS, SITES, grouped_scores, load_proteome, probe_matrix, structural_matrix,
)
from run_qtrp_negative_cleaning_test import COMPILED, all_arm_labels, arm_sites_raw
from run_detection_matched_persulfidation_test import tomato_probe_matrix
from run_cross_protease_detectability_probe import (
    FEATURE_NAMES, HGB_KWARGS, VISIBILITY_ONLY,
)
from v2_stack import FEATURES

MODELS = ROOT / "models_v2"
SCRIPT = ROOT / "scripts" / "run_qtrp_intervals_and_cross_label.py"
DETECTION_MATCHED = RESULTS / "detection_matched_persulfidation_test.csv"
REPLICATES = 5000
SEED = 20260915
MIN_SITES = 100
MIN_CLASS = 20

IN_DOMAIN_DIFFS = (
    ("struct53", "probe_visibility10"),
    ("struct53", "probe_digest25"),
    ("struct53_plus_digest25", "struct53"),
    ("probe_digest25", "probe_visibility10"),
)
ZERO_SHOT_DIFFS = (
    ("v2_chem_zeroshot", "v2_full_zeroshot"),
    ("v2_chem_zeroshot", "probe_visibility10_zeroshot"),
    ("struct53", "v2_chem_zeroshot"),
)
ALL_DIFFS = IN_DOMAIN_DIFFS + ZERO_SHOT_DIFFS
N_FEATURES = {
    "struct53": 53,
    "probe_visibility10": 10,
    "probe_digest25": 25,
    "struct53_plus_digest25": 78,
    "v2_chem_zeroshot": 1046,
    "v2_full_zeroshot": 1107,
    "probe_digest25_zeroshot": 25,
    "probe_visibility10_zeroshot": 10,
}


def fast_auc(y, scores):
    """Mann-Whitney AUC with average ranks for ties."""
    n1 = int(y.sum())
    n0 = int(len(y)) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(scores, method="average")
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def verdict(bounds):
    low, high = bounds
    if low > 0.5:
        return "discriminates"
    if high < 0.5:
        return "inverted"
    return "undetermined_at_this_sample_size"


def cluster_bootstrap(y, groups, score_map, diffs, replicates=REPLICATES, seed=SEED):
    unique = np.unique(groups)
    index_of = {protein: np.flatnonzero(groups == protein) for protein in unique}
    names = list(score_map)
    draws = {name: [] for name in names}
    diff_draws = {f"{a}__minus__{b}": [] for a, b in diffs}
    rng = np.random.default_rng(seed)
    skipped = 0
    for _ in range(replicates):
        picked = rng.integers(0, len(unique), len(unique))
        index = np.concatenate([index_of[unique[i]] for i in picked])
        yy = y[index]
        total = int(yy.sum())
        if total == 0 or total == len(yy):
            skipped += 1
            continue
        values = {name: fast_auc(yy, score_map[name][index]) for name in names}
        for name in names:
            draws[name].append(values[name])
        for left, right in diffs:
            diff_draws[f"{left}__minus__{right}"].append(values[left] - values[right])

    def interval(sample):
        array = np.asarray(sample, dtype=float)
        return [float(np.quantile(array, 0.025)), float(np.quantile(array, 0.975))]

    return {
        "replicates": int(replicates),
        "usable_replicates": int(replicates - skipped),
        "seed": int(seed),
        "unit": "protein",
        "estimand": "conditional on the fixed protein-grouped out-of-fold scores",
        "auc": {
            name: {
                "point": fast_auc(y, score_map[name]),
                "percentile_95_interval": interval(draws[name]),
            }
            for name in names
        },
        "differences": {
            key: {
                "mean_difference": float(np.mean(sample)),
                "percentile_95_interval": interval(sample),
            }
            for key, sample in diff_draws.items()
        },
    }


def zero_shot_scores(keys, sequences):
    """Score the given sites with the four tomato-fitted annotation models."""
    accessions = sorted({accession for accession, _ in keys if accession in sequences})
    records = [{"accession": a, "sequence": sequences[a]} for a in accessions]
    reference_names = np.load(FEATURES / "v2_features.npz", allow_pickle=False)[
        "names"
    ].astype(str)
    encoder = v2_apply.load_source_encoder()
    matrix, matrix_proteins, matrix_positions = v2_apply.build_matrix(
        records, encoder, reference_names
    )
    row_of = {
        (protein, int(position)): i
        for i, (protein, position) in enumerate(zip(matrix_proteins, matrix_positions))
    }
    missing = [key for key in keys if key not in row_of]
    if missing:
        raise RuntimeError(f"{len(missing)} QTRP sites have no v2 feature row, e.g. {missing[:3]}")
    index = np.asarray([row_of[key] for key in keys])

    scores, audit = {}, {}
    for feature_set in ("chem", "full"):
        booster = v2_apply.members.lgb.Booster(
            model_file=str(MODELS / f"v2_lgb_bin_{feature_set}.txt")
        )
        columns = np.asarray(json.loads(
            (MODELS / f"v2_feature_columns_{feature_set}.json").read_text()
        )["column_indices_into_v2_features_npz"])
        scores[f"v2_{feature_set}_zeroshot"] = booster.predict(
            np.ascontiguousarray(matrix[np.ix_(index, columns)])
        )
        audit[f"v2_{feature_set}_zeroshot"] = {
            "source": f"models_v2/v2_lgb_bin_{feature_set}.txt",
            "sha256": sha256(MODELS / f"v2_lgb_bin_{feature_set}.txt"),
            "n_features": int(len(columns)),
            "training": "frozen tomato annotation cohort only; zero-shot here",
            "structural_columns": (
                "the five frozen structural columns are passed as missing, as in "
                "every external application of these models"
            ),
        }

    tomato_x, tomato_y = tomato_probe_matrix()
    human_x, _ = probe_matrix(keys, sequences)
    visibility_columns = [i for i, n in enumerate(FEATURE_NAMES) if n in VISIBILITY_ONLY]
    for label, columns in (
        ("probe_digest25_zeroshot", list(range(len(FEATURE_NAMES)))),
        ("probe_visibility10_zeroshot", visibility_columns),
    ):
        model = HistGradientBoostingClassifier(**HGB_KWARGS)
        model.fit(tomato_x[:, columns], tomato_y)
        scores[label] = model.predict_proba(human_x[:, columns])[:, 1]
        audit[label] = {
            "source": "fitted here on the frozen tomato cohort",
            "n_features": len(columns),
            "classifier": HGB_KWARGS,
            "training": "frozen tomato annotation cohort only; zero-shot here",
        }
    return scores, audit


def main():
    started = time.time()
    rows = read_csv(SITES)
    sequences = load_proteome()
    visibility_columns = [i for i, n in enumerate(FEATURE_NAMES) if n in VISIBILITY_ONLY]
    pooled_ssh, _ = all_arm_labels(rows)
    compiled_ssh = {
        (r["accession"], int(r["site"]))
        for r in read_csv(COMPILED) if r["compiled_ssh"] == "1"
    }

    # ---- one pass over every site that any arm can use
    arm_sets = {}
    for table, ph, description in ARMS:
        positives, negatives_raw = arm_sites_raw(rows, table, ph)
        if not negatives_raw:
            continue
        negatives_arm_clean = negatives_raw - pooled_ssh
        arm_sets[f"{table}_pH{ph}"] = {
            "description": description,
            "positives": positives,
            "variants": {
                "raw": negatives_raw,
                "arm_clean": negatives_arm_clean,
                "full_clean": negatives_arm_clean - compiled_ssh,
            },
        }

    union_keys = sorted({
        key
        for arm in arm_sets.values()
        for key in set(arm["positives"]) | arm["variants"]["raw"]
    })
    residue_ok = np.asarray([
        bool(sequences.get(accession))
        and 1 <= position <= len(sequences[accession])
        and sequences[accession][position - 1] == "C"
        for accession, position in union_keys
    ])
    S_all, _, structure_stats = structural_matrix(union_keys, sequences)
    P_all, probe_ok = probe_matrix(union_keys, sequences)
    covered_all = np.isfinite(S_all).any(axis=1) & residue_ok & probe_ok
    covered_keys = [key for key, keep in zip(union_keys, covered_all) if keep]
    row_of = {key: i for i, key in enumerate(union_keys)}
    print(json.dumps({
        "union_sites": len(union_keys), "usable": int(covered_all.sum()),
        **structure_stats, "seconds": round(time.time() - started, 1),
    }), flush=True)

    zs_all, zero_shot_audit = zero_shot_scores(covered_keys, sequences)
    zs_of = {key: i for i, key in enumerate(covered_keys)}
    print(json.dumps({"zero_shot_models": sorted(zs_all),
                      "seconds": round(time.time() - started, 1)}), flush=True)

    summary, site_rows, intervals, audit_arms = [], [], {}, {}
    for arm, payload in arm_sets.items():
        audit_arms[arm] = {
            "description": payload["description"],
            "positives": len(payload["positives"]),
            "negative_counts": {k: len(v) for k, v in payload["variants"].items()},
            "variants": {},
        }
        for variant, negatives in payload["variants"].items():
            keys = [
                key for key in sorted(payload["positives"]) + sorted(negatives)
                if covered_all[row_of[key]]
            ]
            labels = np.asarray([
                1 if key in payload["positives"] else 0 for key in keys
            ], dtype=int)
            if len(keys) < MIN_SITES or labels.sum() < MIN_CLASS \
                    or (labels == 0).sum() < MIN_CLASS:
                audit_arms[arm]["variants"][variant] = {"skipped": "too few usable sites"}
                print(f"{arm} {variant}: skipped", flush=True)
                continue
            index = np.asarray([row_of[key] for key in keys])
            groups = np.asarray([accession for accession, _ in keys])
            designs = {
                "struct53": S_all[index],
                "probe_visibility10": P_all[np.ix_(index, visibility_columns)],
                "probe_digest25": P_all[index],
                "struct53_plus_digest25": np.hstack([S_all[index], P_all[index]]),
            }
            score_map = {}
            for name, x in designs.items():
                scores = grouped_scores(np.asarray(x, dtype=np.float64), labels, groups)
                if not np.isfinite(scores).all():
                    raise RuntimeError(f"Non-finite out-of-fold scores: {arm} {variant} {name}")
                score_map[name] = scores
            for name, vector in zs_all.items():
                picked = np.asarray([vector[zs_of[key]] for key in keys])
                if not np.isfinite(picked).all():
                    raise RuntimeError(f"Non-finite zero-shot scores: {arm} {variant} {name}")
                score_map[name] = picked

            boot = cluster_bootstrap(labels, groups, score_map, ALL_DIFFS)
            intervals[f"{arm}|{variant}"] = boot
            audit_arms[arm]["variants"][variant] = {
                "n_sites": len(keys),
                "n_positive": int(labels.sum()),
                "n_proteins": int(len(set(groups))),
            }
            for name, vector in score_map.items():
                entry = boot["auc"][name]
                low, high = entry["percentile_95_interval"]
                summary.append({
                    "arm": arm,
                    "arm_description": payload["description"],
                    "negative_definition": variant,
                    "design": name,
                    "fit_regime": "zero_shot_annotation_trained" if name.endswith("_zeroshot")
                                  else "fitted_in_domain_protein_grouped_folds",
                    "n_features": N_FEATURES[name],
                    "n_sites": len(keys),
                    "n_positive": int(labels.sum()),
                    "positive_rate": round(float(labels.mean()), 4),
                    "n_proteins": int(len(set(groups))),
                    "roc_auc": round(entry["point"], 4),
                    "auc_ci_low": round(low, 4),
                    "auc_ci_high": round(high, 4),
                    "verdict": verdict(entry["percentile_95_interval"]),
                    "average_precision": round(
                        float(average_precision_score(labels, vector)), 4),
                })
                print(f"{arm} {variant:10s} {name:28s} AUC {entry['point']:.4f} "
                      f"[{low:.4f}, {high:.4f}] {verdict(entry['percentile_95_interval'])}",
                      flush=True)
            for position, (key, label) in enumerate(zip(keys, labels)):
                row = {
                    "arm": arm, "negative_definition": variant,
                    "accession": key[0], "site": int(key[1]), "label": int(label),
                }
                row.update({
                    name: round(float(vector[position]), 6)
                    for name, vector in score_map.items()
                })
                site_rows.append(row)

    # ---- one table that puts this label next to the annotation label
    merged = []
    for row in summary:
        if row["negative_definition"] == "full_clean":
            continue
        merged.append({
            "experiment": "qtrp_same_run_persulfide_vs_free_thiol",
            "label_semantics": "persulfide form detected in this run vs detected as free thiol only",
            "negative_construction": f"same-run free thiol ({row['negative_definition']})",
            "arm": f"{row['arm']}|{row['negative_definition']}",
            "model_or_design": row["design"],
            "fit_regime": row["fit_regime"],
            "n_sites": row["n_sites"],
            "positive_rate": row["positive_rate"],
            "roc_auc": row["roc_auc"],
            "auc_ci_low": row["auc_ci_low"],
            "auc_ci_high": row["auc_ci_high"],
            "average_precision": row["average_precision"],
            "source": "results/qtrp_interval_summary.csv",
        })
    for row in read_csv(DETECTION_MATCHED):
        merged.append({
            "experiment": "detection_matched_annotation",
            "label_semantics": "persulfidation annotation vs unannotated cysteine",
            "negative_construction": (
                "ABPP-quantified cysteines only" if row["arm"] == "detection_matched_abpp"
                else "every cysteine of the annotated proteins"
            ),
            "arm": row["arm"],
            "model_or_design": row["model"],
            "fit_regime": "zero_shot_annotation_trained",
            "n_sites": int(row["n_sites"]),
            "positive_rate": float(row["annotation_rate"]),
            "roc_auc": float(row["roc_auc"]),
            "auc_ci_low": "",
            "auc_ci_high": "",
            "average_precision": float(row["average_precision"]),
            "source": "results/detection_matched_persulfidation_test.csv",
        })

    write_csv(RESULTS / "qtrp_interval_summary.csv", summary)
    write_csv(RESULTS / "qtrp_site_scores.csv", site_rows)
    write_csv(RESULTS / "qtrp_vs_detection_matched.csv", merged)
    write_json(RESULTS / "qtrp_paired_intervals.json", intervals)
    write_json(RESULTS / "qtrp_intervals_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": (
            "how wide are the QTRP intervals, and do the annotation-trained "
            "models this project deploys predict this persulfidation label"
        ),
        "label": "positive = persulfide form detected in this run; negative = detected as free thiol only, same run",
        "negative_definitions": {
            "raw": "detected as free thiol only in this arm",
            "arm_clean": "raw, minus sites called SSH in any other human QTRP arm",
            "full_clean": "arm_clean, minus sites in the external Free Radic Biol Med 2024 compilation",
        },
        "bootstrap": {
            "unit": "protein",
            "replicates": REPLICATES,
            "seed": SEED,
            "method": "cluster bootstrap over proteins; out-of-fold scores held fixed",
            "reading_rule_predeclared_in_script": (
                "lower bound above 0.5 = discriminates; interval containing 0.5 = "
                "undetermined; upper bound below 0.5 = inverted"
            ),
        },
        "predeclared_differences": {
            "in_domain": [f"{a} - {b}" for a, b in IN_DOMAIN_DIFFS],
            "zero_shot": [f"{a} - {b}" for a, b in ZERO_SHOT_DIFFS],
        },
        "zero_shot_models": zero_shot_audit,
        "union_sites": len(union_keys),
        "usable_sites": int(covered_all.sum()),
        "structure_stats": structure_stats,
        "arms": audit_arms,
        "cross_experiment_table": {
            "path": "results/qtrp_vs_detection_matched.csv",
            "comparability": (
                "the two experiments share the four zero-shot annotation-trained "
                "models and differ in the label; in-domain rows are fitted on the "
                "QTRP label itself and are not comparable with the zero-shot rows"
            ),
        },
        "limits": [
            "the bootstrap conditions on fixed out-of-fold scores, so it does not cover refitting variability",
            "the two arms come from the same paper and cell line, so they are not an independent replication",
            "negatives remain 'no persulfide seen in this run' rather than chemically unmodified",
            "the five frozen structural columns are missing for the zero-shot models, as in every external application",
            "NaHS-treated lysate at pH 5.0 is a strong exogenous sulfide load, not physiological persulfidation",
        ],
        "input_hashes": {
            "results/qtrp_sites_normalised.csv": sha256(SITES),
            "results/frbm_compiled_sites.csv": sha256(COMPILED),
            "features/v2_features.npz": sha256(FEATURES / "v2_features.npz"),
            "scripts/run_qtrp_intervals_and_cross_label.py": sha256(SCRIPT),
        },
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "lightgbm": v2_apply.members.lgb.__version__,
        },
    })
    print(json.dumps(
        {key: value["differences"] for key, value in intervals.items()},
        indent=2)[:2500])


if __name__ == "__main__":
    main()
