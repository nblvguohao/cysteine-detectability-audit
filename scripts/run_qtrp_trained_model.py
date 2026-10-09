"""Train on the persulfidation label instead of scoring it zero-shot.

Everything so far either fitted 53 structural descriptors inside one arm or
applied the annotation models zero-shot.  The open question is whether
training on this label is worth anything: does the project's own chemistry
feature block, fitted on persulfidation, beat both the 53 descriptors and the
deployed zero-shot ranker?  That is the only defensible definition of a "v3"
here, because v2 was fitted on an annotation label the project has since shown
to be detectability-driven.

**The pooled cohort, predeclared.** The two pH 5.0 arms (S1 primary, S2
replication) are pooled so that each site keeps the label from its own run:

* positive - the site's persulfide form was detected in either pH 5.0 arm;
* negative - the site was detected as a free thiol in either pH 5.0 arm and is
  never called SSH in ANY human QTRP arm, including the pH 7.6 arms and the
  S3/S4/S6 tables that have no paired free-thiol sheet.

Pooling is legitimate; mixing experiments is not.  S3 alone carries 3,895 SSH
rows, but it has no paired free-thiol table, so using its positives against
S1/S2 negatives would reintroduce exactly the cross-experiment detection-depth
confound this project keeps finding.  Those rows are therefore used only to
exclude negatives, never to add positives.

`pooled_ph5` is the primary definition.  `pooled_all_arms` additionally admits
the pH 7.6 arms into both classes and is a sensitivity check only: at pH 7.6 the
method barely captures persulfides (21 and 30 positives against about 4,000 free
thiols), so a site seen only there has had little chance to show a persulfide.
`pooled_ph5_full_clean` further removes negatives listed in the external Free
Radic Biol Med 2024 compilation.

**Designs.** All fitted with protein-grouped five folds, the same
HistGradientBoosting configuration as every other probe row:

* struct53                  - the 53 AlphaFold descriptors, the current floor;
* chem1046                  - the v2 chemistry feature block, fitted here on
                              persulfidation rather than on annotation;
* struct53_plus_chem1046    - both;
* struct53_plus_zeroshot    - struct53 plus the deployed zero-shot chemistry
                              score as a single extra column, which involves no
                              QTRP label and so cannot leak;
* probe_digest25, probe_visibility10 - the detectability controls;
* v2_chem_zeroshot          - the deployed model, no training, for reference;
* detection_depth1          - a confound probe: the number of arms in which the
                              site was detected at all.  Pooling gives a site
                              detected in both arms more chances to be called a
                              persulfide, so if this single column discriminates,
                              the pooled design carries a depth confound and the
                              other numbers must be read against it.

Predeclared paired differences:

* chem1046 - struct53                       does the chemistry block add
* struct53_plus_chem1046 - struct53         does the combination add
* struct53_plus_zeroshot - struct53         is the deployed score complementary
* chem1046 - v2_chem_zeroshot               is training on this label worth it
* struct53 - v2_chem_zeroshot               structure against the deployed model
* struct53 - probe_visibility10             the headline contrast
* struct53 - detection_depth1               structure against pooling depth

Predeclared reading rule, identical to the interval round: an AUC interval whose
lower bound stays above 0.5 is discrimination, one containing 0.5 is
undetermined, one whose upper bound is below 0.5 is inverted.
"""
from __future__ import annotations

import collections
import json
import platform
import sys
import time

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import v2_apply
from run_qtrp_persulfidation_test import (
    SITES, grouped_scores, load_proteome, probe_matrix, structural_matrix,
)
from run_qtrp_negative_cleaning_test import COMPILED, all_arm_labels
from run_detection_matched_persulfidation_test import source_probe_matrix
from run_qtrp_intervals_and_cross_label import cluster_bootstrap, fast_auc, verdict
from run_cross_protease_detectability_probe import (
    FEATURE_NAMES, HGB_KWARGS, VISIBILITY_ONLY,
)
from v2_stack import FEATURES

MODELS = ROOT / "models_v2"
SCRIPT = ROOT / "scripts" / "run_qtrp_trained_model.py"
PH5_ARMS = (("S1", "5.0"), ("S2", "5.0"))
PH76_ARMS = (("S1", "7.6"), ("S2", "7.6"))
REPLICATES = 5000
SEED = 20260915

DIFFS = (
    ("chem1046", "struct53"),
    ("struct53_plus_chem1046", "struct53"),
    ("struct53_plus_zeroshot", "struct53"),
    ("chem1046", "v2_chem_zeroshot"),
    ("struct53", "v2_chem_zeroshot"),
    ("struct53", "probe_visibility10"),
    ("struct53", "detection_depth1"),
)


def arm_classes(rows, arms):
    """Sites seen as SSH, and sites seen as SH, in the given arms.

    ``depth`` counts the DISTINCT arms in which the site appears at all, in
    either adduct form.  Counting rows instead would encode the label: every
    persulfide site is also detected as a free thiol in the same arm, so a
    positive contributes at least two rows while a negative contributes one.
    The first version of this probe made that mistake and reported AUC 0.9478
    for a column that was reading the label; see
    reports/QTRP_DEPTH_CONFOUND.md.
    """
    ssh, sh = set(), set()
    seen_arms = collections.defaultdict(set)
    wanted = set(arms)
    for row in rows:
        if row.get("species") != "human":
            continue
        if (row["source_table"], row["ph"]) not in wanted:
            continue
        label = row["class_from_adduct"] or row["class_from_sheet"]
        if label not in ("SH", "SSH"):
            continue
        key = (row["accession"], int(row["site"]))
        (ssh if label == "SSH" else sh).add(key)
        seen_arms[key].add((row["source_table"], row["ph"]))
    depth = collections.Counter({key: len(arms) for key, arms in seen_arms.items()})
    return ssh, sh, depth


def human_matrix(keys, sequences):
    """The v2 feature matrix plus the deployed zero-shot scores for these sites."""
    accessions = sorted({a for a, _ in keys if a in sequences})
    records = [{"accession": a, "sequence": sequences[a]} for a in accessions]
    reference_names = np.load(FEATURES / "v2_features.npz", allow_pickle=False)[
        "names"
    ].astype(str)
    encoder = v2_apply.load_source_encoder()
    matrix, proteins, positions = v2_apply.build_matrix(
        records, encoder, reference_names
    )
    row_of = {(p, int(s)): i for i, (p, s) in enumerate(zip(proteins, positions))}
    missing = [key for key in keys if key not in row_of]
    if missing:
        raise RuntimeError(f"{len(missing)} sites have no v2 feature row, e.g. {missing[:3]}")
    index = np.asarray([row_of[key] for key in keys])

    columns = {}
    scores, audit = {}, {}
    for feature_set in ("chem", "full"):
        columns[feature_set] = np.asarray(json.loads(
            (MODELS / f"v2_feature_columns_{feature_set}.json").read_text()
        )["column_indices_into_v2_features_npz"])
        booster = v2_apply.members.lgb.Booster(
            model_file=str(MODELS / f"v2_lgb_bin_{feature_set}.txt")
        )
        scores[f"v2_{feature_set}_zeroshot"] = booster.predict(
            np.ascontiguousarray(matrix[np.ix_(index, columns[feature_set])])
        )
        audit[f"v2_{feature_set}_zeroshot"] = {
            "source": f"models_v2/v2_lgb_bin_{feature_set}.txt",
            "sha256": sha256(MODELS / f"v2_lgb_bin_{feature_set}.txt"),
            "n_features": int(len(columns[feature_set])),
            "training": "training matrix only; zero-shot here",
        }
    return matrix[index], columns, scores, audit


def main():
    started = time.time()
    rows = read_csv(SITES)
    sequences = load_proteome()
    visibility_columns = [i for i, n in enumerate(FEATURE_NAMES) if n in VISIBILITY_ONLY]
    pooled_ssh_any_arm, _ = all_arm_labels(rows)
    compiled_ssh = {
        (r["accession"], int(r["site"]))
        for r in read_csv(COMPILED) if r["compiled_ssh"] == "1"
    }

    ssh5, sh5, depth5 = arm_classes(rows, PH5_ARMS)
    ssh76, sh76, depth76 = arm_classes(rows, PH76_ARMS)
    depth_all = collections.Counter(depth5)
    depth_all.update(depth76)

    cohorts = {
        "pooled_ph5": {
            "positives": ssh5,
            "negatives": (sh5 - pooled_ssh_any_arm),
            "depth": depth5,
            "note": "primary: both pH 5.0 arms pooled, negatives never SSH in any human arm",
        },
        "pooled_ph5_full_clean": {
            "positives": ssh5,
            "negatives": (sh5 - pooled_ssh_any_arm) - compiled_ssh,
            "depth": depth5,
            "note": "primary minus negatives listed in the external compilation",
        },
        "pooled_all_arms": {
            "positives": ssh5 | ssh76,
            "negatives": ((sh5 | sh76) - pooled_ssh_any_arm),
            "depth": depth_all,
            "note": "sensitivity: pH 7.6 arms admitted, where the method barely captures persulfides",
        },
    }
    print(json.dumps({
        "ssh_ph5": len(ssh5), "sh_ph5": len(sh5),
        "ssh_ph76": len(ssh76), "sh_ph76": len(sh76),
        "pooled_ssh_any_human_arm": len(pooled_ssh_any_arm),
        "cohort_sizes": {k: {"positives": len(v["positives"]), "negatives": len(v["negatives"])}
                         for k, v in cohorts.items()},
    }, indent=2), flush=True)

    # ---- features, once, over the union of every site any cohort can use
    union_keys = sorted({
        key for cohort in cohorts.values()
        for key in set(cohort["positives"]) | set(cohort["negatives"])
    })
    residue_ok = np.asarray([
        bool(sequences.get(a)) and 1 <= p <= len(sequences[a]) and sequences[a][p - 1] == "C"
        for a, p in union_keys
    ])
    S_all, _, structure_stats = structural_matrix(union_keys, sequences)
    P_all, probe_ok = probe_matrix(union_keys, sequences)
    covered = np.isfinite(S_all).any(axis=1) & residue_ok & probe_ok
    keys = [key for key, keep in zip(union_keys, covered) if keep]
    S_all, P_all = S_all[covered], P_all[covered]
    print(json.dumps({"union_sites": len(union_keys), "usable": len(keys),
                      **structure_stats, "seconds": round(time.time() - started, 1)}), flush=True)

    V_all, columns, zero_shot, zero_shot_audit = human_matrix(keys, sequences)
    row_of = {key: i for i, key in enumerate(keys)}
    print(json.dumps({"v2_matrix": list(V_all.shape),
                      "seconds": round(time.time() - started, 1)}), flush=True)

    summary, site_rows, intervals, audit_cohorts = [], [], {}, {}
    for cohort, payload in cohorts.items():
        selected = [
            key for key in keys
            if key in payload["positives"] or key in payload["negatives"]
        ]
        y = np.asarray([1 if key in payload["positives"] else 0 for key in selected])
        index = np.asarray([row_of[key] for key in selected])
        groups = np.asarray([a for a, _ in selected])
        depth = np.asarray([[payload["depth"][key]] for key in selected], dtype=np.float64)

        designs = {
            "struct53": S_all[index],
            "chem1046": V_all[np.ix_(index, columns["chem"])],
            "struct53_plus_chem1046": np.hstack([
                S_all[index], V_all[np.ix_(index, columns["chem"])]]),
            "struct53_plus_zeroshot": np.hstack([
                S_all[index], zero_shot["v2_chem_zeroshot"][index].reshape(-1, 1)]),
            "probe_digest25": P_all[index],
            "probe_visibility10": P_all[np.ix_(index, visibility_columns)],
            "detection_depth1": depth,
        }
        score_map = {}
        for name, x in designs.items():
            score_map[name] = grouped_scores(np.asarray(x, dtype=np.float64), y, groups)
            if not np.isfinite(score_map[name]).all():
                raise RuntimeError(f"Non-finite out-of-fold scores: {cohort} {name}")
        score_map["v2_chem_zeroshot"] = zero_shot["v2_chem_zeroshot"][index]
        score_map["v2_full_zeroshot"] = zero_shot["v2_full_zeroshot"][index]

        boot = cluster_bootstrap(y, groups, score_map, DIFFS,
                                 replicates=REPLICATES, seed=SEED)
        intervals[cohort] = boot
        audit_cohorts[cohort] = {
            "note": payload["note"],
            "n_sites": len(selected),
            "n_positive": int(y.sum()),
            "n_proteins": int(len(set(groups))),
            "sites_detected_in_both_ph5_arms": int((depth == 2).sum()),
        }
        for name in score_map:
            entry = boot["auc"][name]
            low, high = entry["percentile_95_interval"]
            summary.append({
                "cohort": cohort,
                "cohort_note": payload["note"],
                "design": name,
                "fit_regime": "zero_shot_annotation_trained" if name.endswith("_zeroshot")
                              else "fitted_in_domain_protein_grouped_folds",
                "n_features": int(designs[name].shape[1]) if name in designs
                              else int(len(columns["chem" if "chem" in name else "full"])),
                "n_sites": len(selected),
                "n_positive": int(y.sum()),
                "positive_rate": round(float(y.mean()), 4),
                "n_proteins": int(len(set(groups))),
                "roc_auc": round(entry["point"], 4),
                "auc_ci_low": round(low, 4),
                "auc_ci_high": round(high, 4),
                "verdict": verdict(entry["percentile_95_interval"]),
                "average_precision": round(
                    float(average_precision_score(y, score_map[name])), 4),
            })
            print(f"{cohort:22s} {name:24s} AUC {entry['point']:.4f} "
                  f"[{low:.4f}, {high:.4f}] {verdict(entry['percentile_95_interval'])}",
                  flush=True)
        for position, (key, label) in enumerate(zip(selected, y)):
            row = {"cohort": cohort, "accession": key[0], "site": int(key[1]),
                   "label": int(label), "arms_detected": int(payload["depth"][key])}
            row.update({name: round(float(vector[position]), 6)
                        for name, vector in score_map.items()})
            site_rows.append(row)
        for left, right in DIFFS:
            entry = boot["differences"][f"{left}__minus__{right}"]
            print(f"  {left} - {right}: {entry['mean_difference']:+.4f} "
                  f"[{entry['percentile_95_interval'][0]:+.4f}, "
                  f"{entry['percentile_95_interval'][1]:+.4f}]", flush=True)

    write_csv(RESULTS / "qtrp_trained_model_summary.csv", summary)
    write_csv(RESULTS / "qtrp_trained_model_site_scores.csv", site_rows)
    write_json(RESULTS / "qtrp_trained_model_intervals.json", intervals)
    write_json(RESULTS / "qtrp_trained_model_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": (
            "does training on the QTRP persulfidation label beat the 53 structural "
            "descriptors and the deployed zero-shot chemistry ranker"
        ),
        "cohorts": audit_cohorts,
        "cohort_definitions": {
            "positive": "persulfide form detected in the arms admitted by the cohort",
            "negative": "free thiol detected in those arms and never called SSH in any human QTRP arm",
            "why_s3_positives_are_excluded": (
                "S3 carries 3,895 SSH rows but has no paired free-thiol sheet; using its "
                "positives against S1/S2 negatives would reintroduce a cross-experiment "
                "detection-depth confound, so those rows only exclude negatives"
            ),
        },
        "designs": sorted(DIFFS and {d for pair in DIFFS for d in pair} | {"v2_full_zeroshot", "probe_digest25"}),
        "predeclared_differences": [f"{a} - {b}" for a, b in DIFFS],
        "bootstrap": {
            "unit": "protein", "replicates": REPLICATES, "seed": SEED,
            "method": "cluster bootstrap over proteins; out-of-fold scores held fixed",
            "reading_rule_predeclared_in_script": (
                "lower bound above 0.5 = discriminates; interval containing 0.5 = "
                "undetermined; upper bound below 0.5 = inverted"
            ),
        },
        "classifier": HGB_KWARGS,
        "zero_shot_models": zero_shot_audit,
        "structure_stats": structure_stats,
        "limits": [
            "all positives come from one paper and one cell line, so this is not an independent replication",
            "negatives are 'no persulfide seen in any arm of this paper' rather than chemically unmodified",
            "pooling gives sites detected in both arms more chances to be positive; the detection_depth1 probe measures that confound directly",
            "1,046 features against a few hundred positives risks overfitting; the protein-grouped folds measure it but do not remove it",
            "the five frozen structural columns are missing inside chem1046 for human sequences, as in every external application",
            "NaHS-treated lysate at pH 5.0 is a strong exogenous sulfide load, not physiological persulfidation",
        ],
        "input_hashes": {
            "results/qtrp_sites_normalised.csv": sha256(SITES),
            "results/frbm_compiled_sites.csv": sha256(COMPILED),
            "features/v2_features.npz": sha256(FEATURES / "v2_features.npz"),
            "scripts/run_qtrp_trained_model.py": sha256(SCRIPT),
        },
        "versions": {
            "python": sys.version, "numpy": np.__version__,
            "lightgbm": v2_apply.members.lgb.__version__,
        },
    })
    print(json.dumps({k: v["differences"] for k, v in intervals.items()}, indent=2)[:1500])


if __name__ == "__main__":
    main()
