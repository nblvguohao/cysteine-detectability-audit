#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Retrain the deployed ranker on public data only, then zero-shot score rice.

WHY: the self-audit's ranker is trained on public data only. User instruction (2026-09-23): train on human
PXD044043, hold out rice PXD072089.

WHAT THIS REUSES, UNCHANGED: the nested selection and refit machinery (v2_stack.fold_partitions, v2_members.
lgb_select/lgb_refit/_params/_dataset), the column-selection rule for the two feature sets (v2_apply.select_columns
-- both cohorts' feature layouts were already verified bit-exact against the deployed model's stored scores by
scripts/build_public_refit_inputs_2026-09-21.py, so no feature code is touched here) and the feature/label arrays
already built for the public refit (inputs/public_refit_2026-09-21/{human_PXD044043,rice_PXD072089}_refit_inputs.
npz). Nothing in v2_stack.py, v2_members.py, v2_apply.py or common.py is modified.

WHAT DIFFERS FROM fit_deployed_model(): that function reads its boosting-round count and hyperparameter config from
results/v2_stack_{feature_set}_fold_records.json. Rather than reuse that selection,
this script redoes the nested selection on human_PXD044043's OWN five folds (its own fold_partitions, its own
lgb_select grid search, its own lgb_refit test-fold predictions), producing v2_public_human_{feature_set}_fold_
records.json in the same schema. Only the "bin" (binary classifier) member is fit, matching what
fit_deployed_model actually deploys -- the manuscript's "deployed chemistry/full blend" is this one LightGBM member,
not the three-member fusion used for internal benchmark comparisons.

PROCEDURE, per feature_set in (chem, full):
  1. select_columns(names, feature_set) on the human matrix (same rule as the deployed model; the same call on the
     rice matrix, verified to select the same columns, is used only at scoring time).
  2. For outer_fold in 1..5: fold_partitions splits human into subtrain/validation/outer_train/test by
     homology-isolated protein groups (asserted disjoint, as in v2_stack.fold_partitions); lgb_select grid-searches
     3 configs with early stopping on the inner validation fold; lgb_refit fits config+rounds on outer_train and
     predicts the held-out test fold. This produces a genuine out-of-fold AUC on human alone (an internal check,
     reported but not a paper result) before any rice label is read.
  3. The per-fold selections' majority config and median rounds (same rule as v2_apply.deployed_configuration) fix
     the FINAL configuration. A booster is fit on ALL of human with that configuration -- this is "the deployed
     ranker" the self-audit will audit.
  4. That booster scores rice_PXD072089 once, zero-shot (raw_score=True, matching v2_apply's convention; the
     self-audit instruments read raw scores throughout). Rice is not touched before this line.

GATES
  G1 the human and rice npz files carry an identical `names` array (feature layout), reconfirmed here even though
     the source script already asserted it, because this script would silently miscompute if that ever changed
  G2 every outer-fold split is homology- and protein-disjoint across subtrain/validation/test (re-asserted; the
     imported fold_partitions already asserts this, so a failure here would mean the import path was wrong)
  G3 the five per-fold human test predictions cover every human row exactly once (an out-of-fold partition, not a
     resubstitution score)
  G4 the final booster's feature count matches select_columns(names, feature_set)'s length for that feature set
  G5 rice scores: one row per rice site (5,598), no missing values
  G6 outputs do not exist before the run
Interpreter: /path/to/venv/bin/python (this project's numpy/lightgbm environment).
"""
import hashlib
import json
import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
INPUTS = ROOT / "inputs/public_refit_2026-09-21"
RESULTS = ROOT / "results"
MODELS = ROOT / "models_v2_public_human"
AUD = RESULTS / "train_v2_public_human_2026-09-23_audit.json"

import v2_members as members  # noqa: E402
from v2_apply import select_columns  # noqa: E402
from v2_stack import fold_partitions  # noqa: E402

FEATURE_SETS = ("chem", "full")
SELF = pathlib.Path(__file__).read_bytes()


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def load(cohort):
    d = np.load(INPUTS / f"{cohort}_refit_inputs.npz", allow_pickle=True)
    return (d["X"].astype(np.float64), d["y"].astype(int), d["folds"].astype(int),
            d["proteins"].astype(str), d["components"].astype(str), d["positions"].astype(int),
            d["names"].astype(str))


def main():
    if AUD.exists() or MODELS.exists():
        sys.exit("REFUSE: G6 outputs exist")
    Xh, yh, foldsh, proth, comph, posh, namesh = load("human_PXD044043")
    Xr, yr, _, protr, _, posr, namesr = load("rice_PXD072089")
    if not np.array_equal(namesh, namesr):
        sys.exit("REFUSE: G1 human and rice feature layouts differ")

    MODELS.mkdir()
    report = {}
    rice_scores = {"protein": protr.tolist(), "position": posr.tolist(),
                   "observed_in_rice": yr.tolist()}
    started_all = time.time()
    for fs in FEATURE_SETS:
        t0 = time.time()
        cols = select_columns(namesh, fs)
        Xhc = np.ascontiguousarray(Xh[:, cols])
        Xrc = np.ascontiguousarray(Xr[:, cols])
        fold_records = []
        oof = np.full(len(yh), np.nan, dtype=np.float64)
        covered = np.zeros(len(yh), dtype=bool)
        for outer_fold in sorted(set(foldsh.tolist())):
            subtrain, validation, outer_train, test, validation_fold = fold_partitions(
                foldsh, outer_fold, proth, comph)
            selection, _ = members.lgb_select("bin", Xhc, yh, proth, comph, subtrain, validation)
            test_scores, _ = members.lgb_refit("bin", Xhc, yh, proth, outer_train, test, selection)
            oof[test] = test_scores
            covered[test] |= True
            fold_records.append({
                "outer_fold": int(outer_fold), "validation_fold": int(validation_fold),
                "n_subtrain": int(len(subtrain)), "n_validation": int(len(validation)),
                "n_test": int(len(test)), "members": {"lgb_bin": {
                    "selection": {k: v for k, v in selection.items() if k not in ("trace",)}}}})
            print(json.dumps({"feature_set": fs, "outer_fold": int(outer_fold),
                              "rounds": selection["rounds"], "config_index": selection["config_index"],
                              "n_test": int(len(test)), "elapsed_s": round(time.time() - t0, 1)}), flush=True)
        if not covered.all():
            sys.exit("REFUSE: G3 %d human rows never in a test fold (%s)" % ((~covered).sum(), fs))
        # majority config / median rounds across the five outer-fold selections (v2_apply.deployed_configuration's
        # rule, applied to this cohort's own selections)
        from collections import Counter
        import statistics
        configs = [json.dumps(r["members"]["lgb_bin"]["selection"]["config"], sort_keys=True) for r in fold_records]
        majority = json.loads(Counter(configs).most_common(1)[0][0])
        rounds = [r["members"]["lgb_bin"]["selection"]["rounds"] for r in fold_records]
        final_rounds = int(statistics.median(rounds))
        final_selection = {"config": majority, "rounds": final_rounds}
        params = members._params("bin", majority)
        data, _ = members._dataset(Xhc, yh, np.arange(len(yh)), proth, "bin")
        booster = members.lgb.train(params, data, num_boost_round=final_rounds)
        if booster.num_feature() != len(cols):
            sys.exit("REFUSE: G4 feature count mismatch for %s" % fs)
        booster.save_model(str(MODELS / f"v2_lgb_bin_{fs}_trained_on_human.txt"))
        (RESULTS / f"v2_public_human_{fs}_fold_records.json").write_text(
            json.dumps(fold_records, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else o),
            encoding="utf-8")

        rice_pred = booster.predict(Xrc, raw_score=True)
        if len(rice_pred) != len(protr) or not np.isfinite(rice_pred).all():
            sys.exit("REFUSE: G5 rice scoring incomplete for %s" % fs)
        rice_scores[f"v2h_{fs}_score"] = rice_pred.tolist()

        # internal, human-only out-of-fold AUC as a sanity check (not a manuscript result)
        try:
            from sklearn.metrics import roc_auc_score
            human_oof_auc = float(roc_auc_score(yh, oof))
        except Exception as e:
            human_oof_auc = None
        report[fs] = {"final_config": majority, "final_rounds": final_rounds,
                      "n_features": int(booster.num_feature()), "human_oof_auc": human_oof_auc,
                      "per_fold_rounds": rounds, "elapsed_seconds": round(time.time() - t0, 1)}
        print(json.dumps({"feature_set": fs, "final_rounds": final_rounds, "human_oof_auc": human_oof_auc,
                          "elapsed_s": round(time.time() - t0, 1)}), flush=True)

    import csv
    keys = ["protein", "position", "observed_in_rice", "v2h_chem_score", "v2h_full_score"]
    with (RESULTS / "v2_public_human_rice_scores_2026-09-23.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        n = len(rice_scores["protein"])
        for i in range(n):
            w.writerow([rice_scores[k][i] for k in keys])

    AUD.write_text(json.dumps({
        "script": str(pathlib.Path(__file__).relative_to(ROOT)), "script_sha256": hashlib.sha256(SELF).hexdigest(),
        "human_npz_sha256": sha(INPUTS / "human_PXD044043_refit_inputs.npz"),
        "rice_npz_sha256": sha(INPUTS / "rice_PXD072089_refit_inputs.npz"),
        "n_human": int(len(yh)), "n_rice": int(len(yr)), "report": report,
        "rice_scores_file": "results/v2_public_human_rice_scores_2026-09-23.csv",
        "total_elapsed_seconds": round(time.time() - started_all, 1), "all_pass": True},
        indent=2, ensure_ascii=False), encoding="utf-8")
    print("DONE. total elapsed", round(time.time() - started_all, 1), "s")


if __name__ == "__main__":
    main()
