"""Step 2 (POST HOC revision analysis, item D_empirical_background): empirical detectability model on the
PXD063463 global (unenriched) arm, and its aggregation to the cysteines of the trypsin ABE arm.

The primary model repeats the original registered check (scripts/run_phase3_empirical_detectability_2026-09-22.py
in the internal repository): logistic regression (scikit-learn defaults, max_iter 2000) on the 12 standardised
features of the stored candidate table, five folds grouped by protein with the original fold assignment
(numpy default_rng(20260922 + 1) permutation of the sorted protein list, fold = rank mod 5), out-of-fold (OOF)
predictions. It is reused here, post hoc, as a detectability score; nothing about it was pre-specified for
this purpose.

Variants (all with the same folds, OOF):
  lr_full      12 features (length, missed_cleavages, mass, gravy, frac_basic, frac_acidic, dist_n_term,
               dist_c_term, n_cys, log_abundance, abundance_missing, abe_identified) = original M3
  lr_peptide   the 9 peptide features only (no protein abundance, abundance-missing or ABE-identified flags)
  lr_geom      length + missed_cleavages = original M2
  gbm_full     LightGBM on the 12 features (fixed configuration, no tuning, deterministic)
Cysteine-level score = max OOF probability over the candidates that contain the cysteine (0 if none).

Outputs (results/D_empirical_background/)
  candidate_oof_scores.csv.gz   one row per candidate: cand_row, protein, fold, empirical_detected and the OOF scores
  cys_scores_trypsin_arm.csv    cysteine table of step 1 plus the aggregated scores
  empirical_model_summary.json  AUCs (pooled OOF and per fold) against the stored values
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import warnings

import lightgbm as lgb
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from common import IN_CAND, IN_EMP_SUMMARY, MASTER_SEED, N_JOBS, RESULTS, dump_json
from cys_audit.stats import auc as tie_auc  # tie-corrected rank AUC

import json

FEATS = ["length", "missed_cleavages", "mass", "gravy", "frac_basic", "frac_acidic", "dist_n_term", "dist_c_term",
         "n_cys", "log_abundance", "abundance_missing", "abe_identified"]
PEPT = ["length", "missed_cleavages", "mass", "gravy", "frac_basic", "frac_acidic", "dist_n_term", "dist_c_term",
        "n_cys"]
GEOM = ["length", "missed_cleavages"]
ORIG_SEED = 20260922
CV_FOLDS = 5
GBM_PARAMS = dict(objective="binary", n_estimators=400, learning_rate=0.05, num_leaves=31, min_child_samples=100,
                  subsample=1.0, colsample_bytree=1.0, reg_lambda=1.0, random_state=MASTER_SEED, n_jobs=N_JOBS,
                  deterministic=True, force_row_wise=True, verbose=-1)


def fold_ids(proteins):
    uniq = sorted(set(proteins.tolist()))
    rng = np.random.default_rng(ORIG_SEED + 1)
    perm = rng.permutation(len(uniq))
    f = {uniq[perm[i]]: i % CV_FOLDS for i in range(len(uniq))}
    return np.array([f[p] for p in proteins])


def oof_lr(X, y, fold):
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd  # standardised on the full table, as in the original script
    oof = np.zeros(len(y))
    per_fold = []
    for k in range(CV_FOLDS):
        te = fold == k
        clf = LogisticRegression(max_iter=2000)
        clf.fit(Xs[~te], y[~te])
        oof[te] = clf.predict_proba(Xs[te])[:, 1]
        per_fold.append(tie_auc(oof[te], y[te]))
    return oof, per_fold


def oof_gbm(X, y, fold):
    oof = np.zeros(len(y))
    per_fold = []
    for k in range(CV_FOLDS):
        te = fold == k
        clf = lgb.LGBMClassifier(**GBM_PARAMS)
        clf.fit(X[~te], y[~te])
        oof[te] = clf.predict_proba(X[te])[:, 1]
        per_fold.append(tie_auc(oof[te], y[te]))
    return oof, per_fold


def main():
    warnings.filterwarnings("ignore", category=UserWarning)
    C = pd.read_csv(IN_CAND, dtype={"protein": str})
    C["cand_row"] = np.arange(len(C))
    y = C.empirical_detected.to_numpy().astype(int)
    fold = fold_ids(C.protein.to_numpy())
    stored = json.load(open(IN_EMP_SUMMARY, encoding="utf-8"))
    out = {"n_candidates": int(len(C)), "n_positive": int(y.sum()), "n_proteins_with_candidates": int(C.protein.nunique()),
           "stored_summary": {k: stored[k] for k in ("n_candidates", "n_positive", "M1_theoretical_flag_auc",
                                                      "M2_geometry_only_auc", "M3_full_model_auc", "M3_per_fold_auc")}}
    out["M1_theoretical_flag_auc"] = tie_auc(C.theoretical_detectable.to_numpy().astype(float), y)
    scores = {"cand_row": C.cand_row.to_numpy(), "protein": C.protein.to_numpy(), "fold": fold, "empirical_detected": y}
    for name, cols in (("lr_full", FEATS), ("lr_peptide", PEPT), ("lr_geom", GEOM)):
        oof, pf = oof_lr(C[cols].to_numpy(dtype=float), y, fold)
        scores[name] = oof
        out[f"{name}_oof_auc"] = tie_auc(oof, y)
        out[f"{name}_per_fold_auc"] = pf
        print(name, round(out[f"{name}_oof_auc"], 4), [round(a, 4) for a in pf], flush=True)
    oof, pf = oof_gbm(C[FEATS].to_numpy(dtype=float), y, fold)
    scores["gbm_full"] = oof
    out["gbm_full_oof_auc"] = tie_auc(oof, y)
    out["gbm_full_per_fold_auc"] = pf
    out["gbm_params"] = GBM_PARAMS
    print("gbm_full", round(out["gbm_full_oof_auc"], 4), [round(a, 4) for a in pf], flush=True)
    S = pd.DataFrame(scores)
    S.to_csv(f"{RESULTS}/candidate_oof_scores.csv.gz", index=False, compression="gzip", float_format="%.8g")

    # --- aggregate to the cysteines of the trypsin ABE arm (max over containing candidates, 0 if none)
    M = pd.read_csv(f"{RESULTS}/cand_cys_map.csv", dtype={"protein": str})
    M = M.merge(S[["cand_row", "lr_full", "lr_peptide", "lr_geom", "gbm_full"]], on="cand_row", how="left")
    agg = M.groupby(["protein", "position"])[["lr_full", "lr_peptide", "lr_geom", "gbm_full"]].max().reset_index()
    agg.columns = ["protein", "position", "score_lr_full", "score_lr_peptide", "score_lr_geom", "score_gbm_full"]
    T = pd.read_csv(f"{RESULTS}/cys_table_trypsin_arm.csv", dtype={"protein": str})
    T = T.merge(agg, on=["protein", "position"], how="left")
    for c in ("score_lr_full", "score_lr_peptide", "score_lr_geom", "score_gbm_full"):
        T[c] = T[c].fillna(0.0)
    T.to_csv(f"{RESULTS}/cys_scores_trypsin_arm.csv", index=False, float_format="%.10g")

    # how well each cysteine-level score separates positives from the proteome and observed backgrounds (descriptive)
    pos = T.label == 1
    sep = {}
    for c in ("score_lr_full", "score_lr_peptide", "score_lr_geom", "score_gbm_full", "theo_trypsin", "emp_obs"):
        v = T[c].to_numpy(dtype=float)
        sep[c] = {"auc_pos_vs_proteome": tie_auc(v[(T.label == 1) | (T.label == 0)], pos.to_numpy()),
                  "auc_pos_vs_observed": tie_auc(v[pos | ((T.label == 0) & (T.detected == 1))],
                                                  pos[pos | ((T.label == 0) & (T.detected == 1))].to_numpy()),
                  "mean_pos": float(v[pos].mean()), "mean_proteome_bg": float(v[T.label == 0].mean()),
                  "mean_observed_bg": float(v[(T.label == 0) & (T.detected == 1)].mean())}
    out["cysteine_level_separation"] = sep
    # the detection target restricted to the analysis proteins (descriptive)
    ana = set(T.protein)
    sub = S[S.protein.isin(ana)]
    out["analysis_protein_candidates"] = int(len(sub))
    out["analysis_protein_candidates_detected"] = int(sub.empirical_detected.sum())
    out["lr_full_oof_auc_analysis_proteins"] = tie_auc(sub.lr_full.to_numpy(), sub.empirical_detected.to_numpy())
    out["gbm_full_oof_auc_analysis_proteins"] = tie_auc(sub.gbm_full.to_numpy(), sub.empirical_detected.to_numpy())
    dump_json(out, f"{RESULTS}/empirical_model_summary.json")
    print(json.dumps({k: v for k, v in out.items() if k != "gbm_params"}, indent=1, default=float))


if __name__ == "__main__":
    main()
