"""Step 5 (POST HOC revision analysis, item D_empirical_background): diagnostics that explain where the residual
under empirical matching comes from. Descriptive; no interval is used for a decision.

  (a) candidate level (global arm): overall and within-protein AUC of each detectability score for identification
      in the global arm. Within-protein AUC = pooled over all (identified, unidentified) candidate pairs of the same
      protein (per-protein Mann-Whitney AUC weighted by n1*n0).
  (b) cysteine level (trypsin ABE arm proteins): overall and within-protein AUC of the cysteine-level scores for
      (i) observation in the global arm, (ii) identification in the ABE arm, (iii) the site label.
  (c) K/R and D/E band shares by empirical-score bin, positives versus the proteome background.
  (d) K/R distal share of the global-arm-observed background by protein-abundance tertile and by whether the
      protein carries a site.

Output: diagnostics.json, score_bin_table.csv
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

from common import IN_CAND, RESULTS, dump_json
from cys_audit.stats import auc as tie_auc


def within_auc(df, score, target, group="protein"):
    d = df[[group, score, target]].copy()
    d["r"] = d.groupby(group)[score].rank(method="average")
    g = d.groupby(group).agg(n1=(target, "sum"), n=(target, "size"))
    g["n0"] = g.n - g.n1
    s1 = d[d[target] == 1].groupby(group).r.sum()
    g["s1"] = s1.reindex(g.index).fillna(0)
    g = g[(g.n1 > 0) & (g.n0 > 0)]
    u = g.s1 - g.n1 * (g.n1 + 1) / 2
    return {"within_protein_auc": float(u.sum() / (g.n1 * g.n0).sum()), "n_proteins_informative": int(len(g)),
            "n_pairs": float((g.n1 * g.n0).sum())}


def main():
    out = {}
    S = pd.read_csv(f"{RESULTS}/candidate_oof_scores.csv.gz", dtype={"protein": str})
    C = pd.read_csv(IN_CAND, dtype={"protein": str}, usecols=["protein", "theoretical_detectable", "log_abundance",
                                                             "abundance_missing"])
    S["theoretical_detectable"] = C.theoretical_detectable.to_numpy()
    cand = {}
    for sc in ("lr_full", "lr_peptide", "lr_geom", "gbm_full", "theoretical_detectable"):
        cand[sc] = {"overall_auc": tie_auc(S[sc].to_numpy(dtype=float), S.empirical_detected.to_numpy()),
                    **within_auc(S, sc, "empirical_detected")}
    out["candidate_level_global_arm_detection"] = cand

    T = pd.read_csv(f"{RESULTS}/cys_scores_trypsin_arm.csv", dtype={"protein": str})
    ab = C.groupby("protein").agg(log_abundance=("log_abundance", "first"), abundance_missing=("abundance_missing", "first"))
    T = T.merge(ab, left_on="protein", right_index=True, how="left")
    cys = {}
    for target in ("emp_obs", "detected", "label"):
        cys[target] = {}
        for sc in ("score_lr_full", "score_lr_peptide", "score_gbm_full", "theo_trypsin"):
            cys[target][sc] = {"overall_auc": tie_auc(T[sc].to_numpy(dtype=float), T[target].to_numpy()),
                               **within_auc(T, sc, target)}
    out["cysteine_level"] = cys

    # (c) score-bin table
    bins = [-np.inf, 1e-9, 0.01, 0.03, 0.1, 0.2, 0.4, 0.6, np.inf]
    T["score_bin"] = pd.cut(T.score_lr_full, bins, labels=["0", "(0,0.01]", "(0.01,0.03]", "(0.03,0.1]", "(0.1,0.2]",
                                                          "(0.2,0.4]", "(0.4,0.6]", ">0.6"])
    rows = []
    for b, g in T.groupby("score_bin", observed=True):
        for labv, nm in ((1, "positive"), (0, "proteome_background")):
            h = g[g.label == labv]
            rows.append({"score_bin": b, "set": nm, "n": int(len(h)), "kr_dist_share": h.kr_dist.mean() if len(h) else np.nan,
                         "kr_prox_share": h.kr_prox.mean() if len(h) else np.nan,
                         "de_dist_share": h.de_dist.mean() if len(h) else np.nan,
                         "share_identified_in_abe_arm": h.detected.mean() if len(h) else np.nan,
                         "share_emp_obs": h.emp_obs.mean() if len(h) else np.nan})
    B = pd.DataFrame(rows)
    B.to_csv(f"{RESULTS}/score_bin_table.csv", index=False)

    # (d) K/R distal share of backgrounds by abundance tertile and by site-carrying protein
    same = T.protein.isin(set(T.protein[T.label == 1]))
    T["same_protein"] = same.astype(int)
    has_ab = T.abundance_missing == 0
    q = np.quantile(T.loc[has_ab, "log_abundance"].drop_duplicates(), [1 / 3, 2 / 3]) if has_ab.any() else [0, 0]
    T["abundance_tertile"] = np.where(~has_ab, "missing", np.where(T.log_abundance <= q[0], "low",
                                                                    np.where(T.log_abundance <= q[1], "mid", "high")))
    tab = []
    for (tert, sp), g in T.groupby(["abundance_tertile", "same_protein"]):
        for nm, m in (("positive", g.label == 1), ("observed_bg", (g.label == 0) & (g.detected == 1)),
                      ("emp_obs_bg", (g.label == 0) & (g.emp_obs == 1)), ("proteome_bg", g.label == 0)):
            h = g[m]
            tab.append({"abundance_tertile": tert, "same_protein": int(sp), "set": nm, "n": int(len(h)),
                        "kr_dist_share": float(h.kr_dist.mean()) if len(h) else None})
    out["kr_dist_share_by_abundance_and_protein"] = tab
    out["abundance_tertile_cuts_log10"] = [float(x) for x in q]
    out["positives_by_abundance_tertile"] = T[T.label == 1].abundance_tertile.value_counts().to_dict()
    out["positives_score_quantiles_lr_full"] = dict(zip(["min", "p10", "p25", "p50", "p75", "p90", "max"],
                                                        [float(x) for x in np.quantile(T.score_lr_full[T.label == 1],
                                                                                        [0, .1, .25, .5, .75, .9, 1])]))
    out["positives_with_score_le_0.01"] = int((T.score_lr_full[T.label == 1] <= 0.01).sum())
    out["positives_protein_abundance_missing"] = int((T.abundance_missing[T.label == 1] == 1).sum())
    dump_json(out, f"{RESULTS}/diagnostics.json")
    import json
    print(json.dumps({k: v for k, v in out.items() if k != "kr_dist_share_by_abundance_and_protein"}, indent=1, default=float))
    print(pd.DataFrame(tab).to_string())
    print(B.round(3).to_string())


if __name__ == "__main__":
    main()
