# -*- coding: utf-8 -*-
"""B_overadjust_sim step 1f (POST HOC, revision round 2, 2026-09-30).

Per-cysteine K/R-count columns for the universe of step 1, written to a SEPARATE file so that the stored
universe and every stored replicate stay unchanged. They implement the verifier's round-2 alternative
chemistry model: genuine chemistry that acts through the LOCAL K/R COUNT, the feature VIS10 deliberately omits
because "a local basic-residue count can also be a genuine chemical feature (an acid-base motif)"
(Supplemental Note 5).

  kr_count20        number of K and R within +/-20 residues of the cysteine: the repository's own DIG25 column
                    `pep_cleavage_sites_within_20`, computed with phase2_claim_cohorts.site_feature_matrix
                    (subset "DIG25", tryptic rule, imported verbatim)
  F_KRcount20_z     kr_count20 standardised over the universe (mean 0, SD 1; population SD). Chemistry variable
                    of the proxy-chemistry model: P(positive | detected) = expit(a + b ln2 F); b is a log2 odds
                    ratio per SD of the count
  x1_KRcount20_hi   kr_count20 above its universe median: a HYPOTHETICAL claim defined by the local K/R count
                    itself (the verifier's F_count20_hi); binary, analysed like the published attributes
  x2_KR5_ge3        >= 3 of the ten flanking positions within +/-5 are K or R: repository function
                    phase2b_claim_cohorts.flank_count_flag with residues {K, R} (the SNO-016 rule without H);
                    a second hypothetical count-defined claim (the verifier's F_kr5_hi)
  kr_count5, F_KR5_z        number of K/R at the ten positions within +/-5 (centre excluded; same window logic as
                            flank_count_flag), raw and standardised: sensitivity chemistry variable
  akrv_count10, F_AKRV10_z  number of A/K/R/V at the twenty positions within +/-10 (centre excluded), raw and
                            standardised: sensitivity chemistry variable for the A/K/R/V composition claim

Descriptive diagnostics (results/B_overadjust_sim/krcount_universe_summary.json): prevalence, association with
detection, how well VIS10 predicts each column (5-fold protein-grouped CV, the folds of 01_prep_universe.py:
AUC for binary columns, R^2 of a linear fit for counts), correlation with the published attributes, and a
cross-check of F_KRcount20_z against the verifier's independently computed F_count20
(verify_r2/universe_altchem.csv.gz, read only).
Writes results/B_overadjust_sim/universe_krcount.csv.gz (same row order as universe_cysteines.csv.gz).
"""
from __future__ import annotations

import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.linear_model import LinearRegression, LogisticRegression  # noqa: E402
from sklearn.metrics import r2_score, roc_auc_score  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

import phase2_claim_cohorts as p2c  # noqa: E402  (repo, read-only)
import phase2b_claim_cohorts as p2b  # noqa: E402
from run_cross_protease_detectability_probe import FEATURE_NAMES  # noqa: E402

OUT = os.path.join(cb.RESULTS, "universe_krcount.csv.gz")
OUT_JSON = os.path.join(cb.RESULTS, "krcount_universe_summary.json")
VERIFIER_ALT = os.path.join(cb.RESULTS, "verify_r2", "universe_altchem.csv.gz")
VIS10_NAMES = ["pep_len", "pep_mass", "pep_gravy", "pep_detectable_length", "pep_detectable_mass",
               "pep_detectable_both", "pep_log_len", "pep_mc1_detectable", "pep_mc2_detectable",
               "pep_detectable_any_missed_cleavage"]
PUBLISHED = ["a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV", "a3_SNO006_DE3", "a3b_SFE002_E"]


def flank_count(sequence_of, keys, residues, window):
    """Number of flanking positions within +/-window (centre excluded) whose residue is in `residues`; the window
    logic of phase2b_claim_cohorts.flank_count_flag (out-of-sequence positions are not counted)."""
    out = []
    for accession, position in keys:
        s = sequence_of[accession]
        c = position - 1
        n = 0
        for off in range(-window, window + 1):
            if off == 0:
                continue
            i = c + off
            if 0 <= i < len(s) and s[i] in residues:
                n += 1
        out.append(n)
    return np.asarray(out, dtype=float)


def zscore(v):
    v = np.asarray(v, dtype=float)
    return (v - v.mean()) / v.std()


def log2_or(y, a):
    n11 = float(((a == 1) & (y == 1)).sum()); n10 = float(((a == 1) & (y == 0)).sum())
    n01 = float(((a == 0) & (y == 1)).sum()); n00 = float(((a == 0) & (y == 0)).sum())
    return float(np.log2(((n11 + .5) * (n00 + .5)) / ((n10 + .5) * (n01 + .5))))


def main():
    u = pd.read_csv(cb.UNIVERSE, dtype={"protein": str})
    keys = list(zip(u["protein"], u["position"].astype(int)))
    fasta = p2b.read_fasta_gz(cb.FASTA)
    seq_of = {a: fasta[a] for a in u["protein"].unique()}
    assert all(seq_of[a][p - 1] == "C" for a, p in keys)

    # repository DIG25 matrix; the K/R count is the column VIS10 omits
    dig = p2c.site_feature_matrix(seq_of, keys, "DIG25")
    j = FEATURE_NAMES.index("pep_cleavage_sites_within_20")
    kr20 = dig[:, j]
    # the VIS10 columns of the same call must equal the stored VIS10 (same code path)
    vis_idx = [FEATURE_NAMES.index(n) for n in VIS10_NAMES]
    vis_check = float(np.nanmax(np.abs(dig[:, vis_idx] - u[VIS10_NAMES].to_numpy(dtype=float))))

    x1 = (kr20 > np.median(kr20)).astype(int)
    x2 = p2b.flank_count_flag(seq_of, keys, frozenset("KR"))            # >= 3 of +/-5 in {K,R}
    kr5 = flank_count(seq_of, keys, "KR", 5)
    akrv10 = flank_count(seq_of, keys, "AKRV", 10)
    # consistency: x2 is the >= 3 threshold of kr5 (same window rule)
    x2_check = int(((kr5 >= 3).astype(int) != x2).sum())

    out = pd.DataFrame({"protein": u["protein"], "position": u["position"],
                        "kr_count20": kr20, "F_KRcount20_z": zscore(kr20),
                        "x1_KRcount20_hi": x1, "x2_KR5_ge3": x2,
                        "kr_count5": kr5, "F_KR5_z": zscore(kr5),
                        "akrv_count10": akrv10, "F_AKRV10_z": zscore(akrv10)})
    out.to_csv(OUT, index=False, compression="gzip")

    # ---------------- descriptive diagnostics ----------------
    y_det = u["detected"].to_numpy().astype(int)
    X = StandardScaler().fit_transform(u[VIS10_NAMES].to_numpy(dtype=float))
    groups = u["protein"].to_numpy()
    rng = np.random.default_rng(cb.MASTER_SEED)                 # the folds of 01_prep_universe.py
    uprot = np.unique(groups)
    fold_of = dict(zip(uprot, rng.permutation(len(uprot)) % 5))
    fold = np.asarray([fold_of[g] for g in groups])

    def cv_auc(target):
        oof = np.zeros(len(target))
        for k in range(5):
            tr, te = fold != k, fold == k
            m = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs").fit(X[tr], target[tr])
            oof[te] = m.predict_proba(X[te])[:, 1]
        return float(roc_auc_score(target, oof))

    def cv_r2(target):
        oof = np.zeros(len(target))
        for k in range(5):
            tr, te = fold != k, fold == k
            oof[te] = LinearRegression().fit(X[tr], target[tr]).predict(X[te])
        return float(r2_score(target, oof))

    rows = []
    for name, v in (("x1_KRcount20_hi", x1), ("x2_KR5_ge3", x2)):
        rows.append({"column": name, "kind": "binary (hypothetical count-defined claim)",
                     "prevalence_all": float(v.mean()), "prevalence_detected": float(v[y_det == 1].mean()),
                     "prevalence_undetected": float(v[y_det == 0].mean()),
                     "log2_or_vs_detection": log2_or(y_det, v), "vis10_predicts_auc_cv5": cv_auc(v)})
    for name, v in (("kr_count20", kr20), ("kr_count5", kr5), ("akrv_count10", akrv10)):
        rows.append({"column": name, "kind": "count (chemistry variable of the proxy model, standardised)",
                     "mean": float(v.mean()), "sd_population": float(v.std()), "median": float(np.median(v)),
                     "p5": float(np.percentile(v, 5)), "p95": float(np.percentile(v, 95)),
                     "mean_detected": float(v[y_det == 1].mean()), "mean_undetected": float(v[y_det == 0].mean()),
                     "vis10_predicts_r2_cv5": cv_r2(v)})
    corr = {}
    for name, v in (("kr_count20", kr20), ("kr_count5", kr5), ("akrv_count10", akrv10)):
        corr[name] = {a: float(np.corrcoef(v, u[a].to_numpy(dtype=float))[0, 1]) for a in PUBLISHED}
    # cross-check against the verifier's independent computation (read only)
    vcheck = None
    if os.path.exists(VERIFIER_ALT):
        va = pd.read_csv(VERIFIER_ALT, dtype={"protein": str})
        m = out.merge(va[["protein", "position", "F_count20"]], on=["protein", "position"], how="inner")
        vcheck = {"rows_matched": int(len(m)), "rows_universe": int(len(out)),
                  "max_abs_diff_F_KRcount20_z_vs_verifier_F_count20":
                      float(np.max(np.abs(m["F_KRcount20_z"] - m["F_count20"])))}
    summ = {"item": cb.ITEM, "label": "POST HOC (revision round 2, 2026-09-30); not registered",
            "definitions": {
                "kr_count20": "repository DIG25 column pep_cleavage_sites_within_20: K/R in sequence[c-20 : c+21] "
                              "(phase2_claim_cohorts.site_feature_matrix, subset DIG25, tryptic rule)",
                "x1_KRcount20_hi": "kr_count20 > universe median (%g)" % float(np.median(kr20)),
                "x2_KR5_ge3": "phase2b_claim_cohorts.flank_count_flag(residues={K,R}): >= 3 of the ten positions "
                              "within +/-5",
                "kr_count5": "K/R at the ten positions within +/-5 (centre excluded)",
                "akrv_count10": "A/K/R/V at the twenty positions within +/-10 (centre excluded)"},
            "columns": rows, "correlation_with_published_attributes": corr,
            "check_vis10_columns_of_DIG25_call_max_abs_diff_vs_stored": vis_check,
            "check_x2_equals_kr_count5_ge_3_mismatches": x2_check,
            "check_against_verifier": vcheck, "n_rows": int(len(out))}
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
