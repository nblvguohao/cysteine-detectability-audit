# -*- coding: utf-8 -*-
"""B_overadjust_sim step 1g (POST HOC, revision round 3, 2026-09-30).

Round-3 verification proposed a SEQUENCE-ONLY discriminator between the two readings of the SFE-006 attenuation:
in one site table, the association of the label with the local K/R count (x1 = K/R within +/-20 above the universe
median) relative to its association with the claimed SFE-006 offsets (a1). Chemistry that acts through a count
should leave THAT count's association high relative to the offsets; detectability alone should not.

The round-2 sensitivity chemistries followed two other counts (K/R within +/-5; A/K/R/V within +/-10). To test
whether each simulated count chemistry leaves its mark on its own count, two further HYPOTHETICAL count-defined
attributes are defined here in exactly the way x1 was (count above its universe median):

  x3_AKRV10_hi   akrv_count10 (A/K/R/V at the twenty positions within +/-10, centre excluded) > universe median
  x4_KR5_hi      kr_count5 (K/R at the ten positions within +/-5, centre excluded) > universe median

Both counts are read from results/B_overadjust_sim/universe_krcount.csv.gz (01f_krcount_universe.py), which computed
them with the window rule of the repository's phase2b_claim_cohorts.flank_count_flag. Nothing here is registered.

Descriptive diagnostics as in 01f: prevalence, association with detection, how well VIS10 predicts each column
(5-fold protein-grouped CV on the folds of 01_prep_universe.py), and correlation with a1 and with the +/-20 count.
Writes results/B_overadjust_sim/universe_discrim.csv.gz (same row order as universe_cysteines.csv.gz) and
results/B_overadjust_sim/discrim_universe_summary.json.
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
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

KC = os.path.join(cb.RESULTS, "universe_krcount.csv.gz")
OUT = os.path.join(cb.RESULTS, "universe_discrim.csv.gz")
OUT_JSON = os.path.join(cb.RESULTS, "discrim_universe_summary.json")
VIS10_NAMES = ["pep_len", "pep_mass", "pep_gravy", "pep_detectable_length", "pep_detectable_mass",
               "pep_detectable_both", "pep_log_len", "pep_mc1_detectable", "pep_mc2_detectable",
               "pep_detectable_any_missed_cleavage"]


def log2_or(y, a):
    n11 = float(((a == 1) & (y == 1)).sum()); n10 = float(((a == 1) & (y == 0)).sum())
    n01 = float(((a == 0) & (y == 1)).sum()); n00 = float(((a == 0) & (y == 0)).sum())
    return float(np.log2(((n11 + .5) * (n00 + .5)) / ((n10 + .5) * (n01 + .5))))


def main():
    u = pd.read_csv(cb.UNIVERSE, dtype={"protein": str})
    kc = pd.read_csv(KC, dtype={"protein": str})
    assert (kc["protein"].to_numpy() == u["protein"].to_numpy()).all()
    assert (kc["position"].to_numpy() == u["position"].to_numpy()).all()
    akrv10 = kc["akrv_count10"].to_numpy(dtype=float)
    kr5 = kc["kr_count5"].to_numpy(dtype=float)
    kr20 = kc["kr_count20"].to_numpy(dtype=float)
    med10, med5, med20 = float(np.median(akrv10)), float(np.median(kr5)), float(np.median(kr20))
    x3 = (akrv10 > med10).astype(int)
    x4 = (kr5 > med5).astype(int)
    # x1 was defined the same way in 01f; re-derive it as a consistency check
    x1_check = int(((kr20 > med20).astype(int) != kc["x1_KRcount20_hi"].to_numpy()).sum())
    out = pd.DataFrame({"protein": u["protein"], "position": u["position"],
                        "x3_AKRV10_hi": x3, "x4_KR5_hi": x4})
    out.to_csv(OUT, index=False, compression="gzip")

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

    a1 = u["a1_SFE006_KR"].to_numpy(dtype=float)
    x1 = kc["x1_KRcount20_hi"].to_numpy(dtype=float)
    rows = []
    for name, v, rule in (("x1_KRcount20_hi", x1.astype(int), "kr_count20 > %g" % med20),
                          ("x3_AKRV10_hi", x3, "akrv_count10 > %g" % med10),
                          ("x4_KR5_hi", x4, "kr_count5 > %g" % med5)):
        rows.append({"column": name, "rule": rule, "prevalence_all": float(v.mean()),
                     "prevalence_detected": float(v[y_det == 1].mean()),
                     "prevalence_undetected": float(v[y_det == 0].mean()),
                     "log2_or_vs_detection": log2_or(y_det, v), "vis10_predicts_auc_cv5": cv_auc(v),
                     "corr_with_a1_SFE006_KR": float(np.corrcoef(v, a1)[0, 1]),
                     "corr_with_x1_KRcount20_hi": float(np.corrcoef(v, x1)[0, 1]),
                     "corr_with_kr_count20": float(np.corrcoef(v, kr20)[0, 1])})
    summ = {"item": cb.ITEM, "label": "POST HOC (revision round 3, 2026-09-30); not registered",
            "definitions": {"x3_AKRV10_hi": "akrv_count10 > universe median (%g)" % med10,
                            "x4_KR5_hi": "kr_count5 > universe median (%g)" % med5,
                            "x1_KRcount20_hi (01f, for reference)": "kr_count20 > universe median (%g)" % med20},
            "check_x1_rederived_mismatches": x1_check, "columns": rows, "n_rows": int(len(out))}
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
