"""Step 8c (POST HOC; revision after adversarial verification, 2026-09-30): the undetectable-cysteine
contrast in the designs of s04 (unmatched, within iBAQ deciles, 1:1 on iBAQ only, within iBAQ decile x
opportunity tertile), so that the Supplemental Note table can show where the residual sits in each design.
Same rules as s04 (deciles of all 7,692 groups, tertile cut points 3 and 6, strata with <10 per arm dropped,
1:1 nearest log10 iBAQ without replacement). No resampling. Output: s08c_undet_designs.csv
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
from s08_lib import (Design, match_within, rank_biserial, sorted_strata_index, strat_r)

F = pd.read_csv(RES / "s08_group_features_windows.csv")
D = Design(F[["log10_ibaq", "length", "n_cys", "n_det"]])
one = sorted_strata_index(np.zeros(len(F), int), D.x)
rows = []
for lab in ("A", "S", "S1", "S_map"):
    y = F[lab].to_numpy(int)
    P, N = match_within(y, D.x, one)
    for key, v in (("length", D.length), ("n_cys", D.ncys), ("n_undet", D.nundet), ("n_det", D.ndet)):
        rows.append({"label": lab, "feature": key, "r_before": rank_biserial(v[y == 1], v[y == 0]),
                     "r_within_iBAQ_deciles": strat_r(y, v, D.s_dec)[0],
                     "r_1to1_iBAQ_only": rank_biserial(v[P], v[N]),
                     "r_within_decile_x_opp_tertile": strat_r(y, v, D.s_dt)[0]})
R = pd.DataFrame(rows)
# consistency with s04 (length and n_cys must reproduce s04_rank_biserial_strata.csv)
s4 = pd.read_csv(RES / "s04_rank_biserial_strata.csv").set_index(["label", "feature"])
for _, r in R[R.feature.isin(["length", "n_cys"])].iterrows():
    ref = s4.loc[(r["label"], r["feature"])]
    for a, b in (("r_before", "r_before"), ("r_within_iBAQ_deciles", "r_within_iBAQ_deciles"),
                 ("r_1to1_iBAQ_only", "r_1to1_matched"), ("r_within_decile_x_opp_tertile", "r_within_decile_x_opp_tertile")):
        assert abs(r[a] - ref[b]) < 1e-9, (r["label"], r["feature"], a, r[a], ref[b])
R.to_csv(RES / "s08c_undet_designs.csv", index=False)
pd.set_option("display.width", 200)
print(R.round(4).to_string())
