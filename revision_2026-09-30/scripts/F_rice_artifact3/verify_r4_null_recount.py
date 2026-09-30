"""VERIFIER round 4 (adversarial; POST HOC), item F_rice_artifact3, part (a): recount the pooled null results
from the item's raw draws (s10b_null_draws.csv, first 1,000 sets; s12_null_ext_draws.csv.gz, 4,000 sets),
with the verifier's own summary code, and check every null number quoted in the proposed texts.
Output: verify_r4/null_recount_r4.json (+ stdout)."""
import json
import pathlib

import numpy as np
import pandas as pd

W = pathlib.Path("C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30")
RES = W / "results" / "F_rice_artifact3"
OUT = RES / "verify_r4"
STATS14 = ["r_dec_len", "r_dec_cys", "r_dt_len", "r_dt_cys", "r_m_len", "r_m_cys", "r_m_undet", "b_undet_M3u",
           "b_len_M3u", "r_mc_det", "r_mc_len", "dAIC_D_minus_T", "b_det_B", "d_det_minus_undet_L"]
CYS_STATS = ["r_dec_cys", "r_dt_cys", "r_m_cys", "r_m_undet", "b_undet_M3u", "r_mc_det", "dAIC_D_minus_T", "b_det_B",
             "d_det_minus_undet_L"]
LEN_STATS = ["r_dec_len", "r_dt_len", "r_m_len", "b_len_M3u", "r_mc_len"]

old = pd.read_csv(RES / "s10b_null_draws.csv")
ext = pd.read_csv(RES / "s12_null_ext_draws.csv.gz")
obs = json.loads((RES / "s12_null_pooled_summary.json").read_text(encoding="utf-8"))["observed"]
# observed values recomputed independently are checked in verify_r4_null_sim.py; here the item's are used
out = {}
for (lab, null), sub_o in old.groupby(["label", "null"]):
    sub_e = ext[(ext.label == lab) & (ext.null == null)]
    assert len(sub_o) == 1000 and len(sub_e) == 4000
    rec = {}
    for s in STATS14:
        v = np.concatenate([sub_o[s].to_numpy(float), sub_e[s].to_numpy(float)])
        v = v[np.isfinite(v)]
        o = float(obs[lab][s])
        n = len(v)
        pge = (1 + np.sum(v >= o)) / (1 + n); ple = (1 + np.sum(v <= o)) / (1 + n)
        p2 = min(1.0, 2 * min(pge, ple))
        lo, hi = np.percentile(v, [2.5, 97.5])
        # first 1,000 only
        v1 = sub_o[s].to_numpy(float); v1 = v1[np.isfinite(v1)]
        lo1, hi1 = np.percentile(v1, [2.5, 97.5])
        p21 = min(1.0, 2 * min((1 + np.sum(v1 >= o)) / (1 + len(v1)), (1 + np.sum(v1 <= o)) / (1 + len(v1))))
        rec[s] = {"obs": o, "lo": float(lo), "hi": float(hi), "p2": float(p2), "inside": bool(lo <= o <= hi),
                  "n": int(n), "p2_first1000": float(p21), "inside_first1000": bool(lo1 <= o <= hi1)}
    out[f"{lab}|{null[:2]}"] = rec

summary = {}
for key, rec in out.items():
    outside = [s for s in STATS14 if not rec[s]["inside"]]
    edge = [(s, round(rec[s]["p2"], 4)) for s in STATS14 if 0.04 <= rec[s]["p2"] <= 0.06]
    inside_p = [rec[s]["p2"] for s in STATS14 if rec[s]["inside"]]
    summary[key] = {"n_outside": len(outside), "outside": outside, "edge": edge,
                    "min_p_inside": min(inside_p) if inside_p else None,
                    "n_outside_first1000": sum(not rec[s]["inside_first1000"] for s in STATS14),
                    "max_p_outside": max([rec[s]["p2"] for s in outside]) if outside else None,
                    "cys_stats_min_p": min(rec[s]["p2"] for s in CYS_STATS),
                    "len_stats": {s: (round(rec[s]["obs"], 4), round(rec[s]["p2"], 4), rec[s]["inside"]) for s in LEN_STATS}}
    # any statistic outside the range but with P >= 0.05, or inside with P < 0.05 (range/P disagreement)
    summary[key]["range_p_disagreements"] = [(s, round(rec[s]["p2"], 4), rec[s]["inside"]) for s in STATS14
                                            if (rec[s]["inside"] and rec[s]["p2"] < 0.05) or
                                            ((not rec[s]["inside"]) and rec[s]["p2"] >= 0.05)]
(OUT / "null_recount_r4.json").write_text(json.dumps({"per_stat": out, "summary": summary}, indent=1), encoding="utf-8")
for k, v in summary.items():
    print(k, json.dumps(v))
print()
for key in ("A|N4", "A|N2", "S|N4"):
    print(key, {s: (round(out[key][s]["obs"], 4), round(out[key][s]["lo"], 4), round(out[key][s]["hi"], 4),
                    round(out[key][s]["p2"], 4)) for s in STATS14})
