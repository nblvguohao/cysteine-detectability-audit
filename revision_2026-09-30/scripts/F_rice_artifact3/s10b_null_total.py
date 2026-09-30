"""Step 10b (POST HOC; revision round 2, 2026-09-30): total-count nulls beside the detectable-count nulls.

Four calibrated nulls, each fitted to the observed labels on the 6,990 groups with >= 1 detectable
cysteine (as in s08b), with no dependence on length and no dependence on anything but abundance and ONE
cysteine count:
  N1  complementary log-log, free exponent on the DETECTABLE count: 1 - exp(-exp(a + b log10 iBAQ) n_det^g)
  N2  logistic, natural spline in log10 iBAQ (df 5) + DETECTABLE count as a categorical term (capped at 20)
  N3  complementary log-log, free exponent on the TOTAL count:  1 - exp(-exp(a + b log10 iBAQ) n_cys^g)
  N4  logistic, natural spline in log10 iBAQ (df 5) + TOTAL count as a categorical term (capped at 20;
      no category is empty for any label: minimum 7 positives)
Labels S, S_map and A, all four nulls each; 1,000 label sets per (label, null).
Seeds: numpy SeedSequence(entropy, spawn_key = (chunk,)), 4 chunks of 250, 4 worker processes. The five
combinations of s08b keep their s08b entropies (SEED + 810 + k, k = 0..4), so their label sets are the
s08b label sets and every s08b statistic is reproduced exactly (asserted below); the seven new
combinations use SEED + 815 .. SEED + 821.
Statistics per label set: s08_lib.Design.stats (unchanged) plus s10_lib.Design2.stats2 (mirror matching
within exact total count, AIC(D) - AIC(T), each count given the other, symmetric linear per-cysteine
coefficients, length at fixed total count). One-sided P: p_ge = (1 + #null >= observed)/(1 + reps),
p_le likewise; p_two = min(1, 2 min(p_ge, p_le)).
Outputs: s10b_null_draws.csv, s10b_null_summary.csv, s10b_null_summary.json, s10b_contrast_table.csv
"""
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
from s10_lib import Design2, simulate_chunk2

warnings.filterwarnings("ignore")
REPS, CHUNKS, WORKERS = 1000, 4, 4

NULLS = {"N1_cloglog_free_det": ("cloglog", "log10_ibaq + log_det"),
         "N2_logit_spline_categorical_det": ("logit", "cr(log10_ibaq, df=5) + C(opp20)"),
         "N3_cloglog_free_total": ("cloglog", "log10_ibaq + log_cys"),
         "N4_logit_spline_categorical_total": ("logit", "cr(log10_ibaq, df=5) + C(cys20)")}
# (label, null) -> entropy; the first five are the s08b combinations and entropies
S08B_NAME = {"N1_cloglog_free_det": "N1_cloglog_free", "N2_logit_spline_categorical_det": "N2_logit_spline_categorical"}
COMBOS = [("S", "N1_cloglog_free_det", SEED + 810), ("S", "N2_logit_spline_categorical_det", SEED + 811),
          ("S_map", "N1_cloglog_free_det", SEED + 812), ("S_map", "N2_logit_spline_categorical_det", SEED + 813),
          ("A", "N2_logit_spline_categorical_det", SEED + 814),
          ("A", "N1_cloglog_free_det", SEED + 815),
          ("S", "N3_cloglog_free_total", SEED + 816), ("S", "N4_logit_spline_categorical_total", SEED + 817),
          ("S_map", "N3_cloglog_free_total", SEED + 818), ("S_map", "N4_logit_spline_categorical_total", SEED + 819),
          ("A", "N3_cloglog_free_total", SEED + 820), ("A", "N4_logit_spline_categorical_total", SEED + 821)]
STATS = ["n_pos", "r_dec_len", "r_dec_cys", "r_dt_len", "r_dt_cys", "r_m_len", "r_m_cys", "r_m_undet", "r_cell_undet",
         "mean_undet_pos", "mean_undet_neg", "b_len_M3u", "b_undet_M3u",
         "r_mc_det", "r_mc_undet", "r_mc_len", "r_cellc_det", "mean_det_pos_c", "mean_det_neg_c",
         "dAIC_D_minus_T", "LR_det_given_T", "LR_cys_given_D", "b_det_B", "b_cys_B", "b_len_T",
         "b_ndet_L", "b_nundet_L", "d_det_minus_undet_L"]
S08B_STATS = STATS[:13]
# the contrasts that the proposed Results text reports, in the order of the text
CONTRAST = [("r_dec_len", "length, within abundance deciles"), ("r_dec_cys", "cysteine count, within abundance deciles"),
            ("r_dt_len", "length, decile x detectable-count tertile"), ("r_dt_cys", "cysteine count, decile x detectable-count tertile"),
            ("r_m_len", "length, 1:1 on abundance and exact detectable count"),
            ("r_m_cys", "cysteine count, 1:1 on abundance and exact detectable count"),
            ("r_m_undet", "undetectable cysteines, 1:1 on abundance and exact detectable count"),
            ("b_undet_M3u", "log OR per undetectable cysteine (M3u)"), ("b_len_M3u", "log OR per doubling of length (M3u)"),
            ("r_mc_det", "detectable cysteines, 1:1 on abundance and exact total count"),
            ("r_mc_len", "length, 1:1 on abundance and exact total count"),
            ("dAIC_D_minus_T", "AIC(detectable-count model) - AIC(total-count model)"),
            ("b_det_B", "log2(1 + detectable count) given log2(1 + total count)"),
            ("d_det_minus_undet_L", "per-cysteine log OR, detectable minus undetectable")]


def main():
    F = pd.read_csv(RES / "s08_group_features_windows.csv")
    H = F[F.n_det >= 1].reset_index(drop=True)
    H["log_det"] = np.log(H.n_det); H["log_cys"] = np.log(H.n_cys)
    H["opp20"] = np.minimum(H.n_det, 20); H["cys20"] = np.minimum(H.n_cys, 20)
    cols = ["log10_ibaq", "length", "n_cys", "n_det"]
    D = Design2(H[cols])
    tasks, meta = [], {}
    for lab, null, ent in COMBOS:
        link, rhs = NULLS[null]
        fam = sm.families.Binomial() if link == "logit" else sm.families.Binomial(link=sm.families.links.CLogLog())
        m = smf.glm(f"{lab} ~ {rhs}", data=H, family=fam).fit()
        p = m.fittedvalues.to_numpy(float)
        if (lab, "obs") not in meta:
            meta[(lab, "obs")] = D.stats_all(H[lab].to_numpy(int))
        meta[(lab, null)] = {"label": lab, "null": null, "entropy": ent, "link": link, "rhs": rhs,
                             "converged": bool(m.converged), "aic": float(m.aic),
                             "expected_pos": float(p.sum()), "observed_pos": int(H[lab].sum()),
                             "null_params": {t: float(v) for t, v in m.params.items() if not t.startswith("C(")}}
        for c in range(CHUNKS):
            tasks.append(((lab, null), (H[cols], p, REPS // CHUNKS, ent, (c,))))
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        results = list(ex.map(simulate_chunk2, [t[1] for t in tasks]))
    draws = []
    for (key, _), res in zip(tasks, results):
        for i, r in enumerate(res):
            draws.append({"label": key[0], "null": key[1], **r})
    DR = pd.DataFrame(draws)
    DR["draw"] = DR.groupby(["label", "null"]).cumcount()
    DR.to_csv(RES / "s10b_null_draws.csv", index=False)

    # the s08b combinations must reproduce s08b exactly
    old = pd.read_csv(RES / "s08b_null_draws.csv")
    old["draw"] = old.groupby(["label", "null"]).cumcount()
    checked = 0
    for lab, null, ent in COMBOS[:5]:
        a = DR[(DR.label == lab) & (DR.null == null)].sort_values("draw")[S08B_STATS].to_numpy(float)
        b = old[(old.label == lab) & (old.null == S08B_NAME[null])].sort_values("draw")[S08B_STATS].to_numpy(float)
        assert a.shape == b.shape and np.allclose(a, b, rtol=0, atol=1e-12, equal_nan=True), (lab, null)
        checked += 1

    rows = []
    for lab, null, ent in COMBOS:
        sub = DR[(DR.label == lab) & (DR.null == null)]
        obs = meta[(lab, "obs")]
        for s in STATS:
            v = sub[s].to_numpy(float); v = v[np.isfinite(v)]
            o = float(obs[s])
            pge = float((1 + (v >= o).sum()) / (1 + len(v))); ple = float((1 + (v <= o).sum()) / (1 + len(v)))
            rows.append({"label": lab, "null": null, "statistic": s, "observed": o, "null_mean": float(v.mean()),
                         "null_2.5": float(np.percentile(v, 2.5)), "null_97.5": float(np.percentile(v, 97.5)),
                         "p_ge_observed": pge, "p_le_observed": ple, "p_two_sided": float(min(1.0, 2 * min(pge, ple))),
                         "inside_95_range": bool(np.percentile(v, 2.5) <= o <= np.percentile(v, 97.5)),
                         "reps": int(len(v)), "expected_pos_null": meta[(lab, null)]["expected_pos"],
                         "observed_pos": meta[(lab, null)]["observed_pos"]})
    SUM = pd.DataFrame(rows)
    SUM.to_csv(RES / "s10b_null_summary.csv", index=False)
    # compact table of the reported contrasts
    ct = []
    for lab in ("A", "S", "S_map"):
        for s, desc in CONTRAST:
            rec = {"label": lab, "statistic": s, "description": desc, "observed": float(meta[(lab, "obs")][s])}
            for null in NULLS:
                q = SUM[(SUM.label == lab) & (SUM.null == null) & (SUM.statistic == s)].iloc[0]
                rec[f"{null[:2]}_range"] = [float(q["null_2.5"]), float(q["null_97.5"])]
                rec[f"{null[:2]}_p_two"] = float(q["p_two_sided"])
                rec[f"{null[:2]}_inside"] = bool(q["inside_95_range"])
            ct.append(rec)
    CT = pd.DataFrame(ct)
    CT.to_csv(RES / "s10b_contrast_table.csv", index=False)
    counts = {}
    for lab in ("A", "S", "S_map"):
        for null in NULLS:
            sub = CT[CT.label == lab]
            counts[f"{lab}|{null}"] = {"reported_contrasts": int(len(sub)),
                                       "outside_95_range": sub.loc[~sub[f"{null[:2]}_inside"], "statistic"].tolist()}
    write_json(RES / "s10b_null_summary.json", {
        "analysis_label": "POST HOC revision analysis 2026-09-30, round 2 (not registered)",
        "groups": int(len(H)), "reps_per_combo": REPS, "workers": WORKERS,
        "s08b_combinations_reproduced_exactly": checked,
        "combos": [{k: v for k, v in meta[(lab, null)].items()} for lab, null, _ in COMBOS],
        "observed": {lab: meta[(lab, "obs")] for lab in ("S", "S_map", "A")},
        "contrasts_outside_null_range": counts, "summary": rows})
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 40)
    print("s08b combinations reproduced exactly:", checked)
    print(json.dumps({f"{m['label']}|{m['null']}": [m['expected_pos'], m['observed_pos'], m['aic'], m['converged']]
                      for m in [meta[(l, n)] for l, n, _ in COMBOS]}, indent=1))
    print(json.dumps(counts, indent=1))
    print(SUM.round(4).to_string())


if __name__ == "__main__":
    main()
