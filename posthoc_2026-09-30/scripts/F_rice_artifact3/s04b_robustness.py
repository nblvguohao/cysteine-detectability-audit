"""Step 4b (POST HOC): robustness of the opportunity adjustment.

(a) Detectability window: opportunity recomputed under four digestion windows (the manuscript's
    7-30 residues / <=2 missed cleavages, and three wider ones), because 5% of identified sites lie
    on cysteines the primary window calls undetectable, and an imperfect opportunity measure would
    leave a residual cysteine-count effect behind.
(b) One-to-one matching on abundance AND opportunity: each positive group is paired, without
    replacement, with the unmodified group of the SAME detectable-cysteine count (capped at 20)
    nearest in log10 iBAQ -- the direct analogue of the manuscript's one-to-one abundance matching.
(c) Stratification on exact opportunity without the cap at 12 used in s04.
Labels A (stored) and S (>=1 identified Cys-containing peptide) as in s04.
Outputs: s04b_windows.csv, s04b_matched_abundance_opportunity.csv, s04b_summary.json
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
sys.path.insert(0, str(CYS_AUDIT_SRC))
from cys_audit.proteases import digest

F = pd.read_csv(RES / "s04_group_features.csv")
seqs, _ = load_all_sequences()
CY = pd.read_csv(RES / "s04_cysteines.csv", usecols=["gidx", "position", "site"])
sites_by_g = CY[CY.site].groupby("gidx")["position"].apply(set).to_dict()


def rank_biserial(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    rk = stats.rankdata(np.concatenate([x, y]))
    U = rk[:len(x)].sum() - len(x) * (len(x) + 1) / 2.0
    return float(2.0 * U / (len(x) * len(y)) - 1.0)


def strat_r(df, lab, key, cols, min_n=10):
    rs = []
    for _, sub in df.groupby(cols):
        a = sub.loc[sub[lab] == 1, key].values; b = sub.loc[sub[lab] == 0, key].values
        if len(a) >= min_n and len(b) >= min_n:
            rs.append(rank_biserial(a, b))
    return (float(np.mean(rs)) if rs else float("nan")), len(rs)


e = np.percentile(F["log10_ibaq"], np.linspace(0, 100, 11)); e[0] -= 1e-9; e[-1] += 1e-9
F["dec"] = np.clip(np.searchsorted(e, F["log10_ibaq"].values, side="left") - 1, 0, 9)
F["log2_len"] = np.log2(F["length"]); F["cys"] = F["n_cys"].astype(float)

WINDOWS = {"w7_30_mc2": (7, 30, 2), "w6_35_mc2": (6, 35, 2), "w7_40_mc3": (7, 40, 3), "w5_50_mc3": (5, 50, 3)}
rows = []
for wn, (lo, hi, mc) in WINDOWS.items():
    ndet, on_det, tot_sites = [], 0, 0
    for gi, acc in zip(F.gidx, F.seq_accession):
        s = seqs[acc]
        det = np.zeros(len(s), bool)
        for a0, b0 in digest(s, "trypsin", missed=mc, min_len=lo, max_len=hi):
            det[a0 - 1:b0] = True
        cpos = [i for i, ch in enumerate(s) if ch == "C"]
        ndet.append(int(det[cpos].sum()) if cpos else 0)
        st = sites_by_g.get(gi, set())
        tot_sites += len(st); on_det += sum(det[p - 1] for p in st)
    F[f"ndet_{wn}"] = ndet
    F["log2_det1"] = np.log2(1 + F[f"ndet_{wn}"])
    F["opp_exact12"] = np.minimum(F[f"ndet_{wn}"], 12)
    F["opp_exact_nocap"] = F[f"ndet_{wn}"]
    for lab in ("A", "S"):
        m2 = smf.logit(f"{lab} ~ log10_ibaq + log2_len + cys", data=F).fit(disp=0)
        m3 = smf.logit(f"{lab} ~ log10_ibaq + log2_len + cys + log2_det1", data=F).fit(disp=0)
        c3 = m3.conf_int()
        rec = {"window": wn, "label": lab, "fraction_mapped_sites_on_detectable_cys": on_det / tot_sites,
               "fraction_cys_detectable": float(F[f"ndet_{wn}"].sum() / F["n_cys"].sum()),
               "M2_OR_len": float(np.exp(m2.params["log2_len"])), "M2_OR_cys": float(np.exp(m2.params["cys"])),
               "M3_OR_len": float(np.exp(m3.params["log2_len"])),
               "M3_OR_len_ci": [float(np.exp(c3.loc["log2_len", 0])), float(np.exp(c3.loc["log2_len", 1]))],
               "M3_OR_cys": float(np.exp(m3.params["cys"])),
               "M3_OR_cys_ci": [float(np.exp(c3.loc["cys", 0])), float(np.exp(c3.loc["cys", 1]))],
               "M3_p_cys": float(m3.pvalues["cys"]),
               "M3_OR_opp": float(np.exp(m3.params["log2_det1"]))}
        for key, nm in (("length", "len"), ("n_cys", "cys")):
            rec[f"r_{nm}_dec"] = strat_r(F, lab, key, ["dec"])[0]
            rec[f"r_{nm}_dec_x_exact12"], rec["cells_exact12"] = strat_r(F, lab, key, ["dec", "opp_exact12"])
            rec[f"r_{nm}_dec_x_exact_nocap"], rec["cells_exact_nocap"] = strat_r(F, lab, key, ["dec", "opp_exact_nocap"])
        rows.append(rec)
W = pd.DataFrame(rows)
W.to_csv(RES / "s04b_windows.csv", index=False)

# (b) one-to-one matching on abundance AND exact opportunity (primary window)
F["opp20"] = np.minimum(F["ndet_w7_30_mc2"], 20)


def match_ab_opp(df, lab, strata="opp20"):
    """greedy nearest-log10-iBAQ matching without replacement inside exact-opportunity strata"""
    L = [[], [], [], [], [], []]; dropped = 0; dd = []
    for k, sub in df.groupby(strata):
        pos = sub[sub[lab] == 1].sort_values("log10_ibaq", kind="mergesort")
        neg = sub[sub[lab] == 0].sort_values("log10_ibaq", kind="mergesort")
        negv = neg["log10_ibaq"].values.copy(); used = np.zeros(len(neg), bool)
        nl, nc, nd = neg["length"].values, neg["n_cys"].values, neg["ndet_w7_30_mc2"].values
        for t, pl, pc, pd_ in zip(pos["log10_ibaq"].values, pos["length"].values, pos["n_cys"].values,
                                  pos["ndet_w7_30_mc2"].values):
            d = np.abs(negv - t); d[used] = np.inf
            if len(d) == 0 or not np.isfinite(d.min()):
                dropped += 1; continue
            j = int(d.argmin()); used[j] = True; dd.append(d[j])
            for lst, v in zip(L, (pl, nl[j], pc, nc[j], pd_, nd[j])):
                lst.append(v)
    return {"pairs": len(L[0]), "positives_unmatched": dropped,
            "r_length": rank_biserial(L[0], L[1]), "r_cys": rank_biserial(L[2], L[3]),
            "r_n_det_cys_check": rank_biserial(L[4], L[5]),
            "median_abs_dlog10_ibaq": float(np.median(dd)) if dd else None}


mrows = []
for lab in ("A", "S", "S1"):
    rec = {"label": lab, **match_ab_opp(F, lab)}
    rng = np.random.default_rng(SEED + 7)
    bl, bc = [], []
    for _ in range(200):                       # group bootstrap, matching redone in each resample
        b = F.iloc[rng.integers(0, len(F), len(F))]
        r = match_ab_opp(b, lab); bl.append(r["r_length"]); bc.append(r["r_cys"])
    rec["r_length_ci95_boot200"] = [float(np.nanpercentile(bl, 2.5)), float(np.nanpercentile(bl, 97.5))]
    rec["r_cys_ci95_boot200"] = [float(np.nanpercentile(bc, 2.5)), float(np.nanpercentile(bc, 97.5))]
    # abundance-only 1:1 matching on the same data for comparison (strata = one)
    F["_one"] = 0
    rec["abundance_only_1to1"] = match_ab_opp(F, lab, strata="_one")
    mrows.append(rec)
M = pd.DataFrame(mrows)
M.to_csv(RES / "s04b_matched_abundance_opportunity.csv", index=False)

# (d) non-parametric opportunity adjustment in the regression: exact detectable-Cys count as a
# categorical term (capped at 20), abundance as a natural cubic spline
import statsmodels.api as sm
crows = []
# restricted to groups with 1..20 detectable cysteines so that every exact level is its own
# category (no lumped top level in which n_cys would proxy residual opportunity); GLM/IRLS fit
F20 = F[(F.ndet_w7_30_mc2 >= 1) & (F.ndet_w7_30_mc2 <= 20)].copy()
for lab in ("A", "S"):
    for nm, f, dat in (("M3_log2opp_all", f"{lab} ~ log10_ibaq + log2_len + cys + np.log2(1 + ndet_w7_30_mc2)", F),
                       ("M3_log2opp_ndet1to20", f"{lab} ~ log10_ibaq + log2_len + cys + np.log2(1 + ndet_w7_30_mc2)", F20),
                       ("M5_exact_opp_categorical_ndet1to20", f"{lab} ~ cr(log10_ibaq, df=5) + C(ndet_w7_30_mc2) + log2_len + cys", F20),
                       ("M5b_exact_opp_categorical_log2cys_ndet1to20", f"{lab} ~ cr(log10_ibaq, df=5) + C(ndet_w7_30_mc2) + log2_len + np.log2(1 + n_cys)", F20)):
        m = smf.glm(f, data=dat, family=sm.families.Binomial()).fit()
        ci = m.conf_int()
        for term in m.params.index:
            if term in ("log2_len", "cys", "np.log2(1 + n_cys)"):
                crows.append({"label": lab, "model": nm, "term": term, "OR": float(np.exp(m.params[term])),
                              "ci_low": float(np.exp(ci.loc[term, 0])), "ci_high": float(np.exp(ci.loc[term, 1])),
                              "p": float(m.pvalues[term]), "aic": float(m.aic), "n": int(m.nobs),
                              "n_pos": int(dat[lab].sum())})
# same exact-opportunity model under each detectability window (is the residual Cys-count term an
# artefact of the primary window's narrow definition of opportunity?)
for wn in WINDOWS:
    Fw = F[(F[f"ndet_{wn}"] >= 1) & (F[f"ndet_{wn}"] <= 20)].copy()
    Fw["ndet_w"] = Fw[f"ndet_{wn}"]
    Fw["n_undet"] = Fw["n_cys"] - Fw["ndet_w"]
    for lab in ("A", "S"):
        for nm, f in ((f"M5_exact_opp_{wn}", f"{lab} ~ cr(log10_ibaq, df=5) + C(ndet_w) + log2_len + cys"),
                      (f"M6_exact_opp_undetectable_{wn}", f"{lab} ~ cr(log10_ibaq, df=5) + C(ndet_w) + log2_len + n_undet")):
            m = smf.glm(f, data=Fw, family=sm.families.Binomial()).fit()
            ci = m.conf_int()
            for term in ("log2_len", "cys", "n_undet"):
                if term in m.params.index:
                    crows.append({"label": lab, "model": nm, "term": term, "OR": float(np.exp(m.params[term])),
                                  "ci_low": float(np.exp(ci.loc[term, 0])), "ci_high": float(np.exp(ci.loc[term, 1])),
                                  "p": float(m.pvalues[term]), "aic": float(m.aic), "n": int(m.nobs),
                                  "n_pos": int(Fw[lab].sum())})
CR = pd.DataFrame(crows)
CR.to_csv(RES / "s04b_categorical_opportunity.csv", index=False)
print(CR.round(4).to_string())

write_json(RES / "s04b_summary.json", {
    "analysis_label": "POST HOC revision analysis 2026-09-30 (not registered)",
    "windows": W.to_dict("records"), "matched_on_abundance_and_opportunity": M.to_dict("records"),
    "categorical_opportunity_regressions": CR.to_dict("records")})
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 40)
print(W.round(4).to_string())
print(M.round(4).to_string())
