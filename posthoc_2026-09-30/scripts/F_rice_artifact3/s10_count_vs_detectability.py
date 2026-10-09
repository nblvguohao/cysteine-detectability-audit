"""Step 10 (POST HOC; revision round 2, 2026-09-30): detectable versus total cysteine count.

Why this step exists. Verifier round 2 (major 1): the round-1 text credited the abundance-matched contrast
mostly to detection opportunity (the number of theoretically detectable cysteines, n_det) and a smaller
residual to undetectable cysteines. The evidence cannot tell n_det from the total cysteine count n_cys
(Spearman 0.94); the verifier's discriminating tests favoured total count. This step re-does those tests
with this item's own code and features, so that the report and the proposed text rest on them:

(1) protein-level model comparisons, labels A (registered), S, S1, S_map; three samples/families:
      logit, all 7,692 groups; logit, the 6,990 groups with >=1 detectable Cys; complementary log-log
      (free exponents), >=1 detectable Cys; logit with a natural spline in log10 iBAQ (df 5) and either a
      quadratic in log2(1 + count) or the count as a categorical term (capped at 20), >=1 detectable Cys.
    Models: D (iBAQ [+ log2 length] + detectable count), T (same with total count), B (both counts);
    reported: AIC(D) - AIC(T) (equal df), LR of each count given the other, the length term at fixed
    total count, and the fraction-detectable term at fixed total count.
(2) symmetric linear counts: y ~ log10 iBAQ + log2 length + n_det + n_undet; OR per detectable and per
    undetectable cysteine and a Wald test of equality; y ~ log10 iBAQ + log2 length + n_cys for the
    common per-cysteine OR.
(3) the same comparisons under the five detectability windows of s08.
(4) mirror matched design: 1:1 nearest-log10-iBAQ matching without replacement within strata of exact
    total cysteine count (capped at 25): r of n_det, n_undet and length, mean n_det, paired Wilcoxon
    on n_det; 400 group-bootstrap resamples with matching redone (seed SEED + 831 + label index);
    mean r over iBAQ decile x exact n_cys cells (>= 10 per arm).
(5) per-cysteine unit (every detectable cysteine; protein-clustered SEs): n_det and n_undet of the
    protein as separate linear or log2(1 + .) terms, Wald test of equality; total count alone.
(6) fraction of groups with >= 1 identified Cys peptide (S) and in A, by iBAQ decile x tertile of the
    total count and of the detectable count.
Outputs: s10_comparisons.csv, s10_models.csv, s10_linear_symmetric.csv, s10_windows.csv,
         s10_mirror_matched.csv, s10_per_cysteine.csv, s10_tertile_cells.csv, s10_summary.json
"""
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
from s08_lib import decile_codes, match_within, rank_biserial, sorted_strata_index, strat_r, strata_index
from s10_lib import CAP_NCYS

warnings.filterwarnings("ignore")
OUT = RES
LABELS = ("A", "S", "S1", "S_map")
BOOT = 400

F = pd.read_csv(RES / "s08_group_features_windows.csv")
assert len(F) == 7692
F["l2d"] = np.log2(1 + F.n_det)
F["l2c"] = np.log2(1 + F.n_cys)
F["l2u"] = np.log2(1 + F.n_undet)
F["log_det"] = np.log(F.n_det.where(F.n_det > 0))
F["log_cys"] = np.log(F.n_cys.where(F.n_cys > 0))
F["frac_det"] = F.n_det / F.n_cys.where(F.n_cys > 0)
F["ndet20"] = np.minimum(F.n_det, 20)
F["ncys20"] = np.minimum(F.n_cys, 20)
H = F[F.n_det >= 1].reset_index(drop=True)
summary = {"analysis_label": "POST HOC revision analysis 2026-09-30, round 2 (after verifier round 2; not registered)",
           "n_groups": int(len(F)), "n_groups_ndet_ge1": int(len(H)),
           "spearman_ncys_ndet_all": float(stats.spearmanr(F.n_cys, F.n_det)[0]),
           "spearman_ncys_ndet_ndet_ge1": float(stats.spearmanr(H.n_cys, H.n_det)[0]),
           "pearson_log2counts_all": float(np.corrcoef(F.l2c, F.l2d)[0, 1]),
           "frac_cys_detectable": float(F.n_det.sum() / F.n_cys.sum()),
           "frac_groups_no_undetectable": float((F.n_undet == 0).mean()),
           "labels": {k: int(F[k].sum()) for k in LABELS},
           "labels_ndet_ge1": {k: int(H[k].sum()) for k in LABELS}}

LOGIT = sm.families.Binomial()
CLOG = sm.families.Binomial(link=sm.families.links.CLogLog())


def is_conv(m):
    c = getattr(m, "converged", None)
    if c is None:
        c = m.mle_retvals.get("converged", False)
    return bool(c)


def fit(formula, data, fam=LOGIT):
    """GLM by IRLS; if IRLS stops at its iteration cap (it does for the unconstrained cr() spline basis
    with a categorical term, whose columns are collinear with the intercept), refit by Newton-Raphson and
    keep the Newton fit when it converged and is at least as good."""
    m = smf.glm(formula, data=data, family=fam).fit()
    if not is_conv(m):
        m2 = smf.glm(formula, data=data, family=fam).fit(method="newton", maxiter=200)
        if is_conv(m2) and m2.llf >= m.llf - 1e-9:
            m2._fit_method = "newton"
            return m2
    m._fit_method = "irls"
    return m


# ------------------------------------------------------------------ (1) model comparisons -----
FAMS = {
    # name: (data, family, {model: rhs})
    "logit_all": (F, LOGIT, {
        "D0": "log10_ibaq + l2d", "T0": "log10_ibaq + l2c", "B0": "log10_ibaq + l2d + l2c",
        "D": "log10_ibaq + log2_len + l2d", "T": "log10_ibaq + log2_len + l2c",
        "B": "log10_ibaq + log2_len + l2d + l2c"}),
    "logit_ndet_ge1": (H, LOGIT, {
        "D0": "log10_ibaq + l2d", "T0": "log10_ibaq + l2c", "B0": "log10_ibaq + l2d + l2c",
        "D": "log10_ibaq + log2_len + l2d", "T": "log10_ibaq + log2_len + l2c",
        "B": "log10_ibaq + log2_len + l2d + l2c"}),
    "cloglog_ndet_ge1": (H, CLOG, {
        "D0": "log10_ibaq + log_det", "T0": "log10_ibaq + log_cys", "B0": "log10_ibaq + log_det + log_cys",
        "D": "log10_ibaq + log2_len + log_det", "T": "log10_ibaq + log2_len + log_cys",
        "B": "log10_ibaq + log2_len + log_det + log_cys"}),
    "spline_quadratic_ndet_ge1": (H, LOGIT, {
        "D0": "cr(log10_ibaq, df=5) + l2d + I(l2d**2)", "T0": "cr(log10_ibaq, df=5) + l2c + I(l2c**2)",
        "B0": "cr(log10_ibaq, df=5) + l2d + I(l2d**2) + l2c + I(l2c**2)",
        "D": "cr(log10_ibaq, df=5) + log2_len + l2d + I(l2d**2)",
        "T": "cr(log10_ibaq, df=5) + log2_len + l2c + I(l2c**2)",
        "B": "cr(log10_ibaq, df=5) + log2_len + l2d + I(l2d**2) + l2c + I(l2c**2)"}),
    "spline_categorical_ndet_ge1": (H, LOGIT, {
        "D0": "cr(log10_ibaq, df=5) + C(ndet20)", "T0": "cr(log10_ibaq, df=5) + C(ncys20)",
        "D": "cr(log10_ibaq, df=5) + log2_len + C(ndet20)", "T": "cr(log10_ibaq, df=5) + log2_len + C(ncys20)"}),
}
COUNT_TERMS = {"l2d", "l2c", "log_det", "log_cys", "I(l2d ** 2)", "I(l2c ** 2)"}
model_rows, comp_rows = [], []
fits = {}
for fam_name, (dat, fam, models) in FAMS.items():
    for lab in LABELS:
        for mn, rhs in models.items():
            m = fit(f"{lab} ~ {rhs}", dat, fam)
            fits[(fam_name, lab, mn)] = m
            ci = m.conf_int()
            rec = {"family": fam_name, "label": lab, "model": mn, "rhs": rhs, "n": int(m.nobs),
                   "n_pos": int(dat[lab].sum()), "df_model": int(m.df_model), "llf": float(m.llf),
                   "aic": float(m.aic), "converged": is_conv(m), "fit_method": getattr(m, "_fit_method", "irls"),
                   "expected_pos": float(m.fittedvalues.sum())}
            for t in m.params.index:
                if t in COUNT_TERMS or t == "log2_len":
                    rec[f"coef[{t}]"] = float(m.params[t])
                    rec[f"ci[{t}]"] = [float(ci.loc[t, 0]), float(ci.loc[t, 1])]
                    rec[f"p[{t}]"] = float(m.pvalues[t])
            model_rows.append(rec)
        g = lambda mn: fits[(fam_name, lab, mn)]
        row = {"family": fam_name, "label": lab, "n": int(g("D").nobs), "n_pos": int(dat[lab].sum()),
               "AIC_D": float(g("D").aic), "AIC_T": float(g("T").aic), "dAIC_D_minus_T": float(g("D").aic - g("T").aic),
               "AIC_D0": float(g("D0").aic), "AIC_T0": float(g("T0").aic), "dAIC0_D_minus_T": float(g("D0").aic - g("T0").aic)}
        # length at fixed total count and at fixed detectable count
        for base in ("T", "D"):
            lr = 2 * (g(base).llf - g(base + "0").llf)
            row[f"LR_len_given_{base}"] = float(lr)
            row[f"p_len_given_{base}"] = float(stats.chi2.sf(lr, 1))
            c = g(base).conf_int().loc["log2_len"]
            if fam_name.startswith("cloglog"):
                row[f"len_coef_{base}"] = [float(g(base).params["log2_len"]), float(c[0]), float(c[1])]
            else:
                row[f"len_OR_per_doubling_{base}"] = [float(np.exp(g(base).params["log2_len"])), float(np.exp(c[0])), float(np.exp(c[1]))]
        if (fam_name, lab, "B") in fits:
            for base, other, name in (("T", "B", "det_given_total"), ("D", "B", "total_given_det"),
                                      ("T0", "B0", "det_given_total_nolen"), ("D0", "B0", "total_given_det_nolen")):
                lr = 2 * (g(other).llf - g(base).llf); dfd = int(g(other).df_model - g(base).df_model)
                row[f"LR_{name}"] = float(lr); row[f"df_{name}"] = dfd; row[f"p_{name}"] = float(stats.chi2.sf(lr, dfd))
            mb = g("B"); cb = mb.conf_int()
            for t in mb.params.index:
                if t in COUNT_TERMS:
                    row[f"B_coef[{t}]"] = [float(mb.params[t]), float(cb.loc[t, 0]), float(cb.loc[t, 1])]
        if fam_name.startswith("cloglog"):
            for mn in ("D0", "T0", "D", "T"):
                t = "log_det" if mn.startswith("D") else "log_cys"
                c = g(mn).conf_int().loc[t]
                row[f"exponent_{mn}"] = [float(g(mn).params[t]), float(c[0]), float(c[1])]
        comp_rows.append(row)
    # fraction detectable at fixed total count (logit families only, groups with >= 1 Cys)
for lab in LABELS:
    d = F[F.n_cys > 0]
    m = fit(f"{lab} ~ log10_ibaq + log2_len + l2c + frac_det", d)
    c = m.conf_int().loc["frac_det"]
    for r in comp_rows:
        if r["family"] == "logit_all" and r["label"] == lab:
            r["frac_det_given_total_logit"] = [float(m.params["frac_det"]), float(c[0]), float(c[1]), float(m.pvalues["frac_det"])]
MOD = pd.DataFrame(model_rows); MOD.to_csv(OUT / "s10_models.csv", index=False)
COMP = pd.DataFrame(comp_rows); COMP.to_csv(OUT / "s10_comparisons.csv", index=False)


# ------------------------------------------------------------------ (2) symmetric linear counts -----
def symmetric(dat, lab, ndet_col="n_det"):
    d = dat.assign(nd=dat[ndet_col], nu=dat.n_cys - dat[ndet_col])
    mL = fit(f"{lab} ~ log10_ibaq + log2_len + nd + nu", d)
    mT = fit(f"{lab} ~ log10_ibaq + log2_len + n_cys", d)
    V = mL.cov_params(); dd = mL.params["nd"] - mL.params["nu"]
    se = float(np.sqrt(V.loc["nd", "nd"] + V.loc["nu", "nu"] - 2 * V.loc["nd", "nu"]))
    cL = mL.conf_int(); cT = mT.conf_int()
    return {"OR_per_detectable_cys": [float(np.exp(mL.params["nd"])), float(np.exp(cL.loc["nd", 0])), float(np.exp(cL.loc["nd", 1]))],
            "OR_per_undetectable_cys": [float(np.exp(mL.params["nu"])), float(np.exp(cL.loc["nu", 0])), float(np.exp(cL.loc["nu", 1]))],
            "ratio_det_over_undet": float(np.exp(dd)), "z_equal": float(dd / se), "p_equal": float(2 * stats.norm.sf(abs(dd / se))),
            "OR_per_cys_total_only": [float(np.exp(mT.params["n_cys"])), float(np.exp(cT.loc["n_cys", 0])), float(np.exp(cT.loc["n_cys", 1]))],
            "OR_len_per_doubling_total_only": [float(np.exp(mT.params["log2_len"])), float(np.exp(cT.loc["log2_len", 0])), float(np.exp(cT.loc["log2_len", 1]))],
            "p_len_total_only": float(mT.pvalues["log2_len"]),
            "AIC_split": float(mL.aic), "AIC_total": float(mT.aic), "LR_split_vs_total": float(2 * (mL.llf - mT.llf))}


sym_rows = []
for sample, dat in (("all", F), ("ndet_ge1", H)):
    for lab in LABELS:
        sym_rows.append({"sample": sample, "label": lab, **symmetric(dat, lab)})
SYM = pd.DataFrame(sym_rows); SYM.to_csv(OUT / "s10_linear_symmetric.csv", index=False)

# ------------------------------------------------------------------ (3) windows -----
WINDOWS = ["w7_30_mc2", "w6_35_mc2", "w7_40_mc3", "w5_50_mc3", "w7_30_mc2_trypsinP"]
win_rows = []
for wn in WINDOWS:
    d = F.copy()
    d["n_det_w"] = d[f"ndet_{wn}"]; d["l2dw"] = np.log2(1 + d.n_det_w); d["S_map_w"] = d[f"smap_{wn}"]
    for lab in ("A", "S", "S_map"):
        y = "S_map_w" if lab == "S_map" else lab
        mD = fit(f"{y} ~ log10_ibaq + log2_len + l2dw", d); mT = fit(f"{y} ~ log10_ibaq + log2_len + l2c", d)
        mB = fit(f"{y} ~ log10_ibaq + log2_len + l2dw + l2c", d)
        s = symmetric(d.assign(**{lab: d[y]}), lab, ndet_col="n_det_w")
        win_rows.append({"window": wn, "label": lab, "n_pos": int(d[y].sum()),
                         "frac_cys_detectable": float(d.n_det_w.sum() / d.n_cys.sum()),
                         "frac_groups_no_undetectable": float((d.n_cys == d.n_det_w).mean()),
                         "spearman_ncys_ndet": float(stats.spearmanr(d.n_cys, d.n_det_w)[0]),
                         "dAIC_D_minus_T": float(mD.aic - mT.aic),
                         "p_det_given_total": float(stats.chi2.sf(2 * (mB.llf - mT.llf), 1)),
                         "p_total_given_det": float(stats.chi2.sf(2 * (mB.llf - mD.llf), 1)),
                         "OR_per_detectable_cys": s["OR_per_detectable_cys"], "OR_per_undetectable_cys": s["OR_per_undetectable_cys"],
                         "p_equal": s["p_equal"]})
WIN = pd.DataFrame(win_rows); WIN.to_csv(OUT / "s10_windows.csv", index=False)


# ------------------------------------------------------------------ (4) mirror matched design -----
def mirror(df, lab):
    x = df.log10_ibaq.to_numpy(float)
    nc = df.n_cys.to_numpy(float); nd = df.n_det.to_numpy(float); L = df.length.to_numpy(float)
    P, N = match_within(df[lab].to_numpy(int), x, sorted_strata_index(np.minimum(nc, CAP_NCYS).astype(int), x))
    return P, N, {"pairs": int(len(P)), "r_det": rank_biserial(nd[P], nd[N]), "r_undet": rank_biserial(nc[P] - nd[P], nc[N] - nd[N]),
                  "r_len": rank_biserial(L[P], L[N]), "mean_det_pos": float(nd[P].mean()), "mean_det_neg": float(nd[N].mean()),
                  "mean_undet_pos": float((nc[P] - nd[P]).mean()), "mean_undet_neg": float((nc[N] - nd[N]).mean()),
                  "frac_pairs_equal_ncys": float((nc[P] == nc[N]).mean())}


dec_all = decile_codes(F.log10_ibaq.to_numpy(float))
cells_c = strata_index(dec_all * 1000 + np.minimum(F.n_cys.to_numpy(int), CAP_NCYS))
mir_rows = []
for li, lab in enumerate(LABELS):
    P, N, rec = mirror(F, lab)
    nd = F.n_det.to_numpy(float)
    d = nd[P] - nd[N]
    rec["wilcoxon_p_det"] = float(stats.wilcoxon(d[d != 0]).pvalue) if (d != 0).sum() else float("nan")
    y = F[lab].to_numpy(int)
    rec["r_cells_dec_x_ncys_det"], rec["n_cells"] = strat_r(y, nd, cells_c)
    rec["r_cells_dec_x_ncys_len"] = strat_r(y, F.length.to_numpy(float), cells_c)[0]
    rng = np.random.default_rng(SEED + 831 + li)
    slim = F[["log10_ibaq", "length", "n_cys", "n_det", lab]]
    bd, bu, bl = [], [], []
    for _ in range(BOOT):
        b = slim.iloc[rng.integers(0, len(slim), len(slim))].reset_index(drop=True)
        _, _, r = mirror(b, lab)
        bd.append(r["r_det"]); bu.append(r["r_undet"]); bl.append(r["r_len"])
    for k, v in (("r_det", bd), ("r_undet", bu), ("r_len", bl)):
        rec[f"{k}_ci95"] = [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    rec.update({"label": lab, "boot_reps": BOOT, "boot_seed": SEED + 831 + li, "cap_ncys": CAP_NCYS})
    mir_rows.append(rec)
MIR = pd.DataFrame(mir_rows); MIR.to_csv(OUT / "s10_mirror_matched.csv", index=False)

# ------------------------------------------------------------------ (5) per-cysteine unit -----
C = pd.read_csv(RES / "s04_cysteines.csv")
PCD = C[C.detectable].merge(F[["gidx", "log2_len", "n_det", "n_undet", "n_cys"]].rename(columns={"n_cys": "n_cys_p"}), on="gidx")
PCD["y"] = PCD.site.astype(int)
PCD["l2d"] = np.log2(1 + PCD.n_det); PCD["l2u"] = np.log2(1 + PCD.n_undet); PCD["l2c"] = np.log2(1 + PCD.n_cys_p)
pc_rows = []


def pc_fit(name, f, a=None, b=None, data=PCD):
    m = smf.logit(f, data=data).fit(disp=0, cov_type="cluster", cov_kwds={"groups": data["gidx"]}, maxiter=200)
    ci = m.conf_int()
    rec = {"model": name, "formula": f, "n": int(m.nobs), "n_sites": int(data.y.sum()), "n_groups": int(data.gidx.nunique()),
           "aic": float(m.aic)}
    for t in m.params.index:
        if t != "Intercept":
            rec[f"OR[{t}]"] = [float(np.exp(m.params[t])), float(np.exp(ci.loc[t, 0])), float(np.exp(ci.loc[t, 1]))]
    if a is not None:
        V = m.cov_params(); dd = m.params[a] - m.params[b]
        se = float(np.sqrt(V.loc[a, a] + V.loc[b, b] - 2 * V.loc[a, b]))
        rec["equal_test"] = f"{a} == {b}"; rec["z_equal"] = float(dd / se); rec["p_equal"] = float(2 * stats.norm.sf(abs(dd / se)))
    pc_rows.append(rec)


pc_fit("linear_split", "y ~ log10_ibaq + log2_len + n_det + n_undet", "n_det", "n_undet")
pc_fit("linear_split_no_length", "y ~ log10_ibaq + n_det + n_undet", "n_det", "n_undet")
pc_fit("log_split", "y ~ log10_ibaq + log2_len + l2d + l2u", "l2d", "l2u")
pc_fit("log_split_no_length", "y ~ log10_ibaq + l2d + l2u", "l2d", "l2u")
pc_fit("linear_total", "y ~ log10_ibaq + log2_len + n_cys_p")
pc_fit("log_total", "y ~ log10_ibaq + log2_len + l2c")
sub = PCD[(PCD.n_det >= 2) & (PCD.n_cys_p <= 40)]
pc_fit("linear_split_ndet_ge2_ncys_le40", "y ~ log10_ibaq + log2_len + n_det + n_undet", "n_det", "n_undet", data=sub)
PC = pd.DataFrame(pc_rows); PC.to_csv(OUT / "s10_per_cysteine.csv", index=False)

# ------------------------------------------------------------------ (6) decile x tertile tables -----
tt_rows = []
for col in ("n_cys", "n_det"):
    q1, q2 = np.quantile(F[col], [1 / 3, 2 / 3])
    t = np.where(F[col] <= q1, 0, np.where(F[col] <= q2, 1, 2))
    for lab in ("S", "A"):
        tab = F.assign(t=t, dec=dec_all).groupby(["dec", "t"])[lab].agg(["mean", "size"]).reset_index()
        mono = all((np.diff(tab[tab.dec == k].sort_values("t")["mean"].values) >= 0).all() for k in range(10))
        for _, r in tab.iterrows():
            tt_rows.append({"count": col, "cutpoints": [float(q1), float(q2)], "label": lab, "iBAQ_decile": int(r.dec),
                            "tertile": int(r.t), "frac_positive": float(r["mean"]), "n": int(r["size"]),
                            "monotone_in_every_decile": bool(mono)})
TT = pd.DataFrame(tt_rows); TT.to_csv(OUT / "s10_tertile_cells.csv", index=False)

# ------------------------------------------------------------------ summary -----
summary["comparisons"] = COMP.to_dict("records")
summary["linear_symmetric"] = SYM.to_dict("records")
summary["windows"] = WIN.to_dict("records")
summary["mirror_matched"] = MIR.to_dict("records")
summary["per_cysteine"] = PC.to_dict("records")
summary["top_decile_fraction_by_tertile"] = {
    f"{c}|{l}": TT[(TT["count"] == c) & (TT.label == l) & (TT.iBAQ_decile == 9)].sort_values("tertile").frac_positive.tolist()
    for c in ("n_cys", "n_det") for l in ("S", "A")}
write_json(OUT / "s10_summary.json", summary)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 60); pd.set_option("display.max_colwidth", 60)
print(json.dumps({k: v for k, v in summary.items() if not isinstance(v, list)}, indent=1))
print(COMP[["family", "label", "n", "n_pos", "dAIC_D_minus_T", "dAIC0_D_minus_T", "p_det_given_total", "p_total_given_det",
            "p_det_given_total_nolen", "p_total_given_det_nolen", "p_len_given_T"]].round(5).to_string())
print(COMP[[c for c in COMP.columns if c.startswith(("len_", "exponent_", "B_coef", "frac_det"))] + ["family", "label"]].to_string())
print(SYM.to_string())
print(WIN.to_string())
print(MIR.to_string())
print(PC.to_string())
print(json.dumps(summary["top_decile_fraction_by_tertile"], indent=1))
