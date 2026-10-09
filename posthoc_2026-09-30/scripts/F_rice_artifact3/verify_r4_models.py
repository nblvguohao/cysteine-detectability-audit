"""VERIFIER round 4 (adversarial; POST HOC), item F_rice_artifact3: re-derive with own code the numbers quoted in
the proposed Results / Experimental Procedures / Supplemental Note text (detectable versus total count, the
round-3 indicator test, symmetric per-cysteine ORs, length at fixed count, cloglog exponents, windows,
abundance ratios, rank-biserial effect sizes, Spearman correlations).

Input: verify_r4/groups_r4.csv (verify_r4_build.py). Own MLE (verify_r4_lib.fit_binary), not statsmodels GLM.
Output: verify_r4/models_r4.json, models_r4_stdout.txt (by the caller). Bootstrap seed 4040404 (verifier's own).
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from verify_r4_lib import (deciles, dummies, fit_binary, greedy_match, rank_biserial, spline_basis,  # noqa
                           strat_mean_r)

W = pathlib.Path("/path/to/local/_cys_repo_work/public/revision_2026-09-30")
OUT = W / "results" / "F_rice_artifact3" / "verify_r4"
G = pd.read_csv(OUT / "groups_r4.csv")
assert len(G) == 7692
LABELS = ["A", "S", "S1", "S_map"]
res = {}

# ---------------------------------------------------------------- descriptives
res["spearman"] = {
    "logibaq_length": float(stats.spearmanr(G.log10_ibaq, G.length)[0]),
    "logibaq_ncys": float(stats.spearmanr(G.log10_ibaq, G.n_cys)[0]),
    "logibaq_ndet": float(stats.spearmanr(G.log10_ibaq, G.n_det)[0]),
    "ncys_ndet": float(stats.spearmanr(G.n_cys, G.n_det)[0]),
    "median_ncys": float(G.n_cys.median()), "median_ndet": float(G.n_det.median())}

# ---------------------------------------------------------------- abundance ratios
rng = np.random.default_rng(4040404)
neg = G.A == 0


def ratio(pos_mask, col="ibaq", boot=5000):
    a = G.loc[pos_mask, col].to_numpy(float); b = G.loc[neg, col].dropna().to_numpy(float)
    a = a[np.isfinite(a)]
    r = float(np.median(a) / np.median(b))
    bs = []
    for _ in range(boot):
        ia = rng.integers(0, len(a), len(a)); ib = rng.integers(0, len(b), len(b))
        bs.append(np.median(a[ia]) / np.median(b[ib]))
    return [r, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


res["ratios"] = {
    "A_ibaq": ratio(G.A == 1), "S_ibaq": ratio(G.S == 1), "S1_ibaq": ratio(G.S1 == 1),
    "A_not_S_ibaq": ratio((G.A == 1) & (G.S == 0)),
    "A_qty": ratio(G.A == 1, "quantity"), "S_qty": ratio(G.S == 1, "quantity"),
    "note": "bootstrap here resamples the two arms separately (verifier's own choice), so intervals differ "
            "slightly from the item's joint resampling"}

# ---------------------------------------------------------------- rank-biserial, all 7,692 groups, label A
x = G.log10_ibaq.to_numpy(float); yA = G.A.to_numpy(int)
dec_all = deciles(x)
P, N = greedy_match(yA, x, np.zeros(len(x), int))
rb = {}
for f in ("length", "n_cys"):
    v = G[f].to_numpy(float)
    rb[f] = {"unmatched": rank_biserial(v[yA == 1], v[yA == 0]), "deciles": strat_mean_r(yA, v, dec_all),
             "matched_1to1": rank_biserial(v[P], v[N])}
res["rank_biserial_A_all_groups"] = rb

# ---------------------------------------------------------------- model forms
F = G.copy()
F["l2len"] = np.log2(F.length)
F["l2d"] = np.log2(1 + F.n_det); F["l2c"] = np.log2(1 + F.n_cys)
F["anydet"] = (F.n_det >= 1).astype(float)
H = F[F.n_det >= 1].reset_index(drop=True)
one = lambda d: np.ones(len(d))


def design(d, form, which, with_len=True):
    """which in {'D','T','B'}; form in logit_all, logit_h, cloglog_h, splq_h, splc_h, ind_all."""
    cols = []
    if form in ("splq_h", "splc_h"):
        cols.append(spline_basis(d.log10_ibaq.to_numpy(float)))
    else:
        cols += [one(d)[:, None], d.log10_ibaq.to_numpy(float)[:, None]]
    if with_len:
        cols.append(d.l2len.to_numpy(float)[:, None])
    if form == "cloglog_h":
        cd, cc = np.log(d.n_det.to_numpy(float)), np.log(d.n_cys.to_numpy(float))
        terms = {"D": [cd], "T": [cc], "B": [cd, cc]}[which]
    elif form == "splq_h":
        a, c = d.l2d.to_numpy(float), d.l2c.to_numpy(float)
        terms = {"D": [a, a ** 2], "T": [c, c ** 2], "B": [a, a ** 2, c, c ** 2]}[which]
    elif form == "splc_h":
        v = np.minimum(d.n_det if which == "D" else d.n_cys, 20).to_numpy(int)
        Dm, _ = dummies(v)
        terms = [Dm[:, j] for j in range(Dm.shape[1])]
    else:
        a, c = d.l2d.to_numpy(float), d.l2c.to_numpy(float)
        terms = {"D": [a], "T": [c], "B": [a, c]}[which]
    cols += [np.asarray(t, float)[:, None] for t in terms]
    if form == "ind_all":
        cols.append(d.anydet.to_numpy(float)[:, None])
    return np.hstack(cols)


FORMS = {"logit_all": (F, "logit"), "logit_h": (H, "logit"), "cloglog_h": (H, "cloglog"),
         "splq_h": (H, "logit"), "splc_h": (H, "logit"), "ind_all": (F, "logit")}
comp = {}
for lab in LABELS:
    for form, (d, link) in FORMS.items():
        y = d[lab].to_numpy(float)
        fD = fit_binary(y, design(d, form, "D"), link)
        fT = fit_binary(y, design(d, form, "T"), link)
        rec = {"dAIC_D_minus_T": fD["aic"] - fT["aic"], "conv": [fD["converged"], fT["converged"]]}
        if form != "splc_h":
            fB = fit_binary(y, design(d, form, "B"), link)
            ddf_d = fB["k"] - fT["k"]; ddf_c = fB["k"] - fD["k"]
            rec["p_det_given_total"] = float(stats.chi2.sf(2 * (fB["llf"] - fT["llf"]), ddf_d))
            rec["p_total_given_det"] = float(stats.chi2.sf(2 * (fB["llf"] - fD["llf"]), ddf_c))
            rec["conv"].append(fB["converged"])
            if form == "cloglog_h":
                XD0 = np.column_stack([one(d), d.log10_ibaq, np.log(d.n_det)])
                XT0 = np.column_stack([one(d), d.log10_ibaq, np.log(d.n_cys)])
                g0 = fit_binary(y, XD0, "cloglog"); g1 = fit_binary(y, XT0, "cloglog")
                rec["exponent_det_nolen"] = [float(g0["b"][2]), float(g0["b"][2] - 1.96 * g0["se"][2]),
                                             float(g0["b"][2] + 1.96 * g0["se"][2])]
                rec["exponent_total_nolen"] = [float(g1["b"][2]), float(g1["b"][2] - 1.96 * g1["se"][2]),
                                               float(g1["b"][2] + 1.96 * g1["se"][2])]
                rec["length_coef_at_fixed_total"] = [float(fT["b"][2]), float(fT["b"][2] - 1.96 * fT["se"][2]),
                                                     float(fT["b"][2] + 1.96 * fT["se"][2])]
        if form in ("logit_all",):
            b, se = fT["b"][2], fT["se"][2]
            rec["len_OR_per_doubling_at_fixed_l2c"] = [float(np.exp(b)), float(np.exp(b - 1.96 * se)),
                                                       float(np.exp(b + 1.96 * se)), float(2 * stats.norm.sf(abs(b / se)))]
            b, se = fD["b"][2], fD["se"][2]
            rec["len_OR_per_doubling_at_fixed_l2d"] = [float(np.exp(b)), float(np.exp(b - 1.96 * se)),
                                                       float(np.exp(b + 1.96 * se)), float(2 * stats.norm.sf(abs(b / se)))]
            XL = np.column_stack([one(d), d.log10_ibaq, d.l2len, d.n_cys])
            fl = fit_binary(y, XL, "logit")
            b, se = fl["b"][2], fl["se"][2]
            rec["len_OR_per_doubling_at_fixed_linear_ncys"] = [float(np.exp(b)), float(np.exp(b - 1.96 * se)),
                                                               float(np.exp(b + 1.96 * se))]
        if form == "ind_all":
            # without length as well
            fD0 = fit_binary(y, design(d, form, "D", with_len=False), link)
            fT0 = fit_binary(y, design(d, form, "T", with_len=False), link)
            fB0 = fit_binary(y, design(d, form, "B", with_len=False), link)
            rec["p_det_given_total_nolen"] = float(stats.chi2.sf(2 * (fB0["llf"] - fT0["llf"]), 1))
            rec["n_pos_without_det"] = int(d.loc[d.anydet == 0, lab].sum())
            rec["dAIC_DI_minus_TI"] = rec["dAIC_D_minus_T"]
            if rec["n_pos_without_det"] == 0:
                # separation: the indicator sends the n_det = 0 rows to P = 0 (log-lik 0), so the sup of the
                # likelihood equals the fit on the n_det >= 1 rows; LR from those fits (exact limit)
                yh = H[lab].to_numpy(float)
                hD = fit_binary(yh, design(H, "logit_h", "D"), "logit")
                hT = fit_binary(yh, design(H, "logit_h", "T"), "logit")
                hB = fit_binary(yh, design(H, "logit_h", "B"), "logit")
                rec["separation_limit_p_det_given_total"] = float(stats.chi2.sf(2 * (hB["llf"] - hT["llf"]), 1))
                rec["separation_limit_p_total_given_det"] = float(stats.chi2.sf(2 * (hB["llf"] - hD["llf"]), 1))
                rec["separation_limit_dAIC"] = hD["aic"] - hT["aic"]
        comp[f"{lab}|{form}"] = rec
res["comparisons"] = comp

# ---------------------------------------------------------------- symmetric linear counts
lin = {}
for dname, d in (("all", F), ("h", H)):
    for lab in LABELS:
        y = d[lab].to_numpy(float)
        X = np.column_stack([one(d), d.log10_ibaq, d.l2len, d.n_det, d.n_undet])
        f = fit_binary(y, X, "logit")
        b, cov = f["b"], f["cov"]
        diff = b[3] - b[4]; se_diff = np.sqrt(cov[3, 3] + cov[4, 4] - 2 * cov[3, 4])
        XT = np.column_stack([one(d), d.log10_ibaq, d.l2len, d.n_cys])
        ft = fit_binary(y, XT, "logit")
        lin[f"{lab}|{dname}"] = {
            "OR_det": [float(np.exp(b[3])), float(np.exp(b[3] - 1.96 * f["se"][3])), float(np.exp(b[3] + 1.96 * f["se"][3]))],
            "OR_undet": [float(np.exp(b[4])), float(np.exp(b[4] - 1.96 * f["se"][4])), float(np.exp(b[4] + 1.96 * f["se"][4]))],
            "p_equal": float(2 * stats.norm.sf(abs(diff / se_diff))),
            "OR_total_only": [float(np.exp(ft["b"][3])), float(np.exp(ft["b"][3] - 1.96 * ft["se"][3])),
                              float(np.exp(ft["b"][3] + 1.96 * ft["se"][3]))]}
res["linear_symmetric"] = lin

# ---------------------------------------------------------------- windows (logit, all groups)
win = {}
for tag in ("w6_35", "w7_40", "w5_50", "wTP"):
    d = F.copy()
    d["n_det"] = d[f"ndet_{tag}"]; d["n_undet"] = d.n_cys - d.n_det
    d["l2d"] = np.log2(1 + d.n_det)
    d["S_map"] = d[f"smap_{tag}"]
    rec = {"frac_det": float(d.n_det.sum() / d.n_cys.sum()), "spearman": float(stats.spearmanr(d.n_cys, d.n_det)[0]),
           "frac_groups_no_undet": float((d.n_undet == 0).mean())}
    for lab in ("A", "S", "S_map"):
        y = d[lab].to_numpy(float)
        fD = fit_binary(y, design(d, "logit_all", "D"), "logit")
        fT = fit_binary(y, design(d, "logit_all", "T"), "logit")
        X = np.column_stack([one(d), d.log10_ibaq, d.l2len, d.n_det, d.n_undet])
        f = fit_binary(y, X, "logit")
        diff = f["b"][3] - f["b"][4]; se = np.sqrt(f["cov"][3, 3] + f["cov"][4, 4] - 2 * f["cov"][3, 4])
        rec[lab] = {"dAIC": fD["aic"] - fT["aic"], "OR_det": float(np.exp(f["b"][3])), "OR_undet": float(np.exp(f["b"][4])),
                    "p_equal": float(2 * stats.norm.sf(abs(diff / se)))}
    win[tag] = rec
res["windows"] = win

(OUT / "models_r4.json").write_text(json.dumps(res, indent=1), encoding="utf-8")


def r(v, k=3):
    return round(float(v), k)


print("Spearman:", {k: r(v) for k, v in res["spearman"].items()})
print("ratios:", {k: [r(t, 2) for t in v] for k, v in res["ratios"].items() if k != "note"})
print("rank-biserial A (all groups):", {f: {k: r(v, 4) for k, v in d.items()} for f, d in rb.items()})
for k, v in comp.items():
    print(k, {kk: (r(vv, 5) if isinstance(vv, float) else vv) for kk, vv in v.items()})
for k, v in lin.items():
    print("linear", k, {kk: ([r(t) for t in vv] if isinstance(vv, list) else r(vv, 3)) for kk, vv in v.items()})
for k, v in win.items():
    print("window", k, v)
