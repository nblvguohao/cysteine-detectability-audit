"""Step 13 (POST HOC; revision round 3, 2026-09-30): does the detectable count add anything given the total
count once the requirement of at least one detectable cysteine is modelled explicitly?

Why this step exists. Over all 7,692 groups, the round-2 likelihood-ratio test of log2(1 + n_det) given
log2(1 + n_cys) (s10, logistic, with log10 iBAQ and log2 length) gave P = 0.005 for the site-level label S.
The proposed Results text attributes this to the structural requirement that a site-level positive have at
least one detectable cysteine (3 of the 174 groups whose cysteines are all undetectable are S-positive),
supported so far only by the subset analysis (P = 0.10 on the 6,990 groups with n_det >= 1). The round-3
verifier tested the attribution directly by adding an indicator 1[n_det >= 1] to both models (P = 0.149 for
S). This step reproduces that test with this item's code for the four labels, logistic, all 7,692 groups:
    TI  y ~ log10 iBAQ + log2 length + log2(1 + n_cys) + 1[n_det >= 1]
    DI  y ~ log10 iBAQ + log2 length + log2(1 + n_det) + 1[n_det >= 1]
    BI  y ~ log10 iBAQ + log2 length + log2(1 + n_det) + log2(1 + n_cys) + 1[n_det >= 1]
and the same without log2 length (TI0, DI0, BI0). Reported: LR of the detectable count given total count and
indicator (BI vs TI), LR of the total count given detectable count and indicator (BI vs DI), both chi-square
1 df, and AIC(DI) - AIC(TI). Fitting as in s10 (GLM by IRLS, Newton-Raphson if IRLS stops at its cap).
For S_map, groups without a detectable cysteine cannot be positive by construction, so the indicator's
coefficient diverges (quasi-complete separation); the likelihoods and the LR of the count terms are still
defined, and the fits are reported with their convergence flags.
Outputs: s13_indicator_test.csv; stdout saved as s13_stdout.txt by the caller.
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

warnings.filterwarnings("ignore")
LABELS = ("A", "S", "S1", "S_map")
LOGIT = sm.families.Binomial()


def is_conv(m):
    c = getattr(m, "converged", None)
    if c is None:
        c = m.mle_retvals.get("converged", False)
    return bool(c)


def fit(formula, data):
    m = smf.glm(formula, data=data, family=LOGIT).fit()
    if not is_conv(m):
        m2 = smf.glm(formula, data=data, family=LOGIT).fit(method="newton", maxiter=200)
        if is_conv(m2) and m2.llf >= m.llf - 1e-9:
            m2._fit_method = "newton"
            return m2
    m._fit_method = "irls"
    return m


def main():
    F = pd.read_csv(RES / "s08_group_features_windows.csv")
    assert len(F) == 7692
    F["l2d"] = np.log2(1 + F.n_det)
    F["l2c"] = np.log2(1 + F.n_cys)
    F["anydet"] = (F.n_det >= 1).astype(int)
    rows = []
    for lab in LABELS:
        mods = {}
        for name, rhs in {"TI": "log10_ibaq + log2_len + l2c + anydet", "DI": "log10_ibaq + log2_len + l2d + anydet",
                          "BI": "log10_ibaq + log2_len + l2d + l2c + anydet",
                          "TI0": "log10_ibaq + l2c + anydet", "DI0": "log10_ibaq + l2d + anydet",
                          "BI0": "log10_ibaq + l2d + l2c + anydet"}.items():
            mods[name] = fit(f"{lab} ~ {rhs}", F)
        rec = {"label": lab, "n": int(len(F)), "n_pos": int(F[lab].sum()),
               "n_pos_without_detectable_cys": int(F.loc[F.anydet == 0, lab].sum()),
               "n_groups_without_detectable_cys": int((F.anydet == 0).sum())}
        for sfx in ("", "0"):
            lr_d = 2 * (mods["BI" + sfx].llf - mods["TI" + sfx].llf)
            lr_c = 2 * (mods["BI" + sfx].llf - mods["DI" + sfx].llf)
            tag = "" if sfx == "" else "_nolen"
            rec[f"LR_det_given_total_ind{tag}"] = float(lr_d)
            rec[f"p_det_given_total_ind{tag}"] = float(stats.chi2.sf(lr_d, 1))
            rec[f"LR_total_given_det_ind{tag}"] = float(lr_c)
            rec[f"p_total_given_det_ind{tag}"] = float(stats.chi2.sf(lr_c, 1))
            rec[f"dAIC_DI_minus_TI{tag}"] = float(mods["DI" + sfx].aic - mods["TI" + sfx].aic)
        mb = mods["BI"]
        ci = mb.conf_int()
        rec["BI_coef_l2d"] = [float(mb.params["l2d"]), float(ci.loc["l2d", 0]), float(ci.loc["l2d", 1])]
        rec["BI_coef_l2c"] = [float(mb.params["l2c"]), float(ci.loc["l2c", 0]), float(ci.loc["l2c", 1])]
        rec["BI_coef_anydet"] = float(mb.params["anydet"])
        rec["converged"] = {k: is_conv(m) for k, m in mods.items()}
        rec["fit_method"] = {k: getattr(m, "_fit_method", "irls") for k, m in mods.items()}
        rows.append(rec)
    R = pd.DataFrame(rows)
    R.to_csv(RES / "s13_indicator_test.csv", index=False, lineterminator="\n")
    print("POST HOC round 3: detectable count given total count and 1[n_det >= 1], logistic, all 7,692 groups")
    for r in rows:
        print(f"{r['label']}: n_pos {r['n_pos']} ({r['n_pos_without_detectable_cys']} of "
              f"{r['n_groups_without_detectable_cys']} groups without a detectable Cys); "
              f"P(det | total, ind) = {r['p_det_given_total_ind']:.4g} (no length {r['p_det_given_total_ind_nolen']:.4g}); "
              f"P(total | det, ind) = {r['p_total_given_det_ind']:.3g}; dAIC(DI - TI) = {r['dAIC_DI_minus_TI']:.2f}; "
              f"converged {r['converged']}; methods {r['fit_method']}")


if __name__ == "__main__":
    main()
