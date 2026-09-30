"""VERIFIER round 4 (adversarial; POST HOC), item F_rice_artifact3: extra checks behind the verifier's findings.
(a) Do the 447 groups represented only by cysteine-free peptides (A and not S) show the cysteine-count
    dependence that the proposed 'protein-level capture' explanation implies? (logistic at fixed abundance
    [and length], within-decile rank-biserial; S versus negatives for comparison)
(b) How many of the 6,495 'unmodified' groups (the comparison arm described as 'absent from that table') have an
    accession in the peptide table's Proteins column (any peptide; cysteine-containing peptide)?
(c) The item's categorical-spline null fits (IRLS stopped at its cap) against the verifier's converged fits.
Output: verify_r4/extra_r4.json (+ stdout). Deterministic."""
import json
import pathlib
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from verify_r4_lib import deciles, fit_binary, strat_mean_r  # noqa
from verify_r4_null_sim import null_probs  # noqa

W = pathlib.Path("C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30")
OUT = W / "results" / "F_rice_artifact3" / "verify_r4"
G = pd.read_csv(OUT / "groups_r4.csv")
res = {}

# (a)
def or_last(y, X):
    f = fit_binary(y, X, "logit"); b, se = f["b"][-1], f["se"][-1]
    return [float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)), float(2 * stats.norm.sf(abs(b / se)))]

a = {}
for name, pos in (("A_not_S", (G.A == 1) & (G.S == 0)), ("S", G.S == 1)):
    D = G[(G.A == 0) | pos].copy(); y = pos[D.index].astype(float).to_numpy()
    one = np.ones(len(D)); l2len = np.log2(D.length.to_numpy(float))
    a[name] = {
        "n_pos": int(y.sum()), "n_neg": int((y == 0).sum()),
        "OR_per_cys_linear_given_ibaq_len": or_last(y, np.column_stack([one, D.log10_ibaq, l2len, D.n_cys])),
        "OR_per_doubling_l2c_given_ibaq_len": or_last(y, np.column_stack([one, D.log10_ibaq, l2len, np.log2(1 + D.n_cys)])),
        "OR_per_doubling_length_given_ibaq_l2c": or_last(y, np.column_stack([one, D.log10_ibaq, np.log2(1 + D.n_cys), l2len])),
        "within_decile_r_ncys": strat_mean_r(y.astype(int), D.n_cys.to_numpy(float), deciles(D.log10_ibaq.to_numpy(float))),
        "within_decile_r_length": strat_mean_r(y.astype(int), D.length.to_numpy(float), deciles(D.log10_ibaq.to_numpy(float)))}
res["cysteine_dependence_of_positive_subsets"] = a

# (b)
pep = pd.read_csv(W / "external" / "SS-all-peptides.tsv", sep="\t", dtype=str, keep_default_na=False)
lrp = pep["Leading razor protein"]
bad = (pep.Reverse == "+") | lrp.str.startswith("REV__") | (pep["Potential contaminant"] == "+") | lrp.str.startswith("CON__")
real = pep[~bad]
p_any, p_cys = set(), set()
for s, cc in zip(real.Proteins, real["C Count"].astype(int)):
    t = {x.strip() for x in s.split(";") if x.strip()}
    p_any |= t
    if cc > 0:
        p_cys |= t
has = lambda accs, S: any(x in S for x in accs.split(";"))
neg = G[G.A == 0]
res["comparison_arm"] = {"n": int(len(neg)),
                         "with_accession_in_Proteins_column_any_peptide": int(neg.accessions.apply(lambda s: has(s, p_any)).sum()),
                         "with_accession_in_Proteins_column_cys_peptide": int(neg.accessions.apply(lambda s: has(s, p_cys)).sum())}
AnS = G[(G.A == 1) & (G.S == 0)]
res["A_not_S_listed_by_a_cys_peptide_in_Proteins"] = int(AnS.accessions.apply(lambda s: has(s, p_cys)).sum())

# (c)
H = G[G.n_det >= 1].reset_index(drop=True)
H["opp20"] = np.minimum(H.n_det, 20); H["cys20"] = np.minimum(H.n_cys, 20)
fits = {}
for lab in ("A", "S", "S_map"):
    for null, rhs in (("N2", "cr(log10_ibaq, df=5) + C(opp20)"), ("N4", "cr(log10_ibaq, df=5) + C(cys20)")):
        m = smf.glm(f"{lab} ~ {rhs}", data=H, family=sm.families.Binomial()).fit()
        p_mine, f = null_probs(H, lab, null)
        fits[f"{lab}|{null}"] = {"item_irls_converged": bool(m.converged),
                                 "max_abs_diff_fitted_p": float(np.max(np.abs(m.fittedvalues.to_numpy(float) - p_mine))),
                                 "llf_item": float(m.llf), "llf_converged": float(f["llf"])}
res["null_fit_irls_vs_converged"] = fits
(OUT / "extra_r4.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
print(json.dumps(res, indent=1))
