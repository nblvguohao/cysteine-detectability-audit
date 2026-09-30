"""Step 4 (POST HOC): Artifact 3 with an explicit opportunity adjustment.

Question (from the pre-submission review): after abundance matching, persulfidated rice protein
groups are longer and cysteine-richer than unmodified groups (rank-biserial +0.201 and +0.290 within
iBAQ deciles). Is that most parsimoniously an OPPORTUNITY effect -- more theoretically detectable
cysteines give more chances of >=1 identified site?

Unit: protein group of the PXD072035 quantification with a usable iBAQ (7,692; s02). Sequence: the
group's LEADING accession (the accession whose iBAQ is used), from UP000059680 2026_03 or fetched
in s01. Opportunity (Cys-Audit v-released digest, trypsin, cleavage after K/R not before P,
0-2 missed cleavages, 7-30 residues -- the manuscript's detectability model):
  n_det_cys       cysteines covered by >=1 such peptide
  n_det_cys_pep   distinct such peptides containing >=1 C
  n_det_pep       distinct such peptides (any)
Labels:
  A   stored definition: any accession of the group is a leading razor protein of the peptide table
      (any peptide, with or without C)                                   -> 1,197 positives
  S   >=1 identified Cys-containing peptide (decoys/contaminants removed) whose leading razor
      protein is an accession of the group ("has >=1 identified site")
  S1  as S, single-Cys peptides only (the rule behind the manuscript's 479 groups)
  S_map  >=1 identified site mapped onto a DETECTABLE cysteine of the leading sequence (the label
      the per-cysteine simulation regenerates, so simulated and observed labels are defined alike)
Abundance: log10 iBAQ (primary, as stored); raw-intensity (PG.Quantity) deciles and models are a
sensitivity, because iBAQ divides by the number of observable peptides.

Analyses:
  1  descriptives and the stored rank-biserials recomputed on the leading-accession sequences
  2  logistic regressions M1 iBAQ; M2 +log2 length +Cys count; M3 +opportunity; ORs, Wald 95% CI,
     likelihood-ratio tests, VIF; spline-in-iBAQ sensitivity
  3  complementary log-log model with offset log(n_det_cys): the pure-opportunity model for S
  4  iBAQ-decile x opportunity-tertile stratification of the rank-biserials (bootstrap CIs)
  5  per-cysteine analysis: P(site | detectable Cys) vs protein abundance, and whether protein
     length / Cys count change the per-cysteine yield at fixed abundance (cluster-robust SEs)
  6  pure-opportunity simulation: labels regenerated from the per-cysteine abundance-only model,
     rank-biserials recomputed -> the contrast opportunity alone would produce
Seeds: 20260930 for every resampling/simulation here.
"""
import math
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
sys.path.insert(0, str(CYS_AUDIT_SRC))
from cys_audit.proteases import digest  # released tool, same digestion rule as the manuscript

warnings.filterwarnings("ignore", category=RuntimeWarning)
OUT = RES
rng_master = np.random.default_rng(SEED)

# ------------------------------------------------------------------ inputs -----
G = pd.read_csv(RES / "s02_groups.csv")
seqs, src = load_all_sequences()
pep = pd.read_csv(PEP, sep="\t", low_memory=False, dtype=str)
pep["C Count"] = pep["C Count"].astype(int)
pep["lrp"] = pep["Leading razor protein"].fillna("")
bad = (pep["Reverse"].fillna("") == "+") | (pep["Potential contaminant"].fillna("") == "+") | \
      pep.lrp.str.startswith("REV__") | pep.lrp.str.startswith("CON__")
real = pep[~bad]
lrp_any = set(pep.lrp) - {""}                       # stored definition uses every LRP string
cys_peps = real[real["C Count"] > 0]
lrp_cys = set(cys_peps.lrp)
lrp_single = set(cys_peps[cys_peps["C Count"] == 1].lrp)
pep_by_lrp = {}
for s, a, c in zip(cys_peps.Sequence, cys_peps.lrp, cys_peps["C Count"]):
    pep_by_lrp.setdefault(a, []).append((s, c))

# ------------------------------------------------------------------ per-group features -----
feat = []
site_rows = []   # per detectable-or-not cysteine of the leading sequence
for gi, r in G.iterrows():
    accs = r["accessions"].split(";")
    lead = accs[0]
    s = seqs.get(lead)
    seq_acc = lead
    if s is None:                                   # fall back to the first member with a sequence
        for a in accs[1:]:
            if a in seqs:
                s, seq_acc = seqs[a], a
                break
    labA = int(r["persulfidated"])
    labS = int(any(a in lrp_cys for a in accs))
    labS1 = int(any(a in lrp_single for a in accs))
    rec = {"gidx": gi, "group_id": lead, "seq_accession": seq_acc if s else None,
           "lead_is_seq": bool(s is not None and seq_acc == lead),
           "A": labA, "S": labS, "S1": labS1, "log10_ibaq": math.log10(r["ibaq_median"]),
           "ibaq": r["ibaq_median"], "quantity": r["quantity_median"]}
    if s is None:
        feat.append(rec); continue
    L = len(s)
    peps = digest(s, "trypsin", missed=2, min_len=7, max_len=30)
    det = np.zeros(L, bool)
    for a0, b0 in peps:
        det[a0 - 1:b0] = True
    cpos = [i + 1 for i, ch in enumerate(s) if ch == "C"]
    pepseqs = {s[a0 - 1:b0] for a0, b0 in peps}
    cpep = {p for p in pepseqs if "C" in p}
    # single-Cys reach: Cys covered by a detectable peptide that contains no other C
    reach1 = set()
    for a0, b0 in peps:
        sub = s[a0 - 1:b0]
        if sub.count("C") == 1:
            reach1.add(a0 + sub.index("C"))
    # identified sites mapped onto this sequence (peptides whose LRP is any member of the group)
    sites, n_unmapped = set(), 0
    for a in accs:
        for ps, cc in pep_by_lrp.get(a, []):
            k = s.find(ps)
            if k < 0:
                n_unmapped += 1; continue
            while k >= 0:
                for i, ch in enumerate(ps):
                    if ch == "C":
                        sites.add(k + 1 + i)
                k = s.find(ps, k + 1)
    rec.update({"length": L, "n_cys": len(cpos), "n_det_cys": int(sum(det[p - 1] for p in cpos)),
                "n_det_cys_pep": len(cpep), "n_det_pep": len(pepseqs), "n_reach1_cys": len(reach1),
                "n_sites_mapped": len(sites), "n_sites_on_det_cys": int(sum(det[p - 1] for p in sites)),
                "n_cys_peptides_unmapped": n_unmapped,
                "S_map": int(sum(det[p - 1] for p in sites) > 0)})
    feat.append(rec)
    for p in cpos:
        site_rows.append((gi, p, bool(det[p - 1]), p in sites, p in reach1))

F = pd.DataFrame(feat)
F["has_seq"] = F["length"].notna()
F.to_csv(OUT / "s04_group_features.csv", index=False)
C = pd.DataFrame(site_rows, columns=["gidx", "position", "detectable", "site", "single_cys_reach"])
C = C.merge(F[["gidx", "log10_ibaq", "length", "n_cys", "n_det_cys", "A", "S"]], on="gidx")
C.to_csv(OUT / "s04_cysteines.csv", index=False)
H = F[F.has_seq].copy()
H["log2_len"] = np.log2(H["length"])
H["log2_cys1"] = np.log2(1 + H["n_cys"])
H["log2_det1"] = np.log2(1 + H["n_det_cys"])
H["log2_detpep1"] = np.log2(1 + H["n_det_cys_pep"])
H["log2_anypep1"] = np.log2(1 + H["n_det_pep"])
H["cys"] = H["n_cys"].astype(float)
H["log10_q"] = np.log10(H["quantity"].astype(float))
H["det"] = H["n_det_cys"].astype(float)

# ------------------------------------------------------------------ helpers -----
def rank_biserial(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    rk = stats.rankdata(np.concatenate([x, y]))
    U = rk[:len(x)].sum() - len(x) * (len(x) + 1) / 2.0
    return float(2.0 * U / (len(x) * len(y)) - 1.0)


def decile_edges(v):
    e = np.percentile(v, np.linspace(0, 100, 11)); e[0] -= 1e-9; e[-1] += 1e-9
    return e


def strat_r(df, lab, key, strata_cols, min_n=10, weighted=False):
    rs, ws = [], []
    for _, sub in df.groupby(strata_cols, observed=True):
        a = sub.loc[sub[lab] == 1, key].values; b = sub.loc[sub[lab] == 0, key].values
        if len(a) < min_n or len(b) < min_n:
            continue
        rs.append(rank_biserial(a, b)); ws.append(len(a) * len(b))
    if not rs:
        return float("nan"), 0
    return (float(np.average(rs, weights=ws)) if weighted else float(np.mean(rs))), len(rs)


def add_strata(df):
    df = df.copy()
    e = decile_edges(df["log10_ibaq"].values)
    df["dec"] = np.clip(np.searchsorted(e, df["log10_ibaq"].values, side="left") - 1, 0, 9)
    # opportunity tertiles on the pooled distribution of n_det_cys (ties kept together)
    q1, q2 = np.quantile(df["n_det_cys"], [1 / 3, 2 / 3])
    df["opp_t"] = np.where(df["n_det_cys"] <= q1, 0, np.where(df["n_det_cys"] <= q2, 1, 2))
    df["opp_exact"] = np.minimum(df["n_det_cys"], 12)
    # raw-intensity (PG.Quantity) deciles, sensitivity: iBAQ divides by the number of observable
    # peptides, so at fixed iBAQ a longer protein carries more total signal
    if "log10_q" in df:
        ok = df["log10_q"].notna().values
        eq = decile_edges(df.loc[ok, "log10_q"].values)
        dq = np.full(len(df), -1)
        dq[ok] = np.clip(np.searchsorted(eq, df.loc[ok, "log10_q"].values, side="left") - 1, 0, 9)
        df["decq"] = dq
    return df, (float(q1), float(q2))


def matched_11(df, lab, key):
    pos = df[df[lab] == 1].sort_values("log10_ibaq", kind="mergesort")   # stable, as sorted() in the original
    neg = df[df[lab] == 0].sort_values("log10_ibaq", kind="mergesort")
    negv = neg["log10_ibaq"].values.copy(); used = np.zeros(len(neg), bool)
    A, B = [], []
    for t, v in zip(pos["log10_ibaq"].values, pos[key].values):
        d = np.abs(negv - t); d[used] = np.inf
        j = int(d.argmin())
        if not np.isfinite(d[j]):
            break
        used[j] = True; A.append(v); B.append(neg[key].values[j])
    return rank_biserial(A, B), len(A)


results = {"analysis_label": "POST HOC revision analysis 2026-09-30 (not registered)", "seed": SEED}

# ------------------------------------------------------------------ 1 descriptives -----
desc = {
    "n_groups": int(len(F)), "n_with_sequence": int(F.has_seq.sum()),
    "n_lead_sequence_used": int(F.lead_is_seq.sum()),
    "seq_source_of_used_sequence": F.loc[F.has_seq, "seq_accession"].map(lambda a: src.get(a, "none")).value_counts().to_dict(),
    "labels_all_groups": {k: int(F[k].sum()) for k in ("A", "S", "S1")},
    "labels_with_sequence": {k: int(H[k].sum()) for k in ("A", "S", "S1", "S_map")},
    "S_positive_with_no_site_mapped_to_a_detectable_cys_of_the_lead_sequence": int(((H.S == 1) & (H.S_map == 0)).sum()),
    "A_positive_without_any_identified_Cys_peptide": int(((F.A == 1) & (F.S == 0)).sum()),
    "S_positive_not_A": int(((F.S == 1) & (F.A == 0)).sum()),
    "S_positive_with_zero_detectable_cys": int(((H.S == 1) & (H.n_det_cys == 0)).sum()),
    "groups_with_zero_detectable_cys": int((H.n_det_cys == 0).sum()),
    "fraction_cys_detectable": float(H.n_det_cys.sum() / H.n_cys.sum()),
    "identified_sites_mapped": int(H.n_sites_mapped.sum()),
    "identified_sites_on_detectable_cys": int(H.n_sites_on_det_cys.sum()),
    "spearman": {f"{a}~{b}": float(stats.spearmanr(H[a], H[b])[0]) for a, b in
                 [("log10_ibaq", "length"), ("log10_ibaq", "n_cys"), ("log10_ibaq", "n_det_cys"),
                  ("length", "n_cys"), ("length", "n_det_cys"), ("n_cys", "n_det_cys"),
                  ("n_det_cys", "n_det_cys_pep"), ("length", "n_det_pep")]},
}
# in-silico expectation of the Cys-peptide share among identified proteins (enrichment level check)
idprot = sorted(a for a in set(real.lrp) if a in seqs)   # sorted: summation order must not depend on hash seed
fr = []
for a in idprot:
    ps = {seqs[a][x - 1:y] for x, y in digest(seqs[a], "trypsin", missed=0, min_len=7, max_len=30)}
    if ps:
        fr.append(sum("C" in p for p in ps) / len(ps))
desc["observed_share_Cys_peptides_in_table_real"] = float((real["C Count"] > 0).mean())
desc["insilico_share_Cys_peptides_fully_tryptic_7_30_identified_proteins_mean"] = float(np.mean(fr))
results["descriptives"] = desc

# stored rank-biserials recomputed on these sequences, for every label and key
Hs, (q1, q2) = add_strata(H)
results["opportunity_tertile_cutpoints_n_det_cys"] = [q1, q2]
rb_rows = []
for lab in ("A", "S", "S1", "S_map"):
    for key in ("length", "n_cys", "n_det_cys", "n_det_cys_pep", "n_det_pep"):
        a = Hs.loc[Hs[lab] == 1, key]; b = Hs.loc[Hs[lab] == 0, key]
        r_dec, nd = strat_r(Hs, lab, key, ["dec"])
        r_dt, ndt = strat_r(Hs, lab, key, ["dec", "opp_t"])
        r_de, nde = strat_r(Hs, lab, key, ["dec", "opp_exact"])
        r_11, npair = matched_11(Hs, lab, key)
        r_q, nq = strat_r(Hs[Hs.decq >= 0], lab, key, ["decq"])
        rb_rows.append({"label": lab, "feature": key, "n_pos": int(len(a)), "n_neg": int(len(b)),
                        "r_before": rank_biserial(a, b), "r_within_iBAQ_deciles": r_dec, "deciles_used": nd,
                        "r_within_rawintensity_deciles": r_q, "deciles_used_raw": nq,
                        "r_1to1_matched": r_11, "pairs": npair,
                        "r_within_decile_x_opp_tertile": r_dt, "cells_used_tertile": ndt,
                        "r_within_decile_x_exact_n_det_cys": r_de, "cells_used_exact": nde})
RB = pd.DataFrame(rb_rows)
RB.to_csv(OUT / "s04_rank_biserial_strata.csv", index=False)

# bootstrap CIs (group resampling, 1,000 reps) for the key stratified summaries
def boot_strat(df, lab, key, cols, reps=1000, seed=SEED):
    rng = np.random.default_rng(seed)
    n = len(df); out = []
    base = df.reset_index(drop=True)
    for _ in range(reps):
        idx = rng.integers(0, n, n)
        b = base.iloc[idx]
        # strata are re-derived within the resample (deciles and tertiles are data-defined)
        b2, _ = add_strata(b)
        out.append(strat_r(b2, lab, key, cols)[0])
    out = np.array(out, float)
    return [float(np.nanpercentile(out, 2.5)), float(np.nanpercentile(out, 97.5))]

ci_rows = []
combo = 0
for lab in ("A", "S"):
    for key in ("length", "n_cys"):
        for cols, nm in ((["dec"], "deciles"), (["dec", "opp_t"], "decile_x_opp_tertile"),
                         (["dec", "opp_exact"], "decile_x_exact_n_det_cys")):
            combo += 1   # deterministic per-combination seed offset
            point = strat_r(Hs, lab, key, cols)[0]
            ci_rows.append({"label": lab, "feature": key, "strata": nm, "r": point,
                            "ci95": boot_strat(Hs, lab, key, cols, reps=400, seed=SEED + combo),
                            "boot_reps": 400, "boot_seed": SEED + combo})
pd.DataFrame(ci_rows).to_csv(OUT / "s04_rank_biserial_bootstrap.csv", index=False)
results["rank_biserial_bootstrap"] = ci_rows

# ------------------------------------------------------------------ 2 logistic regressions -----
def fit_logit(formula, df):
    m = smf.logit(formula, data=df).fit(disp=0, maxiter=200)
    return m

def or_table(m, model_name, lab):
    rows = []
    ci = m.conf_int()
    for term in m.params.index:
        if term == "Intercept" or term.startswith("cr("):
            continue
        rows.append({"label": lab, "model": model_name, "term": term, "OR": float(np.exp(m.params[term])),
                     "ci_low": float(np.exp(ci.loc[term, 0])), "ci_high": float(np.exp(ci.loc[term, 1])),
                     "p": float(m.pvalues[term]), "n": int(m.nobs), "llf": float(m.llf), "aic": float(m.aic)})
    return rows

MODELS = {
    "M1_iBAQ": "{y} ~ log10_ibaq",
    "M2_iBAQ_len_cys": "{y} ~ log10_ibaq + log2_len + cys",
    "M3_iBAQ_len_cys_opp": "{y} ~ log10_ibaq + log2_len + cys + log2_det1",
    "M3b_iBAQ_len_opp": "{y} ~ log10_ibaq + log2_len + log2_det1",
    "M3c_iBAQ_opp": "{y} ~ log10_ibaq + log2_det1",
    "M3d_iBAQ_len_cys_opppep": "{y} ~ log10_ibaq + log2_len + cys + log2_detpep1",
    "M2L_iBAQ_len_log2cys": "{y} ~ log10_ibaq + log2_len + log2_cys1",
    "M3L_iBAQ_len_log2cys_opp": "{y} ~ log10_ibaq + log2_len + log2_cys1 + log2_det1",
    "M4_iBAQ_len_cys_opp_anypep": "{y} ~ log10_ibaq + log2_len + cys + log2_det1 + log2_anypep1",
    "Q2_raw_len_cys": "{y} ~ log10_q + log2_len + cys",
    "Q3_raw_len_cys_opp": "{y} ~ log10_q + log2_len + cys + log2_det1",
    "S2_spline_iBAQ_len_cys": "{y} ~ cr(log10_ibaq, df=5) + log2_len + cys",
    "S3_spline_iBAQ_len_cys_opp": "{y} ~ cr(log10_ibaq, df=5) + log2_len + cys + log2_det1",
}
or_rows, lr_rows, fits = [], [], {}
for lab in ("A", "S", "S1"):
    for mn, f in MODELS.items():
        m = fit_logit(f.format(y=lab), H.dropna(subset=["log10_q"]) if mn.startswith("Q") else H)
        fits[(lab, mn)] = m
        or_rows += or_table(m, mn, lab)
    for small, big in (("M1_iBAQ", "M2_iBAQ_len_cys"), ("M2_iBAQ_len_cys", "M3_iBAQ_len_cys_opp"),
                       ("M3c_iBAQ_opp", "M3_iBAQ_len_cys_opp"), ("S2_spline_iBAQ_len_cys", "S3_spline_iBAQ_len_cys_opp"),
                       ("M3_iBAQ_len_cys_opp", "M4_iBAQ_len_cys_opp_anypep")):
        ms, mb = fits[(lab, small)], fits[(lab, big)]
        lr = 2 * (mb.llf - ms.llf); dfd = mb.df_model - ms.df_model
        lr_rows.append({"label": lab, "small": small, "big": big, "LR": float(lr), "df": int(dfd),
                        "p": float(stats.chi2.sf(lr, dfd)), "delta_AIC_big_minus_small": float(mb.aic - ms.aic)})
OR = pd.DataFrame(or_rows); OR.to_csv(OUT / "s04_logistic_OR.csv", index=False)
LR = pd.DataFrame(lr_rows); LR.to_csv(OUT / "s04_logistic_LR.csv", index=False)
X = H[["log10_ibaq", "log2_len", "cys", "log2_det1"]].copy(); X.insert(0, "const", 1.0)
results["VIF_M3"] = {c: float(variance_inflation_factor(X.values, i)) for i, c in enumerate(X.columns) if c != "const"}
# standardized ORs (per 1 SD) for M2/M3, label A and S
std_rows = []
Z = H.copy()
for c in ("log10_ibaq", "log2_len", "cys", "log2_det1"):
    Z[c] = (Z[c] - Z[c].mean()) / Z[c].std()
for lab in ("A", "S"):
    for mn in ("M2_iBAQ_len_cys", "M3_iBAQ_len_cys_opp"):
        std_rows += or_table(fit_logit(MODELS[mn].format(y=lab), Z), mn + "_perSD", lab)
pd.DataFrame(std_rows).to_csv(OUT / "s04_logistic_OR_perSD.csv", index=False)

# ------------------------------------------------------------------ 3 cloglog offset model -----
cl_rows = []
Hp = H[H.n_det_cys > 0].copy()
Hp["log_det"] = np.log(Hp["n_det_cys"])
for lab in ("S", "S1", "A"):
    for nm, f, off in (
        ("C0_free_opp", f"{lab} ~ log10_ibaq + log_det", None),
        ("C1_offset", f"{lab} ~ log10_ibaq", "log_det"),
        ("C2_offset_len_cys", f"{lab} ~ log10_ibaq + log2_len + cys", "log_det"),
        ("C3_offset_spline_len_cys", f"{lab} ~ cr(log10_ibaq, df=5) + log2_len + cys", "log_det"),
        ("C4_free_opp_len_cys", f"{lab} ~ log10_ibaq + log_det + log2_len + cys", None)):
        kw = {"offset": Hp[off]} if off else {}
        m = smf.glm(f, data=Hp, family=sm.families.Binomial(link=sm.families.links.CLogLog()), **kw).fit()
        ci = m.conf_int()
        for term in m.params.index:
            if term == "Intercept" or term.startswith("cr("):
                continue
            cl_rows.append({"label": lab, "model": nm, "term": term, "coef": float(m.params[term]),
                            "ci_low": float(ci.loc[term, 0]), "ci_high": float(ci.loc[term, 1]),
                            "exp_coef": float(np.exp(m.params[term])), "p": float(m.pvalues[term]),
                            "n": int(m.nobs), "aic": float(m.aic), "offset": off or ""})
CL = pd.DataFrame(cl_rows); CL.to_csv(OUT / "s04_cloglog_opportunity.csv", index=False)
results["n_groups_with_det_cys_ge1"] = int(len(Hp))
results["positives_excluded_zero_det_cys"] = {lab: int(((H.n_det_cys == 0) & (H[lab] == 1)).sum()) for lab in ("A", "S", "S1")}

# ------------------------------------------------------------------ 4 decile x tertile table -----
cell = []
for (d, t), sub in Hs.groupby(["dec", "opp_t"]):
    cell.append({"iBAQ_decile": int(d) + 1, "opp_tertile": int(t) + 1, "n": int(len(sub)),
                 "n_A": int(sub.A.sum()), "frac_A": float(sub.A.mean()),
                 "n_S": int(sub.S.sum()), "frac_S": float(sub.S.mean()),
                 "median_n_det_cys": float(sub.n_det_cys.median()),
                 "r_len_A": rank_biserial(sub.loc[sub.A == 1, "length"], sub.loc[sub.A == 0, "length"]) if min(sub.A.sum(), (1 - sub.A).sum()) >= 10 else float("nan"),
                 "r_cys_A": rank_biserial(sub.loc[sub.A == 1, "n_cys"], sub.loc[sub.A == 0, "n_cys"]) if min(sub.A.sum(), (1 - sub.A).sum()) >= 10 else float("nan"),
                 "r_len_S": rank_biserial(sub.loc[sub.S == 1, "length"], sub.loc[sub.S == 0, "length"]) if min(sub.S.sum(), (1 - sub.S).sum()) >= 10 else float("nan"),
                 "r_cys_S": rank_biserial(sub.loc[sub.S == 1, "n_cys"], sub.loc[sub.S == 0, "n_cys"]) if min(sub.S.sum(), (1 - sub.S).sum()) >= 10 else float("nan")})
pd.DataFrame(cell).to_csv(OUT / "s04_decile_x_tertile_cells.csv", index=False)

# ------------------------------------------------------------------ 5 per-cysteine analysis -----
D = C[C.detectable & C.gidx.isin(H.gidx)].copy()
D = D.merge(H[["gidx", "log2_len", "cys"]], on="gidx")
D["y"] = D["site"].astype(int)
D["log2_cys1"] = np.log2(1 + D["n_cys"])
e = decile_edges(H["log10_ibaq"].values)
D["dec"] = np.clip(np.searchsorted(e, D["log10_ibaq"].values, side="left") - 1, 0, 9)
pc_dec = D.groupby("dec").agg(n_det_cys=("y", "size"), n_sites=("y", "sum")).reset_index()
pc_dec["p_site"] = pc_dec.n_sites / pc_dec.n_det_cys
pc_dec["iBAQ_decile"] = pc_dec.dec + 1
pc_dec.drop(columns="dec").to_csv(OUT / "s04_per_cysteine_by_decile.csv", index=False)
pcr = []
for nm, f in (("P1_iBAQ", "y ~ log10_ibaq"), ("P2_iBAQ_len_cys", "y ~ log10_ibaq + log2_len + cys"),
              ("P3_spline_len_cys", "y ~ cr(log10_ibaq, df=5) + log2_len + cys"),
              ("P4_iBAQ_len_log2cys", "y ~ log10_ibaq + log2_len + log2_cys1")):
    m = smf.logit(f, data=D).fit(disp=0, cov_type="cluster", cov_kwds={"groups": D["gidx"]}, maxiter=200)
    ci = m.conf_int()
    for term in m.params.index:
        if term == "Intercept" or term.startswith("cr("):
            continue
        pcr.append({"model": nm, "term": term, "OR": float(np.exp(m.params[term])),
                    "ci_low": float(np.exp(ci.loc[term, 0])), "ci_high": float(np.exp(ci.loc[term, 1])),
                    "p": float(m.pvalues[term]), "n_cys": int(m.nobs), "n_groups": int(D.gidx.nunique())})
pd.DataFrame(pcr).to_csv(OUT / "s04_per_cysteine_logit.csv", index=False)
results["per_cysteine"] = {"n_detectable_cys": int(len(D)), "n_sites": int(D.y.sum()),
                           "p_site_overall": float(D.y.mean())}

# ------------------------------------------------------------------ 6 pure-opportunity simulation -----
# per-cysteine abundance-only model (spline), fitted on detectable cysteines; labels S regenerated
m0 = smf.logit("y ~ cr(log10_ibaq, df=5)", data=D).fit(disp=0, maxiter=200)
D["p_hat"] = m0.predict(D)
# expected protein-level probability of >=1 site
exp_g = D.groupby("gidx")["p_hat"].apply(lambda v: 1 - np.prod(1 - v.values))
Hs["p_exp_S"] = Hs["gidx"].map(exp_g).fillna(0.0)
results["expected_vs_observed_S"] = {"expected_n_S_map": float(Hs.p_exp_S.sum()), "observed_n_S_map": int(Hs.S_map.sum()),
                                     "observed_n_S": int(Hs.S.sum())}
# calibration of the opportunity model by length / Cys-count quintile within iBAQ tertile
cal = []
Hs["len_q"] = pd.qcut(Hs["length"], 5, labels=False, duplicates="drop")
Hs["cys_q"] = pd.qcut(Hs["n_cys"].rank(method="first"), 5, labels=False)
Hs["ibaq_t"] = pd.qcut(Hs["log10_ibaq"], 3, labels=False)
for var in ("len_q", "cys_q"):
    for (t, q), sub in Hs.groupby(["ibaq_t", var]):
        cal.append({"stratifier": var, "iBAQ_tertile": int(t) + 1, "quintile": int(q) + 1, "n": int(len(sub)),
                    "observed_S_map": int(sub.S_map.sum()), "expected_S_map_opportunity_only": float(sub.p_exp_S.sum()),
                    "obs_over_exp": float(sub.S_map.sum() / sub.p_exp_S.sum()) if sub.p_exp_S.sum() > 0 else float("nan")})
pd.DataFrame(cal).to_csv(OUT / "s04_opportunity_calibration.csv", index=False)

# simulation
reps = 1000
rng = np.random.default_rng(SEED)
gid = D["gidx"].values
ph = D["p_hat"].values
order_g = np.array(sorted(Hs.gidx.values))
pos_in = {g: i for i, g in enumerate(order_g)}
gi_idx = np.array([pos_in[g] for g in gid])
Hsim = Hs.set_index("gidx").loc[order_g].reset_index()
sim = {"len": [], "cys": [], "len_11": [], "cys_11": [], "len_dt": [], "cys_dt": []}
for k in range(reps):
    hit = rng.random(len(ph)) < ph
    lab = np.zeros(len(order_g), int)
    np.maximum.at(lab, gi_idx, hit.astype(int))
    Hsim["Ssim"] = lab
    sim["len"].append(strat_r(Hsim, "Ssim", "length", ["dec"])[0])
    sim["cys"].append(strat_r(Hsim, "Ssim", "n_cys", ["dec"])[0])
    if k < 200:
        sim["len_dt"].append(strat_r(Hsim, "Ssim", "length", ["dec", "opp_t"])[0])
        sim["cys_dt"].append(strat_r(Hsim, "Ssim", "n_cys", ["dec", "opp_t"])[0])
    if k < 100:
        sim["len_11"].append(matched_11(Hsim, "Ssim", "length")[0])
        sim["cys_11"].append(matched_11(Hsim, "Ssim", "n_cys")[0])
obsS = RB.set_index(["label", "feature"])
simsum = {}
for key, feat, col in (("len", "length", "r_within_iBAQ_deciles"), ("cys", "n_cys", "r_within_iBAQ_deciles"),
                       ("len_dt", "length", "r_within_decile_x_opp_tertile"), ("cys_dt", "n_cys", "r_within_decile_x_opp_tertile"),
                       ("len_11", "length", "r_1to1_matched"), ("cys_11", "n_cys", "r_1to1_matched")):
    v = np.array(sim[key], float)
    simsum[key] = {"observed_S_map": float(obsS.loc[("S_map", feat), col]),
                   "observed_S": float(obsS.loc[("S", feat), col]), "sim_mean": float(np.nanmean(v)),
                   "sim_2.5": float(np.nanpercentile(v, 2.5)), "sim_97.5": float(np.nanpercentile(v, 97.5)),
                   "reps": int(np.isfinite(v).sum())}
results["pure_opportunity_simulation_label_S"] = simsum
pd.DataFrame(sim["len"], columns=["r_len_sim"]).assign(r_cys_sim=sim["cys"]).to_csv(OUT / "s04_simulation_draws.csv", index=False)

write_json(OUT / "s04_opportunity_summary.json", results)
print(json.dumps(results, indent=1)[:6000])
print(RB.to_string())
print(OR[OR.model.isin(["M1_iBAQ", "M2_iBAQ_len_cys", "M3_iBAQ_len_cys_opp"])].to_string())
print(LR.to_string())
print(CL.to_string())
