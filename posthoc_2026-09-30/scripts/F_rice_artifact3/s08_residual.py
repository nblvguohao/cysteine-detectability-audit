"""Step 8 (POST HOC; revision after adversarial verification, 2026-09-30): the cysteine-count residual.

Why this step exists. Verifier round 1 (major 1 and 2) found that the proposed text overstated what the
opportunity adjustment does to the cysteine-count contrast: it quoted the weakest specifications (total
cysteine count next to log2(1 + detectable Cys), 1.02 and 1.01) and read the constant-hazard cloglog
offset model as support for pure opportunity, although the data reject that model (free exponent 0.83).
This step reports the residual in the parametrisation that separates the two kinds of cysteine:
    n_det   = cysteines covered by >=1 peptide of the detectability model (opportunity)
    n_undet = n_cys - n_det (cysteines that cannot yield an identified peptide under the model)
Nothing here is registered; every model below was chosen after the verifier's report.

(1) logistic models, labels A, S, S1, S_map (s04 definitions), primary window (trypsin, 7-30, <=2 mc):
      M3c  y ~ log10 iBAQ + log2(1+n_det)                           (abundance + opportunity only)
      M3   M3c + log2 length + n_cys                                (s04 parametrisation)
      M3u  M3c + log2 length + n_undet
      M3u0 M3c + n_undet
      M3b  M3c + log2 length
      M5u  natural spline in log10 iBAQ (df 5) + exact n_det as categorical (groups with 1..20) + log2 length + n_undet
    likelihood-ratio tests against M3c.
(2) the same under four more detectability windows (6-35/<=2, 7-40/<=3, 5-50/<=3, Trypsin/P 7-30/<=2).
(3) sensitivities of the n_undet term: raw intensity instead of iBAQ; + log2(1 + theoretical iBAQ peptides)
    (Trypsin/P, 0 missed cleavages, 6-30 residues); spline in iBAQ; single-member groups only; sequence AND
    iBAQ of the group member present in UP000059680 2026_03 (instead of the alphabetically first member).
(4) complementary log-log models on groups with >=1 detectable Cys: free exponent on n_det (C0), offset
    = exponent fixed at 1 (C1; LR test of exponent = 1), free + length + n_cys (C4), free + length + n_undet (C5).
(5) 1:1 matching on log10 iBAQ within exact n_det strata (cap 20): r for length, n_cys and n_undet, mean
    n_undet, paired Wilcoxon on n_undet; 400 group-bootstrap resamples with matching redone.
(6) per-cysteine logistic model over detectable cysteines: protein n_det and n_undet as separate terms.
(7) the first-listed accession: PG.ProteinAccessions order and per-member iBAQ segments.
Seeds: bootstrap SEED + 801 + label index. Outputs: s08_*.csv, s08_summary.json.
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
from s08_lib import (Design, features_for_window, match_within, or_rows, rank_biserial,
                     sorted_strata_index)
sys.path.insert(0, str(CYS_AUDIT_SRC))
from cys_audit.proteases import digest

warnings.filterwarnings("ignore")
OUT = RES
summary = {"analysis_label": "POST HOC revision analysis 2026-09-30, after adversarial verification (not registered)",
           "seed_bootstrap_base": SEED + 801}

# ------------------------------------------------------------------ inputs -----
F = pd.read_csv(RES / "s04_group_features.csv")
G = pd.read_csv(RES / "s02_groups.csv")
assert len(F) == len(G) == 7692 and (F.gidx.values == np.arange(len(G))).all()
F["accessions"] = G.accessions.values
F["n_accessions"] = G.n_accessions.values
seqs, src = load_all_sequences()
ref = read_fasta_any(FASTA_REF)
CY = pd.read_csv(RES / "s04_cysteines.csv", usecols=["gidx", "position", "site"])
site_pos = CY[CY.site].groupby("gidx")["position"].apply(lambda v: sorted(v)).to_dict()
F["sites"] = [site_pos.get(g, []) for g in F.gidx]
F["log2_len"] = np.log2(F.length)
F["cys"] = F.n_cys.astype(float)
F["log10_q"] = np.log10(F.quantity.astype(float))

WINDOWS = {"w7_30_mc2": ("trypsin", 7, 30, 2), "w6_35_mc2": ("trypsin", 6, 35, 2),
           "w7_40_mc3": ("trypsin", 7, 40, 3), "w5_50_mc3": ("trypsin", 5, 50, 3),
           "w7_30_mc2_trypsinP": ("trypsin/p", 7, 30, 2)}
for wn, (enz, lo, hi, mc) in WINDOWS.items():
    nd, sm_ = features_for_window(seqs, F.seq_accession, F.sites, digest, enz, lo, hi, mc)
    F[f"ndet_{wn}"] = nd; F[f"smap_{wn}"] = sm_
assert (F["ndet_w7_30_mc2"] == F.n_det_cys).all() and (F["smap_w7_30_mc2"] == F.S_map).all()
F["n_det"] = F.n_det_cys
F["n_undet"] = F.n_cys - F.n_det
F["log2_det1"] = np.log2(1 + F.n_det)
F["n_ibaq_pep"] = [len({seqs[a][x - 1:y] for x, y in digest(seqs[a], "trypsin/p", missed=0, min_len=6, max_len=30)})
                   for a in F.seq_accession]
F["log2_ibaqpep1"] = np.log2(1 + F.n_ibaq_pep)
F["seq_source"] = F.seq_accession.map(lambda a: src.get(a, "none"))

desc = {"n_groups": int(len(F)),
        "spearman_ncys_ndet": float(stats.spearmanr(F.n_cys, F.n_det)[0]),
        "frac_groups_n_undet_0": float((F.n_undet == 0).mean()),
        "median_n_det": float(F.n_det.median()), "mean_n_undet": float(F.n_undet.mean()),
        "labels": {k: int(F[k].sum()) for k in ("A", "S", "S1", "S_map")},
        "seq_source_first_listed_member": F.seq_source.value_counts().to_dict(),
        "single_member_groups": int((F.n_accessions == 1).sum()),
        "single_member_labels": {k: int(F.loc[F.n_accessions == 1, k].sum()) for k in ("A", "S")}}
summary["descriptives"] = desc

# ------------------------------------------------------------------ (7) first-listed member & per-member iBAQ -----
prot = read_tsv(PRO)
IBAQ = [c for c in prot[0] if c.endswith("PG.IBAQ")]


def member_medians(row, n_members):
    vals = [[] for _ in range(n_members)]
    for c in IBAQ:
        raw = (row.get(c) or "").strip()
        if not raw or raw == "NaN":
            continue
        seg = raw.split(";")
        if len(seg) != n_members:
            return None
        for j, t in enumerate(seg):
            try:
                v = float(t)
            except ValueError:
                continue
            if v > 0:
                vals[j].append(v)
    return [float(np.median(v)) if v else np.nan for v in vals]


mm = {}
n_sorted = n_rows = 0
for p in prot:
    accs = [a.strip() for a in (p.get("PG.ProteinAccessions") or "").split(";") if a.strip()]
    if not accs:
        continue
    n_rows += 1; n_sorted += int(accs == sorted(accs))
    mm[";".join(accs)] = (accs, member_medians(p, len(accs)))
ref_member_ibaq, ref_member_acc, max_ratio = [], [], []
for accs_s, ib in zip(F.accessions, F.ibaq):
    accs, med = mm[accs_s]
    assert med is not None and abs(med[0] - ib) <= 1e-9 * max(1.0, ib)
    m = np.array(med, float)
    ok = np.isfinite(m) & (m > 0)
    max_ratio.append(float(np.nanmax(m[ok]) / np.nanmin(m[ok])) if ok.sum() >= 2 else 1.0)
    j = next((k for k, a in enumerate(accs) if a in ref and np.isfinite(m[k]) and m[k] > 0), None)
    ref_member_acc.append(accs[j] if j is not None else None)
    ref_member_ibaq.append(m[j] if j is not None else np.nan)
F["ref_member"] = ref_member_acc
F["ibaq_ref_member"] = ref_member_ibaq
F["member_ibaq_max_over_min"] = max_ratio
multi = F.n_accessions > 1
summary["first_listed_member"] = {
    "report_rows": n_rows, "rows_with_alphabetically_sorted_accessions": n_sorted,
    "groups_multi_member": int(multi.sum()),
    "multi_member_groups_members_iBAQ_differ_gt_2fold": int((F.loc[multi, "member_ibaq_max_over_min"] > 2).sum()),
    "frac_all_groups_members_iBAQ_differ_gt_2fold": float((F.member_ibaq_max_over_min > 2).mean()),
    "groups_with_member_in_UP000059680_2026_03": int(F.ref_member.notna().sum()),
    "first_listed_is_ref_member": int((F.ref_member == F.seq_accession).sum()),
    "first_listed_sequence_source": F.seq_source.value_counts().to_dict()}

# ------------------------------------------------------------------ (1) logistic models -----
LABELS = ("A", "S", "S1", "S_map")
logit_rows, lr_rows = [], []
fits = {}
FORMS = {"M3c": "{y} ~ log10_ibaq + log2_det1",
         "M3": "{y} ~ log10_ibaq + log2_det1 + log2_len + cys",
         "M3u": "{y} ~ log10_ibaq + log2_det1 + log2_len + n_undet",
         "M3u0": "{y} ~ log10_ibaq + log2_det1 + n_undet",
         "M3b": "{y} ~ log10_ibaq + log2_det1 + log2_len"}
for lab in LABELS:
    for mn, f in FORMS.items():
        m = smf.logit(f.format(y=lab), data=F).fit(disp=0, maxiter=200)
        fits[(lab, mn)] = m
        logit_rows += or_rows(m, lab, mn, ["log10_ibaq", "log2_det1", "log2_len", "cys", "n_undet"],
                              {"n_pos": int(F[lab].sum())})
    for big in ("M3", "M3u", "M3u0", "M3b"):
        ms, mb = fits[(lab, "M3c")], fits[(lab, big)]
        lr = 2 * (mb.llf - ms.llf); dfd = int(mb.df_model - ms.df_model)
        lr_rows.append({"label": lab, "small": "M3c", "big": big, "LR": float(lr), "df": dfd,
                        "p": float(stats.chi2.sf(lr, dfd)), "AIC_small": float(ms.aic), "AIC_big": float(mb.aic),
                        "dAIC": float(mb.aic - ms.aic)})
    # M5u: flexible abundance and exact opportunity (categorical), groups with 1..20 detectable Cys
    F20 = F[(F.n_det >= 1) & (F.n_det <= 20)]
    for mn, f in (("M5u", f"{lab} ~ cr(log10_ibaq, df=5) + C(n_det) + log2_len + n_undet"),
                  ("M5c", f"{lab} ~ cr(log10_ibaq, df=5) + C(n_det)")):
        m = smf.glm(f, data=F20, family=sm.families.Binomial()).fit()
        fits[(lab, mn)] = m
        logit_rows += or_rows(m, lab, mn, ["log2_len", "n_undet"], {"n_pos": int(F20[lab].sum())})
    lr = 2 * (fits[(lab, "M5u")].llf - fits[(lab, "M5c")].llf)
    lr_rows.append({"label": lab, "small": "M5c", "big": "M5u", "LR": float(lr), "df": 2,
                    "p": float(stats.chi2.sf(lr, 2)), "AIC_small": float(fits[(lab, "M5c")].aic),
                    "AIC_big": float(fits[(lab, "M5u")].aic), "dAIC": float(fits[(lab, "M5u")].aic - fits[(lab, "M5c")].aic)})
LOG = pd.DataFrame(logit_rows); LOG.to_csv(OUT / "s08_logit_undet.csv", index=False)
LRT = pd.DataFrame(lr_rows); LRT.to_csv(OUT / "s08_logit_LR.csv", index=False)

# ------------------------------------------------------------------ (5) matched pairs + bootstrap -----
def matched(df, lab, key_ndet="n_det"):
    x = df.log10_ibaq.to_numpy(float)
    opp = np.minimum(df[key_ndet].to_numpy(), 20).astype(int)
    P, N = match_within(df[lab].to_numpy(int), x, sorted_strata_index(opp, x))
    L = df.length.to_numpy(float); C = df.n_cys.to_numpy(float); U = C - df[key_ndet].to_numpy(float)
    return P, N, {"pairs": int(len(P)), "r_len": rank_biserial(L[P], L[N]), "r_cys": rank_biserial(C[P], C[N]),
                  "r_undet": rank_biserial(U[P], U[N]), "mean_undet_pos": float(U[P].mean()),
                  "mean_undet_neg": float(U[N].mean()), "frac_pairs_equal_ncys": float((C[P] == C[N]).mean()),
                  "frac_pairs_ndet_equal": float((df[key_ndet].to_numpy()[P] == df[key_ndet].to_numpy()[N]).mean())}


s04b = pd.read_csv(RES / "s04b_matched_abundance_opportunity.csv").set_index("label")
mp_rows = []
for li, lab in enumerate(LABELS):
    P, N, rec = matched(F, lab)
    if lab in s04b.index:     # the new array implementation must reproduce s04b exactly
        assert abs(rec["r_len"] - s04b.loc[lab, "r_length"]) < 1e-12 and abs(rec["r_cys"] - s04b.loc[lab, "r_cys"]) < 1e-12
    U = F.n_undet.to_numpy(float)
    d = U[P] - U[N]
    rec["wilcoxon_p_undet"] = float(stats.wilcoxon(d[d != 0]).pvalue) if (d != 0).sum() > 0 else float("nan")
    rng = np.random.default_rng(SEED + 801 + li)
    bl, bc, bu = [], [], []
    slim = F[["log10_ibaq", "length", "n_cys", "n_det", lab]]
    for _ in range(400):
        b = slim.iloc[rng.integers(0, len(slim), len(slim))].reset_index(drop=True)
        _, _, r = matched(b, lab)
        bl.append(r["r_len"]); bc.append(r["r_cys"]); bu.append(r["r_undet"])
    for k, v in (("r_len", bl), ("r_cys", bc), ("r_undet", bu)):
        rec[f"{k}_ci95"] = [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    rec.update({"label": lab, "boot_reps": 400, "boot_seed": SEED + 801 + li})
    mp_rows.append(rec)
MP = pd.DataFrame(mp_rows); MP.to_csv(OUT / "s08_matched_pairs.csv", index=False)

# ------------------------------------------------------------------ (2) windows -----
w_rows = []
for wn in WINDOWS:
    H = F.copy()
    H["n_det"] = H[f"ndet_{wn}"]; H["n_undet"] = H.n_cys - H.n_det; H["log2_det1"] = np.log2(1 + H.n_det)
    H["S_map"] = H[f"smap_{wn}"]
    for lab in ("A", "S", "S_map"):
        m = smf.logit(FORMS["M3u"].format(y=lab), data=H).fit(disp=0, maxiter=200)
        m0 = smf.logit(FORMS["M3c"].format(y=lab), data=H).fit(disp=0, maxiter=200)
        H20 = H[(H.n_det >= 1) & (H.n_det <= 20)]
        m5 = smf.glm(f"{lab} ~ cr(log10_ibaq, df=5) + C(n_det) + log2_len + n_undet", data=H20,
                     family=sm.families.Binomial()).fit()
        P, N, rec = matched(H, lab)
        c = m.conf_int(); c5 = m5.conf_int()
        w_rows.append({"window": wn, "label": lab, "n_pos": int(H[lab].sum()),
                       "frac_cys_detectable": float(H.n_det.sum() / H.n_cys.sum()),
                       "frac_groups_n_undet_0": float((H.n_undet == 0).mean()),
                       "M3u_OR_undet": float(np.exp(m.params["n_undet"])),
                       "M3u_OR_undet_ci": [float(np.exp(c.loc["n_undet", 0])), float(np.exp(c.loc["n_undet", 1]))],
                       "M3u_p_undet": float(m.pvalues["n_undet"]),
                       "M3u_OR_len": float(np.exp(m.params["log2_len"])),
                       "M3u_OR_len_ci": [float(np.exp(c.loc["log2_len", 0])), float(np.exp(c.loc["log2_len", 1]))],
                       "LR_len_undet_vs_M3c": float(2 * (m.llf - m0.llf)),
                       "p_LR": float(stats.chi2.sf(2 * (m.llf - m0.llf), 2)),
                       "M5u_OR_undet": float(np.exp(m5.params["n_undet"])),
                       "M5u_OR_undet_ci": [float(np.exp(c5.loc["n_undet", 0])), float(np.exp(c5.loc["n_undet", 1]))],
                       **{f"matched_{k}": v for k, v in rec.items()}})
WIN = pd.DataFrame(w_rows); WIN.to_csv(OUT / "s08_windows_undet.csv", index=False)

# ------------------------------------------------------------------ (3) sensitivities -----
sens = []
for lab in ("A", "S", "S_map"):
    for nm, f, dat in (
            ("raw_intensity", f"{lab} ~ log10_q + log2_det1 + log2_len + n_undet", F.dropna(subset=["log10_q"])),
            ("iBAQ_plus_theoretical_peptides", f"{lab} ~ log10_ibaq + log2_ibaqpep1 + log2_det1 + log2_len + n_undet", F),
            ("spline_iBAQ", f"{lab} ~ cr(log10_ibaq, df=5) + log2_det1 + log2_len + n_undet", F),
            ("single_member_groups", FORMS["M3u"].format(y=lab), F[F.n_accessions == 1])):
        m = smf.logit(f, data=dat).fit(disp=0, maxiter=200)
        sens += or_rows(m, lab, nm, ["log2_len", "n_undet", "log2_ibaqpep1"], {"n_pos": int(dat[lab].sum())})
# the member present in UP000059680 2026_03: its sequence AND its own iBAQ segment
R_ = F[F.ref_member.notna()].copy()
nd_r, sm_r = [], []
for acc, st, lead in zip(R_.ref_member, R_.sites, R_.seq_accession):
    s = ref[acc]
    det = np.zeros(len(s), bool)
    for a0, b0 in digest(s, "trypsin", missed=2, min_len=7, max_len=30):
        det[a0 - 1:b0] = True
    cpos = [i for i, ch in enumerate(s) if ch == "C"]
    nd_r.append(int(det[cpos].sum()) if cpos else 0)
R_["length"] = [len(ref[a]) for a in R_.ref_member]
R_["n_cys"] = [ref[a].count("C") for a in R_.ref_member]
R_["n_det"] = nd_r
R_["n_undet"] = R_.n_cys - R_.n_det
R_["log2_len"] = np.log2(R_.length); R_["log2_det1"] = np.log2(1 + R_.n_det)
R_["log10_ibaq"] = np.log10(R_.ibaq_ref_member)
for lab in ("A", "S"):
    m = smf.logit(FORMS["M3u"].format(y=lab), data=R_).fit(disp=0, maxiter=200)
    sens += or_rows(m, lab, "reference_proteome_member_sequence_and_iBAQ", ["log2_len", "n_undet"],
                    {"n_pos": int(R_[lab].sum())})
    P, N, rec = matched(R_, lab)
    sens.append({"label": lab, "model": "reference_proteome_member_matched", "term": "r_undet",
                 "OR": rec["r_undet"], "n": int(len(R_)), "n_pos": int(R_[lab].sum()),
                 "r_len": rec["r_len"], "r_cys": rec["r_cys"]})
SENS = pd.DataFrame(sens); SENS.to_csv(OUT / "s08_sensitivity_undet.csv", index=False)

# ------------------------------------------------------------------ (4) cloglog, free exponent -----
Hp = F[F.n_det >= 1].copy(); Hp["log_det"] = np.log(Hp.n_det)
cl_rows, cl_lr = [], []
for lab in ("S", "S1", "S_map", "A"):
    fam = sm.families.Binomial(link=sm.families.links.CLogLog())
    ms = {}
    for nm, f, off in (("C0_free", f"{lab} ~ log10_ibaq + log_det", None),
                       ("C1_offset_exponent1", f"{lab} ~ log10_ibaq", "log_det"),
                       ("C2_offset_len_cys", f"{lab} ~ log10_ibaq + log2_len + cys", "log_det"),
                       ("C4_free_len_cys", f"{lab} ~ log10_ibaq + log_det + log2_len + cys", None),
                       ("C5_free_len_undet", f"{lab} ~ log10_ibaq + log_det + log2_len + n_undet", None),
                       ("C6_free_undet", f"{lab} ~ log10_ibaq + log_det + n_undet", None)):
        kw = {"offset": Hp[off]} if off else {}
        m = smf.glm(f, data=Hp, family=fam, **kw).fit()
        ms[nm] = m
        ci = m.conf_int()
        for t in m.params.index:
            if t == "Intercept":
                continue
            cl_rows.append({"label": lab, "model": nm, "term": t, "coef": float(m.params[t]),
                            "ci_low": float(ci.loc[t, 0]), "ci_high": float(ci.loc[t, 1]), "p": float(m.pvalues[t]),
                            "aic": float(m.aic), "llf": float(m.llf), "n": int(m.nobs), "n_pos": int(Hp[lab].sum()),
                            "expected_pos": float(m.fittedvalues.sum())})
    for small, big, dfd in (("C1_offset_exponent1", "C0_free", 1), ("C0_free", "C4_free_len_cys", 2),
                            ("C0_free", "C5_free_len_undet", 2), ("C0_free", "C6_free_undet", 1),
                            ("C1_offset_exponent1", "C2_offset_len_cys", 2)):
        lr = 2 * (ms[big].llf - ms[small].llf)
        cl_lr.append({"label": lab, "small": small, "big": big, "LR": float(lr), "df": dfd,
                      "p": float(stats.chi2.sf(lr, dfd)), "dAIC": float(ms[big].aic - ms[small].aic)})
CL = pd.DataFrame(cl_rows); CL.to_csv(OUT / "s08_cloglog_free.csv", index=False)
CLR = pd.DataFrame(cl_lr); CLR.to_csv(OUT / "s08_cloglog_LR.csv", index=False)

# ------------------------------------------------------------------ (6) per-cysteine -----
C = pd.read_csv(RES / "s04_cysteines.csv")
D = C[C.detectable].merge(F[["gidx", "log2_len", "n_det", "n_undet"]], on="gidx")
D["y"] = D.site.astype(int)
pc = []
for nm, f in (("P5_iBAQ_len_ndet_nundet", "y ~ log10_ibaq + log2_len + n_det + n_undet"),
              ("P6_iBAQ_len_log2det_nundet", "y ~ log10_ibaq + log2_len + np.log2(1 + n_det) + n_undet"),
              ("P7_iBAQ_ndet_nundet", "y ~ log10_ibaq + n_det + n_undet")):
    m = smf.logit(f, data=D).fit(disp=0, cov_type="cluster", cov_kwds={"groups": D["gidx"]}, maxiter=200)
    pc += or_rows(m, "per_cysteine", nm, ["log10_ibaq", "log2_len", "n_det", "np.log2(1 + n_det)", "n_undet"],
                  {"n_pos": int(D.y.sum()), "n_groups": int(D.gidx.nunique())})
PC = pd.DataFrame(pc); PC.to_csv(OUT / "s08_per_cysteine_undet.csv", index=False)

# ------------------------------------------------------------------ summary -----
def get(df, **kw):
    q = df
    for k, v in kw.items():
        q = q[q[k] == v]
    assert len(q) == 1, kw
    return q.iloc[0].to_dict()


summary["logit"] = {lab: {mn: {t: {k: get(LOG, label=lab, model=mn, term=t)[k] for k in ("OR", "ci_low", "ci_high", "p")}
                               for t in ("log2_len", "n_undet", "cys") if len(LOG[(LOG.label == lab) & (LOG.model == mn) & (LOG.term == t)])}
                          for mn in ("M3", "M3u", "M5u")} for lab in LABELS}
summary["logit_LR"] = LRT.to_dict("records")
summary["matched_pairs"] = MP.to_dict("records")
summary["windows"] = WIN[["window", "label", "M3u_OR_undet", "M3u_OR_undet_ci", "M5u_OR_undet", "M5u_OR_undet_ci",
                          "M3u_OR_len", "matched_r_undet", "matched_r_cys", "matched_r_len"]].to_dict("records")
summary["sensitivity"] = SENS.to_dict("records")
summary["cloglog_LR"] = CLR.to_dict("records")
summary["cloglog_free_S"] = CL[(CL.label == "S")].to_dict("records")
summary["per_cysteine"] = PC.to_dict("records")
write_json(OUT / "s08_summary.json", summary)
F.drop(columns=["sites"]).to_csv(OUT / "s08_group_features_windows.csv", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(json.dumps(desc, indent=1)); print(json.dumps(summary["first_listed_member"], indent=1))
print(LOG[LOG.term.isin(["log2_len", "cys", "n_undet"])].round(4).to_string())
print(LRT.round(4).to_string())
print(MP.round(4).to_string())
print(WIN.round(4).to_string())
print(SENS.round(4).to_string())
print(CL.round(4).to_string()); print(CLR.round(4).to_string())
print(PC.round(4).to_string())
