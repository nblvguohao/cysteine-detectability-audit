"""Step 3 (POST HOC revision analysis, item D_empirical_background): cleavage-geometry log2 odds ratios of the
PXD063463 trypsin ABE arm under progressively better-matched backgrounds.

Statistic (unchanged from Cys-Audit / Fig. 2): Haldane log2 odds ratio of carrying a feature, positives against
background, protein-clustered percentile bootstrap (5,000 replicates, fixed seed). Features:
  kr_dist / kr_prox   trypsin-competent K/R (not before P) within 6-12 / 1-3 residues (the Fig. 2 statistic)
  de_dist / de_prox   D/E within 6-12 / 1-3 residues (negative-control feature; not a trypsin residue)

Backgrounds (the pool is always the label-0 cysteines of the identified proteins):
  proteome            all label-0 cysteines (Fig. 2a)
  observed            label-0 cysteines identified in the same arm (Fig. 2b)
  theoretical         label-0 cysteines in >= 1 tryptic peptide of 7-30 residues with <= 2 missed cleavages
  empirical           label-0 cysteines in >= 1 candidate peptide identified in the global (unenriched) arm
  match_emp_lr        1:1 nearest-neighbour match (logit of the cysteine-level empirical score, lr_full OOF),
                      greedy, without replacement, random order, ties broken at random
  match_theo          the same matching on the binary theoretical flag (comparator)
  + sensitivity designs (symmetric restriction, other scores, caliper, within-protein matching, matching on the
    empirical-observation flag, Mantel-Haenszel over score strata, HA-specific positive sets, 200 matching seeds)

Intervals:
  per-design intervals use the Cys-Audit convention: clusters = proteins contributing rows to that comparison,
  seeds 20260922 + 101 (kr_prox), + 102 (kr_dist) as in the tool, + 103 (de_prox), + 104 (de_dist); both the 95%
  interval and the Bonferroni 97.5% interval of Fig. 2 are read from the same replicates.
  fraction of the proteome-background signal removed, and differences between backgrounds, use a paired bootstrap:
  one set of protein multiplicities (all 1,991 proteins, seed 20260930 + 7) applied to every design.
  Matching is done once per design (primary matching seed 20260930 + 11); the bootstrap does not re-match.
  Matching variability is shown separately over 200 matching seeds.
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json

import numpy as np
import pandas as pd

from common import BOOT_REPS, MASTER_SEED, RESULTS, TOOL_SEED, dump_json
from cys_audit import stats  # noqa: E402
from cys_audit.stats import auc as tie_auc
from matching import logit, match_design as _match_design

FEATURES = ["kr_prox", "kr_dist", "de_prox", "de_dist"]
FEATURE_SEED = {"kr_prox": TOOL_SEED + 101, "kr_dist": TOOL_SEED + 102, "de_prox": TOOL_SEED + 103,
                "de_dist": TOOL_SEED + 104}
PAIRED_SEED = MASTER_SEED + 7
MATCH_SEED = MASTER_SEED + 11
N_MATCH_SEEDS = 200


def match_design(T, pos_mask, score, rng, caliper=None, within_protein=False):
    return _match_design(T.label.to_numpy(), T.protein.to_numpy(), pos_mask, score, rng, caliper, within_protein)


def balance(T, role, score_logit):
    p, c = role == 1, role == 0
    sd = np.sqrt((score_logit[p].var(ddof=1) + score_logit[c].var(ddof=1)) / 2) if c.sum() > 1 else np.nan
    return {"smd_logit_score_lr_full": float((score_logit[p].mean() - score_logit[c].mean()) / sd) if sd > 0 else None,
            "auc_score_lr_full_pos_vs_bg": tie_auc(score_logit[p | c], p[p | c]),
            "share_theo_trypsin_pos": float(T.theo_trypsin.to_numpy()[p].mean()),
            "share_theo_trypsin_bg": float(T.theo_trypsin.to_numpy()[c].mean()),
            "share_emp_obs_pos": float(T.emp_obs.to_numpy()[p].mean()),
            "share_emp_obs_bg": float(T.emp_obs.to_numpy()[c].mean()),
            "share_identified_in_arm_bg": float(T.detected.to_numpy()[c].mean()),
            "n_proteins_bg": int(T.protein.to_numpy()[c].size and len(set(T.protein.to_numpy()[c])))}


# ----------------------------------------------------------------------------------------------- statistics
def design_ci(T, role, feat, reps=BOOT_REPS):
    """Cys-Audit convention: clusters = proteins with rows in this comparison."""
    keep = role >= 0
    prot = T.protein.to_numpy()
    codes, labels = stats.cluster_index([prot[i] for i in np.flatnonzero(keep)])
    is_pos = role[keep] == 1
    flag = T[feat].to_numpy()[keep].astype(bool)
    counts = stats.two_by_two_counts(codes, len(labels), is_pos, flag)
    point, _, p, vals = stats.boot_log2_or(counts, reps, FEATURE_SEED[feat], 0.95)
    lo95, hi95 = stats.percentile_interval(vals, 0.95)
    lo975, hi975 = stats.percentile_interval(vals, 0.975)
    tot = counts.sum(axis=0)
    return {"estimate": point, "ci95_low": lo95, "ci95_high": hi95, "ci975_low": lo975, "ci975_high": hi975,
            "boot_p": p, "share_flagged_pos": tot[0] / (tot[0] + tot[1]),
            "share_flagged_bg": tot[2] / (tot[2] + tot[3]), "n_clusters": len(labels)}


def paired_boot(T, roles, reps=BOOT_REPS, seed=PAIRED_SEED):
    """One set of protein multiplicities for all designs. Returns {design: {feat: replicate array}} and points."""
    prot = T.protein.to_numpy()
    codes, labels = stats.cluster_index(list(prot))
    K = len(labels)
    names = list(roles)
    blocks = []
    for d in names:
        role = roles[d]
        for f in FEATURES:
            blocks.append(stats.two_by_two_counts(codes, K, role == 1, T[f].to_numpy().astype(bool), mask=role >= 0))
    counts = np.concatenate(blocks, axis=1)  # K x (D*F*4)
    tot = counts.sum(axis=0).reshape(len(names), len(FEATURES), 4)
    point = stats.log2_or(tot[..., 0], tot[..., 1], tot[..., 2], tot[..., 3])
    reps_out = []
    for w in stats.multiplicities(K, reps, seed):
        s = (w @ counts).reshape(w.shape[0], len(names), len(FEATURES), 4)
        reps_out.append(stats.log2_or(s[..., 0], s[..., 1], s[..., 2], s[..., 3]))
    R = np.concatenate(reps_out, axis=0)  # reps x D x F
    return names, point, R


def mh_stratified(T, pos_mask, score, n_strata, reps=BOOT_REPS, seed=MASTER_SEED + 13):
    """Mantel-Haenszel log2 OR over strata of the empirical score. Strata edges = quantiles of the positives'
    scores; background outside the positives' range is dropped (common support). Protein-clustered bootstrap."""
    pool = T.label.to_numpy() == 0
    qs = np.quantile(score[pos_mask], np.linspace(0, 1, n_strata + 1))
    lo, hi = qs[0], qs[-1]
    use = (pos_mask | (pool & (score >= lo) & (score <= hi)))
    strat = np.clip(np.searchsorted(qs[1:-1], score, side="right"), 0, n_strata - 1)
    prot = T.protein.to_numpy()
    codes, labels = stats.cluster_index([prot[i] for i in np.flatnonzero(use)])
    K = len(labels)
    out = {}
    for f in FEATURES:
        flag = T[f].to_numpy().astype(bool)[use]
        isp = pos_mask[use]
        st = strat[use]
        cnt = np.zeros((K, n_strata, 4), dtype=np.int64)
        for j, sel in enumerate((isp & flag, isp & ~flag, ~isp & flag, ~isp & ~flag)):
            np.add.at(cnt, (codes[sel], st[sel], j), 1)
        flat = cnt.reshape(K, n_strata * 4)
        point = stats.mh_log2_or(cnt.sum(axis=0))
        vals = []
        for w in stats.multiplicities(K, reps, seed + FEATURES.index(f)):
            s = (w @ flat).reshape(w.shape[0], n_strata, 4)
            vals.extend(stats.mh_log2_or(s[r]) for r in range(w.shape[0]))
        lo95, hi95 = stats.percentile_interval(vals, 0.95)
        out[f] = {"estimate": point, "ci95_low": lo95, "ci95_high": hi95}
    out["n_positive"] = int(pos_mask.sum())
    out["n_background"] = int((use & pool).sum())
    out["background_dropped_outside_support"] = int((pool & ~use).sum())
    return out


def main():
    T = pd.read_csv(f"{RESULTS}/cys_scores_trypsin_arm.csv", dtype={"protein": str})
    lab = T.label.to_numpy()
    det = T.detected.to_numpy()
    pool = lab == 0
    s_lr = logit(T.score_lr_full.to_numpy())
    s_gbm = logit(T.score_gbm_full.to_numpy())
    s_pep = logit(T.score_lr_peptide.to_numpy())
    theo = T.theo_trypsin.to_numpy().astype(float)
    empf = T.emp_obs.to_numpy().astype(float)
    pos_sets = {
        "all": lab == 1,
        "ha_strict": T.ha_specific.to_numpy() == 1,                              # CAM in HydP, seen without CAM in HydN
        "ha_lenient": (lab == 1) & (T.hydn_label.to_numpy() == 0),               # CAM in HydP, not CAM in HydN
    }
    sd_pooled = float(np.std(s_lr[(lab == 1) | pool], ddof=1))
    caliper = 0.2 * sd_pooled

    def restrict(pos_mask, bg_mask, symmetric_mask=None):
        role = np.full(len(T), -1, dtype=np.int64)
        pm = pos_mask if symmetric_mask is None else (pos_mask & symmetric_mask)
        role[bg_mask] = 0
        role[pm] = 1
        return role

    designs, infos = {}, {}
    for ps_name, pm in pos_sets.items():
        # positives outside the chosen positive set but carrying label 1 are excluded from everything
        tag = "" if ps_name == "all" else f"|{ps_name}"
        designs[f"proteome{tag}"] = restrict(pm, pool)
        designs[f"observed{tag}"] = restrict(pm, pool & (det == 1))
        designs[f"theoretical{tag}"] = restrict(pm, pool & (theo == 1))
        designs[f"empirical{tag}"] = restrict(pm, pool & (empf == 1))
        rng = np.random.default_rng(MATCH_SEED)
        designs[f"match_emp_lr{tag}"], infos[f"match_emp_lr{tag}"] = match_design(T, pm, s_lr, rng)
        rng = np.random.default_rng(MATCH_SEED)
        designs[f"match_theo{tag}"], infos[f"match_theo{tag}"] = match_design(T, pm, theo, rng)
        if ps_name == "all":
            designs["theoretical_symmetric"] = restrict(pm, pool & (theo == 1), theo == 1)
            designs["theoretical_trypsinP"] = restrict(pm, pool & (T.theo_trypsinP.to_numpy() == 1))
            designs["empirical_symmetric"] = restrict(pm, pool & (empf == 1), empf == 1)
            for nm, sc, kw in (("match_emp_gbm", s_gbm, {}), ("match_emp_peptide_only", s_pep, {}),
                               ("match_emp_lr_caliper", s_lr, {"caliper": caliper}),
                               ("match_emp_lr_within_protein", s_lr, {"within_protein": True}),
                               ("match_emp_obs_flag", empf, {})):
                rng = np.random.default_rng(MATCH_SEED)
                designs[nm], infos[nm] = match_design(T, pm, sc, rng, **kw)
            infos["match_emp_lr_caliper"]["caliper_logit"] = caliper
            # decomposition added after the first run (post hoc within this post hoc analysis): the same
            # backgrounds restricted to the proteins that carry at least one positive ("same-protein")
            same = T.protein.isin(set(T.protein[lab == 1])).to_numpy()
            designs["proteome_same_protein"] = restrict(pm, pool & same)
            designs["observed_same_protein"] = restrict(pm, pool & (det == 1) & same)
            designs["theoretical_same_protein"] = restrict(pm, pool & (theo == 1) & same)
            designs["empirical_same_protein"] = restrict(pm, pool & (empf == 1) & same)
            designs["observed_other_protein"] = restrict(pm, pool & (det == 1) & ~same)
            designs["empirical_other_protein"] = restrict(pm, pool & (empf == 1) & ~same)
            rng = np.random.default_rng(MATCH_SEED)
            designs["match_emp_lr_within_protein_caliper"], infos["match_emp_lr_within_protein_caliper"] = \
                match_design(T, pm, s_lr, rng, caliper=caliper, within_protein=True)

    # ---- per-design estimates with Cys-Audit-convention intervals
    rows = []
    for d, role in designs.items():
        r = {"design": d, "n_positive": int((role == 1).sum()), "n_background": int((role == 0).sum()),
             "n_proteins": int(len(set(T.protein.to_numpy()[role >= 0])))}
        for f in FEATURES:
            res = design_ci(T, role, f)
            for k, v in res.items():
                r[f"{f}_{k}"] = v
        b = balance(T, role, s_lr)
        r.update({f"bal_{k}": v for k, v in b.items()})
        rows.append(r)
        print(d, r["n_positive"], r["n_background"], "kr_dist", round(r["kr_dist_estimate"], 4),
              [round(r["kr_dist_ci95_low"], 4), round(r["kr_dist_ci95_high"], 4)], "kr_prox", round(r["kr_prox_estimate"], 4),
              "de_dist", round(r["de_dist_estimate"], 4), flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(f"{RESULTS}/background_log2or_by_design.csv", index=False)

    # ---- paired bootstrap: fraction removed and differences, within each positive set
    paired_rows = []
    for ps_name in pos_sets:
        tag = "" if ps_name == "all" else f"|{ps_name}"
        sub = {d: r for d, r in designs.items() if (tag and d.endswith(tag)) or (not tag and "|" not in d)}
        names, point, R = paired_boot(T, sub)
        i_prot = names.index(f"proteome{tag}")
        i_obs = names.index(f"observed{tag}")
        i_theo = names.index(f"theoretical{tag}")
        i_mth = names.index(f"match_theo{tag}")
        for di, d in enumerate(names):
            for fi, f in enumerate(FEATURES):
                lp, lb, lo_ = R[:, i_prot, fi], R[:, di, fi], R[:, i_obs, fi]
                frac = 1 - lb / lp
                diff_p = lp - lb
                diff_o = lb - lo_
                diff_t = lb - R[:, i_theo, fi]
                diff_mt = lb - R[:, i_mth, fi]
                row = {"positive_set": ps_name, "design": d, "feature": f, "estimate_paired_point": float(point[di, fi]),
                       "proteome_point": float(point[i_prot, fi]), "observed_point": float(point[i_obs, fi]),
                       "fraction_removed_point": float(1 - point[di, fi] / point[i_prot, fi]),
                       "reduction_vs_proteome_point": float(point[i_prot, fi] - point[di, fi]),
                       "excess_over_observed_point": float(point[di, fi] - point[i_obs, fi]),
                       "diff_vs_theoretical_point": float(point[di, fi] - point[i_theo, fi]),
                       "diff_vs_match_theo_primary_seed_point": float(point[di, fi] - point[i_mth, fi])}
                for nm, v in (("fraction_removed", frac), ("reduction_vs_proteome", diff_p),
                              ("excess_over_observed", diff_o), ("diff_vs_theoretical", diff_t),
                              ("diff_vs_match_theo_primary_seed", diff_mt)):
                    lo, hi = stats.percentile_interval(v, 0.95)
                    row[f"{nm}_ci95_low"], row[f"{nm}_ci95_high"] = lo, hi
                paired_rows.append(row)
    Pd = pd.DataFrame(paired_rows)
    Pd.to_csv(f"{RESULTS}/paired_fraction_removed.csv", index=False)

    # ---- Mantel-Haenszel over empirical-score strata (sensitivity)
    mh = {}
    for nst in (5, 10):
        mh[f"lr_full_{nst}_strata"] = mh_stratified(T, lab == 1, T.score_lr_full.to_numpy(), nst)
    dump_json(mh, f"{RESULTS}/mh_stratified.json")
    print("MH", json.dumps(mh, default=float)[:600], flush=True)

    # ---- matching variability over 200 matching seeds (point estimates only)
    var_rows = []
    for nm, sc, kw in (("match_emp_lr", s_lr, {}), ("match_theo", theo, {}), ("match_emp_gbm", s_gbm, {}),
                       ("match_emp_peptide_only", s_pep, {}), ("match_emp_obs_flag", empf, {}),
                       ("match_emp_lr_within_protein", s_lr, {"within_protein": True})):
        for k in range(N_MATCH_SEEDS):
            rng = np.random.default_rng(MATCH_SEED + 1000 + k)
            role, _ = match_design(T, lab == 1, sc, rng, **kw)
            keep = role >= 0
            r = {"design": nm, "match_seed": MATCH_SEED + 1000 + k}
            for f in FEATURES:
                flag = T[f].to_numpy().astype(bool)
                a = int(((role == 1) & flag).sum()); b = int(((role == 1) & ~flag).sum())
                c = int(((role == 0) & flag).sum()); d_ = int(((role == 0) & ~flag).sum())
                r[f] = float(stats.log2_or(a, b, c, d_))
            var_rows.append(r)
    V = pd.DataFrame(var_rows)
    V.to_csv(f"{RESULTS}/matching_seed_variability.csv", index=False)
    vs = V.groupby("design")[FEATURES].describe(percentiles=[0.025, 0.5, 0.975])
    print(vs.T.to_string(), flush=True)

    meta = {"caliper_logit": caliper, "sd_pooled_logit_lr_full": sd_pooled, "match_infos": infos,
            "paired_seed": PAIRED_SEED, "match_seed": MATCH_SEED, "feature_seeds": FEATURE_SEED,
            "n_match_seeds_variability": N_MATCH_SEEDS, "boot_reps": BOOT_REPS,
            "positive_set_sizes": {k: int(v.sum()) for k, v in pos_sets.items()}}
    dump_json(meta, f"{RESULTS}/backgrounds_meta.json")


if __name__ == "__main__":
    main()
