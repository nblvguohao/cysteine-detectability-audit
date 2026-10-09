"""Step 9 (POST HOC revision analysis, item D_empirical_background), written 2026-09-30 AFTER adversarial
verification round 1, to answer the verifier's problems. Nothing here was registered or pre-specified: every quantity
below was decided after all earlier results of this item, and the verifier's, had been seen.

A. Paired contrasts between backgrounds, distal and proximal K/R and D/E, with 95% and 97.5% (Bonferroni over the two
   bands, as in Fig. 2) percentile intervals and a two-sided bootstrap p, from
     - the primary paired seed of steps 3/4b (20260930 + 7; 5,000 replicates; reproduces the stored intervals),
     - 20 alternative seeds (20262930 + k, k = 0..19; 5,000 replicates each): Monte Carlo spread of the bounds,
     - one 50,000-replicate run (seed 20263930): Monte Carlo-precise bounds and p values.
   Restriction designs share one set of protein multiplicities (all 1,991 proteins). Matched designs use the
   matched-sample bootstrap of step 4b: replicate r uses matching (r mod 200) of the 200 matchings of step 3/4b
   (seeds 20260941 + 1000 + k), so the matched samples are identical to the stored ones.
B. A declared family of six distal K/R contrasts used in the text -- (b)-(e), (b')-(e), (c)-(e), (d)-(e), (c)-(b),
   (d)-(b') -- with Holm-adjusted bootstrap p values and Bonferroni-simultaneous 99.17% intervals (1 - 0.05/6).
C. The D/E comparison feature: its association with observation in the unenriched arm (all cysteines, sites only,
   background only), its shift between backgrounds next to the K/R shift, and the symmetric-restriction values.
D. Matching fragility: per-matching Cys-Audit status (97.5% interval, Cys-Audit convention, margin 0.5) over the 200
   matchings of (b'), (d), (d'); and a deterministic flag-standardised (b') ("bpstd"): the background re-weighted to
   the sites' theoretical-flag distribution, which is what 1:1 matching on a binary flag samples from. Its
   like-for-like comparator "empflagstd" re-weights to the sites' mix of observation in the unenriched arm (the
   noise-free counterpart of matching on that flag). Both are carried through the paired bootstrap of A (same protein
   multiplicities), so their contrasts with (a)-(e) are paired. A contrast between designs of different kinds (e.g.
   restriction (c) against bpstd) is written out but is not read on its own.
E. Construction asymmetry (partial check possible with the deposited tables): the observed background (e) restricted to
   cysteines covered by >= 1 phase-3 candidate, and the overlap of (c) and (e).
F. Reproduction deviations against the stored Fig. 2 audits (max |difference| of points and interval bounds).

Outputs (results/D_empirical_background/): r1_paired_contrasts_runs.csv, r1_paired_contrasts_summary.csv,
r1_holm_family.csv, r1_comparison_feature.csv, r1_matching_status.csv, r1_matching_status_summary.csv,
r1_standardised_flags.csv, r1_construction_check.csv, r1_summary.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

import sys

sys.dont_write_bytecode = True
import json

import numpy as np
import pandas as pd

from common import BOOT_REPS, MASTER_SEED, RESULTS, TOOL_SEED, dump_json
from cys_audit import stats  # noqa: E402
from cys_audit.status import status_from_interval  # noqa: E402
from matching import logit, match_design

FEATURES = ["kr_dist", "kr_prox", "de_dist", "de_prox"]
FEATURE_SEED = {"kr_prox": TOOL_SEED + 101, "kr_dist": TOOL_SEED + 102, "de_prox": TOOL_SEED + 103,
                "de_dist": TOOL_SEED + 104}
PAIRED_SEED = MASTER_SEED + 7
ALT_SEEDS = [MASTER_SEED + 2000 + k for k in range(20)]
BIG_SEED = MASTER_SEED + 3000
BIG_REPS = 50000
LINK_SEED = MASTER_SEED + 4000
STD_SEED = PAIRED_SEED
MATCH_SEED = MASTER_SEED + 11
N_MATCH = 200
FAMILY = ["b-e", "bp-e", "c-e", "d-e", "c-b", "d-bp"]
SIMUL_LEVEL = 1 - 0.05 / len(FAMILY)


def l2(x):
    x = np.asarray(x)
    return stats.log2_or(x[..., 0], x[..., 1], x[..., 2], x[..., 3])


def q(v, level):
    return stats.percentile_interval(v, level)


def main():
    T = pd.read_csv(f"{RESULTS}/cys_scores_trypsin_arm.csv", dtype={"protein": str})
    lab, det, prot = T.label.to_numpy(), T.detected.to_numpy(), T.protein.to_numpy()
    pool = lab == 0
    pos = lab == 1
    theo = T.theo_trypsin.to_numpy() == 1
    empf = T.emp_obs.to_numpy() == 1
    cand = T.cand_n.to_numpy() > 0
    same = T.protein.isin(set(T.protein[pos])).to_numpy()
    flags = {f: T[f].to_numpy().astype(bool) for f in FEATURES}
    codes, labels = stats.cluster_index(list(prot))
    K = len(labels)
    NF = len(FEATURES)

    def cc(role):
        return np.stack([stats.two_by_two_counts(codes, K, role == 1, flags[f], mask=role >= 0) for f in FEATURES],
                        axis=1)  # K x F x 4

    def restrict(bg_mask, pos_mask=None):
        role = np.full(len(T), -1, dtype=np.int64)
        role[bg_mask] = 0
        role[pos if pos_mask is None else pos_mask] = 1
        return role

    # ------------------------------------------------------------------ designs
    R = {"a": restrict(pool), "b": restrict(pool & theo), "c": restrict(pool & empf), "e": restrict(pool & (det == 1)),
         "e_cand": restrict(pool & (det == 1) & cand),
         "a_same": restrict(pool & same), "b_same": restrict(pool & theo & same),
         "c_same": restrict(pool & empf & same), "e_same": restrict(pool & (det == 1) & same)}
    rnames = list(R)
    Carr = np.stack([cc(R[d]) for d in rnames], axis=1)  # K x D x F x 4
    flatR = Carr.reshape(K, -1)
    pointR = {d: l2(Carr[:, i].sum(axis=0)) for i, d in enumerate(rnames)}  # F

    scores = {"bp": T.theo_trypsin.to_numpy().astype(float), "d": logit(T.score_lr_full.to_numpy()),
              "dprime": logit(T.score_gbm_full.to_numpy())}
    M, MP = {}, {}
    for nm, sc in scores.items():
        arr = np.zeros((N_MATCH, K, NF, 4), dtype=np.int64)
        for k in range(N_MATCH):
            role, _ = match_design(lab, prot, pos, sc, np.random.default_rng(MATCH_SEED + 1000 + k))
            arr[k] = cc(role)
        M[nm] = arr
        MP[nm] = l2(arr.sum(axis=1))  # N_MATCH x F
        print("matched", nm, "median kr_dist", round(float(np.median(MP[nm][:, 0])), 4), flush=True)
    pointM = {nm: np.median(MP[nm], axis=0) for nm in M}

    # deterministic flag-standardised designs: the proteome background re-weighted to the sites' mix of a binary flag,
    # i.e. the population that 1:1 matching on that flag samples from, without the sampling noise of matching.
    #   bpstd       theoretical flag (noise-free counterpart of (b'))
    #   empflagstd  observation in the unenriched arm (noise-free counterpart of matching on that flag); the
    #               like-for-like comparator of bpstd, added so that no cross-design contrast is read alone
    # Per-cluster columns: sites flagged, sites unflagged, then per stratum s in (0, 1): sites in s, background in s
    # flagged, background in s.
    def std_counts(strata):
        out = []
        for f in FEATURES:
            fl = flags[f]
            cols = [pos & fl, pos & ~fl]
            for s_ in (0, 1):
                cols += [pos & (strata == s_), pool & (strata == s_) & fl, pool & (strata == s_)]
            out.append(np.stack([np.bincount(codes[c], minlength=K) for c in cols], axis=1).astype(np.int64))
        return np.stack(out, axis=1)  # K x F x 8

    STD_DESIGNS = {"bpstd": std_counts(theo.astype(int)), "empflagstd": std_counts(empf.astype(int))}
    CNT_STD = STD_DESIGNS["bpstd"]
    flatS = {nm: c.reshape(K, -1) for nm, c in STD_DESIGNS.items()}

    def std_est(tot):
        a, b = tot[..., 0], tot[..., 1]
        n = a + b
        pstd = 0.0
        for j in range(2):
            ns, bf, bt = tot[..., 2 + 3 * j], tot[..., 3 + 3 * j], tot[..., 4 + 3 * j]
            pstd = pstd + (ns / n) * np.where(bt > 0, bf / np.maximum(bt, 1), 0.0)
        return stats.log2_or(a, b, n * pstd, n * (1 - pstd))

    pointStd = {nm: std_est(c.sum(axis=0)) for nm, c in STD_DESIGNS.items()}  # F each
    pointS = pointStd["bpstd"]

    def run_boot(seed, reps):
        outR, outM, outS = [], {nm: [] for nm in M}, {nm: [] for nm in STD_DESIGNS}
        r0 = 0
        for W in stats.multiplicities(K, reps, seed):
            m = W.shape[0]
            s = (W @ flatR).reshape(m, len(rnames), NF, 4)
            outR.append(l2(s))
            for nm in STD_DESIGNS:
                outS[nm].append(std_est((W @ flatS[nm]).reshape(m, NF, 8)))
            ks = np.arange(r0, r0 + m) % N_MATCH
            for nm, arr in M.items():
                outM[nm].append(l2(np.einsum("rk,rkfc->rfc", W, arr[ks])))
            r0 += m
        Rr = np.concatenate(outR)  # reps x D x F
        reps_by = {d: Rr[:, i, :] for i, d in enumerate(rnames)}
        reps_by.update({nm: np.concatenate(v) for nm, v in outM.items()})
        reps_by.update({nm: np.concatenate(v) for nm, v in outS.items()})
        return reps_by

    points = dict(pointR)
    points.update(pointM)
    points.update(pointStd)
    CONTRASTS = {"a-e": ("a", "e"), "b-e": ("b", "e"), "c-e": ("c", "e"), "c-b": ("c", "b"),
                 "e_cand-e": ("e_cand", "e"), "c-e_cand": ("c", "e_cand"),
                 "c_same-e_same": ("c_same", "e_same"), "b_same-e_same": ("b_same", "e_same"),
                 "c_same-b_same": ("c_same", "b_same"), "a_same-e_same": ("a_same", "e_same"),
                 "d-e": ("d", "e"), "bp-e": ("bp", "e"), "dprime-e": ("dprime", "e"), "d-bp": ("d", "bp"),
                 "dprime-bp": ("dprime", "bp"), "c-a": ("c", "a"), "e-a": ("e", "a"), "b-a": ("b", "a"),
                 "d-a": ("d", "a"), "bp-a": ("bp", "a"),
                 "bpstd-e": ("bpstd", "e"), "c-bpstd": ("c", "bpstd"), "d-bpstd": ("d", "bpstd"),
                 "bpstd-a": ("bpstd", "a"), "bpstd-b": ("bpstd", "b"),
                 "empflagstd-e": ("empflagstd", "e"), "empflagstd-bpstd": ("empflagstd", "bpstd"),
                 "empflagstd-c": ("empflagstd", "c")}
    DESIGNS_SINGLE = ["a", "b", "c", "e", "e_cand", "bp", "bpstd", "empflagstd", "d", "dprime", "c_same", "e_same",
                      "b_same", "a_same"]

    runs = [("primary_5000", PAIRED_SEED, BOOT_REPS)] + [(f"alt_{s}", s, BOOT_REPS) for s in ALT_SEEDS] + \
           [("big_50000", BIG_SEED, BIG_REPS)]
    rows = []
    kr_c_minus_b_big = de_c_minus_b_big = None
    big_reps_store = {}
    for run, seed, reps in runs:
        B = run_boot(seed, reps)
        for cname, (x, y) in CONTRASTS.items():
            for fi, f in enumerate(FEATURES):
                v = B[x][:, fi] - B[y][:, fi]
                lo95, hi95 = q(v, 0.95)
                lo975, hi975 = q(v, 0.975)
                los, his = q(v, SIMUL_LEVEL)
                rows.append({"run": run, "seed": seed, "reps": reps, "kind": "contrast", "name": cname,
                             "feature": f, "point": float(points[x][fi] - points[y][fi]),
                             "ci95_low": lo95, "ci95_high": hi95, "ci975_low": lo975, "ci975_high": hi975,
                             "ci_simul9917_low": los, "ci_simul9917_high": his,
                             "boot_p_two_sided": stats.bootstrap_p(v, 0.0), "share_le0": float(np.mean(v <= 0)),
                             "boot_median": float(np.median(v))})
        for d in DESIGNS_SINGLE:
            for fi, f in enumerate(FEATURES):
                v = B[d][:, fi]
                lo95, hi95 = q(v, 0.95)
                lo975, hi975 = q(v, 0.975)
                frac = 1 - v / B["a"][:, fi]
                flo, fhi = q(frac, 0.95)
                rows.append({"run": run, "seed": seed, "reps": reps, "kind": "design_paired", "name": d, "feature": f,
                             "point": float(points[d][fi]), "ci95_low": lo95, "ci95_high": hi95,
                             "ci975_low": lo975, "ci975_high": hi975, "boot_p_two_sided": stats.bootstrap_p(v, 0.0),
                             "share_le0": float(np.mean(v <= 0)), "boot_median": float(np.median(v)),
                             "fraction_removed_point": float(1 - points[d][fi] / points["a"][fi]),
                             "fraction_removed_ci95_low": flo, "fraction_removed_ci95_high": fhi})
        # difference in the theoretical-to-empirical shift between K/R and D/E (distal band)
        for (x, y) in (("c", "b"), ("d", "bp"), ("c", "a"), ("e", "a")):
            v = (B[x][:, 0] - B[y][:, 0]) - (B[x][:, 2] - B[y][:, 2])
            lo95, hi95 = q(v, 0.95)
            pt = float((points[x][0] - points[y][0]) - (points[x][2] - points[y][2]))
            rows.append({"run": run, "seed": seed, "reps": reps, "kind": "shift_difference_kr_minus_de_distal",
                         "name": f"({x}-{y})_kr_dist - ({x}-{y})_de_dist", "feature": "kr_dist-de_dist",
                         "point": pt, "ci95_low": lo95, "ci95_high": hi95,
                         "boot_p_two_sided": stats.bootstrap_p(v, 0.0), "share_le0": float(np.mean(v <= 0))})
        if run == "big_50000":
            big_reps_store = {k: v[:, 0].copy() for k, v in B.items()}
        print("run", run, "done", flush=True)
    RUNS = pd.DataFrame(rows)
    RUNS.to_csv(f"{RESULTS}/r1_paired_contrasts_runs.csv", index=False)

    # ---- summary across runs (the 21 runs of 5,000 replicates and the 50,000 run)
    summ = []
    for (kind, name, feat), g in RUNS.groupby(["kind", "name", "feature"], sort=False):
        pr = g[g.run == "primary_5000"].iloc[0]
        bg = g[g.run == "big_50000"].iloc[0]
        g5 = g[g.reps == BOOT_REPS]
        r = {"kind": kind, "name": name, "feature": feat, "point": pr.point,
             "primary_ci95_low": pr.ci95_low, "primary_ci95_high": pr.ci95_high,
             "primary_ci975_low": pr.get("ci975_low"), "primary_ci975_high": pr.get("ci975_high"),
             "primary_boot_p": pr.boot_p_two_sided,
             "n_runs_5000": int(len(g5)),
             "runs5000_ci95_low_min": g5.ci95_low.min(), "runs5000_ci95_low_max": g5.ci95_low.max(),
             "runs5000_ci95_high_min": g5.ci95_high.min(), "runs5000_ci95_high_max": g5.ci95_high.max(),
             "runs5000_ci975_low_min": g5.ci975_low.min() if "ci975_low" in g5 else np.nan,
             "runs5000_ci975_low_max": g5.ci975_low.max() if "ci975_low" in g5 else np.nan,
             "runs5000_boot_p_min": g5.boot_p_two_sided.min(), "runs5000_boot_p_max": g5.boot_p_two_sided.max(),
             "big50000_ci95_low": bg.ci95_low, "big50000_ci95_high": bg.ci95_high,
             "big50000_ci975_low": bg.get("ci975_low"), "big50000_ci975_high": bg.get("ci975_high"),
             "big50000_ci_simul9917_low": bg.get("ci_simul9917_low"),
             "big50000_ci_simul9917_high": bg.get("ci_simul9917_high"),
             "big50000_boot_p": bg.boot_p_two_sided}
        if kind == "design_paired":
            r.update({"fraction_removed_point": pr.fraction_removed_point,
                      "primary_fraction_removed_ci95_low": pr.fraction_removed_ci95_low,
                      "primary_fraction_removed_ci95_high": pr.fraction_removed_ci95_high,
                      "big50000_fraction_removed_ci95_low": bg.fraction_removed_ci95_low,
                      "big50000_fraction_removed_ci95_high": bg.fraction_removed_ci95_high})
        summ.append(r)
    SUM = pd.DataFrame(summ)
    SUM.to_csv(f"{RESULTS}/r1_paired_contrasts_summary.csv", index=False)

    # ---- Holm over the declared family (distal K/R), p values from the 50,000 run and from the primary run.
    # "primary6": the six contrasts of the designs in the brief; "extended9": plus the three contrasts of the
    # flag-standardised designs that the revised text also quotes (added after verification).
    families = {"primary6": FAMILY, "extended9": FAMILY + ["bpstd-e", "empflagstd-e", "empflagstd-bpstd"]}
    holm_rows = []
    for fam_name, fam_list in families.items():
        fam = SUM[(SUM.kind == "contrast") & (SUM.feature == "kr_dist") &
                  SUM.name.isin(fam_list)].set_index("name").loc[fam_list]
        for src in ("big50000_boot_p", "primary_boot_p"):
            p = fam[src].to_numpy()
            order = np.argsort(p, kind="stable")
            adj = np.empty(len(p))
            run_max = 0.0
            for rank, i in enumerate(order):
                run_max = max(run_max, min(1.0, (len(p) - rank) * p[i]))
                adj[i] = run_max
            for nm, pv, pa in zip(fam_list, p, adj):
                holm_rows.append({"family": fam_name, "family_size": len(fam_list), "p_source": src, "contrast": nm,
                                  "point": float(fam.loc[nm, "point"]), "boot_p": float(pv),
                                  "holm_adjusted_p": float(pa),
                                  "simul9917_low": float(fam.loc[nm, "big50000_ci_simul9917_low"]),
                                  "simul9917_high": float(fam.loc[nm, "big50000_ci_simul9917_high"])})
    H = pd.DataFrame(holm_rows)
    H.to_csv(f"{RESULTS}/r1_holm_family.csv", index=False)

    # ------------------------------------------------------------------ C. comparison feature D/E
    link_rows = []
    groups = {"all_cysteines": np.ones(len(T), dtype=bool), "sites_only": pos, "background_only": pool}
    for gname, gmask in groups.items():
        for fi, f in enumerate(FEATURES):
            sub = np.flatnonzero(gmask)
            c2, lab2 = stats.cluster_index([prot[i] for i in sub])
            cnt = stats.two_by_two_counts(c2, len(lab2), empf[sub], flags[f][sub])
            point, (lo, hi), p, _ = stats.boot_log2_or(cnt, BOOT_REPS, LINK_SEED + fi, 0.95)
            tot = cnt.sum(axis=0)
            link_rows.append({"group": gname, "feature": f,
                              "contrast": "observed in the unenriched arm vs not (log2 OR of carrying the feature)",
                              "n_observed": int(tot[0] + tot[1]), "n_not_observed": int(tot[2] + tot[3]),
                              "share_feature_observed": tot[0] / (tot[0] + tot[1]),
                              "share_feature_not_observed": tot[2] / (tot[2] + tot[3]),
                              "log2_or": point, "ci95_low": lo, "ci95_high": hi, "boot_p": p})
    LINK = pd.DataFrame(link_rows)
    sens = pd.read_csv(f"{RESULTS}/sensitivity_table.csv")
    sym = sens[(sens.positive_set == "all") & sens.design.isin(["empirical_symmetric", "theoretical_symmetric"])]
    for _, s in sym.iterrows():
        for f in FEATURES:
            LINK = pd.concat([LINK, pd.DataFrame([{
                "group": f"symmetric restriction: {s.design}", "feature": f,
                "contrast": "sites restricted like the background (stored, step 3)",
                "n_observed": int(s.n_positive), "n_not_observed": int(s.n_background),
                "log2_or": s[f"{f}_estimate"], "ci95_low": s[f"{f}_ci95_low"], "ci95_high": s[f"{f}_ci95_high"]}])],
                ignore_index=True)
    LINK.to_csv(f"{RESULTS}/r1_comparison_feature.csv", index=False)

    # ------------------------------------------------------------------ D. matching fragility
    st_rows = []
    for nm, arr in M.items():
        for k in range(N_MATCH):
            for fi, f in enumerate(("kr_dist", "kr_prox")):
                fidx = FEATURES.index(f)
                cnt_all = arr[k][:, fidx, :]
                rows_mask = cnt_all.sum(axis=1) > 0
                cnt = cnt_all[rows_mask]
                point, (lo, hi), p, _ = stats.boot_log2_or(cnt, BOOT_REPS, FEATURE_SEED[f], 0.975)
                status = status_from_interval(point, lo, hi, 0.0, 0.5)[0]
                st_rows.append({"design": nm, "match_seed": MATCH_SEED + 1000 + k, "feature": f, "estimate": point,
                                "ci975_low": lo, "ci975_high": hi, "status": status,
                                "n_clusters": int(rows_mask.sum())})
    ST = pd.DataFrame(st_rows)
    ST.to_csv(f"{RESULTS}/r1_matching_status.csv", index=False)
    sts = []
    for (nm, f), g in ST.groupby(["design", "feature"]):
        r = {"design": nm, "feature": f, "n_matchings": len(g), "estimate_median": g.estimate.median(),
             "estimate_sd": g.estimate.std(ddof=1), "estimate_p2.5": g.estimate.quantile(0.025),
             "estimate_p97.5": g.estimate.quantile(0.975), "estimate_min": g.estimate.min(),
             "estimate_max": g.estimate.max(), "ci975_low_min": g.ci975_low.min(), "ci975_low_max": g.ci975_low.max()}
        for s in ("FAIL", "WARNING", "UNDECIDABLE", "PASS"):
            r[f"share_{s}"] = float((g.status == s).mean())
        sts.append(r)
    STS = pd.DataFrame(sts)
    STS.to_csv(f"{RESULTS}/r1_matching_status_summary.csv", index=False)

    # self-check: the per-matching interval routine reproduces the stored primary-seed interval of step 3
    D3 = pd.read_csv(f"{RESULTS}/background_log2or_by_design.csv").set_index("design")
    chk = {}
    for nm, sc, dname in (("bp", scores["bp"], "match_theo"), ("d", scores["d"], "match_emp_lr")):
        role, _ = match_design(lab, prot, pos, sc, np.random.default_rng(MATCH_SEED))
        cnt_all = cc(role)[:, 0, :]
        cnt = cnt_all[cnt_all.sum(axis=1) > 0]
        point, (lo, hi), _, _ = stats.boot_log2_or(cnt, BOOT_REPS, FEATURE_SEED["kr_dist"], 0.975)
        chk[dname] = {"recomputed": [point, lo, hi],
                      "stored": [float(D3.loc[dname, "kr_dist_estimate"]), float(D3.loc[dname, "kr_dist_ci975_low"]),
                                 float(D3.loc[dname, "kr_dist_ci975_high"])]}

    # deterministic flag-standardised designs (theoretical flag; observation in the unenriched arm), per-design
    # intervals from the primary paired seed (same multiplicities as the primary run of A)
    V3 = pd.read_csv(f"{RESULTS}/matching_seed_variability.csv")
    std_rows = []
    labels_std = {"bpstd": ("background re-weighted to the sites' theoretical-flag mix (noise-free counterpart of b')",
                            "match_theo", "b"),
                  "empflagstd": ("background re-weighted to the sites' mix of observation in the unenriched arm "
                                 "(noise-free counterpart of matching on that flag)", "match_emp_obs_flag", "c")}
    for nm, cnt in STD_DESIGNS.items():
        flat_nm = cnt.reshape(K, -1)
        vals = np.concatenate([std_est((W @ flat_nm).reshape(W.shape[0], NF, 8))
                               for W in stats.multiplicities(K, BOOT_REPS, STD_SEED)])
        desc, matched_name, restr = labels_std[nm]
        for fi, f in enumerate(FEATURES):
            pt = float(pointStd[nm][fi])
            lo95, hi95 = q(vals[:, fi], 0.95)
            lo975, hi975 = q(vals[:, fi], 0.975)
            std_rows.append({"design": nm, "description": desc, "feature": f, "estimate": pt, "ci95_low": lo95,
                             "ci95_high": hi95, "ci975_low": lo975, "ci975_high": hi975,
                             "status_975": status_from_interval(pt, lo975, hi975, 0.0, 0.5)[0],
                             "matched_counterpart": matched_name,
                             "matched_median_200": float(V3[V3.design == matched_name][f].median()),
                             "restriction_counterpart": restr, "restriction_estimate": float(pointR[restr][fi])})
    STD = pd.DataFrame(std_rows)
    STD.to_csv(f"{RESULTS}/r1_standardised_flags.csv", index=False)

    # ------------------------------------------------------------------ E. construction check
    e_mask = pool & (det == 1)
    c_mask = pool & empf
    cons = []
    for nm, m in (("e: observed in the ABE arm", e_mask), ("e with no phase-3 candidate (cand_n = 0)", e_mask & ~cand),
                  ("e covered by >= 1 phase-3 candidate", e_mask & cand), ("c: observed in the unenriched arm", c_mask),
                  ("c and e", c_mask & e_mask), ("c only", c_mask & ~e_mask), ("e only", e_mask & ~c_mask),
                  ("sites", pos), ("sites observed in the unenriched arm", pos & empf),
                  ("sites not observed in the unenriched arm", pos & ~empf)):
        shares = {f"share_{f}": (float(flags[f][m].mean()) if m.any() else float("nan")) for f in FEATURES}
        cons.append({"set": nm, "n": int(m.sum()), **shares})
    CONS = pd.DataFrame(cons)
    CONS.to_csv(f"{RESULTS}/r1_construction_check.csv", index=False)

    # ------------------------------------------------------------------ F. reproduction deviations
    cli = json.load(open(f"{RESULTS}/cli_reproduction.json", encoding="utf-8"))
    dev_pt, dev_b = [], []
    for bgname in ("proteome", "observed"):
        for band in ("proximal_1_3", "distal_6_12"):
            st = cli["runs"][bgname][band]["stored"]
            rp = cli["runs"][bgname][band]["reproduced_fasta_2026_03"]
            dev_pt.append(abs(st["estimate"] - rp["estimate"]))
            dev_b += [abs(st["ci"][0] - rp["ci"][0]), abs(st["ci"][1] - rp["ci"][1])]

    def sel(kind, name, feat="kr_dist"):
        return SUM[(SUM.kind == kind) & (SUM.name == name) & (SUM.feature == feat)].iloc[0].to_dict()

    summary = {
        "label": "POST HOC revision analysis after verification round 1 (2026-09-30); not registered",
        "reproduction_max_abs_dev_point": max(dev_pt), "reproduction_max_abs_dev_bound": max(dev_b),
        "primary_interval_self_check": chk,
        "key_contrasts_kr_dist": {n: sel("contrast", n) for n in
                                  ("b-e", "bp-e", "bpstd-e", "c-e", "d-e", "dprime-e", "c-b", "d-bp", "c-bpstd",
                                   "d-bpstd", "bpstd-b", "empflagstd-e", "empflagstd-bpstd", "empflagstd-c",
                                   "c_same-e_same", "b_same-e_same", "c_same-b_same", "c-e_cand", "e_cand-e")},
        "key_contrasts_de_dist": {n: sel("contrast", n, "de_dist") for n in ("c-b", "d-bp", "c-a", "e-a", "c-e")},
        "shift_difference": RUNS[(RUNS.kind == "shift_difference_kr_minus_de_distal") &
                                 (RUNS.run.isin(["primary_5000", "big_50000"]))].to_dict(orient="records"),
        "holm": H.to_dict(orient="records"),
        "matching_status": STS.to_dict(orient="records"),
        "standardised_flags": STD.to_dict(orient="records"),
        "construction": CONS.to_dict(orient="records"),
        "seeds": {"primary_paired": PAIRED_SEED, "alternative_paired": ALT_SEEDS, "big_run": BIG_SEED,
                  "big_run_reps": BIG_REPS, "detection_link_bootstrap": [LINK_SEED + i for i in range(NF)],
                  "per_matching_status": FEATURE_SEED, "matching_seed_family": "20260941 + 1000 + k, k = 0..199",
                  "standardised_bprime": STD_SEED},
        "family_for_holm": FAMILY, "simultaneous_level": SIMUL_LEVEL,
    }
    dump_json(summary, f"{RESULTS}/r1_summary.json")

    with pd.option_context("display.width", 250, "display.max_columns", 40):
        cols = ["kind", "name", "feature", "point", "primary_ci95_low", "primary_ci95_high", "primary_ci975_low",
                "primary_ci975_high", "primary_boot_p", "runs5000_ci95_low_min", "runs5000_ci95_low_max",
                "big50000_ci95_low", "big50000_ci95_high", "big50000_ci975_low", "big50000_ci975_high",
                "big50000_boot_p"]
        print(SUM[SUM.feature.isin(["kr_dist", "de_dist"])][cols].round(4).to_string())
        print(H.round(4).to_string())
        print(LINK.round(4).to_string())
        print(STS.round(4).to_string())
        print(STD.round(4).to_string())
        print(CONS.round(4).to_string())
        print(json.dumps(chk, indent=1))
        print("reproduction max dev point", max(dev_pt), "bound", max(dev_b))


if __name__ == "__main__":
    main()
