"""Step 10 (POST HOC revision analysis, item D_empirical_background), written 2026-09-30 AFTER adversarial
verification round 2, to answer its two major problems. Nothing here was registered or pre-specified: every quantity
below was chosen after all earlier results of this item, and the round-2 verifier's, had been seen.

The verifier's scripts (verify_r2_*.py) are NOT imported or run. The designs they proposed are rebuilt here with this
item's cysteine table, scores, matching code and seed conventions, and the verifier's stored values are only read
for the comparison in H.

A. Empirical information used WITH the theoretical rule (combined designs):
   bc        background restricted to cysteines that are theoretically detectable AND observed in the unenriched
             arm (asymmetric: all sites, as in Artifact 1); bc_sym restricts the sites the same way.
   cstar_x   whole background standardised to the sites' distribution over theo x emp_obs (4 strata).
   dtheo     1:1 nearest-neighbour matching on logit(empirical score) WITHIN strata of the theoretical flag (exact on
             the flag), greedy, without replacement, random order, random ties; 200 matchings, seeds 20265930 + k.
B. Candidate-aware empirical standardisation cstar_cand: strata {in no phase-3 Trypsin/P candidate; in a candidate
   but not observed; observed}. This replaces cstar (step 9, "empflagstd") as the like-for-like counterpart of the
   theoretical standardisation bstar: cstar pooled the 1,948 background cysteines that lie in no candidate, which the
   unenriched arm's search could not report, with unobserved but observable cysteines. cstar_x3 = theo x emp x
   candidate is a sensitivity.
C. Paired contrasts from one set of protein multiplicities over all 1,991 proteins: primary run (5,000, seed
   20260930 + 7) and a 50,000-replicate run (seed 20260930 + 3000). These are the multiplicities of step 9, so the
   stored step-9 replicates of (b), (c), (e), (b*), (c*), (b') and (d) are reproduced exactly (self-check).
   Matched designs use the matched-sample bootstrap mixed over their 200 matchings (steps 4b and 9).
D. K/R-specific part of each step: (X - Y)_KR - (X - Y)_DE in the distal band, and the proximal K/R band.
E. Decomposition of the (c) residual: counts and shares by theo x emp for sites, background and (e); why the
   background cysteines observed in the unenriched arm fail the theoretical rule (Trypsin/P rule, identified peptide
   length).
F. Search design of the global (unenriched) arm and of the ABE arms, read from the phase-1 probe record of the
   deposit's mqpar_Global.xml and mqpar_ABE.xml (repo/results/phase1_probe_findings_2026-09-22.json; read-only).
G. A Holm family over the nine distal K/R contrasts quoted in the round-2 answer (defined after verification round 2),
   per-design and per-matching Cys-Audit status, balance of dtheo, same-protein and HA-set sensitivities of bc.
H. Comparison with the round-2 verifier's stored values (verify_r2/*.csv, *.json).

Outputs (results/D_empirical_background/): r2_designs.csv, r2_paired_contrasts.csv, r2_holm_round2.csv,
r2_kr_specific_part.csv, r2_standardisation_strata.csv, r2_c_residual_decomposition.csv, r2_c_not_theo_cysteines.csv,
r2_dtheo_matchings.csv, r2_global_arm_search_design.json, r2_verifier_comparison.csv, r2_summary.json
Run after step 9 and before step 07; step 10b (10b_proximal_offsets.py) is a descriptive follow-up.
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

from common import BOOT_REPS, MASTER_SEED, REPO, RESULTS, TOOL_SEED, dump_json, sha256_file
from cys_audit import stats  # noqa: E402
from cys_audit.status import status_from_interval  # noqa: E402
from matching import logit, match_design, nn_match

FEATURES = ["kr_dist", "kr_prox", "de_dist", "de_prox"]
FI = {f: i for i, f in enumerate(FEATURES)}
FEATURE_SEED = {"kr_prox": TOOL_SEED + 101, "kr_dist": TOOL_SEED + 102, "de_prox": TOOL_SEED + 103,
                "de_dist": TOOL_SEED + 104}
PAIRED_SEED = MASTER_SEED + 7
BIG_SEED = MASTER_SEED + 3000
BIG_REPS = 50000
MATCH_SEED = MASTER_SEED + 11          # stored (b') and (d): seeds MATCH_SEED + 1000 + k
DTHEO_SEED = MASTER_SEED + 5000        # new (d_theo): seeds DTHEO_SEED + k
N_MATCH = 200
PROBE_JSON = f"{REPO}/results/phase1_probe_findings_2026-09-22.json"
VR2 = f"{RESULTS}/verify_r2"
ROUND2_FAMILY = ["c-b", "cstar_cand-bstar", "d-bp", "bc-b", "bc-e", "cstar_x-bstar", "cstar_x-e", "dtheo-bp",
                 "dtheo-e"]


def l2(x):
    x = np.asarray(x)
    return stats.log2_or(x[..., 0], x[..., 1], x[..., 2], x[..., 3])


def q(v, level):
    return stats.percentile_interval(v, level)


def main():
    T = pd.read_csv(f"{RESULTS}/cys_scores_trypsin_arm.csv", dtype={"protein": str})
    n_rows = len(T)
    lab, det, prot = T.label.to_numpy(), T.detected.to_numpy() == 1, T.protein.to_numpy()
    pool, pos = lab == 0, lab == 1
    theo = T.theo_trypsin.to_numpy() == 1
    theoP = T.theo_trypsinP.to_numpy() == 1
    empf = T.emp_obs.to_numpy() == 1
    cand = T.cand_n.to_numpy() > 0
    same = T.protein.isin(set(T.protein[pos])).to_numpy()
    ha_strict = T.ha_specific.to_numpy() == 1
    ha_len = pos & (T.hydn_label.to_numpy() == 0)
    flags = {f: T[f].to_numpy().astype(bool) for f in FEATURES}
    codes, labels = stats.cluster_index(list(prot))
    K = len(labels)
    NF = len(FEATURES)
    lg = logit(T.score_lr_full.to_numpy())
    checks = {}

    def cc(role):
        return np.stack([stats.two_by_two_counts(codes, K, role == 1, flags[f], mask=role >= 0) for f in FEATURES],
                        axis=1)  # K x F x 4

    def restrict(bg_mask, pos_mask=None):
        role = np.full(n_rows, -1, dtype=np.int64)
        role[bg_mask] = 0
        role[pos if pos_mask is None else pos_mask] = 1
        return role

    def design_ci(role, feat, level_pair=(0.95, 0.975)):
        """Cys-Audit convention: clusters = proteins with rows in this comparison; the tool's feature seed."""
        keep = role >= 0
        c2, lab2 = stats.cluster_index([prot[i] for i in np.flatnonzero(keep)])
        cnt = stats.two_by_two_counts(c2, len(lab2), role[keep] == 1, flags[feat][keep])
        point, _, p, vals = stats.boot_log2_or(cnt, BOOT_REPS, FEATURE_SEED[feat], 0.95)
        lo95, hi95 = q(vals, level_pair[0])
        lo975, hi975 = q(vals, level_pair[1])
        tot = cnt.sum(axis=0)
        return {"estimate": point, "ci95_low": lo95, "ci95_high": hi95, "ci975_low": lo975, "ci975_high": hi975,
                "boot_p": p, "n_positive": int(tot[0] + tot[1]), "n_background": int(tot[2] + tot[3]),
                "share_flagged_pos": tot[0] / (tot[0] + tot[1]), "share_flagged_bg": tot[2] / (tot[2] + tot[3]),
                "n_clusters": len(lab2)}

    # ------------------------------------------------------------------ A. restriction designs (paired set)
    R = {"a": restrict(pool), "b": restrict(pool & theo), "c": restrict(pool & empf), "e": restrict(pool & det),
         "bc": restrict(pool & theo & empf), "bc_sym": restrict(pool & theo & empf, pos & theo & empf),
         "c_not_theo": restrict(pool & empf & ~theo), "e_theo": restrict(pool & det & theo),
         "a_same": restrict(pool & same), "b_same": restrict(pool & theo & same),
         "c_same": restrict(pool & empf & same), "e_same": restrict(pool & det & same),
         "bc_same": restrict(pool & theo & empf & same)}
    rnames = list(R)
    Carr = np.stack([cc(R[d]) for d in rnames], axis=1)  # K x D x F x 4
    flatR = Carr.reshape(K, -1)
    pointR = {d: l2(Carr[:, i].sum(axis=0)) for i, d in enumerate(rnames)}

    # ------------------------------------------------------------------ B. standardised designs
    STRATA = {"bstar": theo.astype(int), "cstar": empf.astype(int),
              "cstar_cand": np.where(~cand, 2, empf.astype(int)),
              "cstar_x": theo.astype(int) * 2 + empf.astype(int),
              "cstar_x3": theo.astype(int) * 4 + empf.astype(int) * 2 + cand.astype(int)}
    LEVELS = {nm: sorted(set(st.tolist())) for nm, st in STRATA.items()}

    def std_counts(st, levels):
        out = []
        for f in FEATURES:
            fl = flags[f]
            cols = [pos & fl, pos & ~fl]
            for L in levels:
                cols += [pos & (st == L), pool & (st == L) & fl, pool & (st == L)]
            out.append(np.stack([np.bincount(codes[c], minlength=K) for c in cols], axis=1).astype(np.int64))
        return np.stack(out, axis=1)  # K x F x (2 + 3S)

    STD = {nm: std_counts(STRATA[nm], LEVELS[nm]) for nm in STRATA}
    flatS = {nm: c.reshape(K, -1) for nm, c in STD.items()}
    empty_stratum_hits = {nm: 0 for nm in STD}

    def std_est(tot, nm):
        a, b = tot[..., 0], tot[..., 1]
        n = a + b
        pstd = 0.0
        for j in range(len(LEVELS[nm])):
            ns, bf, bt = tot[..., 2 + 3 * j], tot[..., 3 + 3 * j], tot[..., 4 + 3 * j]
            bad = (ns > 0) & (bt == 0)
            if np.any(bad):
                empty_stratum_hits[nm] += int(np.sum(bad))
            pstd = pstd + (ns / n) * np.where(bt > 0, bf / np.maximum(bt, 1), 0.0)
        return stats.log2_or(a, b, n * pstd, n * (1 - pstd))

    pointStd = {nm: std_est(c.sum(axis=0), nm) for nm, c in STD.items()}

    # stratum table (all four standardised constructions): sites, background, background K/R-distal share
    strata_rows = []
    for nm, st in STRATA.items():
        for L in LEVELS[nm]:
            m = st == L
            strata_rows.append({"design": nm, "stratum": int(L), "n_sites": int((pos & m).sum()),
                                "site_weight": float((pos & m).sum() / pos.sum()), "n_background": int((pool & m).sum()),
                                "bg_share_kr_dist": float(flags["kr_dist"][pool & m].mean()) if (pool & m).any() else np.nan,
                                "site_share_kr_dist": float(flags["kr_dist"][pos & m].mean()) if (pos & m).any() else np.nan})
    STRATA_TAB = pd.DataFrame(strata_rows)

    # ------------------------------------------------------------------ A3. matched designs
    def match_within(strata, rng):
        role = np.full(n_rows, -1, dtype=np.int64)
        for s_ in (1, 0):
            pi = np.flatnonzero(pos & (strata == s_))
            ci = np.flatnonzero(pool & (strata == s_))
            m = nn_match(lg[pi], lg[ci], rng)
            ok = m >= 0
            role[pi[ok]] = 1
            role[ci[m[ok]]] = 0
        return role

    def smd(x, y):
        return float((x.mean() - y.mean()) / np.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2))

    M, MP, dtheo_rows = {}, {}, []
    for nm in ("bp", "d", "dtheo"):
        arr = np.zeros((N_MATCH, K, NF, 4), dtype=np.int64)
        for k in range(N_MATCH):
            if nm == "bp":
                role, _ = match_design(lab, prot, pos, T.theo_trypsin.to_numpy().astype(float),
                                       np.random.default_rng(MATCH_SEED + 1000 + k))
            elif nm == "d":
                role, _ = match_design(lab, prot, pos, lg, np.random.default_rng(MATCH_SEED + 1000 + k))
            else:
                role = match_within(theo.astype(int), np.random.default_rng(DTHEO_SEED + k))
                ctrl = role == 0
                assert (role == 1).sum() == pos.sum() == ctrl.sum()
                auc = stats.auc(lg[(role >= 0)], (role == 1)[role >= 0])
                dtheo_rows.append({"match_seed": DTHEO_SEED + k, "n_pairs": int(ctrl.sum()),
                                   "smd_logit_score": smd(lg[role == 1], lg[ctrl]), "auc_score_sites_vs_controls": auc,
                                   "share_theo_controls": float(theo[ctrl].mean()),
                                   "share_emp_obs_controls": float(empf[ctrl].mean()),
                                   "n_controls_in_no_candidate": int((ctrl & ~cand).sum()),
                                   "share_controls_same_protein_as_a_site": float(same[ctrl].mean())})
            arr[k] = cc(role)
        M[nm] = arr
        MP[nm] = l2(arr.sum(axis=1))  # N_MATCH x F
        print("matched", nm, "median kr_dist", round(float(np.median(MP[nm][:, 0])), 4), flush=True)
    pointM = {nm: np.median(MP[nm], axis=0) for nm in M}
    for i, r in enumerate(dtheo_rows):
        r.update({f"estimate_{f}": float(MP["dtheo"][i, FI[f]]) for f in FEATURES})

    # ------------------------------------------------------------------ C. paired bootstrap
    def run_boot(seed, reps):
        outR, outM, outS = [], {nm: [] for nm in M}, {nm: [] for nm in STD}
        r0 = 0
        for Wm in stats.multiplicities(K, reps, seed):
            m = Wm.shape[0]
            outR.append(l2((Wm @ flatR).reshape(m, len(rnames), NF, 4)))
            for nm in STD:
                outS[nm].append(std_est((Wm @ flatS[nm]).reshape(m, NF, -1), nm))
            ks = np.arange(r0, r0 + m) % N_MATCH
            for nm, arr in M.items():
                outM[nm].append(l2(np.einsum("rk,rkfc->rfc", Wm, arr[ks])))
            r0 += m
        Rr = np.concatenate(outR)
        B = {d: Rr[:, i, :] for i, d in enumerate(rnames)}
        B.update({nm: np.concatenate(v) for nm, v in outM.items()})
        B.update({nm: np.concatenate(v) for nm, v in outS.items()})
        return B

    points = dict(pointR)
    points.update(pointM)
    points.update(pointStd)
    CONTRASTS = {
        # reproduced from step 9 (self-check)
        "b-e": ("b", "e"), "c-b": ("c", "b"), "c-e": ("c", "e"), "bp-e": ("bp", "e"), "d-e": ("d", "e"),
        "d-bp": ("d", "bp"), "bstar-e": ("bstar", "e"), "cstar-bstar": ("cstar", "bstar"), "cstar-e": ("cstar", "e"),
        # empirical INSTEAD of the rule, candidate-aware standardisation
        "cstar_cand-bstar": ("cstar_cand", "bstar"), "cstar_cand-e": ("cstar_cand", "e"),
        "cstar_cand-cstar": ("cstar_cand", "cstar"),
        # empirical WITH the rule
        "bc-b": ("bc", "b"), "bc-c": ("bc", "c"), "bc-e": ("bc", "e"), "bc-a": ("bc", "a"),
        "bc_sym-e": ("bc_sym", "e"), "bc-e_theo": ("bc", "e_theo"), "e_theo-e": ("e_theo", "e"),
        "cstar_x-bstar": ("cstar_x", "bstar"), "cstar_x-e": ("cstar_x", "e"), "cstar_x-cstar_cand": ("cstar_x", "cstar_cand"),
        "cstar_x3-bstar": ("cstar_x3", "bstar"), "cstar_x3-e": ("cstar_x3", "e"),
        "dtheo-bp": ("dtheo", "bp"), "dtheo-e": ("dtheo", "e"), "dtheo-d": ("dtheo", "d"), "dtheo-b": ("dtheo", "b"),
        "dtheo-bstar": ("dtheo", "bstar"),
        # within the site-carrying proteins
        "bc_same-e_same": ("bc_same", "e_same"), "bc_same-b_same": ("bc_same", "b_same"),
        "b_same-e_same": ("b_same", "e_same"), "c_same-e_same": ("c_same", "e_same"),
        "c_same-b_same": ("c_same", "b_same")}
    SHIFTS = [("bc", "b"), ("cstar_x", "bstar"), ("dtheo", "bp"), ("c", "b"), ("cstar_cand", "bstar"), ("bc", "e"),
              ("cstar_x", "e"), ("dtheo", "e")]
    DESIGNS_PAIRED = ["a", "b", "c", "e", "bc", "bc_sym", "e_theo", "c_not_theo", "bstar", "cstar", "cstar_cand",
                      "cstar_x", "cstar_x3", "bp", "d", "dtheo", "a_same", "b_same", "c_same", "e_same", "bc_same"]

    rows, shift_rows = [], []
    store = {}
    for run, seed, reps in (("primary_5000", PAIRED_SEED, BOOT_REPS), ("big_50000", BIG_SEED, BIG_REPS)):
        B = run_boot(seed, reps)
        store[run] = B
        for cname, (x, y) in CONTRASTS.items():
            for f in FEATURES:
                v = B[x][:, FI[f]] - B[y][:, FI[f]]
                lo95, hi95 = q(v, 0.95)
                lo975, hi975 = q(v, 0.975)
                rows.append({"run": run, "seed": seed, "reps": reps, "kind": "contrast", "name": cname, "feature": f,
                             "point": float(points[x][FI[f]] - points[y][FI[f]]), "ci95_low": lo95, "ci95_high": hi95,
                             "ci975_low": lo975, "ci975_high": hi975, "boot_p": stats.bootstrap_p(v, 0.0)})
        for d in DESIGNS_PAIRED:
            for f in FEATURES:
                v = B[d][:, FI[f]]
                lo95, hi95 = q(v, 0.95)
                lo975, hi975 = q(v, 0.975)
                frac = 1 - v / B["a"][:, FI[f]]
                flo, fhi = q(frac, 0.95)
                rows.append({"run": run, "seed": seed, "reps": reps, "kind": "design_paired", "name": d, "feature": f,
                             "point": float(points[d][FI[f]]), "ci95_low": lo95, "ci95_high": hi95,
                             "ci975_low": lo975, "ci975_high": hi975, "boot_p": stats.bootstrap_p(v, 0.0),
                             "fraction_removed_point": float(1 - points[d][FI[f]] / points["a"][FI[f]]),
                             "fraction_removed_ci95_low": flo, "fraction_removed_ci95_high": fhi})
        for (x, y) in SHIFTS:
            v = (B[x][:, 0] - B[y][:, 0]) - (B[x][:, 2] - B[y][:, 2])
            lo95, hi95 = q(v, 0.95)
            lo975, hi975 = q(v, 0.975)
            shift_rows.append({"run": run, "seed": seed, "reps": reps, "step": f"{y} -> {x}",
                               "kr_dist_shift": float(points[x][0] - points[y][0]),
                               "de_dist_shift": float(points[x][2] - points[y][2]),
                               "kr_prox_shift": float(points[x][1] - points[y][1]),
                               "kr_minus_de_point": float((points[x][0] - points[y][0]) - (points[x][2] - points[y][2])),
                               "kr_minus_de_ci95_low": lo95, "kr_minus_de_ci95_high": hi95,
                               "kr_minus_de_ci975_low": lo975, "kr_minus_de_ci975_high": hi975,
                               "kr_minus_de_boot_p": stats.bootstrap_p(v, 0.0)})
        print("run", run, "done", flush=True)
    PC = pd.DataFrame(rows)
    PC.to_csv(f"{RESULTS}/r2_paired_contrasts.csv", index=False)
    SH = pd.DataFrame(shift_rows)
    SH.to_csv(f"{RESULTS}/r2_kr_specific_part.csv", index=False)
    checks["std_replicates_with_sites_but_no_background_in_a_stratum"] = dict(empty_stratum_hits)

    def pc(name, feat="kr_dist", run="big_50000", kind="contrast"):
        return PC[(PC.run == run) & (PC.kind == kind) & (PC.name == name) & (PC.feature == feat)].iloc[0]

    # ---- self-check: step-9 stored contrasts (same multiplicities and matchings) are reproduced exactly
    S9 = pd.read_csv(f"{RESULTS}/r1_paired_contrasts_summary.csv")
    map9 = {"b-e": "b-e", "c-b": "c-b", "c-e": "c-e", "bp-e": "bp-e", "d-e": "d-e", "d-bp": "d-bp",
            "bstar-e": "bpstd-e", "cstar-bstar": "empflagstd-bpstd", "cstar-e": "empflagstd-e"}
    sc = {}
    for mine, theirs in map9.items():
        for f in FEATURES:
            s9 = S9[(S9.kind == "contrast") & (S9.name == theirs) & (S9.feature == f)].iloc[0]
            r_b = pc(mine, f)
            r_p = pc(mine, f, run="primary_5000")
            dev = max(abs(r_b.point - s9.point), abs(r_b.ci95_low - s9.big50000_ci95_low),
                      abs(r_b.ci95_high - s9.big50000_ci95_high), abs(r_b.ci975_low - s9.big50000_ci975_low),
                      abs(r_b.ci975_high - s9.big50000_ci975_high), abs(r_b.boot_p - s9.big50000_boot_p),
                      abs(r_p.ci95_low - s9.primary_ci95_low), abs(r_p.ci95_high - s9.primary_ci95_high))
            sc[f"{mine}|{f}"] = float(dev)
    checks["max_abs_dev_vs_step9_stored_contrasts"] = float(max(sc.values()))
    checks["dev_vs_step9_by_contrast"] = sc
    assert checks["max_abs_dev_vs_step9_stored_contrasts"] < 1e-12, sc

    # ------------------------------------------------------------------ per-design table (Cys-Audit convention)
    design_rows = []
    desc = {"a": "(a) proteome: all other cysteines of identified proteins",
            "b": "(b) theoretically detectable (restriction)",
            "c": "(c) observed in the unenriched arm (restriction)",
            "e": "(e) identified in the ABE arm without a site (observed)",
            "bc": "(bc) theoretically detectable AND observed in the unenriched arm (restriction; all sites)",
            "bc_sym": "(bc) symmetric: sites restricted the same way",
            "c_not_theo": "background observed in the unenriched arm but not theoretically detectable (the 290)",
            "e_theo": "(e) restricted to theoretically detectable cysteines",
            "a_same": "proteome, site-carrying proteins", "b_same": "(b), site-carrying proteins",
            "c_same": "(c), site-carrying proteins", "e_same": "(e), site-carrying proteins",
            "bc_same": "(bc), site-carrying proteins"}
    for d in rnames:
        for f in FEATURES:
            r = design_ci(R[d], f)
            r.update({"design": d, "description": desc[d], "feature": f, "interval_type": "Cys-Audit convention",
                      "status_975_margin05": status_from_interval(r["estimate"], r["ci975_low"], r["ci975_high"], 0.0, 0.5)[0]})
            design_rows.append(r)
    for ps_name, ps_mask in (("ha_strict", ha_strict), ("ha_lenient", ha_len)):
        for d, bm in (("a", pool), ("b", pool & theo), ("c", pool & empf), ("bc", pool & theo & empf), ("e", pool & det)):
            for f in ("kr_dist", "kr_prox", "de_dist"):
                r = design_ci(restrict(bm, ps_mask), f)
                r.update({"design": f"{d}|{ps_name}", "description": f"{desc[d]}; positive set {ps_name}", "feature": f,
                          "interval_type": "Cys-Audit convention",
                          "status_975_margin05": status_from_interval(r["estimate"], r["ci975_low"], r["ci975_high"], 0.0, 0.5)[0]})
                design_rows.append(r)
    prim = store["primary_5000"]
    std_desc = {"bstar": "(b*) background standardised to the sites' theoretical-flag mix",
                "cstar": "(c*) step 9: standardised to the sites' mix of observation in the unenriched arm (no-candidate "
                         "cysteines pooled with unobserved ones)",
                "cstar_cand": "(c*cand) candidate-aware: strata no candidate / candidate not observed / observed",
                "cstar_x": "(c*x) standardised on theo x observation in the unenriched arm",
                "cstar_x3": "(c*x3) standardised on theo x observation x candidate"}
    for nm in STD:
        for f in FEATURES:
            v = prim[nm][:, FI[f]]
            lo95, hi95 = q(v, 0.95)
            lo975, hi975 = q(v, 0.975)
            pt = float(pointStd[nm][FI[f]])
            design_rows.append({"design": nm, "description": std_desc[nm], "feature": f, "estimate": pt,
                                "ci95_low": lo95, "ci95_high": hi95, "ci975_low": lo975, "ci975_high": hi975,
                                "boot_p": stats.bootstrap_p(v, 0.0), "n_positive": int(pos.sum()),
                                "n_background": int(pool.sum()),
                                "interval_type": "paired protein bootstrap over all proteins (primary seed), as step 9",
                                "status_975_margin05": status_from_interval(pt, lo975, hi975, 0.0, 0.5)[0]})
    mdesc = {"bp": "(b') 1:1 matched on the theoretical flag", "d": "(d) 1:1 matched on the empirical score",
             "dtheo": "(d_theo) 1:1 matched on the empirical score within theoretical-flag strata"}
    for nm in M:
        for f in FEATURES:
            v = prim[nm][:, FI[f]]
            lo95, hi95 = q(v, 0.95)
            lo975, hi975 = q(v, 0.975)
            design_rows.append({"design": nm, "description": mdesc[nm], "feature": f, "estimate": float(pointM[nm][FI[f]]),
                                "ci95_low": lo95, "ci95_high": hi95, "ci975_low": lo975, "ci975_high": hi975,
                                "boot_p": stats.bootstrap_p(v, 0.0), "n_positive": int(pos.sum()),
                                "n_background": int(pos.sum()),
                                "match_p2.5": float(np.quantile(MP[nm][:, FI[f]], 0.025)),
                                "match_p97.5": float(np.quantile(MP[nm][:, FI[f]], 0.975)),
                                "match_sd": float(np.std(MP[nm][:, FI[f]], ddof=1)),
                                "interval_type": "matched-sample bootstrap mixed over 200 matchings (primary seed); "
                                                 "point = median over matchings"})
    DS = pd.DataFrame(design_rows)
    DS.to_csv(f"{RESULTS}/r2_designs.csv", index=False)
    STRATA_TAB.to_csv(f"{RESULTS}/r2_standardisation_strata.csv", index=False)

    # self-check: stored per-design intervals of the matched designs (step 4b) and stored restriction intervals (step 3)
    CB = pd.read_csv(f"{RESULTS}/conditional_bootstrap_summary.csv")
    MT = pd.read_csv(f"{RESULTS}/main_table_kr_de.csv").set_index("design")
    dev = []
    for nm, stored in (("bp", "match_theo"), ("d", "match_emp_lr")):
        for f in FEATURES:
            cb = CB[(CB.design == stored) & (CB.feature == f)].iloc[0]
            r = DS[(DS.design == nm) & (DS.feature == f)].iloc[0]
            dev += [abs(r.estimate - cb.estimate_median_200_match_seeds), abs(r.ci95_low - cb.cond_ci95_low),
                    abs(r.ci95_high - cb.cond_ci95_high)]
    for nm, stored in (("a", "proteome"), ("b", "theoretical"), ("c", "empirical"), ("e", "observed")):
        for f in FEATURES:
            r = DS[(DS.design == nm) & (DS.feature == f)].iloc[0]
            dev += [abs(r.estimate - MT.loc[stored, f"{f}_estimate"]), abs(r.ci95_low - MT.loc[stored, f"{f}_ci95_low"]),
                    abs(r.ci95_high - MT.loc[stored, f"{f}_ci95_high"])]
    checks["max_abs_dev_vs_stored_per_design_intervals"] = float(max(dev))
    assert checks["max_abs_dev_vs_stored_per_design_intervals"] < 1e-12

    # ------------------------------------------------------------------ dtheo: per-matching status and balance
    st_rows = []
    for k in range(N_MATCH):
        for f in ("kr_dist", "kr_prox"):
            cnt_all = M["dtheo"][k][:, FI[f], :]
            cnt = cnt_all[cnt_all.sum(axis=1) > 0]
            point, (lo, hi), _, _ = stats.boot_log2_or(cnt, BOOT_REPS, FEATURE_SEED[f], 0.975)
            dtheo_rows[k][f"{f}_ci975_low"], dtheo_rows[k][f"{f}_ci975_high"] = lo, hi
            dtheo_rows[k][f"{f}_status_975"] = status_from_interval(point, lo, hi, 0.0, 0.5)[0]
    DT = pd.DataFrame(dtheo_rows)
    DT.to_csv(f"{RESULTS}/r2_dtheo_matchings.csv", index=False)
    dtheo_status = {f: DT[f"{f}_status_975"].value_counts(normalize=True).round(4).to_dict()
                    for f in ("kr_dist", "kr_prox")}

    # ------------------------------------------------------------------ G. Holm, round-2 family (distal K/R)
    hparts = []
    for fam_name, fam_list in (("round2_answer", ROUND2_FAMILY), ("note16_final", ["b-e"] + ROUND2_FAMILY)):
        fam = PC[(PC.run == "big_50000") & (PC.kind == "contrast") & (PC.feature == "kr_dist") &
                 PC.name.isin(fam_list)].set_index("name").loc[fam_list]
        p = fam.boot_p.to_numpy()
        order = np.argsort(p, kind="stable")
        adj = np.empty(len(p))
        run_max = 0.0
        for rank, i in enumerate(order):
            run_max = max(run_max, min(1.0, (len(p) - rank) * p[i]))
            adj[i] = run_max
        hparts.append(pd.DataFrame({"family": fam_name, "family_size": len(p), "contrast": fam_list,
                                    "point": fam.point.to_numpy(), "ci95_low": fam.ci95_low.to_numpy(),
                                    "ci95_high": fam.ci95_high.to_numpy(), "ci975_low": fam.ci975_low.to_numpy(),
                                    "ci975_high": fam.ci975_high.to_numpy(), "boot_p": p, "holm_adjusted_p": adj}))
    H = pd.concat(hparts, ignore_index=True)
    H.to_csv(f"{RESULTS}/r2_holm_round2.csv", index=False)

    # ------------------------------------------------------------------ E. decomposition of the (c) residual
    kr = flags["kr_dist"]
    dec = []

    def add(set_name, m):
        dec.append({"set": set_name, "n": int(m.sum()),
                    "share_kr_dist": float(kr[m].mean()) if m.any() else np.nan,
                    "share_kr_prox": float(flags["kr_prox"][m].mean()) if m.any() else np.nan,
                    "share_de_dist": float(flags["de_dist"][m].mean()) if m.any() else np.nan,
                    "share_not_theo": float((~theo[m]).mean()) if m.any() else np.nan,
                    "share_theoP": float(theoP[m].mean()) if m.any() else np.nan})

    for nm, base in (("sites", pos), ("background (pool)", pool), ("(c) background", pool & empf),
                     ("(e) background", pool & det)):
        add(nm, base)
        for t_ in (True, False):
            add(f"{nm}, theo={int(t_)}", base & (theo == t_))
        if nm in ("sites", "background (pool)"):
            for t_ in (True, False):
                for e_ in (True, False):
                    add(f"{nm}, theo={int(t_)}, emp={int(e_)}", base & (theo == t_) & (empf == e_))
    add("background, emp=0, no candidate", pool & ~empf & ~cand)
    add("background, emp=0, no candidate, theo=0", pool & ~empf & ~cand & ~theo)
    add("background, emp=0, theo=0 (pooled in the unobserved stratum of (c*))", pool & ~empf & ~theo)
    # why the 290 observed-but-not-theoretically-detectable background cysteines fail the rule
    m290 = pool & empf & ~theo
    cmap = pd.read_csv(f"{RESULTS}/cand_cys_map.csv", dtype={"protein": str})
    ids = cmap[cmap.empirical_detected == 1]
    key290 = set(zip(prot[m290], T.position.to_numpy()[m290]))
    ids290 = ids[[(p_, int(x)) in key290 for p_, x in zip(ids.protein, ids.position)]]
    per = ids290.groupby(["protein", "position"]).agg(n_ids=("cand_row", "size"), min_len=("length", "min"),
                                                      max_len=("length", "max")).reset_index()
    assert len(per) == int(m290.sum())
    tp = dict(zip(zip(prot, T.position.to_numpy()), theoP))
    krd = dict(zip(zip(prot, T.position.to_numpy()), kr))
    per["theoP"] = [bool(tp[(a, int(b))]) for a, b in zip(per.protein, per.position)]
    per["kr_dist"] = [bool(krd[(a, int(b))]) for a, b in zip(per.protein, per.position)]
    per["all_ids_longer_than_30"] = per.min_len > 30
    cat = np.where(per.theoP & ~per.all_ids_longer_than_30, "detectable under Trypsin/P (cleavage before proline); "
                   "identified in a peptide of <= 30 residues",
                   np.where(per.all_ids_longer_than_30 & ~per.theoP, "identified only in peptides longer than 30 residues",
                            np.where(per.theoP & per.all_ids_longer_than_30, "both", "neither (other)")))
    per["category"] = cat
    for cname, g in per.groupby("category"):
        dec.append({"set": f"the 290: {cname}", "n": int(len(g)), "share_kr_dist": float(g.kr_dist.mean())})
    DEC = pd.DataFrame(dec)
    DEC.to_csv(f"{RESULTS}/r2_c_residual_decomposition.csv", index=False)
    per.to_csv(f"{RESULTS}/r2_c_not_theo_cysteines.csv", index=False)
    # also: the (e) theo=0 cysteines under Trypsin/P
    checks["e_bg_theo0_share_theoP"] = float(theoP[pool & det & ~theo].mean())
    checks["c_bg_theo0_share_theoP"] = float(theoP[m290].mean())

    # ------------------------------------------------------------------ F. search design of the arms
    probe = json.load(open(PROBE_JSON, encoding="utf-8"))["PXD063463"]
    gl = probe["mqpar_Global"]
    ab = probe["mqpar_ABE"]
    ab_sets = sorted({(tuple(g["enzymes"]), tuple(g["variableModifications"]), tuple(g["fixedModifications"]))
                      for g in ab["parameter_groups"]})
    sd = {"label": "POST HOC (round 2): search design of the PXD063463 arms, read from the phase-1 probe record of the "
                   "deposit's mqpar_Global.xml and mqpar_ABE.xml; the MaxQuant tables themselves are not available "
                   "offline and were not read here",
          "source": PROBE_JSON, "source_sha256": sha256_file(PROBE_JSON),
          "probe_script_recorded_in_source": json.load(open(PROBE_JSON, encoding="utf-8")).get("script"),
          "probe_script_sha256_recorded_in_source": json.load(open(PROBE_JSON, encoding="utf-8")).get("script_sha256"),
          "global_arm": {"maxquant_version": gl["version"], "n_raw_files": gl["n_files"],
                         "experiments": sorted(set(gl["experiments"])),
                         "parameter_groups": gl["parameter_groups"]},
          "abe_arms": {"maxquant_version": ab["version"], "n_raw_files": ab["n_files"],
                       "n_parameter_groups": len(ab["parameter_groups"]),
                       "distinct_enzyme_variable_fixed_sets": [{"enzymes": list(e), "variableModifications": list(v),
                                                                 "fixedModifications": list(fx)} for e, v, fx in ab_sets]},
          "reading": ("The global arm was searched with carbamidomethyl (C) as a FIXED modification and without "
                      "N-ethylmaleimide (C), enzyme Trypsin/P; the ABE arms were searched with carbamidomethyl (C) and "
                      "NEM (C) both VARIABLE and no fixed modification, enzyme Trypsin (or AspN/chymotrypsin/GluC). "
                      "As searched, the global arm's identifications therefore carry no readout of blocking or exchange "
                      "(every cysteine is assumed carbamidomethylated), so it supports an empirical detectability "
                      "background only. Whether the global aliquot was itself exposed to NEM cannot be read from the "
                      "search parameters; the phase-3 pre-registration describes it as standard, unenriched shotgun "
                      "LC-MS/MS of the same cells.")}
    g0 = gl["parameter_groups"][0]
    assert g0["fixedModifications"] == ["Carbamidomethyl (C)"] and "NEM (C)" not in g0["variableModifications"]
    assert all("NEM (C)" in g["variableModifications"] and "Carbamidomethyl (C)" in g["variableModifications"]
               and not g["fixedModifications"] for g in ab["parameter_groups"])
    dump_json(sd, f"{RESULTS}/r2_global_arm_search_design.json")

    # ------------------------------------------------------------------ H. comparison with the round-2 verifier
    cmp_rows = []
    VC = pd.read_csv(f"{VR2}/v2_combined_restriction.csv")
    vmap = {"b_theo": "b", "c_emp": "c", "bc_theo_and_emp": "bc", "bc_symmetric": "bc_sym",
            "c_emp_not_theo": "c_not_theo", "e_observed": "e"}
    for _, r in VC.iterrows():
        if r.quantity in vmap:
            mine = DS[(DS.design == vmap[r.quantity]) & (DS.feature == r.feature)].iloc[0]
            cmp_rows.append({"source": "v2_combined_restriction.csv", "quantity": r.quantity, "feature": r.feature,
                             "mine_point": mine.estimate, "verifier_point": r.point,
                             "mine_ci95": f"[{mine.ci95_low:.4f}, {mine.ci95_high:.4f}]",
                             "verifier_ci95": f"[{r.lo95:.4f}, {r.hi95:.4f}]",
                             "abs_diff_point": abs(mine.estimate - r.point),
                             "abs_diff_bounds": max(abs(mine.ci95_low - r.lo95), abs(mine.ci95_high - r.hi95)),
                             "note": "per-design interval, Cys-Audit convention, same tool seed"})
        elif "(paired)" in str(r.quantity):
            x, y = [s.strip() for s in r.quantity.replace("(paired)", "").split(" - ")]
            nm = f"{vmap[x]}-{vmap[y]}"
            if nm not in CONTRASTS:
                continue
            mine = pc(nm, r.feature)
            cmp_rows.append({"source": "v2_combined_restriction.csv", "quantity": r.quantity, "feature": r.feature,
                             "mine_point": mine.point, "verifier_point": r.point,
                             "mine_ci95": f"[{mine.ci95_low:.4f}, {mine.ci95_high:.4f}]",
                             "verifier_ci95": f"[{r.lo95:.4f}, {r.hi95:.4f}]",
                             "abs_diff_point": abs(mine.point - r.point),
                             "abs_diff_bounds": max(abs(mine.ci95_low - r.lo95), abs(mine.ci95_high - r.hi95)),
                             "note": "paired; mine 50,000 resamples (seed 20263930), verifier 20,000 (own seed)"})
    VS = pd.read_csv(f"{VR2}/v2_cross_std_features.csv")
    smap = {"b_star": "bstar", "c_star": "cstar", "c_star_cand": "cstar_cand", "c_star_x": "cstar_x", "e": "e"}
    for _, r in VS.iterrows():
        qn = str(r.quantity)
        if qn in smap:
            if smap[qn] == "e":
                continue
            mine = DS[(DS.design == smap[qn]) & (DS.feature == r.feature)].iloc[0]
            cmp_rows.append({"source": "v2_cross_std_features.csv", "quantity": qn, "feature": r.feature,
                             "mine_point": mine.estimate, "verifier_point": r.point,
                             "mine_ci95": f"[{mine.ci95_low:.4f}, {mine.ci95_high:.4f}]",
                             "verifier_ci95": f"[{r.lo95:.4f}, {r.hi95:.4f}]",
                             "abs_diff_point": abs(mine.estimate - r.point),
                             "abs_diff_bounds": max(abs(mine.ci95_low - r.lo95), abs(mine.ci95_high - r.hi95)),
                             "note": "standardised; mine primary paired seed 5,000, verifier 20,000 (own seed)"})
        elif " - " in qn and "minus" not in qn:
            x, y = [s.strip() for s in qn.split(" - ")]
            nm = f"{smap.get(x, x)}-{smap.get(y, y)}"
            if nm not in CONTRASTS:
                continue
            mine = pc(nm, r.feature)
            cmp_rows.append({"source": "v2_cross_std_features.csv", "quantity": qn, "feature": r.feature,
                             "mine_point": mine.point, "verifier_point": r.point,
                             "mine_ci95": f"[{mine.ci95_low:.4f}, {mine.ci95_high:.4f}]",
                             "verifier_ci95": f"[{r.lo95:.4f}, {r.hi95:.4f}]",
                             "abs_diff_point": abs(mine.point - r.point),
                             "abs_diff_bounds": max(abs(mine.ci95_low - r.lo95), abs(mine.ci95_high - r.hi95)),
                             "note": "paired; mine 50,000 resamples (seed 20263930), verifier 20,000 (own seed)"})
    for x, y, nm in (("(bc - b) KR minus DE", None, ("bc", "b")),):
        vr = VC[VC.quantity == x].iloc[0]
        mine = SH[(SH.run == "big_50000") & (SH.step == "b -> bc")].iloc[0]
        cmp_rows.append({"source": "v2_combined_restriction.csv", "quantity": x, "feature": "kr_dist minus de_dist",
                         "mine_point": mine.kr_minus_de_point, "verifier_point": vr.point,
                         "mine_ci95": f"[{mine.kr_minus_de_ci95_low:.4f}, {mine.kr_minus_de_ci95_high:.4f}]",
                         "verifier_ci95": f"[{vr.lo95:.4f}, {vr.hi95:.4f}]",
                         "abs_diff_point": abs(mine.kr_minus_de_point - vr.point),
                         "abs_diff_bounds": max(abs(mine.kr_minus_de_ci95_low - vr.lo95), abs(mine.kr_minus_de_ci95_high - vr.hi95)),
                         "note": "paired; different resample sets"})
    VM = json.load(open(f"{VR2}/v2_match_within_theo.json", encoding="utf-8"))
    cmp_rows.append({"source": "v2_match_within_theo.json", "quantity": "match on score within theo strata, median",
                     "feature": "kr_dist", "mine_point": float(pointM["dtheo"][0]),
                     "verifier_point": VM["match_score_within_theo_kr_dist_median"],
                     "mine_ci95": f"over matchings [{np.quantile(MP['dtheo'][:, 0], .025):.4f}, {np.quantile(MP['dtheo'][:, 0], .975):.4f}]",
                     "verifier_ci95": "over matchings [{:.4f}, {:.4f}]".format(*VM["match_score_within_theo_kr_dist_p2.5_p97.5"]),
                     "abs_diff_point": abs(float(pointM["dtheo"][0]) - VM["match_score_within_theo_kr_dist_median"]),
                     "note": "mine 200 matchings, clipped logit (EPS 1e-6) as in (d); verifier 100 matchings, unclipped"})
    for key, mine_v in (("share_theo0_sites", float((~theo[pos]).mean())), ("share_theo0_bg_c", float((~theo[pool & empf]).mean())),
                        ("share_theo0_bg_e", float((~theo[pool & det]).mean())),
                        ("share_theo0_bg_proteome", float((~theo[pool]).mean()))):
        cmp_rows.append({"source": "v2_match_within_theo.json", "quantity": key, "feature": "",
                         "mine_point": mine_v, "verifier_point": VM[key], "abs_diff_point": abs(mine_v - VM[key]),
                         "note": "count share"})
    VX = json.load(open(f"{VR2}/v2_extra.json", encoding="utf-8"))
    for key, mine_v in (("emp0_bg_n", int((pool & ~empf).sum())), ("emp0_bg_nocand_n", int((pool & ~empf & ~cand).sum())),
                        ("emp0_bg_theo0_n", int((pool & ~empf & ~theo).sum())),
                        ("emp0_bg_theo0_kr_share", float(kr[pool & ~empf & ~theo].mean())),
                        ("emp0_bg_theo1_kr_share", float(kr[pool & ~empf & theo].mean())),
                        ("sites_emp0_theo0", int((pos & ~empf & ~theo).sum())),
                        ("sites_emp1_theo0", int((pos & empf & ~theo).sum()))):
        cmp_rows.append({"source": "v2_extra.json", "quantity": key, "feature": "kr_dist" if "share" in key else "",
                         "mine_point": mine_v, "verifier_point": VX[key], "abs_diff_point": abs(mine_v - VX[key]),
                         "note": "count / share"})
    CMP = pd.DataFrame(cmp_rows)
    CMP.to_csv(f"{RESULTS}/r2_verifier_comparison.csv", index=False)

    # ------------------------------------------------------------------ summary
    def dsel(d, f="kr_dist"):
        r = DS[(DS.design == d) & (DS.feature == f)].iloc[0]
        return {k: (float(r[k]) if isinstance(r[k], (int, float, np.floating, np.integer)) and not pd.isna(r[k]) else r[k])
                for k in r.index if not (isinstance(r[k], float) and np.isnan(r[k]))}

    summary = {
        "label": "POST HOC revision analysis after verification round 2 (2026-09-30); not registered",
        "self_checks": checks,
        "designs_kr_dist": {d: dsel(d) for d in ("a", "b", "c", "e", "bc", "bc_sym", "c_not_theo", "e_theo", "bstar", "cstar",
                                                 "cstar_cand", "cstar_x", "cstar_x3", "bp", "d", "dtheo", "bc_same",
                                                 "e_same", "b_same", "c_same")},
        "designs_kr_prox": {d: dsel(d, "kr_prox") for d in ("b", "c", "e", "bc", "bc_sym", "bstar", "cstar_cand", "cstar_x", "dtheo")},
        "designs_de_dist": {d: dsel(d, "de_dist") for d in ("b", "c", "e", "bc", "bc_sym", "bstar", "cstar_cand", "cstar_x", "dtheo")},
        "contrasts_big_run": {f: {nm: {k: float(pc(nm, f)[k]) for k in ("point", "ci95_low", "ci95_high", "ci975_low", "ci975_high", "boot_p")}
                                  for nm in CONTRASTS} for f in FEATURES},
        "kr_specific_part_big_run": SH[SH.run == "big_50000"].to_dict(orient="records"),
        "holm_round2": H.to_dict(orient="records"),
        "dtheo_status_975": dtheo_status,
        "dtheo_balance": {"smd_logit_score_median": float(DT.smd_logit_score.median()),
                          "auc_median": float(DT.auc_score_sites_vs_controls.median()),
                          "share_theo_controls_median": float(DT.share_theo_controls.median()),
                          "share_theo_sites": float(theo[pos].mean()),
                          "share_emp_obs_controls_median": float(DT.share_emp_obs_controls.median()),
                          "share_emp_obs_sites": float(empf[pos].mean()),
                          "n_controls_in_no_candidate_mean": float(DT.n_controls_in_no_candidate.mean())},
        "standardisation_strata": STRATA_TAB.to_dict(orient="records"),
        "decomposition": DEC.to_dict(orient="records"),
        "search_design": sd,
        "verifier_comparison_max_abs_point_diff": float(CMP.abs_diff_point.max()),
        "seeds": {"primary_paired": PAIRED_SEED, "big_run": BIG_SEED, "big_run_reps": BIG_REPS,
                  "stored_matchings_bp_d": "20260941 + 1000 + k, k = 0..199",
                  "dtheo_matchings": f"{DTHEO_SEED} + k, k = 0..199", "tool_feature_seeds": FEATURE_SEED},
        "holm_family_round2": ROUND2_FAMILY,
    }
    dump_json(summary, f"{RESULTS}/r2_summary.json")

    with pd.option_context("display.width", 250, "display.max_columns", 40):
        cols = ["design", "feature", "estimate", "ci95_low", "ci95_high", "ci975_low", "ci975_high", "n_positive",
                "n_background", "share_flagged_pos", "share_flagged_bg", "status_975_margin05"]
        print(DS[DS.feature.isin(["kr_dist"])][cols + ["match_p2.5", "match_p97.5"]].round(4).to_string())
        print(DS[DS.feature.isin(["kr_prox", "de_dist"]) & DS.design.isin(["b", "c", "e", "bc", "bc_sym", "bstar", "cstar_cand",
                                                                           "cstar_x", "dtheo", "bp", "d"])][cols].round(4).to_string())
        big = PC[(PC.run == "big_50000") & (PC.kind == "contrast")]
        print(big[big.feature.isin(["kr_dist", "kr_prox", "de_dist"])][["name", "feature", "point", "ci95_low", "ci95_high",
                                                                          "ci975_low", "ci975_high", "boot_p"]].round(4).to_string())
        print(SH.round(4).to_string())
        print(H.round(4).to_string())
        print(STRATA_TAB.round(4).to_string())
        print(DEC.round(4).to_string())
        print(DT.describe().round(4).to_string())
        print(json.dumps(dtheo_status, indent=1))
        print(CMP.round(4).to_string())
        print(json.dumps({k: v for k, v in checks.items() if k != "dev_vs_step9_by_contrast"}, indent=1))


if __name__ == "__main__":
    main()
