"""Step 11 (POST HOC revision analysis, item D_empirical_background), written 2026-09-30 AFTER adversarial
verification round 3. Nothing here was registered or pre-specified; it was added because the round-3 verifier showed
that the within-strata score matching (f') had not been given the re-matching bootstrap that step 04 gave (d) and
(b'), and that its PASS status and its "removal" reading did not survive that check.

What it does: the step-04 re-matching bootstrap, extended to (f') and to paired differences.
  * Replicates are the step-04 protein multiplicities (Cys-Audit stats.multiplicities over all 1,991 proteins,
    5,000 replicates, seed 20260930 + 7); the replicate rows are built exactly as in step 04 (each drawn copy of a
    protein is its own pseudo-protein).
  * In every replicate the 1:1 greedy nearest-neighbour matching (this item's matching.py) is redone for
      (b')  on the theoretical flag,
      (d)   on logit(empirical score),
      (f')  on logit(empirical score) within strata of the theoretical flag (stratum 1 first, then 0, one generator),
    each with a fresh generator default_rng([20260930 + 11, 7, r]) -- the step-04 convention. The restriction designs
    (a), (b), (e) and (f) are computed on the same rows, so every difference is paired within replicate.
  * Run "A" (primary) reproduces step 04's stored (a), (e), (b') and (d) replicate values (self-check). Run "B" is an
    independent replication (multiplicities seed 20260930 + 6000, matching generators default_rng([20260930 + 6011,
    7, r])) that shows the Monte-Carlo stability of the quantiles.
  * All four features (K/R and D/E, distal and proximal); 95% and 97.5% percentile intervals; Cys-Audit status of the
    97.5% interval with margin 0.5 (the tool rule); two-sided bootstrap p (Cys-Audit stats.bootstrap_p).
  * The conditional (matched-sample) values of step 10 are read back only for side-by-side reporting.

Neither scheme is known to be valid here: the re-matching bootstrap is not valid for nearest-neighbour matching
estimators (Abadie & Imbens, Econometrica 2008), and resampling matched samples ignores the uncertainty of the matching
itself. This step shows how much the matched designs' conclusions depend on the choice; it does not replace one
interval by a "correct" one.

  * Balance diagnostics of the two score-matched designs, (d) and (f'): standardized mean difference and mean
    |pair distance| of the logit score, the share of controls identified in the unenriched arm, and the structure of
    the pairs (same protein; within 12 residues; another protein with an identical +-12-residue context, read from
    the UniProt 2026_03 FASTA; another protein with the identical best candidate peptide, read from the phase-3
    candidate table; identical geometry-only score; K/R-distal agreement and the controls' K/R-distal share given the
    site's flag), in every replicate and in the 200 stored matchings of the original data (which are re-done here and
    must reproduce the stored estimates). They test the explanation, given since step 04, that re-matching shifts
    score-matched designs upward because close controls are scarcer in resamples, and alternatives to it.
  * Retained versus replaced controls: in every replicate, the pairs whose control is one that the site also received
    in any of the 200 original matchings ('kept') are separated from the rest ('new'), with the controls' and sites'
    K/R-distal shares, the mean |pair distance| and the site-control K/R dependence
    P(control K/R | site K/R) - P(control K/R | site not K/R) in each group.
  (The diagnostics were added after the first run of this step had shown the shift; they do not change any estimate.)

Outputs (results/D_empirical_background/): r3_rematch_replicates.csv.gz, r3_rematch_designs.csv,
r3_rematch_contrasts.csv, r3_original_matching_balance.csv, r3_summary.json. Run after step 10b and before step 07.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"  # spawned workers inherit it: no __pycache__ in the Cys-Audit tree

import sys

sys.dont_write_bytecode = True
import json
import time
from multiprocessing import get_context

import numpy as np
import pandas as pd

from common import BOOT_REPS, IN_CAND, IN_FASTA_GZ, MASTER_SEED, RESULTS, dump_json, read_fasta_gz, sha256_file
from cys_audit import stats  # noqa: E402
from cys_audit.status import status_from_interval  # noqa: E402
from matching import logit, match_design, nn_match

FEATURES = ["kr_dist", "kr_prox", "de_dist", "de_prox"]
DESIGNS = ["a", "b", "e", "f", "bp", "d", "fprime"]
RUNS = {"A": {"paired_seed": MASTER_SEED + 7, "match_seed": MASTER_SEED + 11,
              "note": "step-04 multiplicities and matching-generator convention (primary)"},
        "B": {"paired_seed": MASTER_SEED + 6000, "match_seed": MASTER_SEED + 6011,
              "note": "independent replication run (new seeds, same construction)"}}
CONTRASTS = [("fprime", "e"), ("fprime", "bp"), ("fprime", "b"), ("fprime", "d"), ("fprime", "f"),
             ("d", "bp"), ("d", "e"), ("bp", "e"), ("f", "e"), ("f", "b"), ("b", "e")]
STEP10_NAME = {"a": "a", "b": "b", "e": "e", "f": "bc", "bp": "bp", "d": "d", "fprime": "dtheo"}
VR3 = f"{RESULTS}/verify_r3"


def _window_ids(T, half=12):
    """Integer id of each cysteine's +-12-residue sequence context (UniProt 2026_03; '-' beyond the termini)."""
    seqs = read_fasta_gz(IN_FASTA_GZ)
    wins = []
    for p, x in zip(T.protein, T.position):
        s = seqs[p]
        assert s[x - 1] == "C"
        lo, hi = x - 1 - half, x + half
        wins.append("-" * max(0, -lo) + s[max(0, lo):min(len(s), hi)] + "-" * max(0, hi - len(s)))
    return pd.factorize(pd.Series(wins))[0].astype(np.int64)


def _best_peptide_ids(T):
    """Integer id of each cysteine's best candidate (the candidate whose out-of-fold probability is its score): peptide
    sequence plus the cysteine's offset in it; -1 for a cysteine in no candidate. Shared ids across proteins would be
    identical tryptic peptides (paralogs, isoforms)."""
    C = pd.read_csv(IN_CAND, dtype={"protein": str}, usecols=["protein", "start", "sequence"])
    C["cand_row"] = np.arange(len(C))  # as in step 01
    O = pd.read_csv(f"{RESULTS}/candidate_oof_scores.csv.gz", usecols=["cand_row", "lr_full"])
    Mp = pd.read_csv(f"{RESULTS}/cand_cys_map.csv", dtype={"protein": str}, usecols=["cand_row", "protein", "position"])
    Mp = Mp.merge(O, on="cand_row").merge(C[["cand_row", "sequence", "start"]], on="cand_row")
    Mp = Mp.sort_values(["protein", "position", "lr_full", "cand_row"], ascending=[True, True, False, True])
    B = Mp.drop_duplicates(["protein", "position"])
    key = T[["protein", "position"]].merge(B, on=["protein", "position"], how="left")
    assert np.nanmax(np.abs(key.lr_full.fillna(0).to_numpy() - T.score_lr_full.to_numpy())) < 1e-6
    ids = pd.factorize(key.sequence.fillna("") + "|" + (key.position - key.start).fillna(-1).astype(int).astype(str))[0]
    ids = ids.astype(np.int64)
    ids[key.sequence.isna().to_numpy()] = -1
    return ids


def _load():
    T = pd.read_csv(f"{RESULTS}/cys_scores_trypsin_arm.csv", dtype={"protein": str})
    return {"lab": T.label.to_numpy(), "det": T.detected.to_numpy(), "prot": T.protein.to_numpy(),
            "posn": T.position.to_numpy().astype(np.int64), "wid": _window_ids(T), "pid": _best_peptide_ids(T),
            "lg_geom": logit(T.score_lr_geom.to_numpy()),
            "theo": T.theo_trypsin.to_numpy() == 1, "emp": T.emp_obs.to_numpy() == 1,
            "lg": logit(T.score_lr_full.to_numpy()), "theo_f": T.theo_trypsin.to_numpy().astype(float),
            "flags": {f: T[f].to_numpy().astype(bool) for f in FEATURES}}


def _lor(pos_m, bg_m, fl):
    return float(stats.log2_or((pos_m & fl).sum(), (pos_m & ~fl).sum(), (bg_m & fl).sum(), (bg_m & ~fl).sum()))


def _pairs_global(pos_m, pool_m, score, rng):
    """Exactly matching.match_design (no caliper, not within protein), but returning the matched pairs."""
    pool_idx = np.flatnonzero(pool_m)
    pos_idx = np.flatnonzero(pos_m)
    m = nn_match(score[pos_idx], score[pool_idx], rng)
    ok = m >= 0
    return pos_idx[ok], pool_idx[m[ok]]


def _pairs_within(pos_m, pool_m, strata, score, rng):
    """As step 10's match_within: stratum 1 (True) first, then 0, one generator."""
    si, ci = [], []
    for s_ in (True, False):
        pi = np.flatnonzero(pos_m & (strata == s_))
        cj = np.flatnonzero(pool_m & (strata == s_))
        m = nn_match(score[pi], score[cj], rng)
        ok = m >= 0
        si.append(pi[ok])
        ci.append(cj[m[ok]])
    return np.concatenate(si), np.concatenate(ci)


def _role(n, si, ci):
    role = np.full(n, -1, dtype=np.int64)
    role[si] = 1
    role[ci] = 0
    return role


def _balance(si, ci, lg_, emp_, code_, posn_, kr_, wid_, geom_, pid_):
    """Score balance of a matched sample: standardized mean difference and mean |pair distance| of the logit score
    (sites minus controls), the share of controls identified in the unenriched arm, and pair structure: share of
    pairs whose control lies in the same protein as its site (any copy of it), within 12 residues of it, in another
    protein with an identical +-12-residue context, or with an identical geometry-only score (length and missed
    cleavages of the best candidate); mean |difference| of that geometry score; share of pairs that agree on the
    K/R-distal flag; K/R-distal share of the controls, overall and given the site's own K/R-distal flag."""
    xs, xc = lg_[si], lg_[ci]
    same = code_[si] == code_[ci]
    near = same & (np.abs(posn_[si] - posn_[ci]) <= 12)
    ks, kc = kr_[si], kr_[ci]
    return {"smd_logit": float((xs.mean() - xc.mean()) / np.sqrt((xs.var(ddof=1) + xc.var(ddof=1)) / 2)),
            "mean_abs_pair_dist": float(np.mean(np.abs(xs - xc))),
            "share_emp_controls": float(emp_[ci].mean()), "share_emp_sites": float(emp_[si].mean()),
            "share_pairs_same_protein": float(same.mean()), "share_pairs_within_12": float(near.mean()),
            "share_pairs_identical_context_other_protein": float(((wid_[si] == wid_[ci]) & ~same).mean()),
            "share_pairs_shared_best_peptide_other_protein": float(((pid_[si] == pid_[ci]) & (pid_[si] >= 0)
                                                                    & ~same).mean()),
            "share_pairs_identical_geom_score": float(np.isclose(geom_[si], geom_[ci], rtol=0, atol=1e-9).mean()),
            "mean_abs_geom_score_diff": float(np.mean(np.abs(geom_[si] - geom_[ci]))),
            "share_pairs_kr_agree": float((ks == kc).mean()),
            "share_kr_controls_given_site_kr1": float(kc[ks].mean()),
            "share_kr_controls_given_site_kr0": float(kc[~ks].mean()),
            "share_kr_controls": float(kr_[ci].mean()),
            "share_kr_controls_same_protein_pairs": float(kr_[ci][same].mean()) if same.any() else float("nan"),
            "share_kr_controls_other_protein_pairs": float(kr_[ci][~same].mean()) if (~same).any() else float("nan"),
            "share_kr_sites_same_protein_pairs": float(kr_[si][same].mean()) if same.any() else float("nan"),
            "share_kr_sites_other_protein_pairs": float(kr_[si][~same].mean()) if (~same).any() else float("nan")}


def _original_control_sets(X):
    """For (d) and (f'): the set of controls (row indices) that each site received in any of the 200 stored
    matchings of the original data (seeds as in step 10)."""
    lab, lg, theo = X["lab"], X["lg"], X["theo"]
    pos0, pool0 = lab == 1, lab == 0
    sets = {"d": {}, "fprime": {}}
    for k in range(200):
        for nm, (si, ci) in (("d", _pairs_global(pos0, pool0, lg, np.random.default_rng(MASTER_SEED + 11 + 1000 + k))),
                             ("fprime", _pairs_within(pos0, pool0, theo, lg,
                                                      np.random.default_rng(MASTER_SEED + 5000 + k)))):
            for a, b in zip(si.tolist(), ci.tolist()):
                sets[nm].setdefault(a, set()).add(b)
    return sets


def _kept_split(gsi, gci, orig_sets, kr_, lg_):
    """Pairs of a replicate whose control (underlying row) is one the site also received in the original matchings
    ('kept') versus a replacement ('new'): shares and K/R-distal shares of controls and sites, mean |pair distance|,
    and the dependence P(control K/R | site K/R) - P(control K/R | site not K/R)."""
    kept = np.fromiter((b in orig_sets.get(a, ()) for a, b in zip(gsi.tolist(), gci.tolist())), dtype=bool,
                       count=gsi.size)
    ks, kc = kr_[gsi], kr_[gci]
    dist = np.abs(lg_[gsi] - lg_[gci])

    def dep(m):
        return float(kc[m & ks].mean() - kc[m & ~ks].mean()) if (m & ks).any() and (m & ~ks).any() else float("nan")
    out = {"share_pairs_kept": float(kept.mean()), "dependence_all": dep(np.ones_like(kept)),
           "dependence_kept": dep(kept), "dependence_new": dep(~kept)}
    for tag, m in (("kept", kept), ("new", ~kept)):
        out[f"share_kr_controls_{tag}"] = float(kc[m].mean()) if m.any() else float("nan")
        out[f"share_kr_sites_{tag}"] = float(ks[m].mean()) if m.any() else float("nan")
        out[f"mean_abs_pair_dist_{tag}"] = float(dist[m].mean()) if m.any() else float("nan")
    return out


def run_one(run_name):
    """One full re-matching bootstrap (5,000 replicates). Deterministic given the run's seeds."""
    cfg = RUNS[run_name]
    X = _load()
    lab, det, prot, theo, emp, lg, theo_f, flags = (X["lab"], X["det"], X["prot"], X["theo"], X["emp"], X["lg"],
                                                     X["theo_f"], X["flags"])
    orig_sets = _original_control_sets(X)
    codes, labels = stats.cluster_index(list(prot))
    K = len(labels)
    order = np.argsort(codes, kind="stable")
    codes_sorted = codes[order]
    rec = []
    r_global = 0
    t0 = time.time()
    for Wm in stats.multiplicities(K, BOOT_REPS, cfg["paired_seed"]):
        for w in Wm:
            reps = w[codes_sorted]
            rows = np.repeat(order, reps)
            run_start = np.repeat(np.cumsum(reps) - reps, reps)
            copy = np.arange(rows.size) - run_start
            pseudo = codes[rows].astype(np.int64) * 1000 + copy
            lab_r = lab[rows]
            pos_r, pool_r = lab_r == 1, lab_r == 0
            theo_r, emp_r, det_r, lg_r = theo[rows], emp[rows], det[rows] == 1, lg[rows]
            bgs = {"a": pool_r, "b": pool_r & theo_r, "e": pool_r & det_r, "f": pool_r & theo_r & emp_r}
            roles = {}
            rng = np.random.default_rng([cfg["match_seed"], 7, r_global])
            roles["bp"], _ = match_design(lab_r, pseudo, pos_r, theo_f[rows], rng)
            rng = np.random.default_rng([cfg["match_seed"], 7, r_global])
            si_d, ci_d = _pairs_global(pos_r, pool_r, lg_r, rng)  # identical to match_design(..., lg_r, rng)
            roles["d"] = _role(rows.size, si_d, ci_d)
            rng = np.random.default_rng([cfg["match_seed"], 7, r_global])
            si_f, ci_f = _pairs_within(pos_r, pool_r, theo_r, lg_r, rng)
            roles["fprime"] = _role(rows.size, si_f, ci_f)
            out = {"run": run_name, "replicate": r_global,
                   "n_pairs_fprime": int(ci_f.size), "n_sites": int(pos_r.sum())}
            code_r, posn_r, kr_r = codes[rows], X["posn"][rows], flags["kr_dist"][rows]
            wid_r, geom_r, pid_r = X["wid"][rows], X["lg_geom"][rows], X["pid"][rows]
            for nm, (si_, ci_) in (("d", (si_d, ci_d)), ("fprime", (si_f, ci_f))):
                for k_, v_ in _balance(si_, ci_, lg_r, emp_r, code_r, posn_r, kr_r, wid_r, geom_r, pid_r).items():
                    out[f"bal_{nm}_{k_}"] = v_
                for k_, v_ in _kept_split(rows[si_], rows[ci_], orig_sets[nm], flags["kr_dist"], lg).items():
                    out[f"kept_{nm}_{k_}"] = v_
            for f in FEATURES:
                fl = flags[f][rows]
                for d, bm in bgs.items():
                    out[f"{d}|{f}"] = _lor(pos_r, bm, fl)
                for d, ro in roles.items():
                    out[f"{d}|{f}"] = _lor(ro == 1, ro == 0, fl)
            rec.append(out)
            r_global += 1
            if r_global % 500 == 0:
                print(f"run {run_name}: replicates {r_global} ({time.time() - t0:.0f} s)", flush=True)
    return pd.DataFrame(rec)


def q(v, level):
    return stats.percentile_interval(v, level)


def main():
    t0 = time.time()
    ctx = get_context("spawn")
    with ctx.Pool(2) as pool:  # two worker processes (the item's limit is 4)
        parts = pool.map(run_one, list(RUNS))
    RR = pd.concat(parts, ignore_index=True)
    RR.to_csv(f"{RESULTS}/r3_rematch_replicates.csv.gz", index=False, float_format="%.6g",
              compression={"method": "gzip", "mtime": 0})
    checks = {}
    # every site matched in every replicate (1:1, pools are large in both strata)
    checks["fprime_all_sites_matched_every_replicate"] = bool((RR.n_pairs_fprime == RR.n_sites).all())

    # ---- self-check 1: run A reproduces the stored step-04 replicates of (a), (e), (b') and (d)
    S4 = pd.read_csv(f"{RESULTS}/rematch_bootstrap_replicates.csv.gz")
    A = RR[RR.run == "A"].reset_index(drop=True)
    assert len(A) == len(S4) == BOOT_REPS
    dev = {}
    for mine, theirs in (("a", "proteome|all"), ("e", "observed|all"), ("bp", "match_theo"), ("d", "match_emp_lr")):
        for f in FEATURES:
            dev[f"{mine}|{f}"] = float(np.max(np.abs(A[f"{mine}|{f}"].to_numpy() - S4[f"{theirs}|{f}"].to_numpy())))
    checks["max_abs_dev_vs_step04_stored_replicates"] = float(max(dev.values()))
    checks["dev_vs_step04_by_design"] = dev
    assert checks["max_abs_dev_vs_step04_stored_replicates"] < 1e-5, dev  # stored with 6 significant digits

    # ---- self-check 2: run A's restriction designs give step 10's primary-run paired intervals (same multiplicities)
    PC2 = pd.read_csv(f"{RESULTS}/r2_paired_contrasts.csv")
    dev2 = {}
    for mine in ("a", "b", "e", "f"):
        for f in FEATURES:
            s = PC2[(PC2.run == "primary_5000") & (PC2.kind == "design_paired") & (PC2.name == STEP10_NAME[mine]) &
                    (PC2.feature == f)].iloc[0]
            lo, hi = q(A[f"{mine}|{f}"].to_numpy(), 0.95)
            dev2[f"{mine}|{f}"] = float(max(abs(lo - s.ci95_low), abs(hi - s.ci95_high)))
    checks["max_abs_dev_vs_step10_primary_restriction_intervals"] = float(max(dev2.values()))
    assert checks["max_abs_dev_vs_step10_primary_restriction_intervals"] < 1e-5, dev2

    # ---- balance diagnostics: the 200 matchings on the original data (stored seeds) against the replicates.
    # Re-doing the stored matchings also checks that the pair-returning matcher equals the stored designs.
    X = _load()
    lab, lg, theo, emp, kr = X["lab"], X["lg"], X["theo"], X["emp"], X["flags"]["kr_dist"]
    code0, _ = stats.cluster_index(list(X["prot"]))
    posn0 = X["posn"]
    pos0, pool0 = lab == 1, lab == 0
    MS = pd.read_csv(f"{RESULTS}/matching_seed_variability.csv")
    MS_d = MS[MS.design == "match_emp_lr"].reset_index(drop=True)
    DT0 = pd.read_csv(f"{RESULTS}/r2_dtheo_matchings.csv")
    orig_rows = []
    dev3 = []
    for k in range(200):
        si, ci = _pairs_global(pos0, pool0, lg, np.random.default_rng(MASTER_SEED + 11 + 1000 + k))
        est = _lor(_role(lab.size, si, ci) == 1, _role(lab.size, si, ci) == 0, kr)
        dev3.append(abs(est - float(MS_d.loc[MS_d.match_seed == MASTER_SEED + 11 + 1000 + k, "kr_dist"].iloc[0])))
        orig_rows.append({"design": "d", "k": k, "kr_dist": est,
                          **_balance(si, ci, lg, emp, code0, posn0, kr, X["wid"], X["lg_geom"], X["pid"])})
        si, ci = _pairs_within(pos0, pool0, theo, lg, np.random.default_rng(MASTER_SEED + 5000 + k))
        est = _lor(_role(lab.size, si, ci) == 1, _role(lab.size, si, ci) == 0, kr)
        dev3.append(abs(est - float(DT0.loc[DT0.match_seed == MASTER_SEED + 5000 + k, "estimate_kr_dist"].iloc[0])))
        orig_rows.append({"design": "fprime", "k": k, "kr_dist": est,
                          **_balance(si, ci, lg, emp, code0, posn0, kr, X["wid"], X["lg_geom"], X["pid"])})
    checks["max_abs_dev_original_matchings_vs_stored_estimates"] = float(max(dev3))
    assert checks["max_abs_dev_original_matchings_vs_stored_estimates"] < 1e-9
    OB = pd.DataFrame(orig_rows)
    bal = {}
    from scipy.stats import spearmanr
    BAL_COLS = ("smd_logit", "mean_abs_pair_dist", "share_emp_controls", "share_emp_sites", "share_pairs_same_protein",
                "share_pairs_within_12", "share_pairs_identical_context_other_protein",
                "share_pairs_shared_best_peptide_other_protein",
                "share_pairs_identical_geom_score", "mean_abs_geom_score_diff", "share_pairs_kr_agree",
                "share_kr_controls_given_site_kr1", "share_kr_controls_given_site_kr0", "share_kr_controls",
                "share_kr_controls_same_protein_pairs", "share_kr_controls_other_protein_pairs",
                "share_kr_sites_same_protein_pairs", "share_kr_sites_other_protein_pairs")
    KEPT_COLS = ("share_pairs_kept", "dependence_all", "dependence_kept", "dependence_new", "share_kr_controls_kept",
                 "share_kr_controls_new", "share_kr_sites_kept", "share_kr_sites_new", "mean_abs_pair_dist_kept",
                 "mean_abs_pair_dist_new")
    for nm in ("d", "fprime"):
        o = OB[OB.design == nm]
        bal[nm] = {"original_200_matchings_median": {c: float(o[c].median()) for c in ("kr_dist",) + BAL_COLS}}
        bal[nm]["original_200_matchings_median"]["dependence"] = float(
            (o.share_kr_controls_given_site_kr1 - o.share_kr_controls_given_site_kr0).median())
        for run in RUNS:
            Rr = RR[RR.run == run]
            est = Rr[f"{nm}|kr_dist"].to_numpy()
            gap = Rr[f"bal_{nm}_mean_abs_pair_dist"].to_numpy()
            smd_ = Rr[f"bal_{nm}_smd_logit"].to_numpy()
            same_ = Rr[f"bal_{nm}_share_pairs_same_protein"].to_numpy()
            q1, q3 = np.quantile(gap, [0.25, 0.75])
            s1, s3 = np.quantile(same_, [0.25, 0.75])
            rr = {f"median_{c}": float(np.nanmedian(Rr[f"bal_{nm}_{c}"])) for c in BAL_COLS}
            rr.update({
                "share_replicates_pair_dist_above_original_median":
                    float(np.mean(gap > bal[nm]["original_200_matchings_median"]["mean_abs_pair_dist"])),
                "spearman_estimate_vs_mean_abs_pair_dist": float(spearmanr(est, gap)[0]),
                "spearman_estimate_vs_smd_logit": float(spearmanr(est, smd_)[0]),
                "spearman_estimate_vs_share_pairs_same_protein": float(spearmanr(est, same_)[0]),
                "median_estimate_lowest_quartile_pair_dist": float(np.median(est[gap <= q1])),
                "median_estimate_highest_quartile_pair_dist": float(np.median(est[gap >= q3])),
                "median_estimate_lowest_quartile_same_protein": float(np.median(est[same_ <= s1])),
                "median_estimate_highest_quartile_same_protein": float(np.median(est[same_ >= s3]))})
            rr.update({f"median_kept_{c}": float(np.nanmedian(Rr[f"kept_{nm}_{c}"])) for c in KEPT_COLS})
            rr.update({f"mean_kept_{c}": float(np.nanmean(Rr[f"kept_{nm}_{c}"])) for c in KEPT_COLS})
            bal[nm][f"replicates_run_{run}"] = rr
    OB.to_csv(f"{RESULTS}/r3_original_matching_balance.csv", index=False)

    # ---- conditional (step 10) values for side-by-side reporting
    DS2 = pd.read_csv(f"{RESULTS}/r2_designs.csv")
    DT = pd.read_csv(f"{RESULTS}/r2_dtheo_matchings.csv")

    def cond_design(d, f):
        r = DS2[(DS2.design == STEP10_NAME[d]) & (DS2.feature == f)]
        if len(r) == 0:
            return {}
        r = r.iloc[0]
        return {"conditional_point": float(r.estimate), "conditional_ci95_low": float(r.ci95_low),
                "conditional_ci95_high": float(r.ci95_high), "conditional_ci975_low": float(r.ci975_low),
                "conditional_ci975_high": float(r.ci975_high),
                "conditional_interval_type": str(r.interval_type)}

    def cond_contrast(x, y, f):
        nm = f"{STEP10_NAME[x]}-{STEP10_NAME[y]}"
        r = PC2[(PC2.run == "big_50000") & (PC2.kind == "contrast") & (PC2.name == nm) & (PC2.feature == f)]
        if len(r) == 0:
            return {"conditional_name_in_step10": None}
        r = r.iloc[0]
        return {"conditional_name_in_step10": nm, "conditional_point": float(r.point),
                "conditional_ci95_low": float(r.ci95_low), "conditional_ci95_high": float(r.ci95_high),
                "conditional_ci975_low": float(r.ci975_low), "conditional_ci975_high": float(r.ci975_high),
                "conditional_boot_p": float(r.boot_p)}

    drows, crows = [], []
    for run in RUNS:
        Rr = RR[RR.run == run]
        for d in DESIGNS:
            for f in FEATURES:
                v = Rr[f"{d}|{f}"].to_numpy()
                lo95, hi95 = q(v, 0.95)
                lo975, hi975 = q(v, 0.975)
                frac = 1 - v / Rr[f"a|{f}"].to_numpy()
                flo, fhi = q(frac, 0.95)
                row = {"run": run, "design": d, "step10_name": STEP10_NAME[d], "feature": f,
                       "rematch_median": float(np.median(v)), "rematch_ci95_low": lo95, "rematch_ci95_high": hi95,
                       "rematch_ci975_low": lo975, "rematch_ci975_high": hi975,
                       "rematch_boot_p": stats.bootstrap_p(v, 0.0),
                       "fraction_removed_ci95_low": flo, "fraction_removed_ci95_high": fhi}
                row.update(cond_design(d, f))
                pt = row.get("conditional_point", float(np.median(v)))
                row["rematch_status_975_margin05"] = status_from_interval(pt, lo975, hi975, 0.0, 0.5)[0]
                if "conditional_ci975_low" in row:
                    row["conditional_status_975_margin05"] = status_from_interval(
                        pt, row["conditional_ci975_low"], row["conditional_ci975_high"], 0.0, 0.5)[0]
                drows.append(row)
        for x, y in CONTRASTS:
            for f in FEATURES:
                v = Rr[f"{x}|{f}"].to_numpy() - Rr[f"{y}|{f}"].to_numpy()
                lo95, hi95 = q(v, 0.95)
                lo975, hi975 = q(v, 0.975)
                row = {"run": run, "contrast": f"{x}-{y}", "feature": f, "rematch_median": float(np.median(v)),
                       "rematch_ci95_low": lo95, "rematch_ci95_high": hi95, "rematch_ci975_low": lo975,
                       "rematch_ci975_high": hi975, "rematch_boot_p": stats.bootstrap_p(v, 0.0),
                       "share_le_0": float(np.mean(v <= 0))}
                row.update(cond_contrast(x, y, f))
                crows.append(row)
        # K/R-specific part of the matched "with the rule" step (b') -> (f') and of (b') -> (d), distal band
        for x, y in (("fprime", "bp"), ("d", "bp"), ("fprime", "d")):
            v = (Rr[f"{x}|kr_dist"] - Rr[f"{y}|kr_dist"]).to_numpy() - (Rr[f"{x}|de_dist"] - Rr[f"{y}|de_dist"]).to_numpy()
            lo95, hi95 = q(v, 0.95)
            lo975, hi975 = q(v, 0.975)
            crows.append({"run": run, "contrast": f"({x}-{y})_KR minus ({x}-{y})_DE", "feature": "kr_dist minus de_dist",
                          "rematch_median": float(np.median(v)), "rematch_ci95_low": lo95, "rematch_ci95_high": hi95,
                          "rematch_ci975_low": lo975, "rematch_ci975_high": hi975,
                          "rematch_boot_p": stats.bootstrap_p(v, 0.0), "share_le_0": float(np.mean(v <= 0))})
    D = pd.DataFrame(drows)
    C = pd.DataFrame(crows)
    D.to_csv(f"{RESULTS}/r3_rematch_designs.csv", index=False)
    C.to_csv(f"{RESULTS}/r3_rematch_contrasts.csv", index=False)

    # ---- per-matching conditional status of (f') (step 10) for the record
    st = DT.kr_dist_status_975.value_counts(normalize=True).round(4).to_dict()

    # ---- comparison with the round-3 verifier's stored values (read only)
    vcmp = {}
    try:
        V5 = json.load(open(f"{VR3}/v3_rematch_fprime.json", encoding="utf-8"))
        V6 = json.load(open(f"{VR3}/v3_rematch_fprime_itemmatcher.json", encoding="utf-8"))
        a_fp = D[(D.run == "A") & (D.design == "fprime") & (D.feature == "kr_dist")].iloc[0]
        a_fpe = C[(C.run == "A") & (C.contrast == "fprime-e") & (C.feature == "kr_dist")].iloc[0]
        a_d = D[(D.run == "A") & (D.design == "d") & (D.feature == "kr_dist")].iloc[0]
        vcmp = {
            "fprime_rematch": {"mine_runA": [a_fp.rematch_median, a_fp.rematch_ci95_low, a_fp.rematch_ci95_high,
                                             a_fp.rematch_ci975_low, a_fp.rematch_ci975_high],
                               "verifier_own_matcher_2000": [V5["fp"]["median"], *V5["fp"]["ci95"], *V5["fp"]["ci975"]],
                               "verifier_item_matcher_1000": [V6["fprime"]["median"], *V6["fprime"]["ci95"],
                                                              *V6["fprime"]["ci975"]],
                               "order": "median, ci95 low, ci95 high, ci975 low, ci975 high"},
            "fprime_minus_e_rematch": {"mine_runA": [a_fpe.rematch_median, a_fpe.rematch_ci95_low, a_fpe.rematch_ci95_high],
                                       "verifier_own_matcher_2000": [V5["fp-e"]["median"], *V5["fp-e"]["ci95"]],
                                       "verifier_item_matcher_1000": [V6["fprime_minus_e"]["median"],
                                                                      *V6["fprime_minus_e"]["ci95"]]},
            "d_rematch": {"mine_runA": [a_d.rematch_median, a_d.rematch_ci95_low, a_d.rematch_ci95_high],
                          "verifier_own_matcher_2000": [V5["d"]["median"], *V5["d"]["ci95"]]},
            "verifier_status975_fprime": V5["fp"]["status975"],
            "verifier_files_sha256": {"v3_rematch_fprime.json": sha256_file(f"{VR3}/v3_rematch_fprime.json"),
                                      "v3_rematch_fprime_itemmatcher.json":
                                          sha256_file(f"{VR3}/v3_rematch_fprime_itemmatcher.json")}}
    except FileNotFoundError as exc:  # pragma: no cover
        vcmp = {"error": str(exc)}

    def dsel(run, d, f="kr_dist"):
        r = D[(D.run == run) & (D.design == d) & (D.feature == f)].iloc[0]
        return {k: (v if isinstance(v, str) else float(v)) for k, v in r.items()
                if k not in ("run", "design", "feature") and not (isinstance(v, float) and np.isnan(v))}

    def csel(run, c, f="kr_dist"):
        r = C[(C.run == run) & (C.contrast == c) & (C.feature == f)].iloc[0]
        return {k: (v if isinstance(v, (str, type(None))) else float(v)) for k, v in r.items()
                if k not in ("run", "contrast", "feature") and not (isinstance(v, float) and np.isnan(v))}

    summary = {
        "label": "POST HOC revision analysis after verification round 3 (2026-09-30); not registered",
        "what": "re-matching bootstrap of the matched designs, extended to (f') (score matching within theoretical strata)",
        "runs": RUNS, "replicates_per_run": BOOT_REPS,
        "self_checks": checks,
        "fprime_conditional_per_matching_status_975": st,
        "designs_kr_dist": {run: {d: dsel(run, d) for d in DESIGNS} for run in RUNS},
        "contrasts_kr_dist": {run: {f"{x}-{y}": csel(run, f"{x}-{y}") for x, y in CONTRASTS} for run in RUNS},
        "kr_specific_part": {run: {c: csel(run, c, "kr_dist minus de_dist")
                                   for c in ("(fprime-bp)_KR minus (fprime-bp)_DE", "(d-bp)_KR minus (d-bp)_DE",
                                             "(fprime-d)_KR minus (fprime-d)_DE")} for run in RUNS},
        "de_dist_shift_bp_to_fprime": {run: csel(run, "fprime-bp", "de_dist") for run in RUNS},
        "balance_diagnostics_score_matched": bal,
        "verifier_round3_comparison": vcmp,
    }  # no wall-clock time here, so that the file is byte-identical on rerun (the log records it)
    dump_json(summary, f"{RESULTS}/r3_summary.json")
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        cols = ["run", "design", "feature", "conditional_point", "conditional_ci95_low", "conditional_ci95_high",
                "conditional_ci975_low", "conditional_ci975_high", "rematch_median", "rematch_ci95_low",
                "rematch_ci95_high", "rematch_ci975_low", "rematch_ci975_high", "conditional_status_975_margin05",
                "rematch_status_975_margin05"]
        print(D[D.feature.isin(["kr_dist", "de_dist"])][cols].round(4).to_string())
        ccols = ["run", "contrast", "feature", "conditional_point", "conditional_ci95_low", "conditional_ci95_high",
                 "rematch_median", "rematch_ci95_low", "rematch_ci95_high", "rematch_ci975_low", "rematch_ci975_high",
                 "rematch_boot_p", "conditional_boot_p"]
        print(C[C.feature.isin(["kr_dist", "de_dist", "kr_dist minus de_dist"])][ccols].round(4).to_string())
    print(json.dumps({k: v for k, v in checks.items() if k != "dev_vs_step04_by_design"}, indent=1))
    print(json.dumps(st))
    print(json.dumps(bal, indent=1))
    print(json.dumps(vcmp, indent=1, default=float))
    print("elapsed", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
