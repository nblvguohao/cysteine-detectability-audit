# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered, not pre-specified.

Task 3, part 2: claim-matched NULL simulations for every verdict-deciding interval whose bound lies
near 0 (w* between 0.5 and 1.5 in claims_intervals.csv). Each configuration copies the claim's
own cohort shape: number of bootstrap clusters, mean cluster size (singleton clusters for the
protein-level claims), number of sites and positives, attribute prevalence, and (for the
stratified caliber) the propensity AUC, via the covariate shift delta = sqrt(2) * Phi^-1(AUC).
Site-level cohorts carry protein random intercepts (SD 1 on the logit scale for label and
attribute, independent: no confounding), protein-level cohorts are iid. There is no
attribute-label association, so the 95th percentile of |point|/half-width-toward-0 is the factor
w0 by which that interval type must be scaled to have a 5% type-I error at that cohort shape.

Designs (as in run_grid.py): full (baseline), sub (matched-set proxy: all positives + as many random
negatives), mh (stratified caliber on the full cohort), rc (the pipeline's same-size random control).
Writes claims_matched_null_sim.csv (per dataset) and claims_matched_null_summary.csv.
"""
from __future__ import annotations

import os
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats as sps

import simlib as S
from summarize_grid import quantile_ci

OUT = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"
MASTER = 20260931
REPS = 5000
N_DATA = 2000

# cohort shapes from claims_intervals.csv / Supplemental Data 4-6 (clusters from the stored audits)
CONFIGS = [
    dict(claim="PERS-003", K=1089, mean_size=1.853, n=2018, n_pos=724, n_attr=1653, auc=None, designs=["sub", "rc"]),
    dict(claim="PERS-006", K=2032, mean_size=1.0, n=2032, n_pos=79, n_attr=131, auc=None, designs=["sub", "rc"]),
    dict(claim="SFE-006", K=1277, mean_size=7.656, n=9777, n_pos=1745, n_attr=7186, auc="SFE-006", designs=["sub", "mh"]),
    dict(claim="SFE-008", K=3662, mean_size=1.786, n=6542, n_pos=718, n_attr=5284, auc=None, designs=["sub", "rc"]),
    dict(claim="SFE-011", K=27751, mean_size=1.0, n=27751, n_pos=1090, n_attr=215, auc=None, designs=["sub", "rc"]),
    dict(claim="SNO-001", K=164, mean_size=8.951, n=1468, n_pos=238, n_attr=458, auc=None, designs=["sub"]),
    dict(claim="SNO-004", K=164, mean_size=8.951, n=1468, n_pos=238, n_attr=1211, auc=None, designs=["full", "sub"]),
    dict(claim="SNO-012", K=44, mean_size=8.318, n=366, n_pos=55, n_attr=174, auc="SNO-012", designs=["sub", "mh"]),
    dict(claim="SNO-014", K=45, mean_size=8.533, n=384, n_pos=56, n_attr=60, auc=None, designs=["sub", "rc"]),
]
# NEG_A cohorts (all cysteines of site-carrying proteins; every protein has >= 1 positive) are also
# simulated with that selection imposed ('nega'); 're' = same model without the selection.
NEGA = {"SNO-001", "SNO-004", "SNO-012", "SNO-014"}


def propensity_auc(claim):
    ci = pd.read_csv(os.path.join(OUT, "claims_intervals.csv"))
    pkg = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30/inputs/mcp_package"
    for f in ("Supplemental_Data_4_retest_round_b.csv", "Supplemental_Data_5_retest_round_d.csv",
              "Supplemental_Data_6_retest_round_e.csv"):
        t = pd.read_csv(os.path.join(pkg, f))
        t = t[t.claim_id == claim]
        if len(t):
            spec = ci.loc[ci.claim_id == claim, "specification"].iloc[0]
            if "specification" in t:
                t = t[t.specification.astype(str) == str(spec)] if (t.specification.astype(str) == str(spec)).any() else t
            return float(t.propensity_auc.iloc[0])
    return 0.65


def prepared(cfg, variant="re"):
    singleton = cfg["mean_size"] == 1.0
    su = sv = 0.0 if singleton else 1.0
    if variant == "nega":
        par = S.calibrate_selected(cfg["n_pos"] / cfg["n"], cfg["n_attr"] / cfg["n"], cfg["mean_size"], su, sv)
    else:
        par = S.calibrate(cfg["n_pos"] / cfg["n"], cfg["n_attr"] / cfg["n"], su, sv, 0.0)
    auc = propensity_auc(cfg["auc"]) if cfg["auc"] else 0.65
    delta = float(np.sqrt(2) * sps.norm.ppf(auc))
    return dict(cfg, su=su, sv=sv, singleton=singleton, par=par, auc_used=auc, delta=delta,
                variant=("iid" if singleton else variant))


def one(cfg, idx, rep):
    rng = np.random.default_rng(np.random.SeedSequence([MASTER, idx, rep]))
    kind = "singleton" if cfg["singleton"] else "nb"
    if cfg["variant"] == "nega":
        d = S.generate_selected(rng, cfg["K"], kind, cfg["par"], cfg["su"], cfg["sv"], 0.0,
                                cfg["par"]["mean_size_pre"], strat_shift=cfg["delta"])
    else:
        d = S.generate(rng, cfg["K"], kind, cfg["par"], cfg["su"], cfg["sv"], 0.0,
                       mean_size=cfg["mean_size"], strat_shift=cfg["delta"])
    y, x, s, cid, K = d["y"], d["x"], d["s"], d["cid"], d["K"]
    base = 2_000_000_000 + idx * 10_000_000 + rep * 100
    row = {"claim": cfg["claim"], "variant": cfg["variant"], "rep": rep, "n": int(len(y)), "n_pos": int(y.sum()), "n_attr": int(x.sum())}
    if "full" in cfg["designs"]:
        pt, (lo, hi), _, _ = S.tool_stats.boot_log2_or(S.tool_counts(cid, y, x, K), REPS, base + 1, 0.95)
        row.update(full_point=pt, full_lo=lo, full_hi=hi)
    if "mh" in cfg["designs"]:
        m = S.mh_percentile(cid, K, y, x, s, REPS, base + 2)
        row.update(mh_point=m["mh_point"], mh_lo=m["mh_lo"], mh_hi=m["mh_hi"])
    if "sub" in cfg["designs"]:
        sel = S.subsample_1to1(rng, y)
        inv, Ks = S.reindex(cid[sel])
        c = S.tool_counts(inv, y[sel], x[sel], Ks)
        pt, (lo, hi), _, _ = S.tool_stats.boot_log2_or(c, REPS, base + 3, 0.95)
        row.update(sub_point=pt, sub_lo=lo, sub_hi=hi, sub_K=Ks, sub_n_attr=int(c[:, 0].sum() + c[:, 2].sum()))
    if "rc" in cfg["designs"]:
        r = S.random_control_vectorised(cid, y, x, REPS, base + 50)
        row.update(rc_point=r["rc_point"], rc_lo=r["rc_lo"], rc_hi=r["rc_hi"])
    return row


def work(args):
    cfg, idx, start, stop = args
    return [one(cfg, idx, r) for r in range(start, stop)]


def summarise(df):
    out = []
    for (claim, variant), g in df.groupby(["claim", "variant"]):
        for design in ("full", "sub", "mh", "rc"):
            if f"{design}_lo" not in g:
                continue
            gg = g.dropna(subset=[f"{design}_lo", f"{design}_hi"])
            if gg.empty:
                continue
            lo, hi, pt = gg[f"{design}_lo"].to_numpy(), gg[f"{design}_hi"].to_numpy(), gg[f"{design}_point"].to_numpy()
            flag = float(np.mean((lo > 0) | (hi < 0)))
            r = S.ratio_to_zero(pt, lo, hi)
            w0, w0l, w0h = quantile_ci(r)
            out.append(dict(claim=claim, variant=variant, design=design, n_datasets=len(gg), type1=flag,
                            type1_mcse=float(np.sqrt(flag * (1 - flag) / len(gg))),
                            w0=w0, w0_ci_lo=w0l, w0_ci_hi=w0h, mean_width=float(np.mean(hi - lo)),
                            mean_point=float(np.mean(pt))))
    return pd.DataFrame(out)


def main():
    cfgs = [prepared(c, "re") for c in CONFIGS] + [prepared(c, "nega") for c in CONFIGS if c["claim"] in NEGA]
    tasks = []
    for idx, c in enumerate(cfgs):
        for st in range(0, N_DATA, 100):
            tasks.append((c, idx, st, min(N_DATA, st + 100)))
    t0 = time.time()
    rows = []
    from multiprocessing import Pool
    with Pool(4) as pool:
        for res in pool.imap_unordered(work, tasks):
            rows.extend(res)
    df = pd.DataFrame(rows).sort_values(["claim", "variant", "rep"])
    df.to_csv(os.path.join(OUT, "claims_matched_null_sim.csv.gz"), index=False)
    summ = summarise(df)
    meta = pd.DataFrame([{k: v for k, v in c.items() if k not in ("par", "designs")} |
                         {"designs": ",".join(c["designs"])} | {f"par_{k}": v for k, v in c["par"].items()}
                         for c in cfgs])
    meta.to_csv(os.path.join(OUT, "claims_matched_null_configs.csv"), index=False)
    summ.to_csv(os.path.join(OUT, "claims_matched_null_summary.csv"), index=False)
    pd.set_option("display.width", 200)
    print(summ.round(4).to_string())
    print("seconds", round(time.time() - t0))


if __name__ == "__main__":
    main()
