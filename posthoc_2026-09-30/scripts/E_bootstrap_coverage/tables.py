# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered.

Renders the Markdown tables used in reports/E_bootstrap_coverage.md from the stored result files,
so that every number in the report is traceable to a CSV. Writes results/E_bootstrap_coverage/report_tables.md
"""
from __future__ import annotations

import json
import os
import sys

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

OUT = r"/path/to/local/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"


def f(x, d=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    return f"{x:.{d}f}"


def md(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def revision_tables():
    """Tables added in the revision after verification (round 1, 2026-09-30; POST HOC)."""
    parts = ["## Revision after verification (round 1), POST HOC"]
    p = os.path.join(OUT, "r1_rc_decomposition_by_K.csv")
    if os.path.exists(p):
        b = pd.read_csv(p)
        rows = []
        for _, r in b.iterrows():
            rows.append({"K": int(r.K),
                         "mixing factor width(rc)/width(one draw)": f"{r.mixing_factor_mean:.3f} [{r.mixing_factor_min:.3f}-{r.mixing_factor_max:.3f}]",
                         "SD(rc point)/SD(one-draw point)": f"{r.point_sd_ratio_mean:.3f} [{r.point_sd_ratio_min:.3f}-{r.point_sd_ratio_max:.3f}]",
                         "SD(whole-cohort point)/SD(one-draw point)": f"{r.sd_full_over_sub_mean:.3f}",
                         "w0_rc (nominal)": f"{r.w0_rc_mean:.3f}",
                         "over-width vs own point (1/w0_rc)": f"{r.overwidth_vs_own_point_mean:.2f}",
                         "w_single (single-draw yardstick)": f"{r.w_single_mean:.3f} [{r.w_single_min:.3f}-{r.w_single_max:.3f}]",
                         "over-width vs one calibrated draw": f"{r.overwidth_vs_single_draw_mean:.2f} [{r.overwidth_vs_single_draw_min:.2f}-{r.overwidth_vs_single_draw_max:.2f}]",
                         "w_uncal": f"{r.w_uncal_mean:.3f}", "type-I rc": f"{r.type1_rc_mean:.4f}"})
        parts.append("### R1-a. Random control: decomposition of the over-width (8 null cells per K; mean [range])\n\n" + md(pd.DataFrame(rows)))
    p = os.path.join(OUT, "r1_rc_yardsticks_claims.csv")
    if os.path.exists(p):
        c = pd.read_csv(p)
        rows = [{"claim": r.claim, "run": f"first claim-matched run, {r.variant}", "datasets": int(r.n_datasets),
                 "w0_rc (nominal)": f"{r.w0_rc:.3f} [{r.w0_rc_lo:.3f}, {r.w0_rc_hi:.3f}]",
                 "w_single": f"{r.w_single:.3f} [{r.w_single_lo:.3f}, {r.w_single_hi:.3f}]", "w_uncal": f"{r.w_uncal:.3f}",
                 "mixing": f"{r.mixing_factor:.3f}", "point SD ratio": f"{r.point_sd_ratio:.3f}"} for _, r in c.iterrows()]
        q = os.path.join(OUT, "r1_sfe008_yardsticks.json")
        if os.path.exists(q):
            s = json.load(open(q, encoding="utf-8"))
            for name, r in s["variants"].items():
                rows.append({"claim": "SFE-008", "run": f"r1_sfe008_sim {name} (singletons {r['mean_singleton_fraction']:.3f})",
                             "datasets": int(r["n_datasets"]),
                             "w0_rc (nominal)": f"{r['w0_rc']:.3f} [{r['w0_rc_lo']:.3f}, {r['w0_rc_hi']:.3f}]",
                             "w_single": f"{r['w_single']:.3f} [{r['w_single_lo']:.3f}, {r['w_single_hi']:.3f}]",
                             "w_uncal": f"{r['w_uncal']:.3f}", "mixing": f"{r['mixing_factor']:.3f}",
                             "point SD ratio": f"{r['point_sd_ratio']:.3f}"})
        parts.append("### R1-b. Random-control yardsticks at claim-matched cohort shapes (95% CI from 2,000 dataset resamples)\n\n" + md(pd.DataFrame(rows)))
    p = os.path.join(OUT, "r1_rc_intervals.csv")
    if os.path.exists(p):
        d = pd.read_csv(p)
        d = d[d.decisive]
        rows = [{"claim": r.claim_id, "stored random control": f"{r.point:.4f} [{r.lo:.4f}, {r.hi:.4f}]", "w*": f"{r.wstar:.3f}",
                 "w_single (range)": f"{r.w_single:.3f} ({r.w_single_lo:.3f}-{r.w_single_hi:.3f})",
                 "bound at w_single (range)": f"{r.near_bound_at_w_single:.3f} ({min(r.near_bound_at_w_single_lo, r.near_bound_at_w_single_hi):.3f} to {max(r.near_bound_at_w_single_lo, r.near_bound_at_w_single_hi):.3f})",
                 "bound at w_uncal": f"{r.near_bound_at_w_uncal:.3f}", "w0_rc (nominal)": f"{r.w0_rc:.3f}",
                 "flips at w_single": "yes" if r.flips_at_w_single else "no",
                 "could flip in range": "yes" if r.could_flip_in_range else "no", "source": r.source} for _, r in d.iterrows()]
        parts.append("### R1-c. Decisive random-control intervals under the single-draw yardstick (primary)\n\n" + md(pd.DataFrame(rows)))
    p = os.path.join(OUT, "r1_calibrated_verdicts.csv")
    if os.path.exists(p):
        v = pd.read_csv(p)
        cols = ["verdict_primary", "verdict_primary_low_end", "verdict_primary_high_end",
                "verdict_rc_uncalibrated_single_draw_only", "verdict_secondary_rc_nominal"]
        ch = v[(v[cols].ne(v.stored_verdict, axis=0)).any(axis=1)]
        parts.append("### R1-d. Verdicts that change under any calibration reading (primary = percentile intervals at w0, random control at the single-draw yardstick)\n\n"
                     + md(ch[["claim_id", "is_transfer", "stored_verdict", "reading"] + cols]))
    p = os.path.join(OUT, "r1_sno004_sensitivity.json")
    if os.path.exists(p):
        s = json.load(open(p, encoding="utf-8"))
        rows = []
        for fam, dd in s["adjusted_normal_approx"].items():
            for t, a in dd.items():
                rows.append({"adjustment (post hoc, normal approximation)": fam, "interval": t,
                             "adjusted interval": f"[{a['interval'][0]:.3f}, {a['interval'][1]:.3f}]",
                             "excludes 0": "yes" if a["excludes0"] else "no",
                             "null resolution (half-width, log2)": f"{a['half_width_null_resolution']:.3f}",
                             "largest OR inside": f"{a['max_odds_ratio_in_interval']:.2f}"})
        parts.append("### R1-e. SNO-004: registered reading and post hoc sensitivity\n\n"
                     f"stored verdict `{s['stored_verdict']}`; registered bucket `{s['registered_bucket_of_stored_verdict']}`; "
                     f"released tool: `{s['released_tool_claim_verdict']['verdict']}` ({s['released_tool_status']}); "
                     f"stored matched null resolution {s['stored_null_resolution_matched_half_width']:.3f}; "
                     f"P (normal approx.) matched {s['p_normal_approx']['matched']:.4f}, stratified {s['p_normal_approx']['stratified']:.4f}\n\n"
                     + md(pd.DataFrame(rows)))
    p = os.path.join(OUT, "r1_b1_null.json")
    if os.path.exists(p):
        b1 = json.load(open(p, encoding="utf-8"))
        rows = []
        for key in ("registered_first50", "registered_300", "fresh_10000"):
            r = b1[key]
            rows.append({"set": key, "datasets": r["n"], "flagged": r["flagged"],
                         "rate [Wilson 95%]": f"{r['rate']:.4f} [{r['wilson_lo']:.4f}, {r['wilson_hi']:.4f}]",
                         "below 0 / above 0": f"{r['below0']:.4f} / {r['above0']:.4f}",
                         "w0 [95% CI]": f"{r['w0']:.3f} [{r['w0_ci'][0]:.3f}, {r['w0_ci'][1]:.3f}]",
                         "clusters": f"{r['mean_clusters']:.1f}", "positives": f"{r['mean_positives']:.1f}",
                         "near-cut prevalence": f"{r['mean_near_prevalence']:.3f}"})
        parts.append("### R1-f. Registered benchmark B1 (cleavage) at strength 0, regenerated with its own generator\n\n"
                     "validation: " + json.dumps(b1["validation_against_registered"]) + "\n\n" + md(pd.DataFrame(rows)))
    p = os.path.join(OUT, "r1_confounding_anchored_cells.csv")
    if os.path.exists(p):
        c = pd.read_csv(p)
        rows = []
        for att, g in c.groupby("attribute", sort=False):
            r0 = g.iloc[0]
            rec = {"attribute (anchor)": att, "label SD": f"{r0.su:.2f}", "attribute SD": f"{r0.sv:.2f}", "rho": f"{r0.rho:.3f}",
                   "prevalence": f"{r0.prevalence:.3f}", "pooled log2 OR (no within effect)": f"{r0.theta_pooled_log2:.4f}"}
            for _, r in g.iterrows():
                rec[f"K={int(r.K)}"] = f"{r.p_excludes0:.3f} / {r.coverage_pooled_target:.3f}"
            rows.append(rec)
        parts.append("### R1-g. SUPERSEDED IN ROUND 2 (parameters from the inaccurate 20 x 20 Gauss-Hermite fit; see R2-d). "
                     "Anchored between-protein confounding: P(excludes 0) / coverage of the pooled target (2,000 datasets per cell)\n\n" + md(pd.DataFrame(rows)))
        t = pd.read_csv(os.path.join(OUT, "r1_confounding_anchored_targets.csv"))
        t["pooled log2 OR"] = t.pooled_log2or_no_within_association.map(lambda v: f"{v:.4f}")
        parts.append("### R1-h. SUPERSEDED IN ROUND 2 (see R2-a, R2-f). Pooled log2 OR with no within-protein association (Gauss-Hermite), anchors and rho +/- 1.96 SE\n\n"
                     + md(t[["attribute", "rho_setting", "rho", "su", "sv", "prevalence", "pooled log2 OR"]].round(3)))
    return parts


def revision2_tables():
    """Tables added in revision round 2 (2026-09-30; POST HOC): accurate refit of the anchoring model and the
    anchored confounding rerun (r2_anchor_refit.py, r2_confounding_anchored.py)."""
    parts = ["## Revision round 2, POST HOC"]
    lab = {"A_nbcys": "neighbouring Cys", "A_acidic": "acidic", "A_distal": "distal cleavage"}
    p = os.path.join(OUT, "r2_anchor_refit.json")
    if not os.path.exists(p):
        return parts
    r = json.load(open(p, encoding="utf-8"))
    rows, comp = [], []
    for attr, a in r["attributes"].items():
        m, w, pr, ag, old = a["mle_dense"], a["wald"], a["profile_rho"], a["mle_aghq_independent"], a["round1_gh20_fit"]
        tg = a["pooled_log2or_no_within_association"]
        rows.append({"attribute": lab[attr], "prevalence": f"{a['prevalence']:.3f}", "label SD": f"{m['sigma_u']:.3f}",
                     "attribute SD": f"{m['sigma_v']:.3f}",
                     "rho [profile 95% CI]": f"{m['rho']:.4f} [{pr['ci95'][0]:.4f}, {pr['ci95'][1]:.4f}]",
                     "rho Wald SE": f"{w['rho_se_delta']:.4f}", "LR test rho = 0: P": f"{pr['lr_test_rho0']['p']:.3g}",
                     "within log2 OR (model)": f"{m['gamma_log2']:.3f}", "nll (dense)": f"{m['nll']:.4f}",
                     "rho, independent AGHQ-25 fit": f"{ag['rho']:.4f}",
                     "pooled log2 OR, no within effect: MLE [CI ends]": f"{tg['mle']:.4f} [{tg['profile_ci_lo']:.4f}, {tg['profile_ci_hi']:.4f}]",
                     "round-1 GH20 rho (SE)": f"{old['rho']:.4f} ({old['rho_se_delta']:.4f})"})
        for where, rec in a["nll_integrator_comparison"].items():
            base = rec["dense_z8_h0.04 (primary)"]
            comp.append({"attribute": lab[attr], "theta": where, "dense primary nll": f"{base:.6f}",
                         **{k.split(" ")[0]: f"{v - base:+.2e}" for k, v in rec.items() if not k.startswith("dense_z8")}})
    parts.append("### R2-a. Anchoring model refitted with accurate integration (dense trapezoid grid; independent adaptive Gauss-Hermite fit)\n\n"
                 + md(pd.DataFrame(rows)))
    parts.append("### R2-b. Negative log-likelihood under other integrators, minus the primary dense grid (same theta)\n\n" + md(pd.DataFrame(comp)))
    pf = pd.read_csv(os.path.join(OUT, "r2_anchor_profile.csv"))
    fx = pf[pf.grid == "fixed_round"]
    keep = np.round(np.arange(-0.40, 0.2001, 0.04), 2)
    rows = []
    for attr, g in fx.groupby("attribute", sort=False):
        rec = {"attribute": lab[attr]}
        for rr in keep:
            h = g[np.isclose(g.rho, rr)]
            rec[f"{rr + 0.0:+.2f}"] = f"{h.lr_stat.iloc[0]:.2f}" if len(h) else ""
        rows.append(rec)
    parts.append("### R2-c. Profile likelihood of rho: 2 x (nll_profile - nll_min) on a fixed grid (95% cut-off 3.84)\n\n" + md(pd.DataFrame(rows)))
    p = os.path.join(OUT, "r2_confounding_anchored_cells.csv")
    if os.path.exists(p):
        c = pd.read_csv(p)
        old = pd.read_csv(os.path.join(OUT, "r1_confounding_anchored_cells.csv"))
        rows = []
        for attr in ("A_nbcys", "A_acidic", "A_distal"):
            for st in ("mle", "profile_lo", "profile_hi"):
                g = c[(c.attribute == attr) & (c.setting == st) & (c.sizes == "lognormal_mean8.5")].sort_values("K")
                if g.empty:
                    continue
                r0 = g.iloc[0]
                rec = {"attribute": lab[attr], "setting": st, "SDs (label, attribute)": f"{r0.su:.3f}, {r0.sv:.3f}",
                       "rho": f"{r0.rho:.4f}", "pooled target log2 OR": f"{r0.theta_pooled_log2:.4f}",
                       "datasets per cell": int(r0.n_done)}
                for _, x in g.iterrows():
                    rec[f"K={int(x.K)}"] = f"{x.p_excludes0:.3f} ({x.p_excludes0_mcse:.3f}) / {x.coverage_pooled_target:.3f}"
                rows.append(rec)
            go = old[old.attribute == attr].sort_values("K")
            if len(go):
                r0 = go.iloc[0]
                rec = {"attribute": lab[attr], "setting": "round 1, GH20 parameters (superseded)",
                       "SDs (label, attribute)": f"{r0.su:.3f}, {r0.sv:.3f}", "rho": f"{r0.rho:.4f}",
                       "pooled target log2 OR": f"{r0.theta_pooled_log2:.4f}", "datasets per cell": int(r0.n_datasets)}
                for _, x in go.iterrows():
                    rec[f"K={int(x.K)}"] = f"{x.p_excludes0:.3f} ({x.p_excludes0_mcse:.3f}) / {x.coverage_pooled_target:.3f}"
                rows.append(rec)
        parts.append("### R2-d. Anchored between-protein confounding at the accurate estimates and across the profile-likelihood CI of rho: "
                     "P(excludes 0) (MC SE) / coverage of the pooled target; whole cohort, log-normal cluster sizes mean 8.5\n\n"
                     + md(pd.DataFrame(rows)))
        s = c[c.sizes != "lognormal_mean8.5"]
        rows = [{"attribute": lab[x.attribute], "setting": x.setting, "rho": f"{x.rho:.4f}", "K": int(x.K),
                 "mean cysteines": f"{x.mean_n:.0f}", "pooled target log2 OR": f"{x.theta_pooled_log2:.4f}",
                 "datasets": int(x.n_done), "P(excludes 0) (MC SE)": f"{x.p_excludes0:.3f} ({x.p_excludes0_mcse:.3f})",
                 "below / above 0": f"{x.p_below0:.3f} / {x.p_above0:.3f}", "coverage": f"{x.coverage_pooled_target:.3f}"}
                for _, x in s.iterrows()]
        parts.append("### R2-e. Sensitivity: the anchoring cohort's own 766 cluster sizes (fixed design, mean 12.4, maximum 332)\n\n"
                     + md(pd.DataFrame(rows)))
        t = pd.read_csv(os.path.join(OUT, "r2_confounding_anchored_targets.csv"))
        t["pooled log2 OR"] = t.pooled_log2or_no_within_association.map(lambda v: f"{v:.4f}")
        parts.append("### R2-f. Pooled log2 OR with no within-protein association (simlib.calibrate, 80-node Gauss-Hermite on the population integrand)\n\n"
                     + md(t[["attribute", "setting", "rho", "su", "sv", "prevalence", "pooled log2 OR"]].round(4)))
    return parts


def main():
    sm = pd.read_csv(os.path.join(OUT, "sim_summary_long.csv"))
    parts = []
    names = {"pct": "percentile (paper)", "bca": "BCa", "cr0": "sandwich CR0, z", "cr1t": "sandwich CR1, t(G-1)",
             "mdt": "sandwich Mancl-DeRouen, t(G-1)"}
    # T1: type-I error (null_re), by K x method, full and sub designs, pooled over pi and size family
    for design, label in (("full", "whole cohort (baseline interval)"), ("sub", "1:1 subsample (matched-set proxy)")):
        rows = []
        for m in ("pct", "bca", "cr0", "cr1t", "mdt"):
            rec = {"method": names[m]}
            for K in (44, 164, 776, 1475):
                for scen in ("null_iid", "null_re"):
                    g = sm[(sm.design == design) & (sm.method == m) & (sm.K == K) & (sm.scenario == scen)]
                    if g.empty:
                        rec[f"K={K} {scen}"] = "n/a"
                        continue
                    n = g.n_finite.sum()
                    p = float((g.flag0 * g.n_finite).sum() / n)
                    se = np.sqrt(p * (1 - p) / n)
                    rec[f"K={K} {scen}"] = f"{p:.3f} ({se:.3f})"
            rows.append(rec)
        parts.append(f"### Type-I error, {label}: P(95% interval excludes 0) under no association (MC SE)\n\n" + md(pd.DataFrame(rows)))
    # T2: coverage at theta = 0.5 and 1.0
    for design, label in (("full", "whole cohort"), ("sub", "1:1 subsample")):
        rows = []
        for m in ("pct", "bca", "cr0", "cr1t", "mdt"):
            rec = {"method": names[m]}
            for scen in ("alt05", "alt10"):
                for K in (44, 164, 776, 1475):
                    g = sm[(sm.design == design) & (sm.method == m) & (sm.K == K) & (sm.scenario == scen)]
                    if g.empty:
                        rec[f"{scen} K={K}"] = "n/a"
                        continue
                    n = g.n_finite.sum()
                    p = float((g.cover * g.n_finite).sum() / n)
                    rec[f"{scen} K={K}"] = f"{p:.3f}"
            rows.append(rec)
        parts.append(f"### Coverage of the pooled log2 OR 0.5 (alt05) and 1.0 (alt10), {label}\n\n" + md(pd.DataFrame(rows)))
    # T3: calibration multiplier w0 of the percentile interval by design and K (null_re; range over cells)
    rows = []
    for design in ("full", "full_mh", "sub", "rc"):
        rec = {"design": design}
        for K in (44, 164, 776, 1475):
            g = sm[(sm.design == design) & (sm.method == "pct") & (sm.K == K) & (sm.scenario.isin(["null_iid", "null_re"]))]
            if g.empty:
                rec[f"K={K}"] = "n/a"
                continue
            t1 = float((g.flag0 * g.n_finite).sum() / g.n_finite.sum())
            rec[f"K={K}"] = f"type-I {t1:.3f}; w0 {g.w0.mean():.3f} [{g.w0.min():.3f}-{g.w0.max():.3f}]"
        rows.append(rec)
    parts.append("### Percentile interval by design: pooled type-I error and calibration multiplier w0 (mean [range] over the 8 null_iid/null_re cells)\n\n" + md(pd.DataFrame(rows)))
    # T4: between-protein confounding
    rows = []
    for scen in ("null_re", "null_conf15", "null_conf50"):
        for design in ("full", "sub"):
            rec = {"scenario": scen, "design": design}
            for K in (44, 164, 776, 1475):
                g = sm[(sm.design == design) & (sm.method == "pct") & (sm.K == K) & (sm.scenario == scen)]
                if g.empty:
                    rec[f"K={K}"] = "n/a"
                    continue
                rec[f"K={K}"] = (f"{float((g.flag0 * g.n_finite).sum() / g.n_finite.sum()):.3f} / "
                                 f"{float((g.cover * g.n_finite).sum() / g.n_finite.sum()):.3f}")
            th = sm[(sm.scenario == scen)].theta.unique()
            rec["pooled target log2 OR"] = ", ".join(f"{t:.3f}" for t in sorted(th))
            rows.append(rec)
    parts.append("### Between-protein confounding (no within-protein association): P(excludes 0) / coverage of the pooled target, percentile\n\n" + md(pd.DataFrame(rows)))
    # T5: size family and prevalence detail for K=44 percentile
    rows = []
    for scen in ("null_iid", "null_re"):
        for K in (44, 164):
            for pix in (0.15, 0.35):
                for size in ("nb", "lognormal"):
                    g = sm[(sm.design == "full") & (sm.method == "pct") & (sm.K == K) & (sm.scenario == scen)
                           & (sm.pi_x == pix) & (sm["size"] == size)]
                    if g.empty:
                        continue
                    r = g.iloc[0]
                    rows.append({"scenario": scen, "K": K, "pi_x": pix, "size": size, "datasets": int(r.n_datasets),
                                 "type-I": f"{r.flag0:.3f} ({r.flag0_mcse:.3f})",
                                 "miss below / above": f"{r.miss_lo:.3f} / {r.miss_hi:.3f}",
                                 "w0 [95% CI]": f"{r.w0:.3f} [{r.w0_ci_lo:.3f}, {r.w0_ci_hi:.3f}]",
                                 "mean width": f"{r.mean_width:.2f}"})
    parts.append("### Percentile, whole cohort, per cell (few-cluster detail)\n\n" + md(pd.DataFrame(rows)))
    # claim-matched
    p = os.path.join(OUT, "claims_matched_null_summary.csv")
    if os.path.exists(p):
        cm = pd.read_csv(p)
        cm["type-I (MC SE)"] = [f"{a:.3f} ({b:.3f})" for a, b in zip(cm.type1, cm.type1_mcse)]
        cm["w0 [95% CI]"] = [f"{a:.3f} [{b:.3f}, {c:.3f}]" for a, b, c in zip(cm.w0, cm.w0_ci_lo, cm.w0_ci_hi)]
        parts.append("### Claim-matched null simulations (2,000 datasets each)\n\n" +
                     md(cm[["claim", "variant", "design", "n_datasets", "type-I (MC SE)", "w0 [95% CI]"]]))
    p = os.path.join(OUT, "calibrated_intervals.csv")
    if os.path.exists(p):
        ci = pd.read_csv(p)
        d = ci[ci.decisive & (ci.wstar > 0.5) & (ci.wstar < 1.6)].copy()
        d["interval"] = [f"{a:.3f} [{b:.3f}, {c:.3f}]" for a, b, c in zip(d.point, d.lo, d.hi)]
        d["w*"] = d.wstar.map(lambda v: f"{v:.3f}")
        d["w0 (range)"] = [f"{a:.3f} ({b:.3f}-{c:.3f})" for a, b, c in zip(d.w0, d.w0_lo, d.w0_hi)]
        d["flips at w0"] = d.flips_at_w0.map({True: "yes", False: "no"})
        d["could flip in range"] = d.could_flip_in_range.map({True: "yes", False: "no"})
        d["z_MC"] = d.z_mc.map(lambda v: f"{v:.1f}")
        parts.append("### Decisive intervals with a bound near 0 (0.5 < w* < 1.6)\n\n" +
                     md(d[["claim_id", "type", "interval", "w*", "w0 (range)", "flips at w0", "could flip in range",
                           "verdict_if_toggled", "z_MC", "w0_source"]]))
        vd = pd.read_csv(os.path.join(OUT, "calibrated_verdicts.csv"))
        parts.append("### Verdicts re-derived with calibrated intervals\n\n" +
                     md(vd[vd.could_flip_in_range][["claim_id", "unit", "clusters", "stored_verdict", "calibrated_verdict_w0",
                                                     "verdict_w0_low_end", "verdict_w0_high_end"]]))
    p = os.path.join(OUT, "multiplicity_verdicts.csv")
    if os.path.exists(p):
        mv = pd.read_csv(p)
        ch = mv[(mv.drop(columns=["claim_id", "unit", "stored_verdict", "stored_95"]).ne(mv.stored_verdict, axis=0)).any(axis=1)]
        parts.append("### Verdicts under multiplicity-adjusted intervals (normal approximation)\n\n" + md(ch))
        ms = json.load(open(os.path.join(OUT, "multiplicity_summary.json")))
        parts.append("```\n" + json.dumps({k: ms[k] for k in ("levels", "expected_chance_exclusions_global_null", "observed")},
                                          indent=1) + "\n```")
    parts += revision_tables()
    parts += revision2_tables()
    with open(os.path.join(OUT, "report_tables.md"), "w", encoding="utf-8") as fh:
        fh.write("\n\n".join(parts) + "\n")
    print("\n\n".join(parts))


if __name__ == "__main__":
    main()
