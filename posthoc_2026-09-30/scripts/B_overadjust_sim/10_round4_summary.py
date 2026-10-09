# -*- coding: utf-8 -*-
"""B_overadjust_sim step 10 (POST HOC, revision round 4, 2026-09-30): summaries of the round-4 runs.

Round 4 answers the fourth adversarial verification, which showed that a DETECTABILITY-ONLY process (selection of
positives on tryptic-peptide length among detected cysteines, no chemistry; verify_r4_08_lensel_calibrated.py)
reproduces the observed SFE-006 baseline and propensity AUC, and that VIS10 matching then removes MORE than the
observed loss. That model is re-implemented in this item's engine (detection mode 'lensel7') and run here with the
item's own seeds and the full re-test analysis (all attributes, observed-negative and detection-model estimators,
verdicts with intervals).

Writes (results/B_overadjust_sim/):
  lensel_grid_summary.csv          mean baselines along the gamma grid (L0grid), every analysed attribute
  lensel_calibration.json          the gamma reproducing the observed SFE-006 baseline (rule of 02_run_simulation)
  lensel_summary.csv               the calibrated cells L0cal (rate 0.20), L0cal16 (rate 0.16), L0ci (intervals):
                                   baseline, matched, loss, retention band, propensity AUC, observed-negative and
                                   detection-model estimators, count/offset ratios, verdict frequencies (L0ci)
  lensel_observed_comparison.csv   the observed re-tests against the calibrated length-selection model
  lensel_verifier_check.csv        this item's length-selection numbers against the round-4 verifier's
  kr5_calibration.csv              chemistry on the K/R count within +/-5: extended grid, exact b, calibrated cell,
                                   and the verifier's values
  negb_detmatch_comparison.csv     observed-negative background and detection-model matching: detectability-only
                                   models (real; length selection at the observed baseline) against K/R-count
                                   chemistry at the observed baseline
  discriminator_overlap.csv        EXPLORATORY: the count/offset association ratio for every generating model,
                                   including the detectability-only model at the observed baseline
  discriminator_overlap_groups.csv the union ranges per group of models and which ranges overlap
"""
from __future__ import annotations

import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

R = cb.RESULTS
V4 = os.path.join(R, "verify_r4")
OBS_A1 = 1.3840344206722544          # Source_Data_text_phase2c_sfe006_reproduction.json (checked below)
COUNTS = {"x1_KRcount20_hi": "K/R within +/-20", "x3_AKRV10_hi": "A/K/R/V within +/-10", "x4_KR5_hi": "K/R within +/-5"}
OBS_ATTR = {"a1_SFE006_KR": "SFE-006 (Yang 2014, own data)", "a1b_SNO021_K": "SNO-021", "a2b_PERS009_AKRV": "PERS-009",
            "a2_SNO016_KRH": "SNO-016 (transfer, HeLa)", "a3_SNO006_DE3": "SNO-006", "a3b_SFE002_E": "SFE-002"}


def band(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return np.nan, np.nan, np.nan
    return float(np.median(v)), float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))


def observed():
    j = json.load(open(cb.SFE006_JSON, encoding="utf-8"))
    assert abs(float(j["baseline_log2_or"]) - OBS_A1) < 1e-12
    o = pd.read_csv(os.path.join(R, "observed_claims.csv"))
    return o.set_index("claim")


def first_crossing(x, y, target):
    for i in range(len(x) - 1):
        if (y[i] - target) * (y[i + 1] - target) <= 0 and y[i + 1] != y[i]:
            return float(x[i] + (target - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))
    return None


def per_table(d, scen):
    """Wide per-table frame for one scenario: one row per replicate, columns '<col>|<attribute>'."""
    q = d[d.scenario == scen]
    cols = [c for c in ("baseline", "matched", "neg_b_oracle", "detmatch_crude", "pair_conditional", "matched_lo",
                        "matched_hi", "baseline_lo", "baseline_hi", "verdict_final_guarded") if c in q.columns]
    w = q.pivot(index="rep", columns="attribute", values=cols)
    w.columns = ["%s|%s" % (a, b) for a, b in w.columns]
    t = q.groupby("rep")[["n_pos", "prop_auc_before", "prop_auc_after", "gamma", "rate"]].first()
    return w.join(t)


def summarise_cell(w, attrs, obs, label):
    rows = []
    for a in attrs:
        if "baseline|" + a not in w.columns:
            continue
        b = w["baseline|" + a].astype(float)
        m = w["matched|" + a].astype(float)
        n = int(b.notna().sum())
        rec = {"cell": label, "attribute": a, "n_tables": n, "gamma": float(w.gamma.iloc[0]),
               "rate": float(w.rate.iloc[0]), "n_pos_mean": float(w.n_pos.mean()),
               "prop_auc_before_mean": float(w.prop_auc_before.mean()),
               "prop_auc_after_mean": float(w.prop_auc_after.mean()),
               "baseline_mean": float(b.mean()), "baseline_mcse": float(b.std() / np.sqrt(n)),
               "matched_mean": float(m.mean()), "matched_mcse": float(m.std() / np.sqrt(n)),
               "loss_mean": float((b - m).mean())}
        rec["retention_p50"], rec["retention_p2_5"], rec["retention_p97_5"] = band(m / b)
        for col, key in (("neg_b_oracle", "neg_b"), ("detmatch_crude", "detmatch"), ("pair_conditional", "pairc")):
            if col + "|" + a in w.columns:
                v = w[col + "|" + a].astype(float)
                rec[key + "_mean"] = float(v.mean())
                rec[key + "_over_baseline"] = float(v.mean() / b.mean())
        claim = OBS_ATTR.get(a)
        if claim is not None and claim in obs.index:
            ob = obs.loc[claim]
            rec["observed_baseline"] = float(ob.baseline)
            rec["observed_matched"] = float(ob.matched)
            rec["observed_retention"] = float(ob.retention_matched_over_baseline)
            rec["frac_tables_retention_ge_observed"] = float((m / b >= ob.retention_matched_over_baseline).mean())
            rec["frac_tables_matched_ge_observed"] = float((m >= ob.matched).mean())
        if a == "a1_SFE006_KR":
            for x in COUNTS:
                if "baseline|" + x in w.columns:
                    rec[x + "_baseline_mean"] = float(w["baseline|" + x].astype(float).mean())
                    rec[x + "_ratio_p50"], rec[x + "_ratio_p2_5"], rec[x + "_ratio_p97_5"] = band(
                        w["baseline|" + x].astype(float) / b)
        if "verdict_final_guarded|" + a in w.columns:
            vv = w["verdict_final_guarded|" + a].astype(str)
            for k, f in vv.value_counts(normalize=True).items():
                rec["verdict_" + k] = float(f)
            lo, hi = w["matched_lo|" + a].astype(float), w["matched_hi|" + a].astype(float)
            rec["frac_matched_ci_excludes_0"] = float(((lo > 0) | (hi < 0)).mean())
            blo, bhi = w["baseline_lo|" + a].astype(float), w["baseline_hi|" + a].astype(float)
            rec["frac_baseline_ci_excludes_0"] = float(((blo > 0) | (bhi < 0)).mean())
        rows.append(rec)
    return rows


def main():
    obs = observed()
    # ---------------- length-selection grid and calibration ----------------
    g = pd.read_csv(os.path.join(R, "reps_lenselgrid.csv.gz"), low_memory=False)
    gs = (g.groupby(["gamma", "attribute"])
          .agg(n_tables=("rep", "size"), baseline_mean=("baseline", "mean"), baseline_sd=("baseline", "std"),
               matched_mean=("matched", "mean"), prop_auc_before_mean=("prop_auc_before", "mean"))
          .reset_index())
    gs["baseline_mcse"] = gs.baseline_sd / np.sqrt(gs.n_tables)
    gs.drop(columns="baseline_sd").to_csv(os.path.join(R, "lensel_grid_summary.csv"), index=False)
    q = gs[gs.attribute == "a1_SFE006_KR"].sort_values("gamma")
    gstar = first_crossing(q.gamma.to_numpy(), q.baseline_mean.to_numpy(), OBS_A1)
    cal = pd.read_csv(os.path.join(R, "reps_lenselcal.csv.gz"), low_memory=False)
    ci = pd.read_csv(os.path.join(R, "reps_lenselci.csv.gz"), low_memory=False)
    used = sorted(set(np.round(cal.gamma.unique(), 12)) | set(np.round(ci.gamma.unique(), 12)))
    assert len(used) == 1 and abs(used[0] - gstar) < 1e-9, (used, gstar)
    vcal = pd.read_csv(os.path.join(V4, "lensel_calib_summary.csv"))
    json.dump({"label": "POST HOC revision round 4; rule fixed before the lenselcal and lenselci runs "
                        "(02_run_simulation.exact_lensel_gamma)",
               "model": "detectability only: among detected cysteines P(positive) = expit(a + gamma * -ln max(pep_len,"
                        " 7)), no chemistry (re-implementation of verify_r4_08_lensel_calibrated.py)",
               "observed_sfe006_baseline": OBS_A1, "gamma_star": gstar,
               "grid": q[["gamma", "baseline_mean", "baseline_mcse"]].round(6).values.tolist(),
               "gamma_used_in_lenselcal_and_lenselci": used[0],
               "verifier_gamma_star": float(vcal.gamma.iloc[0])},
              open(os.path.join(R, "lensel_calibration.json"), "w", encoding="utf-8"), indent=1)

    # ---------------- calibrated cells ----------------
    attrs = ["a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV", "a3_SNO006_DE3", "a3b_SFE002_E",
             "a4_LATENT30", "x1_KRcount20_hi", "x2_KR5_ge3", "x3_AKRV10_hi", "x4_KR5_hi"]
    rows = []
    for d, scen in ((cal, "L0cal"), (cal, "L0cal16"), (ci, "L0ci")):
        rows += summarise_cell(per_table(d, scen), attrs, obs, scen)
    ls = pd.DataFrame(rows)
    ls.to_csv(os.path.join(R, "lensel_summary.csv"), index=False)

    # ---------------- observed re-tests against the calibrated length-selection model ----------------
    rows = []
    l0 = ls[ls.cell == "L0cal"].set_index("attribute")
    for a, claim in OBS_ATTR.items():
        if claim not in obs.index or a not in l0.index:
            continue
        ob, r = obs.loc[claim], l0.loc[a]
        rows.append({"claim": claim, "attribute": a, "observed_baseline": ob.baseline, "observed_matched": ob.matched,
                     "observed_retention": ob.retention_matched_over_baseline,
                     "observed_propensity_auc": ob.propensity_auc,
                     "lensel_baseline_mean": r.baseline_mean, "lensel_matched_mean": r.matched_mean,
                     "lensel_retention_p50": r.retention_p50, "lensel_retention_p2_5": r.retention_p2_5,
                     "lensel_retention_p97_5": r.retention_p97_5, "lensel_prop_auc_mean": r.prop_auc_before_mean,
                     "observed_retention_position": ("above the band" if ob.retention_matched_over_baseline >
                                                     r.retention_p97_5 else "below the band"
                                                     if ob.retention_matched_over_baseline < r.retention_p2_5
                                                     else "within the band")})
    pd.DataFrame(rows).to_csv(os.path.join(R, "lensel_observed_comparison.csv"), index=False)

    # ---------------- verifier cross-check (length selection) ----------------
    rows = []
    for vcell, cell in (("lensel7_calib_r20", "L0cal"), ("lensel7_calib_r16", "L0cal16")):
        v = vcal[vcal.cell == vcell].iloc[0]
        r = ls[(ls.cell == cell) & (ls.attribute == "a1_SFE006_KR")].iloc[0]
        rx = ls[(ls.cell == cell) & (ls.attribute == "a3b_SFE002_E")].iloc[0]
        for name, vv, ii in (("gamma", v.gamma, r.gamma), ("n tables", v.n, r.n_tables),
                             ("positives per table", v.n_pos_mean, r.n_pos_mean),
                             ("SFE-006 baseline (mean)", v.a1_base_mean, r.baseline_mean),
                             ("SFE-006 matched (mean)", v.a1_match_mean, r.matched_mean),
                             ("SFE-006 loss (mean)", v.a1_loss_mean, r.loss_mean),
                             ("propensity AUC (mean)", v.prop_auc_mean, r.prop_auc_before_mean),
                             ("SFE-006 retention p50", v.a1_ret_p50, r.retention_p50),
                             ("SFE-006 retention p2.5", v.a1_ret_lo, r.retention_p2_5),
                             ("SFE-006 retention p97.5", v.a1_ret_hi, r.retention_p97_5),
                             ("tables with retention >= 0.43", v["frac_ret_ge_0.43"],
                              r.frac_tables_retention_ge_observed),
                             ("K/R +/-20 attribute baseline (mean)", v.x1_base_mean, r.x1_KRcount20_hi_baseline_mean),
                             ("K/R +/-20 / offsets p50", v.x1_ratio_p50, r.x1_KRcount20_hi_ratio_p50),
                             ("K/R +/-20 / offsets p2.5", v.x1_ratio_lo, r.x1_KRcount20_hi_ratio_p2_5),
                             ("K/R +/-20 / offsets p97.5", v.x1_ratio_hi, r.x1_KRcount20_hi_ratio_p97_5),
                             ("SFE-002 Glu attribute baseline (mean)", v.a3b_base_mean, rx.baseline_mean)):
            rows.append({"cell": cell, "verifier_cell": vcell, "quantity": name, "verifier": float(vv),
                         "this_item": float(ii)})
    pd.DataFrame(rows).to_csv(os.path.join(R, "lensel_verifier_check.csv"), index=False)

    # ---------------- observed-negative background and detection-model matching ----------------
    disc = pd.concat([pd.read_csv(os.path.join(R, f), low_memory=False)
                      for f in ("reps_discrim.csv.gz", "reps_discrim2.csv.gz")], ignore_index=True)
    vnd = pd.read_csv(os.path.join(V4, "negb_detmatch_summary.csv"))
    vmap = {"D0real": "real null (implementer's D0real model)", "L0cal": "lensel7_calib (artifact, no chemistry)",
            "DM2real": "M2 real, b 2.711"}
    rows = []
    for scen, lab, src in (("D0real", "detectability only, observed identifications (no selection)", disc),
                           ("D0steep", "detectability only, steeper selection on peptide detectability", disc),
                           ("L0cal", "detectability only, selection on tryptic-peptide length (observed baseline)",
                            cal),
                           ("DM2real", "chemistry on the K/R count within +/-20 (observed baseline)", disc)):
        qq = src[(src.scenario == scen) & (src.attribute == "a1_SFE006_KR")]
        rec = {"cell": scen, "model": lab, "n_tables": len(qq), "baseline_mean": qq.baseline.mean(),
               "matched_mean": qq.matched.mean(), "neg_b_mean": qq.neg_b_oracle.mean(),
               "detmatch_mean": qq.detmatch_crude.mean()}
        rec["matched_over_baseline"] = rec["matched_mean"] / rec["baseline_mean"]
        rec["neg_b_over_baseline"] = rec["neg_b_mean"] / rec["baseline_mean"]
        rec["detmatch_over_baseline"] = rec["detmatch_mean"] / rec["baseline_mean"]
        if scen in vmap:
            v = vnd[vnd.cell == vmap[scen]].iloc[0]
            rec.update(verifier_baseline=v.base, verifier_neg_b=v.negb, verifier_detmatch=v.detmatch,
                       verifier_neg_b_over_baseline=v.negb_over_base,
                       verifier_detmatch_over_baseline=v.detmatch_over_base)
        rows.append(rec)
    pd.DataFrame(rows).to_csv(os.path.join(R, "negb_detmatch_comparison.csv"), index=False)

    # ---------------- K/R +/-5 chemistry: extended grid and calibrated cell ----------------
    parts = []
    for f in ("reps_proxysens.csv.gz", "reps_proxysensext.csv.gz", "reps_kr5ext.csv.gz"):
        x = pd.read_csv(os.path.join(R, f), low_memory=False,
                        usecols=["scenario", "attribute", "effect_log2_or", "baseline"])
        parts.append(x[(x.scenario == "P1kr5grid") & (x.attribute == "a1_SFE006_KR")])
    kg = pd.concat(parts).groupby("effect_log2_or")["baseline"].agg(["mean", "std", "size"]).reset_index()
    bstar = first_crossing(kg.effect_log2_or.to_numpy(), kg["mean"].to_numpy(), OBS_A1)
    kc = pd.read_csv(os.path.join(R, "reps_kr5cal.csv.gz"), low_memory=False)
    assert abs(float(kc.effect_log2_or.iloc[0]) - bstar) < 1e-9
    w = per_table(kc, "DM2kr5")
    r = summarise_cell(w, ["a1_SFE006_KR"], obs, "DM2kr5")[0]
    vcp = json.load(open(os.path.join(V4, "calib_points.json"), encoding="utf-8"))["points"]["M2kr5|real"]
    vcs = pd.read_csv(os.path.join(V4, "cells_summary.csv")).set_index("cell").loc["DM2kr5_6.06"]
    rows = [{"row": "grid", "b": float(e), "sfe006_baseline_mean": float(m), "mcse": float(s / np.sqrt(n)),
             "n_tables": int(n)} for e, m, s, n in kg[["effect_log2_or", "mean", "std", "size"]].values]
    rows.append({"row": "calibrated (this item)", "b": bstar, "sfe006_baseline_mean": r["baseline_mean"],
                 "n_tables": r["n_tables"], "retention_p50": r["retention_p50"],
                 "retention_p2_5": r["retention_p2_5"], "retention_p97_5": r["retention_p97_5"],
                 "prop_auc_mean": r["prop_auc_before_mean"], "x1_ratio_p50": r["x1_KRcount20_hi_ratio_p50"],
                 "x1_ratio_p2_5": r["x1_KRcount20_hi_ratio_p2_5"], "x1_ratio_p97_5": r["x1_KRcount20_hi_ratio_p97_5"],
                 "x4_ratio_p50": r["x4_KR5_hi_ratio_p50"], "x4_ratio_p2_5": r["x4_KR5_hi_ratio_p2_5"],
                 "x4_ratio_p97_5": r["x4_KR5_hi_ratio_p97_5"],
                 "frac_tables_retention_ge_observed": r["frac_tables_retention_ge_observed"]})
    rows.append({"row": "calibrated (round-4 verifier)", "b": float(vcp["exact_strength"]),
                 "sfe006_baseline_mean": float(vcs.a1_base_mean), "n_tables": int(vcs.n),
                 "retention_p50": float(vcs.a1_ret_p50), "retention_p2_5": float(vcs.a1_ret_lo),
                 "retention_p97_5": float(vcs.a1_ret_hi), "prop_auc_mean": float(vcs.prop_auc_mean),
                 "x1_ratio_p50": float(vcs.x1_ratio_p50), "x1_ratio_p2_5": float(vcs.x1_ratio_lo),
                 "x1_ratio_p97_5": float(vcs.x1_ratio_hi), "x4_ratio_p50": float(vcs.x4_ratio_p50),
                 "x4_ratio_p2_5": float(vcs.x4_ratio_lo), "x4_ratio_p97_5": float(vcs.x4_ratio_hi)})
    pd.DataFrame(rows).to_csv(os.path.join(R, "kr5_calibration.csv"), index=False)

    # ---------------- EXPLORATORY count/offset ratio, every generating model, overlaps stated ----------------
    ds = pd.read_csv(os.path.join(R, "discriminator_summary.csv"))
    dm = pd.read_csv(os.path.join(R, "discriminator_mixture.csv"))
    group_of = {"D0real": "detectability only, below the observed baseline",
                "D0steep": "detectability only, below the observed baseline",
                "D0steep4": "detectability only, below the observed baseline",
                "D0steepvis": "detectability only, below the observed baseline",
                "D0steepnc": "detectability only, below the observed baseline",
                "D0r16": "detectability only, below the observed baseline",
                "DM1real": "chemistry on the SFE-006 offsets", "DM1steep": "chemistry on the SFE-006 offsets",
                "DM1perm": "chemistry on the SFE-006 offsets", "DM1r16": "chemistry on the SFE-006 offsets",
                "DM2real": "chemistry on the K/R count within +/-20", "DM2steep": "chemistry on the K/R count within +/-20",
                "DM2perm": "chemistry on the K/R count within +/-20", "DM2r16": "chemistry on the K/R count within +/-20",
                "DM2akrv10": "chemistry on the A/K/R/V count within +/-10"}
    rows = []
    for _, r in ds.iterrows():
        if r.scenario not in group_of:          # D0perm: SFE-006 baseline 0.09, ratio not interpretable
            continue
        rec = {"cell": r.scenario, "source": "this item (round 3)", "group": group_of[r.scenario],
               "detection": r.detection, "n_tables": int(r.n_tables), "a1_baseline_mean": r.a1_baseline_mean,
               "a1_retention_p50": r.a1_retention_p50}
        for x in COUNTS:
            for s in ("p50", "p2_5", "p97_5"):
                rec["%s_ratio_%s" % (x, s)] = r["%s_ratio_%s" % (x, s)]
        rows.append(rec)
    for _, r in dm[dm.scenario == "DMIX"].iterrows():
        rec = {"cell": "DMIX b=%.1f" % r.effect_log2_or, "source": "this item (round 3)",
               "group": "mixture: chemistry on the K/R +/-20 count and on the offsets", "detection": "real",
               "n_tables": int(r.n_tables), "a1_baseline_mean": r.a1_baseline_mean,
               "a1_retention_p50": r.a1_retention_p50}
        for x in COUNTS:
            for s in ("p50", "p2_5", "p97_5"):
                rec["%s_ratio_%s" % (x, s)] = r["%s_ratio_%s" % (x, s)]
        rows.append(rec)
    for cell, grp, det in (("L0cal", "detectability only, at the observed baseline (length selection)", "real, rate 0.20"),
                           ("L0cal16", "detectability only, at the observed baseline (length selection)", "real, rate 0.16")):
        r = ls[(ls.cell == cell) & (ls.attribute == "a1_SFE006_KR")].iloc[0]
        rec = {"cell": cell, "source": "this item (round 4)", "group": grp, "detection": det,
               "n_tables": int(r.n_tables), "a1_baseline_mean": r.baseline_mean, "a1_retention_p50": r.retention_p50}
        for x in COUNTS:
            for s in ("p50", "p2_5", "p97_5"):
                rec["%s_ratio_%s" % (x, s)] = r["%s_ratio_%s" % (x, s)]
        rows.append(rec)
    kr = pd.read_csv(os.path.join(R, "kr5_calibration.csv")).set_index("row").loc["calibrated (this item)"]
    wk = per_table(kc, "DM2kr5")
    rec = {"cell": "DM2kr5", "source": "this item (round 4)", "group": "chemistry on the K/R count within +/-5",
           "detection": "real", "n_tables": int(kr.n_tables), "a1_baseline_mean": kr.sfe006_baseline_mean,
           "a1_retention_p50": kr.retention_p50}
    for x in COUNTS:
        rec["%s_ratio_p50" % x], rec["%s_ratio_p2_5" % x], rec["%s_ratio_p97_5" % x] = band(
            wk["baseline|" + x].astype(float) / wk["baseline|a1_SFE006_KR"].astype(float))
    rows.append(rec)
    for _, v in vcal.iterrows():
        rec = {"cell": "verifier " + v.cell, "source": "round-4 verifier (verify_r4)",
               "group": "detectability only, at the observed baseline (length selection)",
               "detection": "real, rate %s" % ("0.20" if v.cell.endswith("r20") else "0.16"), "n_tables": int(v.n),
               "a1_baseline_mean": v.a1_base_mean, "a1_retention_p50": v.a1_ret_p50}
        for x, k in (("x1_KRcount20_hi", "x1"), ("x3_AKRV10_hi", "x3"), ("x4_KR5_hi", "x4")):
            rec["%s_ratio_p50" % x], rec["%s_ratio_p2_5" % x], rec["%s_ratio_p97_5" % x] = (
                v[k + "_ratio_p50"], v[k + "_ratio_lo"], v[k + "_ratio_hi"])
        rows.append(rec)
    ov = pd.DataFrame(rows)
    ov.to_csv(os.path.join(R, "discriminator_overlap.csv"), index=False)
    # union ranges per group (K/R +/-20 ratio) and pairwise overlaps
    grp = (ov.groupby("group").agg(n_cells=("cell", "size"), a1_baseline_min=("a1_baseline_mean", "min"),
                                   a1_baseline_max=("a1_baseline_mean", "max"),
                                   x1_ratio_p50_min=("x1_KRcount20_hi_ratio_p50", "min"),
                                   x1_ratio_p50_max=("x1_KRcount20_hi_ratio_p50", "max"),
                                   x1_union_lo=("x1_KRcount20_hi_ratio_p2_5", "min"),
                                   x1_union_hi=("x1_KRcount20_hi_ratio_p97_5", "max")).reset_index())
    names = grp.group.tolist()
    for gname in names:
        lo, hi = grp.set_index("group").loc[gname, ["x1_union_lo", "x1_union_hi"]]
        over = []
        for other in names:
            if other == gname:
                continue
            lo2, hi2 = grp.set_index("group").loc[other, ["x1_union_lo", "x1_union_hi"]]
            if lo <= hi2 and lo2 <= hi:
                over.append(other)
        grp.loc[grp.group == gname, "x1_range_overlaps"] = "; ".join(over)
    grp.to_csv(os.path.join(R, "discriminator_overlap_groups.csv"), index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 60)
    print("gamma* =", gstar, "(verifier", float(vcal.gamma.iloc[0]), ")")
    print(ls[ls.attribute.isin(["a1_SFE006_KR", "a1b_SNO021_K", "a2b_PERS009_AKRV", "a3b_SFE002_E",
                                "x1_KRcount20_hi"])].round(3).T.to_string())
    print(pd.read_csv(os.path.join(R, "lensel_observed_comparison.csv")).round(3).to_string(index=False))
    print(pd.read_csv(os.path.join(R, "lensel_verifier_check.csv")).round(3).to_string(index=False))
    print(pd.read_csv(os.path.join(R, "negb_detmatch_comparison.csv")).round(3).to_string(index=False))
    print(pd.read_csv(os.path.join(R, "kr5_calibration.csv")).round(3).to_string(index=False))
    print(ov[["cell", "group", "a1_baseline_mean", "x1_KRcount20_hi_ratio_p50", "x1_KRcount20_hi_ratio_p2_5",
              "x1_KRcount20_hi_ratio_p97_5"]].round(2).to_string(index=False))
    print(grp.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
