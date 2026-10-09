# -*- coding: utf-8 -*-
"""B_overadjust_sim step 6 (POST HOC): print the report's markdown tables from the stored CSVs,
so that every number in reports/B_overadjust_sim.md is read from a file rather than typed.
Writes results/B_overadjust_sim/report_tables.txt."""
from __future__ import annotations

import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ATTR = ["a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV", "a3_SNO006_DE3",
        "a3b_SFE002_E", "a4_LATENT30"]
SHORT = {"a1_SFE006_KR": "a1 K/R SFE-006 offsets", "a1b_SNO021_K": "a1b K SNO-021 offsets",
         "a2_SNO016_KRH": "a2 >=3 KRH (SNO-016)", "a2b_PERS009_AKRV": "a2b >=3 AKRV (PERS-009)",
         "a3_SNO006_DE3": "a3 D/E at +3 (SNO-006)", "a3b_SFE002_E": "a3b E SFE-002 offsets",
         "a4_LATENT30": "a4 latent 30%", "a1c_KRH_p5p8": "a1c K/R/H +5..+8 (Artifact-1 analogue)",
         "x1_KRcount20_hi": "x1 K/R count +/-20 > median (count claim)",
         "x2_KR5_ge3": "x2 >=3 K/R within +/-5 (count claim)"}
PROXY_ATTR = ["a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV", "a3_SNO006_DE3",
              "a3b_SFE002_E", "a4_LATENT30", "x1_KRcount20_hi", "x2_KR5_ge3"]


def round2_tables():
    """Revision round 2 (post hoc): proxy-chemistry model and count-defined claims."""
    out = []
    kc = os.path.join(cb.RESULTS, "krcount_universe_summary.json")
    if os.path.exists(kc):
        k = json.load(open(kc, encoding="utf-8"))
        rows = []
        for c in k["columns"]:
            if c["kind"].startswith("binary"):
                rows.append({"column": c["column"], "prevalence": f(c["prevalence_all"]),
                             "detected / undetected": "%s / %s" % (f(c["prevalence_detected"]), f(c["prevalence_undetected"])),
                             "log2 OR vs detection": f(c["log2_or_vs_detection"]),
                             "VIS10 predicts (5-fold CV)": "AUC " + f(c["vis10_predicts_auc_cv5"], 3)})
            else:
                rows.append({"column": c["column"], "prevalence": "mean %s, SD %s, median %s" % (
                    f(c["mean"]), f(c["sd_population"]), f(c["median"], 0)),
                             "detected / undetected": "%s / %s" % (f(c["mean_detected"]), f(c["mean_undetected"])),
                             "log2 OR vs detection": "",
                             "VIS10 predicts (5-fold CV)": "R2 " + f(c["vis10_predicts_r2_cv5"], 3)})
        out.append("### T1b K/R-count columns (01f; the count is the repository's pep_cleavage_sites_within_20)\n\n"
                   + md(pd.DataFrame(rows)) + "\n\nchecks: VIS10 columns of the DIG25 call vs stored max diff %s; "
                   "x2 vs kr_count5>=3 mismatches %s; verifier F_count20 %s" % (
                       k["check_vis10_columns_of_DIG25_call_max_abs_diff_vs_stored"],
                       k["check_x2_equals_kr_count5_ge_3_mismatches"], json.dumps(k["check_against_verifier"])))
    pg_path = os.path.join(cb.RESULTS, "proxy_grid_summary.csv")
    if os.path.exists(pg_path):
        pg = pd.read_csv(pg_path)
        for scen, g in pg.groupby("scenario", sort=False):
            rows = []
            for b, gb in g.groupby("b_log2_or_per_sd"):
                row = {"b (log2 OR per SD)": f(b, 1), "prop AUC": f(gb["prop_auc_mean"].iloc[0], 3)}
                for a in PROXY_ATTR:
                    r = gb[gb.attribute == a]
                    if r.empty:
                        continue
                    r = r.iloc[0]
                    row[SHORT[a].split(" ")[0]] = "%s->%s; %s [%s, %s]" % (
                        f(r["baseline_mean"]), f(r["matched_mean"]), f(r["ratio_p50"]), f(r["ratio_p2_5"]),
                        f(r["ratio_p97_5"]))
                x1 = gb[gb.attribute == "x1_KRcount20_hi"]
                if len(x1):
                    row["x1 within detected"] = f(x1["realised_within_detected_mean"].iloc[0])
                rows.append(row)
            out.append("### T13 proxy chemistry %s (%s): per attribute baseline->matched; per-table matched/baseline "
                       "median [2.5, 97.5]\n\n" % (scen, g["detection_model"].iloc[0]) + md(pd.DataFrame(rows)))
    ovp_path = os.path.join(cb.RESULTS, "observed_vs_proxy.csv")
    if os.path.exists(ovp_path):
        o = pd.read_csv(ovp_path)
        rows = []
        for _, r in o.iterrows():
            band = ("%s [%s, %s] at b %s" % (f(r["band_p50"]), f(r["band_p2_5"]), f(r["band_p97_5"]),
                                             f(r["nearest_grid_b"], 1))) if np.isfinite(r.get("band_p50", np.nan)) else \
                "not reached (max baseline %s at b %s)" % (f(r["max_sim_baseline"]), f(r["b_at_max_sim_baseline"], 1))
            pos = ""
            if np.isfinite(r.get("band_p50", np.nan)):
                pos = "below" if r["observed_below_band"] else ("above" if r["observed_above_band"] else "within")
            rows.append({"claim": r["claim"], "proxy grid": r["proxy_scenario"], "model": r["detection_model"],
                         "obs baseline -> matched (retention)": "%s -> %s (%s)" % (
                             f(r["observed_baseline"]), f(r["observed_matched"]), f(r["observed_retention"])),
                         "b reproducing obs baseline": f(r["b_reproducing_observed_baseline"]),
                         "proxy band": band, "obs vs band": pos,
                         "prop AUC sim / obs": "%s / %s" % (f(r.get("sim_prop_auc_at_nearest", np.nan), 3),
                                                            f(r["observed_propensity_auc"], 3)),
                         "x1 baseline (implied count assoc.)": f(r.get("x1_baseline_at_nearest", np.nan)),
                         "excess diag. proxy / obs": "%s / %s" % (f(r.get("proxy_excess_baseline", np.nan)),
                                                                  f(r.get("observed_excess_baseline", np.nan))),
                         "retention medians over b>=1": "%s-%s" % (f(r["retention_median_min_b_ge_1"]),
                                                                  f(r["retention_median_max_b_ge_1"])),
                         "attribute-planted band": r.get("attribute_planted_band", "")})
        out.append("### T14 observed re-tests vs the proxy-chemistry calibration (chemistry on the local K/R count, "
                   "not on the claim's attribute)\n\n" + md(pd.DataFrame(rows)))
    pj = os.path.join(cb.RESULTS, "proxy_joint_at_calibration.csv")
    if os.path.exists(pj):
        j = pd.read_csv(pj)
        j = j[j.calibrated_to == "SFE-006 (Yang 2014, own data)"]
        rows = [{"proxy grid": r["proxy_scenario"], "b": f(r["b"], 1), "attribute": SHORT.get(r["attribute"], r["attribute"]),
                 "baseline -> matched": "%s -> %s" % (f(r["baseline_mean"]), f(r["matched_mean"])),
                 "matched/baseline median [2.5, 97.5]": "%s [%s, %s]" % (f(r["ratio_p50"]), f(r["ratio_p2_5"]),
                                                                        f(r["ratio_p97_5"]))}
                for _, r in j.iterrows()]
        out.append("### T14b all attributes in the proxy tables calibrated to the SFE-006 baseline (same cohort as "
                   "SFE-002)\n\n" + md(pd.DataFrame(rows)))
    pv = os.path.join(cb.RESULTS, "proxy_verdicts.csv")
    if os.path.exists(pv) and os.path.getsize(pv) > 5:
        v = pd.read_csv(pv)
        rows = []
        for _, r in v.iterrows():
            rows.append({"scenario": r["scenario"], "b": f(r["b"], 1), "attribute": SHORT.get(r["attribute"], r["attribute"]),
                         "n": int(r["n_tables"]), "prop AUC": f(r["prop_auc_mean"], 3),
                         "baseline -> matched": "%s -> %s" % (f(r["baseline_mean"]), f(r["matched_mean"])),
                         "ratio median [2.5, 97.5]": "%s [%s, %s]" % (f(r["ratio_p50"]), f(r["ratio_p2_5"]), f(r["ratio_p97_5"])),
                         "P(matched CI excl. 0)": f(r["p_matched_ci_excludes_0"]),
                         "final: attenuated (MC SE)": "%s (%s)" % (f(r["final_attenuated"]), f(r["final_attenuated_mcse"], 3)),
                         "final: survives (MC SE)": "%s (%s)" % (f(r["final_survives"]), f(r["final_survives_mcse"], 3)),
                         "final: vanishes / undecidable": "%s / %s" % (f(r["final_vanishes"]), f(r["final_undecidable"])),
                         "final: bcc / reverses": "%s / %s" % (f(r["final_baseline_contradicts_claim"]), f(r["final_reverses"]))})
        out.append("### T16 pipeline verdicts under proxy chemistry at the SFE-006-calibrated b (intervals, B = 1000)\n\n"
                   + md(pd.DataFrame(rows)))
    pa = os.path.join(cb.RESULTS, "proxy_artifact1.csv")
    if os.path.exists(pa) and os.path.getsize(pa) > 5:
        a1 = pd.read_csv(pa)
        rows = []
        for _, r in a1.iterrows():
            if r["b"] not in (0.0, 1.0, 2.0, 3.0, 4.0):
                continue
            rows.append({"grid": r["scenario"], "attribute": SHORT.get(r["attribute"], r["attribute"]), "b": f(r["b"], 1),
                         "MH all -> restricted": "%s -> %s" % (f(r["mh_all"]), f(r["mh_restricted"])),
                         "shift (MC SE)": "%s (%s)" % (f(r["mh_shift"]), f(r["mh_shift_mcse"], 3)),
                         "z all -> restricted": "%s -> %s" % (f(r["z_all"]), f(r["z_restricted"])),
                         "share of proxy increment passing": f(r.get("share_of_increment_passing_restriction", np.nan))})
        out.append("### T17 Artifact-1 restriction under proxy chemistry (protein-stratified MH log2 OR)\n\n"
                   + md(pd.DataFrame(rows)))
    cc = os.path.join(cb.RESULTS, "count_claim_summary.csv")
    if os.path.exists(cc) and os.path.getsize(cc) > 5:
        c = pd.read_csv(cc)
        rows = []
        for _, r in c.iterrows():
            rows.append({"model": r["detection_model"], "attribute": SHORT.get(r["attribute"], r["attribute"]),
                         "planted": f(r["effect_log2_or"], 1), "n": int(r["n_tables"]),
                         "S0 baseline -> matched": "%s -> %s" % (f(r["S0_baseline"]), f(r["S0_matched"])),
                         "S1 baseline -> matched": "%s -> %s" % (f(r["S1_baseline"]), f(r["S1_matched"])),
                         "A rel. latent (MC SE)": "%s (%s)" % (f(r["A_rel_latent"]), f(r["A_rel_latent_mcse"], 3)),
                         "B rel. latent (MC SE)": "%s (%s)" % (f(r["B_rel_latent"]), f(r["B_rel_latent_mcse"], 3)),
                         "loss beyond latent A / B": "%s / %s" % (f(r["loss_beyond_latent_A"]), f(r["loss_beyond_latent_B"])),
                         "matched/baseline median [2.5, 97.5]": "%s [%s, %s]" % (f(r["ratio_p50"]), f(r["ratio_p2_5"]),
                                                                                f(r["ratio_p97_5"]))})
        out.append("### T15 count-defined claims next to the published definitions: loss beyond the latent attribute "
                   "(planted on the attribute itself)\n\n" + md(pd.DataFrame(rows)))
        # compact ranges over 0.5-1.5 and both normalisations
        rows = []
        for (mode, att), g in c.groupby(["detection_model", "attribute"], sort=False):
            lo = min(g["loss_beyond_latent_A"].min(), g["loss_beyond_latent_B"].min())
            hi = max(g["loss_beyond_latent_A"].max(), g["loss_beyond_latent_B"].max())
            g1 = g[np.isclose(g.effect_log2_or, 1.0)]
            rows.append({"model": mode, "attribute": SHORT.get(att, att), "loss range 0.5-1.5, A and B": "%s-%s" % (f(lo), f(hi)),
                         "loss at 1.0, A / B": "%s / %s" % (f(g1["loss_beyond_latent_A"].iloc[0]),
                                                           f(g1["loss_beyond_latent_B"].iloc[0])) if len(g1) else ""})
        out.append("### T15b loss-beyond-latent ranges (planted 0.5-1.5, both normalisations)\n\n" + md(pd.DataFrame(rows)))
    return out


def round3_tables():
    """Revision round 3 (post hoc): the sequence-only discriminator (count vs offsets association)."""
    out = []
    p = os.path.join(cb.RESULTS, "discriminator_summary.csv")
    if not os.path.exists(p):
        return out
    ds = pd.read_csv(p)
    lab = {"x1_KRcount20_hi": "K/R +/-20", "x3_AKRV10_hi": "A/K/R/V +/-10", "x4_KR5_hi": "K/R +/-5"}
    rows = []
    for _, r in ds.iterrows():
        row = {"cell": r.scenario, "generating model": r.model, "detection": r.detection, "tables": int(r.n_tables),
               "planted": f(r.effect_log2_or), "SFE-006 baseline": f(r.a1_baseline_mean),
               "SFE-006 retention": "%s [%s, %s]" % (f(r.a1_retention_p50), f(r.a1_retention_p2_5),
                                                     f(r.a1_retention_p97_5)),
               "prop AUC": f(r.prop_auc_mean, 3)}
        for k, nm in lab.items():
            row["%s / offsets" % nm] = "%s [%s, %s]" % (f(r[k + "_ratio_p50"]), f(r[k + "_ratio_p2_5"]),
                                                        f(r[k + "_ratio_p97_5"]))
        rows.append(row)
    out.append("### T18 sequence-only discriminator: the association of the label with a local count divided by its "
               "association with the SFE-006 offsets, in the same table (baseline log2 OR; median [2.5, 97.5%])\n\n"
               + md(pd.DataFrame(rows)))
    p = os.path.join(cb.RESULTS, "discriminator_mixture.csv")
    if os.path.exists(p):
        dm = pd.read_csv(p)
        rows = []
        for _, r in dm.iterrows():
            rows.append({"cell": r.scenario, "b (K/R +/-20, per SD)": f(r.effect_log2_or),
                         "c (SFE-006 offsets)": f(r.effect2_log2_or), "tables": int(r.n_tables),
                         "SFE-006 baseline": f(r.a1_baseline_mean),
                         "SFE-006 retention": "%s [%s, %s]" % (f(r.a1_retention_p50), f(r.a1_retention_p2_5),
                                                              f(r.a1_retention_p97_5)),
                         "K/R +/-20 / offsets": "%s [%s, %s]" % (f(r.x1_KRcount20_hi_ratio_p50),
                                                                 f(r.x1_KRcount20_hi_ratio_p2_5),
                                                                 f(r.x1_KRcount20_hi_ratio_p97_5))})
        out.append("### T19 mixtures of chemistry on the K/R +/-20 count and on the SFE-006 offsets (DMIX = the c "
                   "that restores the observed baseline at each b)\n\n" + md(pd.DataFrame(rows)))
    p = os.path.join(cb.RESULTS, "discriminator_verifier_check.csv")
    if os.path.exists(p):
        c = pd.read_csv(p)
        rows = [{"cell": r.cell,
                 "verifier ratio": "%s [%s, %s]" % (f(r.verifier_ratio_p50), f(r.verifier_lo), f(r.verifier_hi)),
                 "this item": "%s [%s, %s]" % (f(r.item_ratio_p50), f(r.item_lo), f(r.item_hi)),
                 "SFE-006 baseline, verifier / item": "%s / %s" % (f(r.verifier_a1_baseline), f(r.item_a1_baseline)),
                 "count baseline, verifier / item": "%s / %s" % (f(r.verifier_x1_baseline), f(r.item_x1_baseline))}
                for _, r in c.iterrows()]
        out.append("### T20 round-3 verifier's discriminator numbers against this item's rerun (K/R +/-20 ratio)\n\n"
                   + md(pd.DataFrame(rows)))
    p = os.path.join(cb.RESULTS, "sfe006_count_test.json")
    if os.path.exists(p):
        j = json.load(open(p, encoding="utf-8"))
        # revision round 4: the reading rule is withdrawn; the script is descriptive only
        sha = j.get("reading_rule_sha256_withdrawn", j.get("reading_rule_sha256"))
        rows = [{"field": "status", "value": j.get("status")},
                {"field": "label", "value": j.get("label")},
                {"field": "round-3 reading rule sha256 (first 16; withdrawn in round 4)", "value": str(sha)[:16]},
                {"field": "reading", "value": j.get("reading") if isinstance(j.get("reading"), str) else "see json"}]
        if j.get("status") == "not_computed_inputs_missing":
            for m in j.get("missing", []):
                rows.append({"field": "missing input", "value": m})
        else:
            res = j.get("result", {})
            for k, v in res.get("baseline_log2_or", {}).items():
                rows.append({"field": "baseline log2 OR " + k, "value": f(v, 3)})
            for k, v in res.get("ratio_to_a1", {}).items():
                ci = ((res.get("ratio_to_a1_bootstrap", {}) or {}).get(k, {}) or {}).get("ci", [np.nan, np.nan])
                rows.append({"field": "ratio " + k, "value": "%s [%s, %s]" % (f(v), f(ci[0]), f(ci[1]))})
        out.append("### T21 the count/offset ratio on the authors' SFE-006 cohort (09_sfe006_count_test.py; "
                   "descriptive only since round 4)\n\n" + md(pd.DataFrame(rows)))
    return out


def round4_tables():
    """Revision round 4 (post hoc): the detectability-only length-selection model at the observed baseline, the
    K/R +/-5 calibration, the observed-negative and detection-model estimators, and the exploratory ratio overlaps."""
    out = []
    R = cb.RESULTS
    p = os.path.join(R, "lensel_calibration.json")
    if not os.path.exists(p):
        return out
    cal = json.load(open(p, encoding="utf-8"))
    gs = pd.read_csv(os.path.join(R, "lensel_grid_summary.csv"))
    g1 = gs[gs.attribute == "a1_SFE006_KR"].sort_values("gamma")
    rows = [{"gamma": f(r.gamma, 1), "tables": int(r.n_tables),
             "SFE-006 baseline, mean (MC SE)": "%s (%s)" % (f(r.baseline_mean, 3), f(r.baseline_mcse, 3)),
             "SFE-006 matched, mean": f(r.matched_mean, 3), "propensity AUC": f(r.prop_auc_before_mean, 3)}
            for _, r in g1.iterrows()]
    out.append("### T22 detectability only, selection on tryptic-peptide length (no chemistry): calibration grid "
               "(L0grid). gamma reproducing the observed SFE-006 baseline 1.3840: %s (round-4 verifier %s)\n\n"
               % (f(cal["gamma_star"], 4), f(cal["verifier_gamma_star"], 4)) + md(pd.DataFrame(rows)))
    ls = pd.read_csv(os.path.join(R, "lensel_summary.csv"))
    rows = []
    for _, r in ls.iterrows():
        if r.cell == "L0ci" or (r.cell == "L0cal16" and r.attribute != "a1_SFE006_KR"):
            continue
        rows.append({"cell": r.cell, "attribute": SHORT.get(r.attribute, r.attribute), "tables": int(r.n_tables),
                     "positives": f(r.n_pos_mean, 0), "prop AUC": f(r.prop_auc_before_mean, 3),
                     "baseline": f(r.baseline_mean, 3), "matched": f(r.matched_mean, 3), "loss": f(r.loss_mean, 3),
                     "retention [2.5, 97.5%]": "%s [%s, %s]" % (f(r.retention_p50), f(r.retention_p2_5),
                                                             f(r.retention_p97_5)),
                     "NEG_B / baseline": f(r.get("neg_b_over_baseline", np.nan)),
                     "detection-model matched / baseline": f(r.get("detmatch_over_baseline", np.nan)),
                     "observed retention": f(r.get("observed_retention", np.nan)),
                     "tables with retention >= observed": f(r.get("frac_tables_retention_ge_observed", np.nan), 3)})
    out.append("### T23 the length-selection model at the observed SFE-006 baseline (L0cal, rate 0.20, 300 tables; "
               "L0cal16, rate 0.16, SFE-006 row only)\n\n" + md(pd.DataFrame(rows)))
    rows = []
    for _, r in ls[ls.cell == "L0ci"].iterrows():
        vc = {k.replace("verdict_", ""): r[k] for k in ls.columns if k.startswith("verdict_") and np.isfinite(r[k])}
        rows.append({"attribute": SHORT.get(r.attribute, r.attribute), "tables": int(r.n_tables),
                     "baseline CI excludes 0": f(r.frac_baseline_ci_excludes_0, 3),
                     "matched CI excludes 0": f(r.frac_matched_ci_excludes_0, 3),
                     "final verdict (guarded)": ", ".join("%s %s" % (k, f(v, 3)) for k, v in
                                                          sorted(vc.items(), key=lambda kv: -kv[1]))})
    out.append("### T24 verdicts under the length-selection model (L0ci, 200 tables with intervals, B = 1,000)\n\n"
               + md(pd.DataFrame(rows)))
    v = pd.read_csv(os.path.join(R, "lensel_verifier_check.csv"))
    rows = [{"cell": r.cell, "quantity": r.quantity, "round-4 verifier": f(r.verifier, 3), "this item": f(r.this_item, 3)}
            for _, r in v.iterrows()]
    out.append("### T25 the length-selection model: round-4 verifier (verify_r4_08) against this item's rerun\n\n"
               + md(pd.DataFrame(rows)))
    n = pd.read_csv(os.path.join(R, "negb_detmatch_comparison.csv"))
    rows = [{"cell": r.cell, "model": r.model, "tables": int(r.n_tables), "baseline": f(r.baseline_mean, 3),
             "VIS10-matched / baseline": f(r.matched_over_baseline, 3),
             "observed-negative (NEG_B) / baseline": f(r.neg_b_over_baseline, 3),
             "detection-model matched / baseline": f(r.detmatch_over_baseline, 3),
             "verifier NEG_B / baseline": f(r.get("verifier_neg_b_over_baseline", np.nan), 3),
             "verifier detection-model / baseline": f(r.get("verifier_detmatch_over_baseline", np.nan), 3)}
            for _, r in n.iterrows()]
    out.append("### T26 the SFE-006 attribute against an observed-negative background and under detection-model "
               "matching (exploratory estimators)\n\n" + md(pd.DataFrame(rows)))
    k = pd.read_csv(os.path.join(R, "kr5_calibration.csv"))
    rows = []
    for _, r in k.iterrows():
        row = {"row": r["row"], "b (log2 OR per SD)": f(r.b, 3), "SFE-006 baseline": f(r.sfe006_baseline_mean, 3),
               "tables": int(r.n_tables)}
        if np.isfinite(r.get("retention_p50", np.nan)):
            row["retention [2.5, 97.5%]"] = "%s [%s, %s]" % (f(r.retention_p50), f(r.retention_p2_5),
                                                             f(r.retention_p97_5))
            row["prop AUC"] = f(r.prop_auc_mean, 3)
            row["K/R +/-20 / offsets"] = "%s [%s, %s]" % (f(r.x1_ratio_p50), f(r.x1_ratio_p2_5), f(r.x1_ratio_p97_5))
        rows.append(row)
    out.append("### T27 chemistry on the K/R count within +/-5 (real detection): grid extended to b = 8 and the "
               "calibrated cell\n\n" + md(pd.DataFrame(rows).fillna("")))
    g = pd.read_csv(os.path.join(R, "discriminator_overlap_groups.csv"))
    rows = [{"group of generating models": r.group, "cells": int(r.n_cells),
             "SFE-006 baselines": "%s-%s" % (f(r.a1_baseline_min), f(r.a1_baseline_max)),
             "K/R +/-20 / offsets, medians": "%s-%s" % (f(r.x1_ratio_p50_min), f(r.x1_ratio_p50_max)),
             "union of 2.5-97.5% ranges": "[%s, %s]" % (f(r.x1_union_lo), f(r.x1_union_hi)),
             "overlaps": r.x1_range_overlaps} for _, r in g.iterrows()]
    out.append("### T28 EXPLORATORY count/offset ratio (K/R +/-20): union ranges per group of models and which "
               "overlap (discriminator_overlap_groups.csv)\n\n" + md(pd.DataFrame(rows)))
    return out


def f(x, d=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return ("%." + str(d) + "f") % x


def ms(r, e, d=2):
    """mean (MC SE)"""
    return "%s (%s)" % (f(r[e + "_mean"], d), f(r[e + "_mcse"], 3))


def rng(r, e, d=2):
    return "[%s, %s]" % (f(r[e + "_p2_5"], d), f(r[e + "_p97_5"], d))


def md(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out)


def pick(s, scen, attr, eff=None):
    q = s[(s.scenario == scen) & (s.attribute == attr)]
    if eff is not None:
        q = q[np.isclose(q.effect_log2_or, eff)]
    return q.iloc[0] if len(q) else None


def main():
    s = pd.read_csv(os.path.join(cb.RESULTS, "sim_summary.csv"))
    prep = json.load(open(cb.PREP_SUMMARY, encoding="utf-8"))
    out = []
    # T1 descriptive
    rows = []
    for a in prep["attributes"]:
        rows.append({"attribute": SHORT[a["attribute"]], "prevalence (all Cys)": f(a["prevalence_all"]),
                     "detected": f(a["prevalence_detected"]), "undetected": f(a["prevalence_undetected"]),
                     "log2 OR attribute vs detection": f(a["log2_or_attribute_vs_detection"]),
                     "VIS10 predicts attribute, AUC (5-fold CV)": f(a["vis10_predicts_attribute_auc_cv5"], 3)})
    out.append("### T1 universe and attributes\n\n" + md(pd.DataFrame(rows)))

    # T2 null chemistry
    for scen, title in (("S0", "real detection"), ("S0steep", "steep selection model (gamma 2)"),
                        ("S0perm", "detection permuted within protein"),
                        ("S0steepvis", "revision: selection on the VIS10 detection logit (gamma 0.25)"),
                        ("S0steepnc", "revision: selection on the composition-free peptide model (gamma 3)")):
        rows = []
        for a in ATTR:
            r = pick(s, scen, a)
            if r is None:
                continue
            row = {"attribute": SHORT[a], "n": int(r["n_reps"]), "baseline": ms(r, "baseline"),
                   "matched": ms(r, "matched"), "matched 2.5-97.5%": rng(r, "matched"),
                   "stratified": ms(r, "stratified"), "random control": ms(r, "random_control"),
                   "NEG_B oracle": ms(r, "neg_b_oracle"), "pair-conditional": ms(r, "pair_conditional"),
                   "share of baseline removed": f(r["artifact_removed_fraction"])}
            if r["n_reps_with_ci"] > 0:
                row["P(matched CI excl. 0)"] = "%s (%s)" % (f(r["frac_matched_ci_excludes_0"]),
                                                          f(r["frac_matched_ci_excludes_0_mcse"], 3))
                row["P(baseline CI excl. 0)"] = f(r["frac_baseline_ci_excludes_0"])
            rows.append(row)
        if rows:
            out.append("### T2 null chemistry, %s (mean (MC SE))\n\n" % title + md(pd.DataFrame(rows)))

    # T3 verdicts under null (interval cells)
    for scen in ("S0", "S0steep", "S0steepvis", "S0steepnc"):
        rows = []
        for a in ATTR:
            r = pick(s, scen, a)
            if r is None or r["n_reps_with_ci"] == 0:
                continue
            row = {"attribute": SHORT[a], "n": int(r["n_reps_with_ci"])}
            for v in ("survives", "attenuated", "vanishes", "undecidable", "reverses", "baseline_contradicts_claim"):
                row["final: " + v] = f(r["verdict_final_guarded_" + v])
            row["primary: vanishes"] = f(r["verdict_primary_vanishes"])
            rows.append(row)
        if rows:
            out.append("### T3 verdict frequencies under null chemistry, %s (final verdict with the phase-2b "
                       "direction guard)\n\n" % scen + md(pd.DataFrame(rows)))

    # T4 planted chemistry
    for scen, title in (("S1", "real detection"), ("S1steep", "steep selection model (gamma 2)")):
        rows = []
        for a in ATTR:
            for e in (0.5, 1.0, 1.5):
                r = pick(s, scen, a, e)
                if r is None:
                    continue
                row = {"attribute": SHORT[a], "true log2 OR": f(e, 1), "n": int(r["n_reps"]),
                       "realised within detected": f(r["realised_within_detected_mean"]),
                       "baseline": ms(r, "baseline"), "matched": ms(r, "matched"),
                       "matched / truth": "%s (%s)" % (f(r["matched_retention_vs_truth"]),
                                                        f(r["matched_retention_vs_truth_mcse"], 3)),
                       "matched / baseline": "%s (%s)" % (f(r["matched_over_baseline_ratio_of_means"]),
                                                           f(r.get("matched_over_baseline_ratio_mcse", np.nan), 3)),
                       "P(matched/baseline < 0.5)": f(r.get("frac_matched_over_baseline_lt_0_5", np.nan)),
                       "NEG_B oracle / truth": f(r["neg_b_oracle_retention_vs_truth"])}
                if r["n_reps_with_ci"] > 0:
                    row["P(attenuated)"] = f(r["verdict_final_guarded_attenuated"])
                    row["P(survives)"] = f(r["verdict_final_guarded_survives"])
                    row["P(vanishes)"] = f(r["verdict_final_guarded_vanishes"])
                    row["P(undecidable)"] = f(r["verdict_final_guarded_undecidable"])
                    row["P(bcc)"] = f(r["verdict_final_guarded_baseline_contradicts_claim"])
                rows.append(row)
        if rows:
            out.append("### T4 planted chemistry, %s\n\n" % title + md(pd.DataFrame(rows)))

    # T5 over-adjustment decomposition
    oa_path = os.path.join(cb.RESULTS, "overadjustment_summary.csv")
    if os.path.exists(oa_path):
        oa = pd.read_csv(oa_path)
        for (mode, scen), g in oa.groupby(["detection_model", "planted_scenario"]):
            if scen not in ("S1", "S1steep", "S1perm", "S1art1", "S1art1steep", "S1steepvis", "S1steepnc"):
                continue
            rows = []
            for a in ATTR:
                for _, r in g[g.attribute == a].sort_values("effect_log2_or").iterrows():
                    row = {"attribute": SHORT[a], "true": f(r["effect_log2_or"], 1)}
                    for est in ("matched", "pair_conditional", "stratified", "detmatch_crude",
                                "a1_wp_mh_restricted"):
                        c = est + "_genuine_part_retention"
                        if c in r and np.isfinite(r[c]):
                            row[est + " genuine-part"] = "%s (%s)" % (f(r[c]), f(r[est + "_genuine_part_retention_mcse"], 3))
                    c = "matched_genuine_part_retention_relative_to_latent"
                    if c in r and np.isfinite(r.get(c, np.nan)):
                        row["matched rel. latent (A: /truth)"] = "%s (%s)" % (
                            f(r[c]), f(r.get(c + "_mcse", np.nan), 3))
                    c = "matched_increment_normalised_retention"
                    if c in r and np.isfinite(r.get(c, np.nan)):
                        row["matched increment-normalised (B)"] = "%s (%s)" % (f(r[c]), f(r.get(c + "_mcse", np.nan), 3))
                        row["B rel. latent"] = "%s (%s)" % (f(r.get(c + "_relative_to_latent", np.nan)),
                                                            f(r.get(c + "_relative_to_latent_mcse", np.nan), 3))
                    rows.append(row)
            out.append("### T5 genuine-part retention: (A) (S1 - S0 matched)/truth, (B) (S1 - S0 matched)/"
                       "(S1 - S0 baseline): %s / %s\n\n" % (mode, scen) + md(pd.DataFrame(rows)))
        # compact ranges across effects for the report text
        rows = []
        for (mode, scen), g in oa.groupby(["detection_model", "planted_scenario"]):
            if scen not in ("S1", "S1steep", "S1perm", "S1steepvis", "S1steepnc"):
                continue
            for a in ATTR[:-1]:
                q = g[g.attribute == a]
                if q.empty:
                    continue
                A = q["matched_genuine_part_retention_relative_to_latent"]
                B = q["matched_increment_normalised_retention_relative_to_latent"]
                rows.append({"model": mode, "scenario": scen, "attribute": SHORT[a],
                             "effects": ",".join(f(e, 1) for e in q["effect_log2_or"]),
                             "A rel. latent range": "%s-%s" % (f(A.min()), f(A.max())),
                             "B rel. latent range": "%s-%s" % (f(B.min()), f(B.max())),
                             "loss range over both (1 - rel.)": "%s-%s" % (f(1 - max(A.max(), B.max())),
                                                                          f(1 - min(A.min(), B.min())))})
        out.append("### T5b over-adjustment ranges, both normalisations\n\n" + md(pd.DataFrame(rows)))

    # T6 rate sensitivity
    rows = []
    for rate in (0.1, 0.2, 0.4):
        for a in ATTR:
            r0 = s[(s.scenario.isin(["S0rate", "S0"])) & (s.attribute == a) & np.isclose(s.rate, rate)
                   & (s.detection_mode == "real")]
            r1 = s[(s.scenario.isin(["S1rate", "S1"])) & (s.attribute == a) & np.isclose(s.rate, rate)
                   & np.isclose(s.effect_log2_or, 1.0) & (s.detection_mode == "real")]
            if r0.empty or r1.empty:
                continue
            r0, r1 = r0.iloc[0], r1.iloc[0]
            rows.append({"rate": f(rate, 2), "attribute": SHORT[a], "cohort positive fraction": f(r1["mean_pos_rate"], 3),
                         "S0 baseline": f(r0["baseline_mean"]), "S0 matched": f(r0["matched_mean"]),
                         "S1 baseline": f(r1["baseline_mean"]), "S1 matched": f(r1["matched_mean"]),
                         "S1 matched/truth": f(r1["matched_retention_vs_truth"]),
                         "S1 matched/baseline": f(r1["matched_over_baseline_ratio_of_means"])})
    if rows:
        out.append("### T6 positive-rate sensitivity (planted 1.0 for S1)\n\n" + md(pd.DataFrame(rows)))

    # T7 observed vs simulated
    ovs_path = os.path.join(cb.RESULTS, "observed_vs_simulated.csv")
    if os.path.exists(ovs_path):
        o = pd.read_csv(ovs_path)
        rows = []
        for _, r in o.iterrows():
            rows.append({"claim": r["claim"], "model": r["detection_model"],
                         "obs baseline": f(r["observed_baseline"]), "obs matched": f(r["observed_matched"]),
                         "obs retention": f(r["observed_retention"]),
                         "sim null baseline": f(r["sim_null_baseline_mean"]),
                         "sim null matched": f(r["sim_null_matched_mean"]),
                         "effect reproducing obs baseline": f(r["effect_reproducing_observed_baseline"]),
                         "sim retention at that effect, median [2.5, 97.5]": "%s [%s, %s]" % (
                             f(r.get("sim_retention_p50_at_nearest", np.nan)),
                             f(r.get("sim_retention_p2_5_at_nearest", np.nan)),
                             f(r.get("sim_retention_p97_5_at_nearest", np.nan))),
                         "nearest grid": f(r.get("nearest_grid_effect", np.nan)),
                         "obs below 2.5%": str(r.get("observed_retention_below_sim_p2_5", "")),
                         "effect reproducing obs matched": f(r["effect_reproducing_observed_matched"]),
                         "sim baseline at that effect": f(r.get("sim_baseline_at_matched_effect", np.nan)),
                         "excess obs baseline": f(r.get("excess_baseline_over_sim", np.nan))})
        out.append("### T7 observed re-tests vs simulated calibration curves (real, steep: planted effect + the "
                   "model's own detectability association; permuted: whole baseline planted)\n\n" + md(pd.DataFrame(rows)))

    # T8 Artifact-1 control (revision: artifact1_shift_summary.csv, all selection models)
    a1p = os.path.join(cb.RESULTS, "artifact1_shift_summary.csv")
    if os.path.exists(a1p):
        a1s = pd.read_csv(a1p)
        rows = []
        for _, r in a1s.iterrows():
            if r["attribute"] not in ("a1_SFE006_KR", "a1c_KRH_p5p8", "a1b_SNO021_K", "a2_SNO016_KRH",
                                      "a3_SNO006_DE3", "a3b_SFE002_E", "a4_LATENT30"):
                continue
            rows.append({"model": "%s (gamma %s)" % (r["detection_model"], f(r["gamma"], 2)),
                         "attribute": SHORT.get(r["attribute"], r["attribute"]), "planted": f(r["effect_log2_or"], 1),
                         "n": int(r["n_tables"]), "prop AUC": f(r["mean_prop_auc_before"], 3),
                         "z all -> restricted": "%s -> %s" % (f(r["z_all"]), f(r["z_restricted"])),
                         "z ratio": f(r["z_ratio_of_means"]),
                         "MH all -> restricted": "%s -> %s" % (f(r["mh_all"], 3), f(r["mh_restricted"], 3)),
                         "MH shift (MC SE)": "%s (%s)" % (f(r["mh_shift"], 3), f(r["mh_shift_mcse"], 3)),
                         "share of planted MH increment passing": f(r.get("share_of_planted_increment_passing_restriction", np.nan))})
        out.append("### T8 Artifact-1 control (within-protein z / protein-stratified MH log2 OR), all selection "
                   "models\n\n" + md(pd.DataFrame(rows)))

    # T10 released tool vs pipeline (revision)
    tvp = os.path.join(cb.RESULTS, "tool_vs_pipeline.csv")
    if os.path.exists(tvp):
        t = pd.read_csv(tvp)
        rows = []
        for _, r in t.iterrows():
            rows.append({"scenario": r["scenario"], "planted": f(r["effect_log2_or"], 1),
                         "attribute": SHORT.get(r["attribute"], r["attribute"]), "n": int(r["n_tables"]),
                         "tool survives/attenuated": "%d (%s)" % (r["tool_survives_or_attenuated_n"],
                                                                 f(r["tool_survives_or_attenuated_frac"])),
                         "pipeline final survives/attenuated": "%d (%s)" % (r["pipeline_final_survives_or_attenuated_n"],
                                                                           f(r["pipeline_final_survives_or_attenuated_frac"])),
                         "tool CI excl. 0": f(r["tool_ci_excludes_0"]),
                         "stratified CI excl. 0": f(r["pipeline_stratified_ci_excludes_0"]),
                         "matched CI excl. 0": f(r["pipeline_matched_ci_excludes_0"]),
                         "tool = pipeline final": f(r["tool_equals_pipeline_final"]),
                         "tool = pipeline secondary": f(r["tool_equals_pipeline_secondary_guarded"]),
                         "tool verdicts": r["tool_verdict_counts"], "pipeline final verdicts": r["pipeline_final_counts"]})
        out.append("### T10 released tool (v0.2.2) vs pipeline final verdict, same tables\n\n" + md(pd.DataFrame(rows)))

    # T11 null interval exclusion vs latent (revision)
    nel = os.path.join(cb.RESULTS, "null_exclusion_vs_latent.csv")
    if os.path.exists(nel):
        t = pd.read_csv(nel)
        rows = [{"scenario": r["scenario"], "attribute": SHORT.get(r["attribute"], r["attribute"]), "n": int(r["n_tables"]),
                 "P(matched CI excl. 0) (MC SE)": "%s (%s)" % (f(r["p_matched_ci_excludes_0"], 3), f(r["mcse"], 3)),
                 "paired diff vs latent (MC SE)": "%s (%s)" % (f(r["paired_diff_vs_latent"], 3), f(r["paired_diff_mcse"], 3))}
                for _, r in t.iterrows()]
        out.append("### T11 null chemistry: matched-interval exclusion of 0 and paired difference from the latent "
                   "attribute\n\n" + md(pd.DataFrame(rows)))

    # T12 alternative selection models (revision): planted 1.0 and verdicts at null
    rows = []
    for s0, s1, mode in (("S0steepvis", "S1steepvis", "steep_vis"), ("S0steepnc", "S1steepnc", "steep_nocomp"),
                         ("S0steep", "S1steep", "steep"), ("S0", "S1", "real")):
        for a in ATTR:
            r0 = pick(s, s0, a)
            r1 = pick(s, s1, a, 1.0)
            if r0 is None:
                continue
            row = {"model": mode, "attribute": SHORT[a], "prop AUC": f(r0["mean_prop_auc_before"], 3),
                   "null baseline -> matched": "%s -> %s" % (f(r0["baseline_mean"]), f(r0["matched_mean"])),
                   "null P(matched CI excl. 0)": f(r0.get("frac_matched_ci_excludes_0", np.nan)),
                   "null final survives": f(r0.get("verdict_final_guarded_survives", np.nan)),
                   "null final attenuated": f(r0.get("verdict_final_guarded_attenuated", np.nan)),
                   "null final bcc": f(r0.get("verdict_final_guarded_baseline_contradicts_claim", np.nan))}
            if r1 is not None:
                row["planted 1.0 matched/truth"] = f(r1["matched_retention_vs_truth"])
            rows.append(row)
    if rows:
        out.append("### T12 selection models compared (null verdicts with intervals; planted 1.0 point estimates)\n\n"
                   + md(pd.DataFrame(rows)))

    # revision round 2 (post hoc)
    out.extend(round2_tables())
    # revision round 3 (post hoc)
    out.extend(round3_tables())
    # revision round 4 (post hoc)
    out.extend(round4_tables())

    # T9 checks
    for name in ("boot_check.csv", "tool_check.csv"):
        p = os.path.join(cb.RESULTS, name)
        if os.path.exists(p) and os.path.getsize(p) > 5:
            d = pd.read_csv(p)
            out.append("### %s\n\n" % name + md(d.round(4)))
    text = "\n\n".join(out) + "\n"
    with open(os.path.join(cb.RESULTS, "report_tables.txt"), "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
