# -*- coding: utf-8 -*-
"""B_overadjust_sim step 3 (POST HOC): summarise the replicate tables.

Outputs (results/B_overadjust_sim/):
  sim_summary.csv              one row per cell (scenario x detection model x rate x attribute x
                               planted effect): mean / sd / Monte Carlo SE / percentiles of every
                               estimator; retention = mean(estimator) / truth (S1 cells); matched /
                               baseline ratio of means; interval-exclusion and verdict frequencies
                               with Monte Carlo SE where intervals were computed
  verdict_summary.csv          verdict frequencies (pipeline primary, final, final with the
                               phase-2b direction guard; released-tool verdict) per interval cell
  overadjustment_summary.csv   per attribute and detection model: null residual (S0 matched),
                               retention of a planted effect, and the genuine-part retention
                               (S1 matched - S0 matched) / effect, relative to the latent attribute
  observed_claims.csv          observed site-level VIS10 re-tests read from Supplemental Data 3-6
                               and the SFE-006 reproduction Source Data
  observed_vs_simulated.csv    each observed K/R or acidic claim against the simulated calibration
                               curves of its own attribute
  boot_check.csv               verdict agreement, B = 1000 vs B = 5000 on the same 50 S0 tables
  tool_check.csv               released tool (v0.2.2) vs pipeline secondary specification

Added in revision round 2 (post hoc):
  proxy_grid_summary.csv       proxy-chemistry grids (chemistry on the local K/R count): per model, b and
                               analysed attribute, mean baseline / matched (MC SE), per-table matched/baseline
                               percentiles, propensity AUC and the implied association of the label with the
                               count-defined attribute x1
  observed_vs_proxy.csv        each observed claim against the proxy-chemistry grids: the b whose mean baseline
                               of the claim's attribute equals the observed baseline, the per-table retention band
                               there, and the implementer's 'excess baseline' diagnostic computed on the proxy
                               tables
  proxy_joint_at_calibration.csv  every analysed attribute in the proxy tables at the b calibrated to one claim
                               (e.g. the SFE-002 Glu attribute in the tables calibrated to SFE-006, same cohort)
  proxy_verdicts.csv           verdict frequencies of the pipeline under proxy chemistry (interval tables)
  proxy_artifact1.csv          Artifact-1 restriction under proxy chemistry: MH all -> restricted, and the share of
                               the proxy-induced MH increment that passes the restriction
  count_claim_summary.csv      the two hypothetical count-defined claims: null residual, planted retention and
                               loss beyond the latent attribute (both normalisations, MC SE), per detection model

Added in the revision after adversarial verification (post hoc):
  overadjustment_summary.csv   now also carries the INCREMENT-normalised genuine-part retention,
                               (matched_S1 - matched_S0) / (baseline_S1 - baseline_S0), i.e. the share
                               of the chemistry present in the baseline that survives matching, with
                               a delta-method MC SE that uses the within-table covariance of matched
                               and baseline, both absolute and relative to the latent attribute
  observed_vs_simulated.csv    third calibration grid: detection permuted within protein
                               (S1permgrid), where the whole baseline is planted chemistry; columns
                               renamed (sim_retention_at_calibrated_effect) because under the real and
                               steep models the calibrated baseline also contains the model's own
                               detectability association
  tool_vs_pipeline.csv         released tool (stratified specification, permuted-stratum random
                               control) against the pipeline's final verdict on the same tables
  null_exclusion_vs_latent.csv P(matched interval excludes 0) under null chemistry per attribute and
                               its paired difference from the latent attribute on the same tables
  artifact1_shift_summary.csv  Artifact-1 restriction: z and MH before/after, z ratio, MH shift and the
                               share of a planted MH increment that passes the restriction, per
                               selection model
"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

EST = ["baseline", "matched", "stratified", "random_control", "neg_b_oracle", "pair_conditional",
       "detmatch_crude", "detmatch_pair", "realised_within_detected", "tool_controlled",
       "a1_wp_z_all", "a1_wp_z_restricted", "a1_wp_mh_all", "a1_wp_mh_restricted", "a1_restricted_crude",
       "a1_mh_shift", "a1_z_shift"]
CELL = ["scenario", "detection_mode", "gamma", "rate", "attribute", "effect_log2_or"]
ATTR_ORDER = ["a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV",
              "a3_SNO006_DE3", "a3b_SFE002_E", "a4_LATENT30", "a1c_KRH_p5p8",
              "x1_KRcount20_hi", "x2_KR5_ge3"]
ATTR_CLASS = {"a1_SFE006_KR": "basic", "a1b_SNO021_K": "basic", "a2_SNO016_KRH": "basic",
              "a2b_PERS009_AKRV": "basic", "a3_SNO006_DE3": "acidic", "a3b_SFE002_E": "acidic",
              "a4_LATENT30": "latent", "a1c_KRH_p5p8": "basic",
              "x1_KRcount20_hi": "basic_count", "x2_KR5_ge3": "basic_count"}
# revision round 2: proxy-chemistry grids (chemistry on a continuous K/R-count feature)
PROXY_GRIDS = [("P1grid", "real", "F_KRcount20_z"), ("P1steepgrid", "steep", "F_KRcount20_z"),
               ("P1permgrid", "permuted_within_protein", "F_KRcount20_z"),
               ("P1kr5grid", "real", "F_KR5_z"), ("P1akrv10grid", "real", "F_AKRV10_z")]
ATTR_GRID_OF_MODE = {"real": "S1grid", "steep": "S1steepgrid", "permuted_within_protein": "S1permgrid"}
VERDICTS = ["survives", "attenuated", "vanishes", "undecidable", "reverses", "baseline_contradicts_claim"]


def load_reps():
    frames = []
    for f in sorted(glob.glob(os.path.join(cb.RESULTS, "reps_*.csv.gz"))):
        d = pd.read_csv(f)
        d["source_file"] = os.path.basename(f)
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)
    if "gamma" not in df.columns:
        df["gamma"] = 0.0
    df["gamma"] = df["gamma"].fillna(0.0)
    if "a1_wp_mh_all" in df.columns:  # per-table shift produced by the Artifact-1 restriction
        df["a1_mh_shift"] = df["a1_wp_mh_restricted"] - df["a1_wp_mh_all"]
        df["a1_z_shift"] = df["a1_wp_z_restricted"] - df["a1_wp_z_all"]
    return df


def prop_se(p, n):
    return float(np.sqrt(p * (1 - p) / n)) if n > 0 else np.nan


def summarise_cell(g):
    row = {"n_reps": int(len(g)), "n_reps_with_ci": int(g["with_ci"].sum()),
           "truth_log2_or": float(g["truth_log2_or"].iloc[0]),
           "mean_n_obs": g["n_obs"].mean(), "mean_n_pos": g["n_pos"].mean(),
           "mean_n_proteins": g["n_proteins"].mean(), "mean_pos_rate": g["pos_rate"].mean(),
           "mean_n_pairs": g["n_pairs"].mean(), "mean_prop_auc_before": g["prop_auc_before"].mean(),
           "mean_prop_auc_after": g["prop_auc_after"].mean(),
           "frac_precheck_blocked": g["precheck_blocked"].astype(bool).mean(),
           "mean_matched_neg_detected_frac": g["matched_neg_detected_frac"].mean(),
           "mean_all_neg_detected_frac": g["all_neg_detected_frac"].mean(),
           "mean_attr_prev_pos": g["attr_prev_pos"].mean(), "mean_attr_prev_neg": g["attr_prev_neg"].mean(),
           "mean_attr_prev_matched_neg": g["attr_prev_matched_neg"].mean()}
    n = len(g)
    truth = row["truth_log2_or"]
    for e in EST:
        if e not in g.columns:
            continue
        v = g[e].astype(float).dropna()
        if len(v) == 0:
            continue
        row[e + "_n"] = int(len(v))
        row[e + "_mean"] = v.mean()
        row[e + "_sd"] = v.std(ddof=1)
        row[e + "_mcse"] = v.std(ddof=1) / np.sqrt(len(v))
        for q in (2.5, 25, 50, 75, 97.5):
            row["%s_p%s" % (e, str(q).replace(".", "_"))] = float(np.percentile(v, q))
        if truth != 0:
            row[e + "_retention_vs_truth"] = v.mean() / truth
            row[e + "_retention_vs_truth_mcse"] = v.std(ddof=1) / np.sqrt(len(v)) / abs(truth)
    # matched / baseline, the quantity the manuscript calls retention
    b, m = g["baseline"].astype(float), g["matched"].astype(float)
    row["matched_over_baseline_ratio_of_means"] = m.mean() / b.mean() if b.mean() != 0 else np.nan
    # within-cell covariance of the matched and baseline estimates (revision: delta-method MC SEs)
    ok = np.isfinite(m) & np.isfinite(b)
    if ok.sum() > 2:
        row["matched_baseline_cov"] = float(np.cov(np.vstack([m[ok], b[ok]]), ddof=1)[0, 1])
    # delta-method MC SE of a ratio of means (paired replicates)
    if b.mean() != 0 and n > 2:
        rb, rm = b.mean(), m.mean()
        cov = np.cov(np.vstack([m, b]), ddof=1)
        var = (cov[0, 0] / rb ** 2 - 2 * rm * cov[0, 1] / rb ** 3 + rm ** 2 * cov[1, 1] / rb ** 4) / n
        row["matched_over_baseline_ratio_mcse"] = float(np.sqrt(max(var, 0)))
    ratio = (m / b).replace([np.inf, -np.inf], np.nan).dropna()
    if len(ratio):
        for q in (2.5, 50, 97.5):
            row["matched_over_baseline_p%s" % str(q).replace(".", "_")] = float(np.percentile(ratio, q))
        row["frac_matched_over_baseline_lt_0_5"] = float((ratio.abs() < 0.5).mean())
    row["artifact_removed_fraction"] = (1 - m.mean() / b.mean()) if (truth == 0 and b.mean() != 0) else np.nan
    # interval-based quantities
    ci = g[g["with_ci"].astype(bool)]
    k = len(ci)
    if k:
        for nm, lo, hi in (("baseline", "baseline_lo", "baseline_hi"), ("matched", "matched_lo", "matched_hi"),
                           ("stratified", "stratified_lo", "stratified_hi"), ("random", "random_lo", "random_hi")):
            exc = ((ci[lo] > 0) | (ci[hi] < 0)).astype(float)
            pos = (ci[lo] > 0).astype(float)
            neg = (ci[hi] < 0).astype(float)
            row["frac_%s_ci_excludes_0" % nm] = exc.mean()
            row["frac_%s_ci_excludes_0_mcse" % nm] = prop_se(exc.mean(), k)
            row["frac_%s_ci_above_0" % nm] = pos.mean()
            row["frac_%s_ci_below_0" % nm] = neg.mean()
            if truth != 0:
                cover = ((ci[lo] <= truth) & (ci[hi] >= truth)).astype(float)
                row["frac_%s_ci_covers_truth" % nm] = cover.mean()
        for col in ("verdict_primary", "verdict_final", "verdict_final_guarded"):
            for v in VERDICTS:
                p = float((ci[col] == v).mean())
                row["%s_%s" % (col, v)] = p
                row["%s_%s_mcse" % (col, v)] = prop_se(p, k)
        if "tool_verdict" in ci.columns:
            t = ci["tool_verdict"].dropna()
            row["n_tool"] = int(len(t))
            for v in VERDICTS:
                row["tool_verdict_%s" % v] = float((t == v).mean()) if len(t) else np.nan
    return row


def build_summary(df):
    rows = []
    for key, g in df.groupby(CELL, dropna=False):
        r = dict(zip(CELL, key))
        r["attribute_class"] = ATTR_CLASS.get(r["attribute"], "")
        r.update(summarise_cell(g))
        rows.append(r)
    out = pd.DataFrame(rows)
    out["attr_order"] = out["attribute"].map({a: i for i, a in enumerate(ATTR_ORDER)})
    out = out.sort_values(["scenario", "detection_mode", "gamma", "rate", "attr_order", "effect_log2_or"])
    return out.drop(columns="attr_order").reset_index(drop=True)


def verdict_table(df):
    ci = df[df["with_ci"].astype(bool)]
    rows = []
    for key, g in ci.groupby(CELL, dropna=False):
        base = dict(zip(CELL, key))
        k = len(g)
        for col in ("verdict_primary", "verdict_final", "verdict_final_guarded", "tool_verdict"):
            if col not in g.columns:
                continue
            gg = g[col].dropna()
            if not len(gg):
                continue
            for v in VERDICTS:
                p = float((gg == v).mean())
                rows.append(dict(base, verdict_column=col, verdict=v, n=int(len(gg)), fraction=p,
                                 mcse=prop_se(p, len(gg))))
    return pd.DataFrame(rows)


def _increment_ratio(r0, r1):
    """(m1 - m0) / (b1 - b0) with a delta-method MC SE; S0 and S1 are independent tables, while m and b
    of the same table covary (matched_baseline_cov / n)."""
    n0, n1 = r0["matched_n"], r1["matched_n"]
    N = r1["matched_mean"] - r0["matched_mean"]
    D = r1["baseline_mean"] - r0["baseline_mean"]
    if not np.isfinite(D) or abs(D) < 1e-9:
        return np.nan, np.nan
    vN = r1["matched_mcse"] ** 2 + r0["matched_mcse"] ** 2
    vD = r1["baseline_mcse"] ** 2 + r0["baseline_mcse"] ** 2
    cND = r1.get("matched_baseline_cov", np.nan) / n1 + r0.get("matched_baseline_cov", np.nan) / n0
    if not np.isfinite(cND):
        cND = 0.0
    var = vN / D ** 2 + N ** 2 * vD / D ** 4 - 2 * N * cND / D ** 3
    return N / D, float(np.sqrt(max(var, 0.0)))


PAIRS = [("real", "S0", "S1", 0.0), ("real", "S0", "S1grid", 0.0), ("steep", "S0steep", "S1steep", 2.0),
         ("steep", "S0steep", "S1steepgrid", 2.0), ("permuted_within_protein", "S0perm", "S1perm", 0.0),
         ("permuted_within_protein", "S0perm", "S1permgrid", 0.0),
         ("steep_vis", "S0steepvis", "S1steepvis", None), ("steep_nocomp", "S0steepnc", "S1steepnc", None),
         ("real", "S0art1", "S1art1", 0.0), ("steep", "S0art1steep", "S1art1steep", 2.0),
         ("steep", "S0art1steep4", "S1art1steep4", 4.0)]


def overadjustment(summary):
    """Genuine-part retention, two normalisations (revision):
       truth-normalised      (mean matched S1 at effect e - mean matched S0) / e
       increment-normalised  (mean matched S1 - mean matched S0) / (mean baseline S1 - mean baseline S0)
    each also relative to the latent attribute in the same detection model / scenario / effect."""
    rows = []
    for mode, s0, s1, gam in PAIRS:
        q0 = (summary.scenario == s0) & (summary.detection_mode == mode) & (summary.rate == 0.2)
        q1 = (summary.scenario == s1) & (summary.detection_mode == mode) & (summary.rate == 0.2)
        if gam is not None:
            q0 &= np.isclose(summary.gamma, gam)
            q1 &= np.isclose(summary.gamma, gam)
        s0rows, s1rows = summary[q0], summary[q1]
        if s0rows.empty or s1rows.empty:
            continue
        for a in ATTR_ORDER:
            r0 = s0rows[s0rows.attribute == a]
            if r0.empty:
                continue
            r0 = r0.iloc[0]
            for _, r1 in s1rows[(s1rows.attribute == a) & (s1rows.effect_log2_or > 0)].iterrows():
                e = r1["effect_log2_or"]
                out = {"detection_model": mode, "gamma": r1["gamma"], "null_scenario": s0, "planted_scenario": s1,
                       "attribute": a, "attribute_class": ATTR_CLASS[a], "effect_log2_or": e}
                for est in ("baseline", "matched", "stratified", "pair_conditional", "detmatch_crude",
                            "detmatch_pair", "neg_b_oracle", "a1_wp_mh_all", "a1_wp_mh_restricted",
                            "a1_restricted_crude"):
                    if est + "_mean" not in r1 or pd.isna(r1.get(est + "_mean")) or pd.isna(r0.get(est + "_mean")):
                        continue
                    out[est + "_S0_mean"] = r0[est + "_mean"]
                    out[est + "_S1_mean"] = r1[est + "_mean"]
                    out[est + "_retention_vs_truth"] = r1[est + "_mean"] / e
                    out[est + "_genuine_part_retention"] = (r1[est + "_mean"] - r0[est + "_mean"]) / e
                    out[est + "_genuine_part_retention_mcse"] = np.sqrt(r1[est + "_mcse"] ** 2 + r0[est + "_mcse"] ** 2) / e
                inc, inc_se = _increment_ratio(r0, r1)
                out["baseline_increment"] = r1["baseline_mean"] - r0["baseline_mean"]
                out["baseline_increment_over_truth"] = (r1["baseline_mean"] - r0["baseline_mean"]) / e
                out["matched_increment_normalised_retention"] = inc
                out["matched_increment_normalised_retention_mcse"] = inc_se
                out["matched_over_baseline_S1"] = r1["matched_over_baseline_ratio_of_means"]
                out["matched_over_baseline_S0"] = r0["matched_over_baseline_ratio_of_means"]
                rows.append(out)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # relative to the latent attribute in the same detection model / scenario / effect
    key = ["detection_model", "planted_scenario", "gamma", "effect_log2_or"]
    lat = df[df.attribute == "a4_LATENT30"].set_index(key)
    for col in [c for c in df.columns if c.endswith("_genuine_part_retention")] + ["matched_increment_normalised_retention"]:
        ref = [lat.loc[(r.detection_model, r.planted_scenario, r.gamma, r.effect_log2_or), col]
               if (r.detection_model, r.planted_scenario, r.gamma, r.effect_log2_or) in lat.index else np.nan
               for r in df.itertuples()]
        df[col + "_relative_to_latent"] = df[col] / np.asarray(ref, dtype=float)
    # MC SE of the relative figures (independent attribute and latent cells: S1 cells differ, S0 shared;
    # the S0 covariance is ignored, so the SE is approximate)
    for col, se in (("matched_genuine_part_retention", "matched_genuine_part_retention_mcse"),
                    ("matched_increment_normalised_retention", "matched_increment_normalised_retention_mcse")):
        ref_v = np.asarray([lat.loc[(r.detection_model, r.planted_scenario, r.gamma, r.effect_log2_or), col]
                            if (r.detection_model, r.planted_scenario, r.gamma, r.effect_log2_or) in lat.index else np.nan
                            for r in df.itertuples()], dtype=float)
        ref_s = np.asarray([lat.loc[(r.detection_model, r.planted_scenario, r.gamma, r.effect_log2_or), se]
                            if (r.detection_model, r.planted_scenario, r.gamma, r.effect_log2_or) in lat.index else np.nan
                            for r in df.itertuples()], dtype=float)
        rel = df[col] / ref_v
        df[col + "_relative_to_latent_mcse"] = np.abs(rel) * np.sqrt((df[se] / df[col]) ** 2 + (ref_s / ref_v) ** 2)
    return df


OBS_MAP = {  # observed claim -> simulated attribute with the same definition
    "SFE-006 (Yang 2014, own data)": "a1_SFE006_KR",
    "SFE-006 (transfer, Arabidopsis)": "a1_SFE006_KR",
    "SNO-021": "a1b_SNO021_K",
    "SNO-016 (transfer, HeLa)": "a2_SNO016_KRH",
    "PERS-009": "a2b_PERS009_AKRV",
    "SNO-006": "a3_SNO006_DE3",
    "SFE-002": "a3b_SFE002_E",
    "SFE-001 (RSA > 25%)": "a4_LATENT30",
}


def observed_claims():
    rows = []
    j = json.load(open(cb.SFE006_JSON, encoding="utf-8"))
    rows.append({"claim": "SFE-006 (Yang 2014, own data)", "source": "Source_Data_text_phase2c_sfe006_reproduction.json",
                 "own_data": True, "baseline": j["baseline_log2_or"], "baseline_lo": j["baseline_interval"][0],
                 "baseline_hi": j["baseline_interval"][1], "matched": j["matched_log2_or"],
                 "matched_lo": j["matched_interval"][0], "matched_hi": j["matched_interval"][1],
                 "random_control": j["random_control_log2_or"], "stratified": j["stratified_log2_or"],
                 "propensity_auc": j["propensity_auc"], "n_obs": j["n_observations"], "n_pos": j["n_positive"],
                 "verdict": j["verdict"]})
    sd4 = pd.read_csv(cb.SD4)
    names = {"SFE-006": "SFE-006 (transfer, Arabidopsis)", "SNO-021": "SNO-021", "SNO-016": "SNO-016 (transfer, HeLa)",
             "PERS-009": "PERS-009", "PERS-002": "PERS-002", "PERS-003": "PERS-003", "SFE-007": "SFE-007 (transfer)",
             "SFE-008": "SFE-008 (transfer)"}
    for _, r in sd4[sd4.unit == "site"].iterrows():
        rows.append({"claim": names.get(r.claim_id, r.claim_id), "source": "Supplemental_Data_4_retest_round_b.csv",
                     "own_data": not bool(r.is_transfer), "baseline": r.baseline_log2_or, "baseline_lo": r.baseline_ci_low,
                     "baseline_hi": r.baseline_ci_high, "matched": r.matched_log2_or, "matched_lo": r.matched_ci_low,
                     "matched_hi": r.matched_ci_high, "random_control": r.random_control_log2_or,
                     "stratified": r.stratified_log2_or, "propensity_auc": r.propensity_auc,
                     "n_obs": r.n_observations, "n_pos": r.n_positive, "verdict": r.verdict})
    sd5 = pd.read_csv(cb.SD5)
    for _, r in sd5[(sd5.unit == "site") & (sd5.specification == "primary")].iterrows():
        rows.append({"claim": r.claim_id, "source": "Supplemental_Data_5_retest_round_d.csv", "own_data": True,
                     "baseline": r.baseline_log2_or, "baseline_lo": r.baseline_ci_low, "baseline_hi": r.baseline_ci_high,
                     "matched": r.matched_log2_or, "matched_lo": r.matched_ci_low, "matched_hi": r.matched_ci_high,
                     "random_control": r.random_control_log2_or, "stratified": r.stratified_log2_or,
                     "propensity_auc": r.propensity_auc, "n_obs": r.n_observations, "n_pos": r.n_positive,
                     "verdict": r.verdict})
    sd6 = pd.read_csv(cb.SD6)
    for _, r in sd6[(sd6.unit == "site") & (sd6.specification == "primary")].iterrows():
        nm = "SFE-001 (RSA > 25%)" if r.claim_id == "SFE-001" else r.claim_id + " (structural)"
        rows.append({"claim": nm, "source": "Supplemental_Data_6_retest_round_e.csv", "own_data": True,
                     "baseline": r.baseline_log2_or, "baseline_lo": r.baseline_ci_low, "baseline_hi": r.baseline_ci_high,
                     "matched": r.matched_log2_or, "matched_lo": r.matched_ci_low, "matched_hi": r.matched_ci_high,
                     "random_control": r.random_control_log2_or, "stratified": r.stratified_log2_or,
                     "propensity_auc": r.propensity_auc, "n_obs": r.n_observations, "n_pos": r.n_positive,
                     "verdict": r.verdict})
    ob = pd.DataFrame(rows)
    ob["retention_matched_over_baseline"] = ob["matched"] / ob["baseline"]
    ob["pos_rate"] = ob["n_pos"] / ob["n_obs"]
    ob["simulated_attribute"] = ob["claim"].map(OBS_MAP)
    return ob


def interp_inverse(x_grid, y_grid, y_target):
    """Smallest x on a piecewise-linear curve with y(x) == y_target (curve assumed increasing
    where crossed); NaN if outside the simulated range."""
    for i in range(len(x_grid) - 1):
        y0, y1 = y_grid[i], y_grid[i + 1]
        if (y0 - y_target) * (y1 - y_target) <= 0 and y1 != y0:
            return x_grid[i] + (y_target - y0) * (x_grid[i + 1] - x_grid[i]) / (y1 - y0)
    return np.nan


def observed_vs_simulated(ob, reps, summary):
    """For each observed claim and detection model: the grid effect whose mean simulated baseline equals the
    observed baseline (e_b), and the per-table matched/baseline distribution at the nearest grid effect.
    Under the real and steep models the simulated baseline at e_b is the planted effect PLUS the model's own
    detectability association (the band is 'a genuine effect that, together with the simulated detectability
    association, produces the observed baseline'); under detection permuted within protein (revision) there is
    no detectability association, so the whole baseline is planted chemistry and the band is the conservative,
    pure over-adjustment bound."""
    rows = []
    grids = [("real", "S1grid"), ("steep", "S1steepgrid"), ("permuted_within_protein", "S1permgrid")]
    for _, o in ob.dropna(subset=["simulated_attribute"]).iterrows():
        a = o["simulated_attribute"]
        for mode, scen in grids:
            s = summary[(summary.scenario == scen) & (summary.attribute == a)].sort_values("effect_log2_or")
            if s.empty:
                continue
            x = s["effect_log2_or"].to_numpy()
            yb = s["baseline_mean"].to_numpy()
            ym = s["matched_mean"].to_numpy()
            e_b = interp_inverse(x, yb, o["baseline"])
            e_m = interp_inverse(x, ym, o["matched"])
            out = {"claim": o["claim"], "simulated_attribute": a, "detection_model": mode,
                   "observed_baseline": o["baseline"], "observed_matched": o["matched"],
                   "observed_retention": o["retention_matched_over_baseline"],
                   "observed_propensity_auc": o["propensity_auc"],
                   "sim_null_baseline_mean": yb[0], "sim_null_matched_mean": ym[0],
                   "sim_null_propensity_auc": s["mean_prop_auc_before"].iloc[0],
                   "effect_reproducing_observed_baseline": e_b,
                   "effect_reproducing_observed_matched": e_m}
            if np.isfinite(e_b):
                mb = np.interp(e_b, x, ym)
                out["sim_matched_at_that_effect"] = mb
                out["sim_retention_at_calibrated_effect"] = mb / o["baseline"]
                out["sim_null_detectability_share_of_observed_baseline"] = yb[0] / o["baseline"]
                # per-replicate ratio distribution at the nearest grid effect
                near = x[np.argmin(np.abs(x - e_b))]
                rr = reps[(reps.scenario == scen) & (reps.attribute == a) & (reps.effect_log2_or == near)]
                ratio = (rr["matched"] / rr["baseline"]).replace([np.inf, -np.inf], np.nan).dropna()
                out["nearest_grid_effect"] = near
                out["sim_retention_p2_5_at_nearest"] = float(np.percentile(ratio, 2.5))
                out["sim_retention_p50_at_nearest"] = float(np.percentile(ratio, 50))
                out["sim_retention_p97_5_at_nearest"] = float(np.percentile(ratio, 97.5))
                out["observed_retention_below_sim_p2_5"] = bool(o["retention_matched_over_baseline"]
                                                               < out["sim_retention_p2_5_at_nearest"])
            if np.isfinite(e_m):
                bm = np.interp(e_m, x, yb)
                out["sim_baseline_at_matched_effect"] = bm
                out["excess_baseline_over_sim"] = o["baseline"] - bm
                out["implied_retention_vs_truth"] = o["matched"] / e_m if e_m > 0 else np.nan
            rows.append(out)
    return pd.DataFrame(rows)


def boot_check(reps):
    a = reps[reps.scenario == "S0"]
    b = reps[reps.scenario == "S0boot5000"]
    if a.empty or b.empty:
        return pd.DataFrame()
    m = b.merge(a, on=["rep", "attribute"], suffixes=("_5000", "_1000"))
    rows = []
    for col in ("verdict_primary", "verdict_final", "verdict_final_guarded"):
        rows.append({"verdict_column": col, "n_tables": int(len(m)),
                     "agreement": float((m[col + "_5000"] == m[col + "_1000"]).mean())})
    for nm in ("matched", "baseline", "stratified", "random"):
        for s in ("lo", "hi"):
            c = "%s_%s" % (nm, s)
            rows.append({"verdict_column": "endpoint_" + c, "n_tables": int(len(m)),
                         "mean_abs_diff": float((m[c + "_5000"] - m[c + "_1000"]).abs().mean())})
        exc5 = (m[nm + "_lo_5000"] > 0) | (m[nm + "_hi_5000"] < 0)
        exc1 = (m[nm + "_lo_1000"] > 0) | (m[nm + "_hi_1000"] < 0)
        rows.append({"verdict_column": "excludes_0_" + nm, "n_tables": int(len(m)),
                     "agreement": float((exc5 == exc1).mean())})
    return pd.DataFrame(rows)


def tool_check(reps):
    t = reps.dropna(subset=["tool_controlled"]) if "tool_controlled" in reps.columns else pd.DataFrame()
    if t.empty:
        return t
    rows = []
    for (scen, att), g in t.groupby(["scenario", "attribute"]):
        rows.append({"scenario": scen, "attribute": att, "n": int(len(g)),
                     "mean_abs_diff_tool_vs_pipeline_stratified": float((g["tool_controlled"] - g["stratified"]).abs().mean()),
                     "max_abs_diff": float((g["tool_controlled"] - g["stratified"]).abs().max()),
                     "tool_ci_excludes_0": float(((g["tool_controlled_lo"] > 0) | (g["tool_controlled_hi"] < 0)).mean()),
                     "pipeline_stratified_ci_excludes_0": float(((g["stratified_lo"] > 0) | (g["stratified_hi"] < 0)).mean()),
                     "tool_verdict_equals_pipeline_secondary_guarded": float(
                         (g["tool_verdict"] == [direction(v, b, lo, hi) for v, b, lo, hi in
                                               zip(g["verdict_secondary"], g["baseline"], g["baseline_lo"], g["baseline_hi"])]).mean())})
    return pd.DataFrame(rows)


def direction(v, b, lo, hi):
    if np.isfinite(lo) and np.isfinite(hi) and not (lo <= 0 <= hi) and b < 0:
        return "baseline_contradicts_claim"
    return v


def tool_vs_pipeline(reps):
    """Revision: the released tool (v0.2.2 claim_retest: propensity-quintile Mantel-Haenszel only, with a
    permuted-stratum random control when <= 5% of rows are dropped, and its own verdict rule, which also
    requires the baseline interval to exclude zero) against the pipeline's final verdict (matched primary,
    stratified secondary, downgrade on disagreement, phase-2b direction guard) on the same tables."""
    t = reps.dropna(subset=["tool_verdict"]) if "tool_verdict" in reps.columns else pd.DataFrame()
    if t.empty:
        return t
    rows = []
    for (scen, eff, att), g in t.groupby(["scenario", "effect_log2_or", "attribute"]):
        n = len(g)
        tool_sa = g["tool_verdict"].isin(["survives", "attenuated"])
        pipe_sa = g["verdict_final_guarded"].isin(["survives", "attenuated"])
        sec_guarded = [direction(v, b, lo, hi) for v, b, lo, hi in
                       zip(g["verdict_secondary"], g["baseline"], g["baseline_lo"], g["baseline_hi"])]
        rows.append({
            "scenario": scen, "effect_log2_or": eff, "attribute": att, "n_tables": n,
            "tool_survives_or_attenuated_n": int(tool_sa.sum()),
            "tool_survives_or_attenuated_frac": float(tool_sa.mean()),
            "tool_survives_or_attenuated_mcse": prop_se(float(tool_sa.mean()), n),
            "pipeline_final_survives_or_attenuated_n": int(pipe_sa.sum()),
            "pipeline_final_survives_or_attenuated_frac": float(pipe_sa.mean()),
            "pipeline_final_survives_or_attenuated_mcse": prop_se(float(pipe_sa.mean()), n),
            "tool_ci_excludes_0": float(((g["tool_controlled_lo"] > 0) | (g["tool_controlled_hi"] < 0)).mean()),
            "pipeline_stratified_ci_excludes_0": float(((g["stratified_lo"] > 0) | (g["stratified_hi"] < 0)).mean()),
            "pipeline_matched_ci_excludes_0": float(((g["matched_lo"] > 0) | (g["matched_hi"] < 0)).mean()),
            "tool_equals_pipeline_final": float((g["tool_verdict"] == g["verdict_final_guarded"]).mean()),
            "tool_equals_pipeline_secondary_guarded": float((g["tool_verdict"] == np.asarray(sec_guarded)).mean()),
            "tool_verdict_counts": json.dumps(g["tool_verdict"].value_counts().to_dict(), sort_keys=True),
            "pipeline_final_counts": json.dumps(g["verdict_final_guarded"].value_counts().to_dict(), sort_keys=True),
            "tool_random_control_kind": json.dumps(g["tool_random_kind"].value_counts().to_dict(), sort_keys=True),
        })
    return pd.DataFrame(rows)


def null_exclusion_vs_latent(reps):
    """Revision: P(matched interval excludes 0) under null chemistry, and its paired difference from the
    latent attribute analysed on the same tables (MC SE from the per-table differences)."""
    rows = []
    for scen in ("S0", "S0steep", "S0steepvis", "S0steepnc"):
        d = reps[(reps.scenario == scen) & reps["with_ci"].astype(bool)]
        if d.empty:
            continue
        exc = ((d["matched_lo"] > 0) | (d["matched_hi"] < 0)).astype(float)
        d = d.assign(exc=exc.values)
        wide = d.pivot_table(index="rep", columns="attribute", values="exc")
        if "a4_LATENT30" not in wide.columns:
            continue
        for a in [c for c in ATTR_ORDER if c in wide.columns]:
            diff = (wide[a] - wide["a4_LATENT30"]).dropna()
            p = float(wide[a].mean())
            rows.append({"scenario": scen, "attribute": a, "n_tables": int(wide[a].notna().sum()),
                         "p_matched_ci_excludes_0": p, "mcse": prop_se(p, int(wide[a].notna().sum())),
                         "paired_diff_vs_latent": float(diff.mean()),
                         "paired_diff_mcse": float(diff.std(ddof=1) / np.sqrt(len(diff))) if len(diff) > 2 else np.nan})
    return pd.DataFrame(rows)


ART1_MODELS = [("real", "S0art1", "S1art1", 0.0), ("steep", "S0art1steep", "S1art1steep", 2.0),
               ("steep", "S0art1steep4", "S1art1steep4", 4.0),
               ("steep_vis", "S0steepvis", "S1steepvis", None), ("steep_nocomp", "S0steepnc", "S1steepnc", None)]


def artifact1_shift(summary):
    """Revision: behaviour of the Artifact-1 restriction (background = positives + theoretically detectable
    cysteines) with and without planted chemistry, per selection model and attribute."""
    rows = []
    for mode, s0, s1, gam in ART1_MODELS:
        for a in ATTR_ORDER:
            cells = []
            for scen in (s0, s1):
                q = summary[(summary.scenario == scen) & (summary.detection_mode == mode) & (summary.attribute == a)]
                if gam is not None:
                    q = q[np.isclose(q.gamma, gam)]
                q = q[q["a1_wp_z_all_mean"].notna()] if "a1_wp_z_all_mean" in q.columns else q.iloc[0:0]
                cells.append(q)
            r0 = cells[0].iloc[0] if len(cells[0]) else None
            for _, r in pd.concat(cells).sort_values("effect_log2_or").iterrows():
                out = {"detection_model": mode, "gamma": r["gamma"], "scenario": r["scenario"], "attribute": a,
                       "effect_log2_or": r["effect_log2_or"], "n_tables": int(r["n_reps"]),
                       "mean_prop_auc_before": r["mean_prop_auc_before"], "mean_n_pos": r["mean_n_pos"],
                       "z_all": r["a1_wp_z_all_mean"], "z_restricted": r["a1_wp_z_restricted_mean"],
                       "z_all_mcse": r["a1_wp_z_all_mcse"], "z_restricted_mcse": r["a1_wp_z_restricted_mcse"],
                       "z_ratio_of_means": r["a1_wp_z_restricted_mean"] / r["a1_wp_z_all_mean"],
                       "mh_all": r["a1_wp_mh_all_mean"], "mh_restricted": r["a1_wp_mh_restricted_mean"],
                       "mh_shift": r["a1_mh_shift_mean"], "mh_shift_mcse": r["a1_mh_shift_mcse"],
                       "z_shift": r["a1_z_shift_mean"], "z_shift_mcse": r["a1_z_shift_mcse"]}
                if r0 is not None and r["effect_log2_or"] > 0:
                    d_all = r["a1_wp_mh_all_mean"] - r0["a1_wp_mh_all_mean"]
                    d_res = r["a1_wp_mh_restricted_mean"] - r0["a1_wp_mh_restricted_mean"]
                    out["planted_mh_increment_all"] = d_all
                    out["planted_mh_increment_restricted"] = d_res
                    out["share_of_planted_increment_passing_restriction"] = d_res / d_all if d_all else np.nan
                rows.append(out)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------------
# revision round 2 (post hoc): proxy-chemistry model and count-defined claims
# ---------------------------------------------------------------------------------------------------------------
def _ratio(g):
    return (g["matched"] / g["baseline"]).replace([np.inf, -np.inf], np.nan).dropna()


def _mcse(v):
    v = pd.Series(v).dropna()
    return float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else np.nan


def proxy_grid_summary(reps):
    """Per proxy grid, planted strength b and analysed attribute."""
    rows = []
    for scen, mode, feat in PROXY_GRIDS:
        d = reps[reps.scenario == scen]
        for (eff, att), g in d.groupby(["effect_log2_or", "attribute"]):
            r = _ratio(g)
            row = {"scenario": scen, "detection_model": mode, "gamma": float(g["gamma"].iloc[0]),
                   "planted_feature": feat, "b_log2_or_per_sd": eff, "attribute": att,
                   "attribute_class": ATTR_CLASS.get(att, ""), "n_tables": int(len(g)),
                   "prop_auc_mean": g["prop_auc_before"].mean(), "n_pos_mean": g["n_pos"].mean(),
                   "baseline_mean": g["baseline"].mean(), "baseline_mcse": _mcse(g["baseline"]),
                   "matched_mean": g["matched"].mean(), "matched_mcse": _mcse(g["matched"]),
                   "stratified_mean": g["stratified"].mean(),
                   "matched_over_baseline_ratio_of_means": g["matched"].mean() / g["baseline"].mean(),
                   "ratio_p2_5": float(np.percentile(r, 2.5)) if len(r) else np.nan,
                   "ratio_p50": float(np.percentile(r, 50)) if len(r) else np.nan,
                   "ratio_p97_5": float(np.percentile(r, 97.5)) if len(r) else np.nan,
                   "realised_within_detected_mean": g["realised_within_detected"].mean()}
            if "a1_wp_mh_all" in g.columns and g["a1_wp_mh_all"].notna().any():
                row.update(mh_all_mean=g["a1_wp_mh_all"].mean(), mh_restricted_mean=g["a1_wp_mh_restricted"].mean(),
                           z_all_mean=g["a1_wp_z_all"].mean(), z_restricted_mean=g["a1_wp_z_restricted"].mean())
            rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["attr_order"] = out["attribute"].map({a: i for i, a in enumerate(ATTR_ORDER)})
    return out.sort_values(["scenario", "b_log2_or_per_sd", "attr_order"]).drop(columns="attr_order").reset_index(drop=True)


def observed_vs_proxy(ob, reps, summary, ovs):
    """Each observed sequence-attribute claim against the proxy-chemistry grids (chemistry on a K/R-count feature,
    not on the claim's attribute). Band = per-table matched/baseline of the CLAIM'S attribute at the grid b nearest
    to the b whose mean simulated baseline equals the observed baseline. Where no b on the grid reproduces the
    observed baseline, the band is not defined and the largest simulated baseline is reported instead.
    'proxy_excess_baseline' applies the implementer's excess diagnostic to the proxy tables: the proxy cell's mean
    baseline minus the baseline that the attribute-planted grid of the same detection model gives at the effect
    reproducing the proxy cell's mean matched estimate (the same computation as for the observed claims)."""
    rows = []
    for _, o in ob.dropna(subset=["simulated_attribute"]).iterrows():
        a = o["simulated_attribute"]
        if a == "a4_LATENT30":   # SFE-001 stand-in is sequence-independent; the proxy model does not apply
            continue
        for scen, mode, feat in PROXY_GRIDS:
            d = reps[(reps.scenario == scen) & (reps.attribute == a)]
            if d.empty:
                continue
            grp = d.groupby("effect_log2_or")
            x = np.asarray(sorted(d["effect_log2_or"].unique()))
            yb = grp["baseline"].mean().reindex(x).to_numpy()
            ym = grp["matched"].mean().reindex(x).to_numpy()
            med = np.asarray([np.percentile(_ratio(grp.get_group(e)), 50) for e in x])
            big = x >= 1.0
            out = {"claim": o["claim"], "simulated_attribute": a, "proxy_scenario": scen, "detection_model": mode,
                   "planted_feature": feat, "observed_baseline": o["baseline"], "observed_matched": o["matched"],
                   "observed_retention": o["retention_matched_over_baseline"],
                   "observed_propensity_auc": o["propensity_auc"],
                   "grid_b": ";".join("%g" % e for e in x),
                   "max_sim_baseline": float(np.nanmax(yb)), "b_at_max_sim_baseline": float(x[np.nanargmax(yb)]),
                   "retention_median_min_b_ge_1": float(med[big].min()) if big.any() else np.nan,
                   "retention_median_max_b_ge_1": float(med[big].max()) if big.any() else np.nan}
            e_b = interp_inverse(x, yb, o["baseline"])
            out["b_reproducing_observed_baseline"] = e_b
            # band at the b with the largest mean baseline (used when the observed baseline is not reached)
            gmax = grp.get_group(float(x[np.nanargmax(yb)]))
            rmax = _ratio(gmax)
            out.update(band_at_max_p2_5=float(np.percentile(rmax, 2.5)), band_at_max_p50=float(np.percentile(rmax, 50)),
                       band_at_max_p97_5=float(np.percentile(rmax, 97.5)))
            if np.isfinite(e_b):
                near = float(x[np.argmin(np.abs(x - e_b))])
                g = grp.get_group(near)
                r = _ratio(g)
                out.update(nearest_grid_b=near, n_tables=int(len(g)),
                           band_p2_5=float(np.percentile(r, 2.5)), band_p50=float(np.percentile(r, 50)),
                           band_p97_5=float(np.percentile(r, 97.5)),
                           sim_baseline_at_nearest=float(g["baseline"].mean()),
                           sim_matched_at_nearest=float(g["matched"].mean()),
                           sim_prop_auc_at_nearest=float(g["prop_auc_before"].mean()),
                           sim_retention_of_means_at_b=float(np.interp(e_b, x, ym)) / o["baseline"])
                out["observed_below_band"] = bool(o["retention_matched_over_baseline"] < out["band_p2_5"])
                out["observed_above_band"] = bool(o["retention_matched_over_baseline"] > out["band_p97_5"])
                out["observed_within_band"] = not (out["observed_below_band"] or out["observed_above_band"])
                # implied association of the label with the count-defined attribute in the same tables
                cx = reps[(reps.scenario == scen) & (reps.attribute == "x1_KRcount20_hi")
                          & np.isclose(reps.effect_log2_or, near)]
                if len(cx):
                    out["x1_baseline_at_nearest"] = float(cx["baseline"].mean())
                    out["x1_realised_within_detected_at_nearest"] = float(cx["realised_within_detected"].mean())
                # the excess-baseline diagnostic on the proxy tables
                if mode in ATTR_GRID_OF_MODE:
                    sg = summary[(summary.scenario == ATTR_GRID_OF_MODE[mode]) & (summary.attribute == a)] \
                        .sort_values("effect_log2_or")
                    if len(sg):
                        xs = sg["effect_log2_or"].to_numpy()
                        e_m = interp_inverse(xs, sg["matched_mean"].to_numpy(), out["sim_matched_at_nearest"])
                        if np.isfinite(e_m):
                            bg = float(np.interp(e_m, xs, sg["baseline_mean"].to_numpy()))
                            out["proxy_effect_reproducing_matched_on_attribute_grid"] = e_m
                            out["proxy_excess_baseline"] = out["sim_baseline_at_nearest"] - bg
            q = ovs[(ovs.claim == o["claim"]) & (ovs.detection_model == mode)]
            if len(q):
                out["observed_excess_baseline"] = q["excess_baseline_over_sim"].iloc[0]
                out["attribute_planted_band"] = "%.2f [%.2f, %.2f]" % (
                    q["sim_retention_p50_at_nearest"].iloc[0], q["sim_retention_p2_5_at_nearest"].iloc[0],
                    q["sim_retention_p97_5_at_nearest"].iloc[0]) \
                    if np.isfinite(q["sim_retention_p50_at_nearest"].iloc[0]) else ""
            rows.append(out)
    return pd.DataFrame(rows)


def proxy_joint_at_calibration(ovp, reps):
    """Every analysed attribute in the proxy tables at the grid b calibrated to one observed claim."""
    rows = []
    for _, c in ovp.dropna(subset=["nearest_grid_b"]).iterrows():
        d = reps[(reps.scenario == c["proxy_scenario"]) & np.isclose(reps.effect_log2_or, c["nearest_grid_b"])]
        for att, g in d.groupby("attribute"):
            r = _ratio(g)
            rows.append({"calibrated_to": c["claim"], "proxy_scenario": c["proxy_scenario"],
                         "detection_model": c["detection_model"], "b": c["nearest_grid_b"], "attribute": att,
                         "n_tables": int(len(g)), "baseline_mean": g["baseline"].mean(),
                         "matched_mean": g["matched"].mean(),
                         "ratio_p2_5": float(np.percentile(r, 2.5)) if len(r) else np.nan,
                         "ratio_p50": float(np.percentile(r, 50)) if len(r) else np.nan,
                         "ratio_p97_5": float(np.percentile(r, 97.5)) if len(r) else np.nan,
                         "prop_auc_mean": g["prop_auc_before"].mean()})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["attr_order"] = out["attribute"].map({a: i for i, a in enumerate(ATTR_ORDER)})
    return out.sort_values(["calibrated_to", "proxy_scenario", "attr_order"]).drop(columns="attr_order")


def proxy_verdicts(reps):
    """Pipeline verdicts (final, with the specification downgrade and the phase-2b direction guard) under proxy
    chemistry, interval tables."""
    rows = []
    for scen in ("P1ci", "P1permci"):
        d = reps[(reps.scenario == scen) & reps["with_ci"].astype(bool)]
        for att, g in d.groupby("attribute"):
            k = len(g)
            r = _ratio(g)
            exc = ((g["matched_lo"] > 0) | (g["matched_hi"] < 0)).astype(float)
            row = {"scenario": scen, "detection_model": g["detection_mode"].iloc[0],
                   "b": float(g["effect_log2_or"].iloc[0]), "attribute": att, "n_tables": int(k),
                   "prop_auc_mean": g["prop_auc_before"].mean(), "baseline_mean": g["baseline"].mean(),
                   "matched_mean": g["matched"].mean(), "ratio_p50": float(np.percentile(r, 50)),
                   "ratio_p2_5": float(np.percentile(r, 2.5)), "ratio_p97_5": float(np.percentile(r, 97.5)),
                   "p_matched_ci_excludes_0": float(exc.mean()), "p_matched_ci_excludes_0_mcse": prop_se(exc.mean(), k)}
            for v in VERDICTS:
                p = float((g["verdict_final_guarded"] == v).mean())
                row["final_" + v] = p
                row["final_" + v + "_mcse"] = prop_se(p, k)
                row["primary_" + v] = float((g["verdict_primary"] == v).mean())
            rows.append(row)
    return pd.DataFrame(rows)


def proxy_artifact1(reps):
    """Artifact-1 restriction (positives + theoretically detectable cysteines) under proxy chemistry: per grid
    and b, the MH log2 OR against all cysteines and after the restriction, the shift, and the share of the
    proxy-induced increment (relative to b = 0 of the same grid) that passes the restriction."""
    rows = []
    for scen, mode in (("P1gridA1", "real"), ("P1steepgridA1", "steep")):   # tables identical to P1grid / P1steepgrid
        d = reps[(reps.scenario == scen) & reps["a1_wp_mh_all"].notna()] if "a1_wp_mh_all" in reps.columns else reps.iloc[0:0]
        for att in ("a1_SFE006_KR", "a1c_KRH_p5p8", "x1_KRcount20_hi", "a3_SNO006_DE3", "a4_LATENT30"):
            q = d[d.attribute == att]
            if q.empty:
                continue
            g = q.groupby("effect_log2_or")
            all_m, res_m = g["a1_wp_mh_all"].mean(), g["a1_wp_mh_restricted"].mean()
            sh = g["a1_mh_shift"].agg(["mean", "std", "size"])
            for e in sorted(q["effect_log2_or"].unique()):
                row = {"scenario": scen, "detection_model": mode, "attribute": att, "b": e,
                       "n_tables": int(sh.loc[e, "size"]), "mh_all": all_m[e], "mh_restricted": res_m[e],
                       "mh_shift": sh.loc[e, "mean"], "mh_shift_mcse": sh.loc[e, "std"] / np.sqrt(sh.loc[e, "size"]),
                       "z_all": g["a1_wp_z_all"].mean()[e], "z_restricted": g["a1_wp_z_restricted"].mean()[e]}
                if e > 0 and 0.0 in all_m.index:
                    inc_all = all_m[e] - all_m[0.0]
                    inc_res = res_m[e] - res_m[0.0]
                    row["proxy_increment_all"] = inc_all
                    row["proxy_increment_restricted"] = inc_res
                    row["share_of_increment_passing_restriction"] = inc_res / inc_all if inc_all else np.nan
                rows.append(row)
    return pd.DataFrame(rows)


def count_claim_summary(oa, reps):
    """The two hypothetical count-defined claims next to the published definitions: loss beyond the latent
    attribute under both normalisations, per detection model and planted effect, with MC SEs, and the per-table
    matched/baseline distribution of the planted cells."""
    pairs = [("real", "S1"), ("steep", "S1steep"), ("permuted_within_protein", "S1permgrid")]
    rows = []
    for mode, scen in pairs:
        q = oa[(oa.detection_model == mode) & (oa.planted_scenario == scen) & oa.effect_log2_or.isin([0.5, 1.0, 1.5])]
        for _, r in q.iterrows():
            if r["attribute"] in ("a4_LATENT30", "a1c_KRH_p5p8"):
                continue
            A = r.get("matched_genuine_part_retention_relative_to_latent", np.nan)
            B = r.get("matched_increment_normalised_retention_relative_to_latent", np.nan)
            g = reps[(reps.scenario == scen) & (reps.detection_mode == mode) & (reps.attribute == r["attribute"])
                     & np.isclose(reps.effect_log2_or, r["effect_log2_or"])]
            rr = _ratio(g)
            rows.append({"detection_model": mode, "planted_scenario": scen, "attribute": r["attribute"],
                         "attribute_class": r["attribute_class"], "effect_log2_or": r["effect_log2_or"],
                         "n_tables": int(len(g)),
                         "S0_baseline": r.get("baseline_S0_mean"), "S0_matched": r.get("matched_S0_mean"),
                         "S1_baseline": r.get("baseline_S1_mean"), "S1_matched": r.get("matched_S1_mean"),
                         "matched_over_truth": r.get("matched_retention_vs_truth"),
                         "A_rel_latent": A, "A_rel_latent_mcse": r.get("matched_genuine_part_retention_relative_to_latent_mcse"),
                         "B_rel_latent": B,
                         "B_rel_latent_mcse": r.get("matched_increment_normalised_retention_relative_to_latent_mcse"),
                         "loss_beyond_latent_A": 1 - A, "loss_beyond_latent_B": 1 - B,
                         "ratio_p2_5": float(np.percentile(rr, 2.5)) if len(rr) else np.nan,
                         "ratio_p50": float(np.percentile(rr, 50)) if len(rr) else np.nan,
                         "ratio_p97_5": float(np.percentile(rr, 97.5)) if len(rr) else np.nan})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["attr_order"] = out["attribute"].map({a: i for i, a in enumerate(ATTR_ORDER)})
    return out.sort_values(["detection_model", "attr_order", "effect_log2_or"]).drop(columns="attr_order")


def main():
    reps = load_reps()
    summary = build_summary(reps)
    summary.to_csv(os.path.join(cb.RESULTS, "sim_summary.csv"), index=False)
    verdict_table(reps).to_csv(os.path.join(cb.RESULTS, "verdict_summary.csv"), index=False)
    oa = overadjustment(summary)
    oa.to_csv(os.path.join(cb.RESULTS, "overadjustment_summary.csv"), index=False)
    ob = observed_claims()
    ob.to_csv(os.path.join(cb.RESULTS, "observed_claims.csv"), index=False)
    ovs = observed_vs_simulated(ob, reps, summary)
    ovs.to_csv(os.path.join(cb.RESULTS, "observed_vs_simulated.csv"), index=False)
    bc = boot_check(reps)
    bc.to_csv(os.path.join(cb.RESULTS, "boot_check.csv"), index=False)
    tc = tool_check(reps)
    tc.to_csv(os.path.join(cb.RESULTS, "tool_check.csv"), index=False)
    # revision after adversarial verification
    tool_vs_pipeline(reps).to_csv(os.path.join(cb.RESULTS, "tool_vs_pipeline.csv"), index=False)
    null_exclusion_vs_latent(reps).to_csv(os.path.join(cb.RESULTS, "null_exclusion_vs_latent.csv"), index=False)
    artifact1_shift(summary).to_csv(os.path.join(cb.RESULTS, "artifact1_shift_summary.csv"), index=False)
    # revision round 2
    proxy_grid_summary(reps).to_csv(os.path.join(cb.RESULTS, "proxy_grid_summary.csv"), index=False)
    ovp = observed_vs_proxy(ob, reps, summary, ovs)
    ovp.to_csv(os.path.join(cb.RESULTS, "observed_vs_proxy.csv"), index=False)
    if len(ovp):
        proxy_joint_at_calibration(ovp, reps).to_csv(os.path.join(cb.RESULTS, "proxy_joint_at_calibration.csv"),
                                                     index=False)
    proxy_verdicts(reps).to_csv(os.path.join(cb.RESULTS, "proxy_verdicts.csv"), index=False)
    proxy_artifact1(reps).to_csv(os.path.join(cb.RESULTS, "proxy_artifact1.csv"), index=False)
    count_claim_summary(oa, reps).to_csv(os.path.join(cb.RESULTS, "count_claim_summary.csv"), index=False)
    print("cells", len(summary), "replicate rows", len(reps))


if __name__ == "__main__":
    main()
