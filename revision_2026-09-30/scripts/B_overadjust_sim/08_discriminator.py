# -*- coding: utf-8 -*-
"""B_overadjust_sim step 8 (POST HOC, revision round 3, 2026-09-30): the sequence-only discriminator.

Round-3 verification observed that the item's proposed text offered only one way to separate the two readings of
the SFE-006 attenuation (a background of detected but unmodified cysteines, which the published cohort lacks),
while the simulation itself contains a cheaper one: in ONE site table, the association of the label with a local
residue count RELATIVE to its association with the claimed SFE-006 offsets. This script measures that ratio in
every simulated generating model, at the strength that reproduces the observed SFE-006 baseline, and turns the
simulated bands into a reading rule for the authors' own cohort (09_sfe006_count_test.py).

Discriminators (all defined as "count above its universe median", like the round-2 count-defined claim x1):
  x1_KRcount20_hi   K/R within +/-20   (the count VIS10 omits; the M2 chemistry variable)
  x3_AKRV10_hi      A/K/R/V within +/-10 (the M2' sensitivity count for PERS-009)
  x4_KR5_hi         K/R within +/-5    (the other M2' sensitivity count)

Generating models compared (all analysed in the SAME tables, 200 per cell unless stated):
  detectability only, no chemistry: real (500 tables), steep gamma 2 (300), steep gamma 4 (200), steep_vis (200),
      steep_nocomp (200), permuted (200; reported but NOT used in the rule, because its SFE-006 baseline is 0.09
      and a ratio on it is not interpretable);
  M1  chemistry on the SFE-006 offsets, at the strength reproducing the observed baseline (real/steep/permuted);
  M2  chemistry on the K/R count within +/-20, likewise;
  M2' chemistry on the A/K/R/V count within +/-10, likewise (real detection; the K/R +/-5 count does not reach the
      observed baseline for b <= 6 and has no calibrated cell);
  mixtures of M2 and M1 (count chemistry at b = 0.5-2.0 per SD plus offsets chemistry at the c that restores the
      observed baseline), real detection;
  cohort size: the real-detection cells re-run at rate 0.16 (about 1,100 positives, as in the SFE-006 cohort).

Outputs (results/B_overadjust_sim/):
  discriminator_summary.csv      one row per cell: baselines, retention, every ratio with its 2.5/50/97.5 band
  discriminator_mixture.csv      the mixture grid and the calibrated mixtures
  discriminator_reading_rule.json  the rule applied to the authors' cohort, written BEFORE any cohort value was
                                 computed; 09_sfe006_count_test.py hashes it into its own output
  discriminator_verifier_check.csv  this item's numbers against the round-3 verifier's (verify_r3, read only)
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

DISCRIM = ["x1_KRcount20_hi", "x3_AKRV10_hi", "x4_KR5_hi"]
A1 = "a1_SFE006_KR"
OBSERVED = {"baseline": 1.3840344206722544, "matched": 0.594987270476461, "retention": 0.4299, "prop_auc": 0.7616}
# cells: scenario -> (model family, detection model, note)
CELLS = {
    "D0real": ("detectability only", "real", "no chemistry"),
    "D0steep": ("detectability only", "steep (gamma 2)", "no chemistry"),
    "D0steep4": ("detectability only", "steep (gamma 4)", "no chemistry; steepest simulated selection"),
    "D0steepvis": ("detectability only", "steep_vis (gamma 0.25)", "no chemistry"),
    "D0steepnc": ("detectability only", "steep_nocomp (gamma 3)", "no chemistry"),
    "D0perm": ("detectability only", "permuted", "no chemistry; SFE-006 baseline ~0, ratio not interpretable"),
    "DM1real": ("M1 chemistry on the SFE-006 offsets", "real", "calibrated to the observed baseline"),
    "DM1steep": ("M1 chemistry on the SFE-006 offsets", "steep (gamma 2)", "calibrated"),
    "DM1perm": ("M1 chemistry on the SFE-006 offsets", "permuted", "calibrated"),
    "DM2real": ("M2 chemistry on the K/R count within +/-20", "real", "calibrated"),
    "DM2steep": ("M2 chemistry on the K/R count within +/-20", "steep (gamma 2)", "calibrated"),
    "DM2perm": ("M2 chemistry on the K/R count within +/-20", "permuted", "calibrated"),
    "DM2akrv10": ("M2' chemistry on the A/K/R/V count within +/-10", "real", "calibrated"),
    "D0r16": ("detectability only", "real, rate 0.16", "cohort-size check"),
    "DM1r16": ("M1 chemistry on the SFE-006 offsets", "real, rate 0.16", "cohort-size check"),
    "DM2r16": ("M2 chemistry on the K/R count within +/-20", "real, rate 0.16", "cohort-size check"),
}
RULE_NULL = ["D0real", "D0steep", "D0steep4", "D0steepvis", "D0steepnc"]   # permuted excluded (baseline ~0)
RULE_CHEM = {"x1_KRcount20_hi": ["DM2real", "DM2steep", "DM2perm"],
             "x3_AKRV10_hi": ["DM2akrv10"], "x4_KR5_hi": []}
RULE_M1 = ["DM1real", "DM1steep", "DM1perm"]


def load():
    parts = []
    for f in ("reps_discrim.csv.gz", "reps_discrim2.csv.gz"):
        p = os.path.join(cb.RESULTS, f)
        if os.path.exists(p):
            parts.append(pd.read_csv(p, low_memory=False))
    d = pd.concat(parts, ignore_index=True)
    if "effect2_log2_or" not in d.columns:
        d["effect2_log2_or"] = 0.0
    d["effect2_log2_or"] = d["effect2_log2_or"].fillna(0.0)
    return d


def wide(d, keys):
    b = d.pivot_table(index=keys, columns="attribute", values="baseline")
    m = d.pivot_table(index=keys, columns="attribute", values="matched")
    nb = d.pivot_table(index=keys, columns="attribute", values="neg_b_oracle")
    auc = d[d.attribute == A1].set_index(keys)["prop_auc_before"]
    w = pd.DataFrame(index=b.index)
    w["a1_baseline"] = b[A1]
    w["a1_matched"] = m[A1]
    w["a1_retention"] = m[A1] / b[A1]
    w["a1_neg_b"] = nb[A1]
    w["prop_auc"] = auc
    for k in DISCRIM:
        w[k + "_baseline"] = b[k]
        w[k + "_ratio"] = b[k] / b[A1]
    return w


def band(s):
    return float(s.median()), float(s.quantile(0.025)), float(s.quantile(0.975))


def summarise(w, group_keys):
    rows = []
    for key, g in w.groupby(level=group_keys):
        rec = {"n_tables": int(len(g)), "a1_baseline_mean": float(g.a1_baseline.mean()),
               "a1_retention_p50": float(g.a1_retention.median()),
               "a1_retention_p2_5": float(g.a1_retention.quantile(.025)),
               "a1_retention_p97_5": float(g.a1_retention.quantile(.975)),
               "a1_neg_b_mean": float(g.a1_neg_b.mean()), "prop_auc_mean": float(g.prop_auc.mean())}
        if isinstance(key, tuple):
            rec.update(dict(zip(group_keys, key)))
        else:
            rec[group_keys[0]] = key
        for k in DISCRIM:
            rec[k + "_baseline_mean"] = float(g[k + "_baseline"].mean())
            p50, lo, hi = band(g[k + "_ratio"])
            rec[k + "_ratio_p50"], rec[k + "_ratio_p2_5"], rec[k + "_ratio_p97_5"] = p50, lo, hi
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    d = load()
    keys = ["scenario", "effect_log2_or", "effect2_log2_or", "rep"]
    w = wide(d, keys)

    # ---- one row per simulation cell (no mixtures) ----
    main_idx = w.index.get_level_values("scenario").isin(CELLS)
    s = summarise(w[main_idx], ["scenario"])
    s["model"] = s.scenario.map(lambda x: CELLS[x][0])
    s["detection"] = s.scenario.map(lambda x: CELLS[x][1])
    s["note"] = s.scenario.map(lambda x: CELLS[x][2])
    strength = (d[d.attribute == A1].groupby("scenario")[["effect_log2_or", "effect2_log2_or"]].mean())
    s = s.merge(strength, left_on="scenario", right_index=True, how="left")
    s = s.sort_values(["model", "scenario"])
    cols = (["scenario", "model", "detection", "note", "n_tables", "effect_log2_or", "effect2_log2_or",
             "a1_baseline_mean", "a1_retention_p50", "a1_retention_p2_5", "a1_retention_p97_5", "prop_auc_mean",
             "a1_neg_b_mean"]
            + [k + suff for k in DISCRIM for suff in ("_baseline_mean", "_ratio_p50", "_ratio_p2_5", "_ratio_p97_5")])
    s[cols].to_csv(os.path.join(cb.RESULTS, "discriminator_summary.csv"), index=False)

    # ---- mixtures ----
    mix_idx = w.index.get_level_values("scenario").isin(["DMIXgrid", "DMIX"])
    mx = summarise(w[mix_idx], ["scenario", "effect_log2_or", "effect2_log2_or"])
    mx = mx.sort_values(["scenario", "effect_log2_or", "effect2_log2_or"])
    mx.to_csv(os.path.join(cb.RESULTS, "discriminator_mixture.csv"), index=False)

    # ---- reading rule (written BEFORE any value is computed on the authors' cohort) ----
    idx = s.set_index("scenario")
    rule = {"item": cb.ITEM,
            "label": "POST HOC (revision round 3, 2026-09-30); NOT registered. Written from the simulation only, "
                     "before any value was computed on the authors' SFE-006 cohort.",
            "statistic": "baseline log2 odds ratio of a count-defined attribute divided by the baseline log2 odds "
                         "ratio of the SFE-006 offsets, in the same site table (NEG_A background, no matching)",
            "thresholds_of_the_count_attributes": "count > the median of the simulation universe "
                                                  "(K/R within +/-20 > 4; A/K/R/V within +/-10 > 5; K/R within +/-5 > 1)",
            "observed_sfe006": OBSERVED,
            "bands_used": {}, "discriminators": {}}
    for k in DISCRIM:
        nulls = idx.loc[RULE_NULL]
        det_hi = float(nulls[k + "_ratio_p97_5"].max())
        det_lo = float(nulls[k + "_ratio_p2_5"].min())
        chem = RULE_CHEM[k]
        chem_lo = float(idx.loc[chem][k + "_ratio_p2_5"].min()) if chem else None
        chem_hi = float(idx.loc[chem][k + "_ratio_p97_5"].max()) if chem else None
        m1 = idx.loc[RULE_M1]
        rule["bands_used"][k] = {
            "detectability_only_models": RULE_NULL, "detectability_only_range": [det_lo, det_hi],
            "chemistry_on_this_count_cells": chem, "chemistry_on_this_count_range": [chem_lo, chem_hi],
            "chemistry_on_the_offsets_cells": RULE_M1,
            "chemistry_on_the_offsets_range": [float(m1[k + "_ratio_p2_5"].min()), float(m1[k + "_ratio_p97_5"].max())],
            "excluded_from_the_rule": {"D0perm": CELLS["D0perm"][2]}}
        rule["discriminators"][k] = {
            "chemistry_lower": chem_lo, "detectability_upper": det_hi,
            "reading_high": ("at or above every simulated detectability-only band and inside the band for genuine "
                             "chemistry following this count: consistent with chemistry on this count, not with a "
                             "detectability artifact of the simulated kinds"),
            "reading_low": ("inside the simulated detectability-only band and below the band for chemistry following "
                            "this count: not evidence for chemistry on this count; chemistry on the claimed offsets "
                            "themselves gives an even lower value"),
            "reading_between": ("between the simulated bands: neither reading is excluded; a mixture of chemistry on "
                                "this count and on the claimed offsets covers this range")}
    rule["caveats"] = [
        "The ratio identifies WHICH feature an effect follows, not whether it is chemistry: chemistry on a count "
        "other than the one measured can leave this ratio inside the detectability-only band (A/K/R/V within +/-10 "
        "chemistry gives %.2f for the K/R +/-20 ratio), so the three counts are read as a panel."
        % float(idx.loc["DM2akrv10", "x1_KRcount20_hi_ratio_p50"]),
        "Mixtures of count chemistry and offsets chemistry that reproduce the observed baseline span %.2f to %.2f "
        "for the K/R +/-20 ratio, so an intermediate value bounds the mixture rather than deciding it."
        % (float(mx[mx.scenario == "DMIX"]["x1_KRcount20_hi_ratio_p2_5"].min()),
           float(mx[mx.scenario == "DMIX"]["x1_KRcount20_hi_ratio_p97_5"].max())),
        "The simulation is one mouse detection process; the cohort is human with an unknown one. The bands are "
        "simulated bands, not confidence intervals for the cohort.",
        "The cohort ratio carries its own sampling error; 09_sfe006_count_test.py reports a protein-clustered "
        "bootstrap interval for it and that interval must be read against these bands.",
    ]
    with open(os.path.join(cb.RESULTS, "discriminator_reading_rule.json"), "w", encoding="utf-8") as fh:
        json.dump(rule, fh, indent=1)

    # ---- cross-check against the round-3 verifier (read only) ----
    vpath = os.path.join(cb.RESULTS, "verify_r3", "discriminator_r3_summary.csv")
    rows = []
    if os.path.exists(vpath):
        v = pd.read_csv(vpath)
        pair = {("none", "real"): "D0real", ("none", "steep"): "D0steep", ("none", "perm"): "D0perm",
                ("M1", "real"): "DM1real", ("M1", "steep"): "DM1steep", ("M1", "perm"): "DM1perm",
                ("M2", "real"): "DM2real", ("M2", "steep"): "DM2steep", ("M2", "perm"): "DM2perm"}
        for _, r in v.iterrows():
            scen = pair.get((r["cell"], r["mode"]))
            if scen is None or scen not in idx.index:
                continue
            rows.append({"cell": "%s/%s" % (r["cell"], r["mode"]), "scenario": scen,
                         "verifier_ratio_p50": r["ratio_p50"], "item_ratio_p50": idx.loc[scen, "x1_KRcount20_hi_ratio_p50"],
                         "verifier_lo": r["ratio_p2_5"], "item_lo": idx.loc[scen, "x1_KRcount20_hi_ratio_p2_5"],
                         "verifier_hi": r["ratio_p97_5"], "item_hi": idx.loc[scen, "x1_KRcount20_hi_ratio_p97_5"],
                         "verifier_a1_baseline": r["a1_baseline"], "item_a1_baseline": idx.loc[scen, "a1_baseline_mean"],
                         "verifier_x1_baseline": r["x1_baseline"], "item_x1_baseline": idx.loc[scen, "x1_KRcount20_hi_baseline_mean"]})
    g4 = os.path.join(cb.RESULTS, "verify_r3", "r3_G4x.csv.gz")
    if os.path.exists(g4):
        vg = pd.read_csv(g4, low_memory=False)
        vw = vg.pivot_table(index="rep", columns="attribute", values="baseline")
        vr = vw["x1_KRcount20_hi"] / vw[A1]
        rows.append({"cell": "none/steep gamma 4", "scenario": "D0steep4",
                     "verifier_ratio_p50": float(vr.median()), "item_ratio_p50": idx.loc["D0steep4", "x1_KRcount20_hi_ratio_p50"],
                     "verifier_lo": float(vr.quantile(.025)), "item_lo": idx.loc["D0steep4", "x1_KRcount20_hi_ratio_p2_5"],
                     "verifier_hi": float(vr.quantile(.975)), "item_hi": idx.loc["D0steep4", "x1_KRcount20_hi_ratio_p97_5"],
                     "verifier_a1_baseline": float(vw[A1].mean()), "item_a1_baseline": idx.loc["D0steep4", "a1_baseline_mean"],
                     "verifier_x1_baseline": float(vw["x1_KRcount20_hi"].mean()),
                     "item_x1_baseline": idx.loc["D0steep4", "x1_KRcount20_hi_baseline_mean"]})
    if rows:
        c = pd.DataFrame(rows)
        c["abs_diff_p50"] = (c.verifier_ratio_p50 - c.item_ratio_p50).abs()
        c.to_csv(os.path.join(cb.RESULTS, "discriminator_verifier_check.csv"), index=False)

    pd.set_option("display.width", 250)
    print(s[["scenario", "model", "detection", "n_tables", "a1_baseline_mean", "a1_retention_p50",
             "x1_KRcount20_hi_ratio_p50", "x1_KRcount20_hi_ratio_p2_5", "x1_KRcount20_hi_ratio_p97_5",
             "x3_AKRV10_hi_ratio_p50", "x4_KR5_hi_ratio_p50"]].round(3).to_string(index=False))
    print()
    print(mx[["scenario", "effect_log2_or", "effect2_log2_or", "a1_baseline_mean", "a1_retention_p50",
              "x1_KRcount20_hi_ratio_p50", "x1_KRcount20_hi_ratio_p2_5",
              "x1_KRcount20_hi_ratio_p97_5"]].round(3).to_string(index=False))
    if rows:
        print()
        print(c.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
