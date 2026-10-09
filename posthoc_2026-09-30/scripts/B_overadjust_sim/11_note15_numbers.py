# -*- coding: utf-8 -*-
"""B_overadjust_sim step 11 (POST HOC, revision round 4, 2026-09-30): number traceability for Supplemental Note 15.

Every number printed in results/B_overadjust_sim/Supplemental_Note_15.md is computed here from a results file (or is a
design constant of the scripts, labelled as such) and written, with its source and the selection used, to
results/B_overadjust_sim/note15_numbers.csv. With --check, the script then extracts every numeric token from the Note
(after removing identifiers such as claim codes, accessions and note numbers) and reports any token that is not the
formatted text of an entry of that table; it exits 1 if one is found.

Formatting: 'f2' two decimals, 'f3' three, 'f1' one, 'p0' percent without decimals, 'p1' percent with one decimal,
'i' integer with thousands separators, 'r10' / 'r100' integer rounded to the nearest 10 / 100. Negative values are
printed with a Unicode minus in the Note; the check normalises it.
"""
from __future__ import annotations

import json
import os
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

R = cb.RESULTS
NOTE = os.path.join(R, "Supplemental_Note_15.md")
ROWS = []


def fmt(v, kind):
    v = float(v)
    if kind == "f1":
        return "%.1f" % v
    if kind == "f2":
        return "%.2f" % v
    if kind == "f3":
        return "%.3f" % v
    if kind == "p0":
        return "%.0f" % (100 * v)
    if kind == "p1":
        return "%.1f" % (100 * v)
    if kind == "i":
        return "{:,}".format(int(round(v)))
    if kind == "r10":
        return "{:,}".format(int(round(v / 10.0) * 10))
    if kind == "r100":
        return "{:,}".format(int(round(v / 100.0) * 100))
    raise ValueError(kind)


def add(key, value, kind, source, how):
    ROWS.append({"id": key, "value": float(value), "text": fmt(value, kind), "format": kind, "source": source,
                 "selection": how})


def csv(name, **kw):
    return pd.read_csv(os.path.join(R, name), low_memory=False, **kw)


def main():
    prep = json.load(open(os.path.join(R, "prep_summary.json"), encoding="utf-8"))
    ver = prep["verification"]
    add("universe proteins", prep["n_universe_proteins"], "i", "prep_summary.json", "n_universe_proteins")
    add("universe cysteines", prep["n_universe_cysteines"], "i", "prep_summary.json", "n_universe_cysteines")
    add("detected cysteines", prep["n_detected_cysteines"], "i", "prep_summary.json", "n_detected_cysteines")
    add("detected fraction", prep["detected_fraction"], "p0", "prep_summary.json", "detected_fraction")
    add("candidate peptides", ver["candidate_rows"], "i", "prep_summary.json", "verification.candidate_rows")
    add("candidate rows matching 2026_03", ver["row_match_rate"], "p1", "prep_summary.json",
        "verification.row_match_rate")
    add("proteins removed", ver["proteins_absent_from_fasta"] + ver["proteins_with_any_mismatch"], "i",
        "prep_summary.json", "verification.proteins_absent_from_fasta + proteins_with_any_mismatch")
    auc = {a["attribute"]: a["vis10_predicts_attribute_auc_cv5"] for a in prep["attributes"]}
    basic = [auc[a] for a in ("a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV")]
    acidic = [auc[a] for a in ("a3_SNO006_DE3", "a3b_SFE002_E")]
    add("VIS10 AUC basic min", min(basic), "f2", "prep_summary.json", "attributes[a1,a1b,a2,a2b].vis10_predicts_attribute_auc_cv5 min")
    add("VIS10 AUC basic max", max(basic), "f2", "prep_summary.json", "same, max")
    add("VIS10 AUC acidic min", min(acidic), "f2", "prep_summary.json", "attributes[a3,a3b] min")
    add("VIS10 AUC acidic max", max(acidic), "f2", "prep_summary.json", "attributes[a3,a3b] max")
    kc = json.load(open(os.path.join(R, "krcount_universe_summary.json"), encoding="utf-8"))
    r2 = [c["vis10_predicts_r2_cv5"] for c in kc["columns"] if c["column"] == "kr_count20"][0]
    add("VIS10 R2 K/R count +/-20", r2, "f2", "krcount_universe_summary.json", "columns[kr_count20].vis10_predicts_r2_cv5")
    uk = pd.read_csv(os.path.join(R, "universe_krcount.csv.gz"), usecols=["kr_count20"])
    add("SD of K/R count +/-20", float(np.std(uk.kr_count20.to_numpy(float))), "f1", "universe_krcount.csv.gz",
        "population SD of kr_count20 (the standardisation of F_KRcount20_z)")

    # ---- design constants (scripts) ----
    for key, v, k, src in (("positive rate", 0.2, "p0", "02_run_simulation.py: rate=0.2"),
                           ("cohort-size rate", 0.16, "p0", "02_run_simulation.py: rate=0.16"),
                           ("caliper SD", 0.2, "f1", "repo run_phase2_claims_under_detectability_control.CALIPER_SD"),
                           ("random-control draws", 20, "i", "repo ...RANDOM_DRAWS"),
                           ("bootstrap replicates", 1000, "i", "sim_engine.N_BOOT"),
                           ("bootstrap replicates check", 5000, "i", "02_run_simulation.py group boot5000"),
                           ("planted min", 0.25, "f2", "02_run_simulation.EFFECTS_GRID"),
                           ("planted max", 3.0, "f1", "02_run_simulation.EFFECTS_GRID"),
                           ("count chemistry min", 0.5, "f1", "02_run_simulation.PROXY_GRID"),
                           ("count chemistry max", 6, "i", "02_run_simulation.PROXY_GRID"),
                           ("K/R +/-5 chemistry max", 8, "i", "02_run_simulation.KR5_EXT"),
                           ("count-claim planted min", 0.5, "f1", "02_run_simulation.EFFECTS_MAIN"),
                           ("count-claim planted max", 1.5, "f1", "02_run_simulation.EFFECTS_MAIN"),
                           ("mixture b min", 0.5, "f1", "02_run_simulation.MIX_B"),
                           ("mixture b max", 2.0, "f1", "02_run_simulation.MIX_B"),
                           ("length floor", 7, "i", "sim_engine lensel7_s = -ln max(pep_len, 7)"),
                           ("tables no chemistry", 500, "i", "02_run_simulation group S0"),
                           ("tables length selection", 300, "i", "02_run_simulation group lenselcal"),
                           ("tables per cell", 200, "i", "02_run_simulation (every other cell)"),
                           ("latent prevalence", 0.30, "p0", "sim_engine.LATENT_PREVALENCE"),
                           ("attribute-tables in bootstrap check", 350, "i", "boot_check.csv n_tables")):
        add(key, v, k, "design constant: " + src, "")

    s = csv("sim_summary.csv")
    sq = lambda scen, att: s[(s.scenario == scen) & (s.attribute == att)].iloc[0]
    r0 = sq("S0", "a1_SFE006_KR")
    add("positives per table", r0.mean_n_pos, "r10", "sim_summary.csv", "S0, a1: mean_n_pos")
    add("cysteines per table", r0.mean_n_obs, "r100", "sim_summary.csv", "S0, a1: mean_n_obs")
    add("proteins per table", r0.mean_n_proteins, "r10", "sim_summary.csv", "S0, a1: mean_n_proteins")
    add("matched negatives detected", r0.mean_matched_neg_detected_frac, "p0", "sim_summary.csv",
        "S0: mean_matched_neg_detected_frac")
    add("all negatives detected", r0.mean_all_neg_detected_frac, "p0", "sim_summary.csv", "S0: mean_all_neg_detected_frac")
    add("null a1 baseline", r0.baseline_mean, "f2", "sim_summary.csv", "S0, a1: baseline_mean")
    add("null a1 matched", r0.matched_mean, "f2", "sim_summary.csv", "S0, a1: matched_mean")
    r2_ = sq("S0", "a2_SNO016_KRH")
    add("null a2 baseline", r2_.baseline_mean, "f2", "sim_summary.csv", "S0, a2: baseline_mean")
    add("null a2 matched", r2_.matched_mean, "f2", "sim_summary.csv", "S0, a2: matched_mean")
    add("null a2 baseline CI below 0", r2_.frac_baseline_ci_below_0, "p0", "sim_summary.csv",
        "S0, a2: frac_baseline_ci_below_0 (= every table)")
    pub = ["a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV", "a3_SNO006_DE3", "a3b_SFE002_E"]
    ex = [sq("S0", a).frac_matched_ci_excludes_0 for a in pub]
    add("null CI excludes 0 min", min(ex), "p0", "sim_summary.csv", "S0, published attributes: frac_matched_ci_excludes_0 min")
    add("null CI excludes 0 max", max(ex), "p0", "sim_summary.csv", "same, max")
    add("null CI excludes 0 latent", sq("S0", "a4_LATENT30").frac_matched_ci_excludes_0, "p1", "sim_summary.csv",
        "S0, a4: frac_matched_ci_excludes_0")
    fs = [sq("S0", a).verdict_final_guarded_survives + sq("S0", a).verdict_final_guarded_attenuated
          for a in ("a1_SFE006_KR", "a1b_SNO021_K", "a3b_SFE002_E")]
    add("null false survives/attenuated min", min(fs), "p0", "sim_summary.csv",
        "S0, a1/a1b/a3b: verdict_final_guarded_survives + attenuated, min")
    add("null false survives/attenuated max", max(fs), "p0", "sim_summary.csv", "same, max")
    for e in (0.5, 1.0):
        q = s[(s.scenario == "S1") & (s.attribute == "a2_SNO016_KRH") & np.isclose(s.effect_log2_or, e)].iloc[0]
        add("SNO-016 planted %.1f bcc" % e, q.verdict_final_guarded_baseline_contradicts_claim, "p0", "sim_summary.csv",
            "S1, a2, planted %.1f: verdict_final_guarded_baseline_contradicts_claim" % e)
    lat = s[(s.scenario == "S1") & (s.attribute == "a4_LATENT30")].matched_retention_vs_truth
    add("latent kept min", lat.min(), "f2", "sim_summary.csv", "S1, a4, planted 0.5-1.5: matched_retention_vs_truth min")
    add("latent kept max", lat.max(), "f2", "sim_summary.csv", "same, max")
    for e in (0.5, 1.0):
        q = s[(s.scenario == "S1") & (s.attribute == "a1_SFE006_KR") & np.isclose(s.effect_log2_or, e)].iloc[0]
        add("SFE-006 planted %.1f attenuated" % e, q.verdict_final_guarded_attenuated, "p1" if e == 1.0 else "p0",
            "sim_summary.csv", "S1, a1, planted %.1f: verdict_final_guarded_attenuated" % e)

    # ---- over-adjustment of chemistry on the attribute (loss beyond the sequence-independent attribute) ----
    o = csv("overadjustment_summary.csv")
    o = o[o.effect_log2_or.isin([0.5, 1.0, 1.5])]
    cls = {"a1_SFE006_KR": "offsets", "a1b_SNO021_K": "offsets", "a2_SNO016_KRH": "compositions",
           "a2b_PERS009_AKRV": "compositions", "a3_SNO006_DE3": "acidic", "a3b_SFE002_E": "acidic"}
    for scen, lab in (("S1", "observed detection"), ("S1permgrid", "permuted detection")):
        for c in ("offsets", "compositions", "acidic"):
            q = o[(o.planted_scenario == scen) & o.attribute.map(cls).eq(c)]
            loss = np.concatenate([1 - q.matched_genuine_part_retention_relative_to_latent.to_numpy(),
                                   1 - q.matched_increment_normalised_retention_relative_to_latent.to_numpy()])
            add("loss beyond latent, %s, %s, min" % (c, lab), loss.min(), "p0", "overadjustment_summary.csv",
                "planted_scenario %s, planted 0.5-1.5, both normalisations: 1 - *_relative_to_latent, min" % scen)
            add("loss beyond latent, %s, %s, max" % (c, lab), loss.max(), "p0", "overadjustment_summary.csv", "same, max")
    cc = csv("count_claim_summary.csv")
    cc = cc[cc.effect_log2_or.isin([0.5, 1.0, 1.5])]
    for att, lab in (("x1_KRcount20_hi", "count +/-20 claim"), ("x2_KR5_ge3", ">=3 K/R +/-5 claim")):
        for dm in ("real", "steep", "permuted_within_protein"):
            q = cc[(cc.attribute == att) & (cc.detection_model == dm)]
            v = np.concatenate([q.loss_beyond_latent_A.to_numpy(), q.loss_beyond_latent_B.to_numpy()])
            add("%s loss, %s, min" % (lab, dm), v.min(), "p0", "count_claim_summary.csv",
                "attribute %s, detection %s, planted 0.5-1.5: loss_beyond_latent_A/B min" % (att, dm))
            add("%s loss, %s, max" % (lab, dm), v.max(), "p0", "count_claim_summary.csv", "same, max")

    # ---- chemistry on the K/R count (exact calibration cells) ----
    d = csv("discriminator_summary.csv").set_index("scenario")
    for scen, lab in (("DM2real", "observed"), ("DM2steep", "steeper"), ("DM2perm", "permuted")):
        r = d.loc[scen]
        add("count chemistry b, %s" % lab, r.effect_log2_or, "f2", "discriminator_summary.csv", scen + ": effect_log2_or")
        add("count chemistry retention p50, %s" % lab, r.a1_retention_p50, "f2", "discriminator_summary.csv",
            scen + ": a1_retention_p50")
        add("count chemistry retention p2.5, %s" % lab, r.a1_retention_p2_5, "f2", "discriminator_summary.csv",
            scen + ": a1_retention_p2_5")
        add("count chemistry retention p97.5, %s" % lab, r.a1_retention_p97_5, "f2", "discriminator_summary.csv",
            scen + ": a1_retention_p97_5")
        add("count chemistry propensity AUC, %s" % lab, r.prop_auc_mean, "f2", "discriminator_summary.csv",
            scen + ": prop_auc_mean")
    r = d.loc["DM2akrv10"]
    add("A/K/R/V chemistry b", r.effect_log2_or, "f2", "discriminator_summary.csv", "DM2akrv10: effect_log2_or")
    add("A/K/R/V chemistry retention p50", r.a1_retention_p50, "f2", "discriminator_summary.csv", "DM2akrv10")
    add("A/K/R/V chemistry retention p2.5", r.a1_retention_p2_5, "f2", "discriminator_summary.csv", "DM2akrv10")
    add("A/K/R/V chemistry retention p97.5", r.a1_retention_p97_5, "f2", "discriminator_summary.csv", "DM2akrv10")
    pv = csv("proxy_verdicts.csv")
    for scen, k in (("P1ci", "p0"), ("P1permci", "p1")):
        q = pv[(pv.scenario == scen) & (pv.attribute == "a1_SFE006_KR")].iloc[0]
        add("count chemistry attenuated, %s" % scen, q.final_attenuated, k, "proxy_verdicts.csv",
            "%s (b %.1f), a1: final_attenuated" % (scen, q.b))
        add("count chemistry verdict b, %s" % scen, q.b, "f1", "proxy_verdicts.csv", scen + ": b")
    pg = csv("proxy_grid_summary.csv")
    inc = []
    for sc in ("P1grid", "P1steepgrid", "P1permgrid"):
        for a in ("a3_SNO006_DE3", "a3b_SFE002_E"):
            q = pg[(pg.scenario == sc) & (pg.attribute == a)].sort_values("b_log2_or_per_sd")
            inc.append(float((q.baseline_mean - q.baseline_mean.iloc[0]).max()))
    add("acidic increment under count chemistry, max", max(inc), "f2", "proxy_grid_summary.csv",
        "P1grid/P1steepgrid/P1permgrid, a3/a3b: max(baseline_mean(b) - baseline_mean(0))")
    k5 = csv("kr5_calibration.csv").set_index("row")
    r = k5.loc["calibrated (this item)"]
    add("K/R +/-5 chemistry b", r.b, "f2", "kr5_calibration.csv", "calibrated (this item): b")
    add("K/R +/-5 chemistry retention p50", r.retention_p50, "f2", "kr5_calibration.csv", "calibrated (this item)")
    add("K/R +/-5 chemistry retention p2.5", r.retention_p2_5, "f2", "kr5_calibration.csv", "calibrated (this item)")
    add("K/R +/-5 chemistry retention p97.5", r.retention_p97_5, "f2", "kr5_calibration.csv", "calibrated (this item)")
    rv = k5.loc["calibrated (round-4 verifier)"]
    add("K/R +/-5 chemistry b, independent", rv.b, "f2", "kr5_calibration.csv", "calibrated (round-4 verifier): b")
    add("K/R +/-5 chemistry retention p50, independent", rv.retention_p50, "f2", "kr5_calibration.csv", "verifier row")
    add("K/R +/-5 chemistry retention p2.5, independent", rv.retention_p2_5, "f2", "kr5_calibration.csv", "verifier row")
    add("K/R +/-5 chemistry retention p97.5, independent", rv.retention_p97_5, "f2", "kr5_calibration.csv",
        "verifier row")
    dm = csv("discriminator_mixture.csv")
    dm = dm[dm.scenario == "DMIX"].sort_values("effect_log2_or")
    for _, q in dm.iterrows():
        add("mixture retention p50, b %.1f" % q.effect_log2_or, q.a1_retention_p50, "f2", "discriminator_mixture.csv",
            "DMIX b %.1f: a1_retention_p50" % q.effect_log2_or)

    # ---- detectability only: selection on tryptic-peptide length at the observed baseline ----
    lc = json.load(open(os.path.join(R, "lensel_calibration.json"), encoding="utf-8"))
    add("length selection gamma", lc["gamma_star"], "f2", "lensel_calibration.json", "gamma_star")
    add("length selection gamma, independent", lc["verifier_gamma_star"], "f2", "lensel_calibration.json",
        "verifier_gamma_star (verify_r4/lensel_calib_summary.csv)")
    ls = csv("lensel_summary.csv")
    lq = lambda cell, att: ls[(ls.cell == cell) & (ls.attribute == att)].iloc[0]
    r = lq("L0cal", "a1_SFE006_KR")
    add("length selection baseline", r.baseline_mean, "f2", "lensel_summary.csv", "L0cal, a1: baseline_mean")
    add("length selection matched", r.matched_mean, "f2", "lensel_summary.csv", "L0cal, a1: matched_mean")
    add("length selection loss", r.loss_mean, "f2", "lensel_summary.csv", "L0cal, a1: loss_mean")
    add("length selection propensity AUC", r.prop_auc_before_mean, "f3", "lensel_summary.csv", "L0cal: prop_auc_before_mean")
    add("length selection retention p50", r.retention_p50, "f2", "lensel_summary.csv", "L0cal, a1: retention_p50")
    add("length selection retention p2.5", r.retention_p2_5, "f2", "lensel_summary.csv", "L0cal, a1: retention_p2_5")
    add("length selection retention p97.5", r.retention_p97_5, "f2", "lensel_summary.csv", "L0cal, a1: retention_p97_5")
    add("length selection tables with retention >= observed", r.frac_tables_retention_ge_observed * r.n_tables, "i",
        "lensel_summary.csv", "L0cal, a1: frac_tables_retention_ge_observed x n_tables (= 0)")
    add("length selection K/R count attribute baseline", lq("L0cal", "x1_KRcount20_hi").baseline_mean, "f2",
        "lensel_summary.csv", "L0cal, x1: baseline_mean")
    for att, lab in (("a1b_SNO021_K