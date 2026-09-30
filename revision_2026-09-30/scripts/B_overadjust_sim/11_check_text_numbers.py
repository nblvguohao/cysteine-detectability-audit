# -*- coding: utf-8 -*-
"""B_overadjust_sim step 11 (POST HOC, revision round 4): traceability check of the numbers used in the proposed
manuscript text and in Supplemental Note 15.

Every claim below is recomputed from a stored file in results/B_overadjust_sim/ and compared with the value written in
the Note (results/B_overadjust_sim/Supplemental_Note_15_matched_background_genuine_chemistry.md) and in the report's
proposed text. The check also greps both texts for phrases the round-4 verification refuted, so that they cannot come
back unnoticed. Writes results/B_overadjust_sim/text_number_check.csv and exits non-zero on any failure.
"""
from __future__ import annotations

import io
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

R = cb.RESULTS
NOTE = os.path.join(R, "Supplemental_Note_15_matched_background_genuine_chemistry.md")
REPORT = cb.REPORT
# phrases withdrawn in round 4; they may appear in the report's history sections, never in the Note or in the
# proposed main text (section 7 of the report)
FORBIDDEN = ["Neither reading is cheap", "detectability alone gives that count", "5- to 10-fold", "5--10x",
             "do not overlap", "non-overlapping", "would separate artifact from chemistry directly",
             "three simulated 95"]


def rd(p):
    with io.open(p, encoding="utf-8") as fh:
        return fh.read()


def band(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.median(v)), float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))


def main():
    rows = []

    def chk(name, value, expected, tol, source):
        ok = abs(float(value) - float(expected)) <= tol
        rows.append({"quantity": name, "recomputed": round(float(value), 6), "in_text": expected, "tol": tol,
                     "ok": ok, "source_file": source})

    prep = json.load(open(cb.PREP_SUMMARY, encoding="utf-8"))
    auc = {a["attribute"]: a["vis10_predicts_attribute_auc_cv5"] for a in prep["attributes"]}
    basic = [auc[a] for a in ("a1_SFE006_KR", "a1b_SNO021_K", "a2_SNO016_KRH", "a2b_PERS009_AKRV")]
    acid = [auc[a] for a in ("a3_SNO006_DE3", "a3b_SFE002_E")]
    chk("VIS10 AUC for basic attributes, min", min(basic), 0.66, 0.005, "prep_summary.json")
    chk("VIS10 AUC for basic attributes, max", max(basic), 0.84, 0.005, "prep_summary.json")
    chk("VIS10 AUC for acidic attributes, min", min(acid), 0.59, 0.005, "prep_summary.json")
    chk("VIS10 AUC for acidic attributes, max", max(acid), 0.61, 0.005, "prep_summary.json")
    kc = json.load(open(os.path.join(R, "krcount_universe_summary.json"), encoding="utf-8"))
    r2 = {c["column"]: c.get("vis10_predicts_r2_cv5") for c in kc["columns"]}
    chk("VIS10 R2 for the K/R count within +/-20", r2["kr_count20"], 0.44, 0.005, "krcount_universe_summary.json")

    o = pd.read_csv(os.path.join(R, "overadjustment_summary.csv"))
    o = o[(o.detection_model == "real") & (o.planted_scenario == "S1") & o.effect_log2_or.isin([0.5, 1.0, 1.5])]
    lo = pd.DataFrame({"attribute": pd.concat([o.attribute, o.attribute], ignore_index=True),
                       "loss": pd.concat([1 - o.matched_genuine_part_retention_relative_to_latent,
                                          1 - o.matched_increment_normalised_retention_relative_to_latent],
                                         ignore_index=True)})
    off = lo.loss[lo.attribute.isin(["a1_SFE006_KR", "a1b_SNO021_K"])]
    comp = lo.loss[lo.attribute.isin(["a2_SNO016_KRH", "a2b_PERS009_AKRV"])]
    chk("loss beyond latent, named offsets, min (real)", 100 * off.min(), 3, 1.0, "overadjustment_summary.csv")
    chk("loss beyond latent, named offsets, max (real)", 100 * off.max(), 11, 1.0, "overadjustment_summary.csv")
    chk("loss beyond latent, composition, min (real)", 100 * comp.min(), 12, 1.0, "overadjustment_summary.csv")
    chk("loss beyond latent, composition, max (real)", 100 * comp.max(), 23, 1.0, "overadjustment_summary.csv")

    d = pd.read_csv(os.path.join(R, "discriminator_summary.csv")).set_index("scenario")
    ret = [d.loc[s, "a1_retention_p50"] for s in ("DM2real", "DM2steep", "DM2perm")]
    chk("K/R-count chemistry, SFE-006 retention, min", min(ret), 0.37, 0.006, "discriminator_summary.csv")
    chk("K/R-count chemistry, SFE-006 retention, max", max(ret), 0.47, 0.006, "discriminator_summary.csv")

    obs = pd.read_csv(os.path.join(R, "observed_claims.csv")).set_index("claim")
    ob = obs.loc["SFE-006 (Yang 2014, own data)"]
    chk("observed SFE-006 baseline", ob.baseline, 1.3840, 5e-4, "observed_claims.csv")
    chk("observed SFE-006 matched", ob.matched, 0.5950, 5e-4, "observed_claims.csv")
    chk("observed SFE-006 retention", ob.retention_matched_over_baseline, 0.43, 0.005, "observed_claims.csv")
    chk("observed SFE-006 loss", ob.baseline - ob.matched, 0.789, 0.001, "observed_claims.csv")
    chk("observed SFE-006 propensity AUC", ob.propensity_auc, 0.7616, 5e-4, "observed_claims.csv")

    ls = pd.read_csv(os.path.join(R, "lensel_summary.csv"))
    l0 = ls[(ls.cell == "L0cal") & (ls.attribute == "a1_SFE006_KR")].iloc[0]
    chk("length selection, SFE-006 baseline", l0.baseline_mean, 1.381, 0.002, "lensel_summary.csv")
    chk("length selection, SFE-006 matched", l0.matched_mean, 0.143, 0.002, "lensel_summary.csv")
    chk("length selection, SFE-006 loss", l0.loss_mean, 1.237, 0.002, "lensel_summary.csv")
    chk("length selection, propensity AUC", l0.prop_auc_before_mean, 0.765, 0.002, "lensel_summary.csv")
    chk("length selection, retention median", l0.retention_p50, 0.11, 0.005, "lensel_summary.csv")
    chk("length selection, tables reaching the observed retention", l0.frac_tables_retention_ge_observed, 0.0, 0.0,
        "lensel_summary.csv")
    chk("length selection, K/R +/-20 baseline (no chemistry)", l0.x1_KRcount20_hi_baseline_mean, 1.179, 0.002,
        "lensel_summary.csv")
    l16 = ls[(ls.cell == "L0cal16") & (ls.attribute == "a1_SFE006_KR")].iloc[0]
    chk("length selection at rate 0.16, baseline", l16.baseline_mean, 1.421, 0.002, "lensel_summary.csv")
    chk("length selection at rate 0.16, loss", l16.loss_mean, 1.281, 0.002, "lensel_summary.csv")
    for a, v in (("a1b_SNO021_K", 0.555), ("a2b_PERS009_AKRV", 0.479), ("a3b_SFE002_E", 0.153)):
        q = ls[(ls.cell == "L0cal") & (ls.attribute == a)].iloc[0]
        chk("length selection, baseline of " + a, q.baseline_mean, v, 0.002, "lensel_summary.csv")
    cal = json.load(open(os.path.join(R, "lensel_calibration.json"), encoding="utf-8"))
    chk("length selection gamma", cal["gamma_star"], 3.4978, 5e-4, "lensel_calibration.json")
    chk("length selection gamma, verifier", cal["verifier_gamma_star"], 3.4804, 5e-4, "lensel_calibration.json")

    nb = pd.read_csv(os.path.join(R, "negb_detmatch_comparison.csv")).set_index("cell")
    chk("observed-negative / baseline, length-selection artifact", nb.loc["L0cal", "neg_b_over_baseline"], 0.885,
        0.002, "negb_detmatch_comparison.csv")
    chk("observed-negative / baseline, K/R-count chemistry", nb.loc["DM2real", "neg_b_over_baseline"], 0.920, 0.002,
        "negb_detmatch_comparison.csv")
    chk("observed-negative / baseline, real-detection null", nb.loc["D0real", "neg_b_over_baseline"], 0.135, 0.002,
        "negb_detmatch_comparison.csv")
    chk("detection-model matched / baseline, length selection", nb.loc["L0cal", "detmatch_over_baseline"], 0.733,
        0.002, "negb_detmatch_comparison.csv")
    chk("detection-model matched / baseline, K/R-count chemistry", nb.loc["DM2real", "detmatch_over_baseline"], 0.796,
        0.002, "negb_detmatch_comparison.csv")

    k = pd.read_csv(os.path.join(R, "kr5_calibration.csv")).set_index("row")
    chk("K/R +/-5 chemistry, calibrated b", k.loc["calibrated (this item)", "b"], 6.09, 0.005, "kr5_calibration.csv")
    chk("K/R +/-5 chemistry, calibrated b (verifier)", k.loc["calibrated (round-4 verifier)", "b"], 6.06, 0.005,
        "kr5_calibration.csv")
    chk("K/R +/-5 chemistry, retention median", k.loc["calibrated (this item)", "retention_p50"], 0.32, 0.005,
        "kr5_calibration.csv")
    chk("K/R +/-5 chemistry, retention 2.5%", k.loc["calibrated (this item)", "retention_p2_5"], 0.18, 0.005,
        "kr5_calibration.csv")
    chk("K/R +/-5 chemistry, retention 97.5%", k.loc["calibrated (this item)", "retention_p97_5"], 0.43, 0.005,
        "kr5_calibration.csv")

    g = pd.read_csv(os.path.join(R, "discriminator_overlap_groups.csv")).set_index("group")
    pairs = [("detectability only, at the observed baseline (length selection)", 0.70, 1.02),
             ("detectability only, below the observed baseline", 0.45, 1.95),
             ("chemistry on the K/R count within +/-20", 1.82, 2.72),
             ("chemistry on the SFE-006 offsets", 0.20, 0.53),
             ("chemistry on the K/R count within +/-5", 0.96, 1.17),
             ("chemistry on the A/K/R/V count within +/-10", 0.71, 0.95),
             ("mixture: chemistry on the K/R +/-20 count and on the offsets", 0.58, 2.06)]
    for grp, lo_, hi_ in pairs:
        chk("ratio range low, " + grp, g.loc[grp, "x1_union_lo"], lo_, 0.006, "discriminator_overlap_groups.csv")
        chk("ratio range high, " + grp, g.loc[grp, "x1_union_hi"], hi_, 0.006, "discriminator_overlap_groups.csv")
    # the two ranges the round-3 text called separate must overlap
    a_lo, a_hi = g.loc["detectability only, at the observed baseline (length selection)", ["x1_union_lo", "x1_union_hi"]]
    for other in ("chemistry on the A/K/R/V count within +/-10", "chemistry on the K/R count within +/-5"):
        b_lo, b_hi = g.loc[other, ["x1_union_lo", "x1_union_hi"]]
        rows.append({"quantity": "ranges overlap: length selection vs " + other, "recomputed": "yes",
                     "in_text": "yes", "tol": 0, "ok": bool(a_lo <= b_hi and b_lo <= a_hi),
                     "source_file": "discriminator_overlap_groups.csv"})

    pv = pd.read_csv(os.path.join(R, "proxy_verdicts.csv"))
    for scen, v in (("P1ci", 0.84), ("P1permci", 0.925)):
        q = pv[(pv.scenario == scen) & (pv.attribute == "a1_SFE006_KR")].iloc[0]
        chk("attenuated share at grid b 2.5, " + scen, q.final_attenuated, v, 0.001, "proxy_verdicts.csv")
        chk("b of the interval run, " + scen, q.b, 2.5, 1e-9, "proxy_verdicts.csv")

    note, report = rd(NOTE), rd(REPORT)
    i7 = report.index("## 7. Proposed manuscript text")
    i8 = report.index("## 8. Supplemental Note 15")
    proposed = report[i7:i8]
    for phrase in FORBIDDEN:
        for label, text in (("Note", note), ("proposed main text", proposed)):
            rows.append({"quantity": "withdrawn phrase absent from the %s: %r" % (label, phrase), "recomputed": "-",
                         "in_text": "-", "tol": 0, "ok": phrase.lower() not in text.lower(), "source_file": "-"})
    rows.append({"quantity": "Note title", "recomputed": "-", "in_text": "-", "tol": 0,
                 "ok": note.startswith("# Supplemental Note 15. Does the matched background remove genuine "
                                       "chemistry? (post hoc)"), "source_file": "-"})
    for bad in ("M1", "M2", "round 3", "round 4", "S0", "NEG_A", "TODO", "TBD", "XXX"):
        rows.append({"quantity": "internal jargon absent from the Note: %r" % bad, "recomputed": "-", "in_text": "-",
                     "tol": 0, "ok": bad not in note, "source_file": "-"})

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(R, "text_number_check.csv"), index=False)
    pd.set_option("display.width", 220)
    pd.set_option("display.max_rows", 200)
    print(out.to_string(index=False))
    bad = out[~out.ok.astype(bool)]
    print("checks:", len(out), "failed:", len(bad))
    sys.exit(1 if len(bad) else 0)


if __name__ == "__main__":
    main()
