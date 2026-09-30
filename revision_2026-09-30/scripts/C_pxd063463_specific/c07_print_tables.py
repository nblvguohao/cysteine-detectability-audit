"""Print the report tables of item C_pxd063463_specific as Markdown to stdout (nothing is written to disk).

POST HOC revision analysis (2026-09-30). Values are rounded for display only; the CSVs keep full precision.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402

R = C.RES
ARMN = {"Trypsin": "trypsin", "AspN": "AspN", "CT": "chymotrypsin", "GluC": "GluC"}
BAND = {"proximal_1_3": "1-3", "distal_6_12": "6-12"}


def f(x, k=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "-"
    return f"{x:.{k}f}"


def iv(e, lo, hi, k=2):
    if e is None or (isinstance(e, float) and not np.isfinite(e)):
        return "n/a"
    return f"{e:.{k}f} [{lo:.{k}f}, {hi:.{k}f}]"


def md(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out)


def t1():
    d = pd.read_csv(f"{R}/t1_fig2_reproduction.csv")
    rows = []
    for _, r in d.iterrows():
        q = r.quantity
        rows.append({"arm": ARMN[r.arm], "rule": r.rule, "background": r.background,
                     "quantity": BAND.get(q, q),
                     "stored (2025 FASTA)": iv(r.stored_estimate, r.stored_ci_low, r.stored_ci_high, 3)
                     if r.in_stored_figure else "not in figure",
                     "re-run (2026_03 FASTA)": iv(r.new_estimate, r.new_ci_low, r.new_ci_high, 3),
                     "delta": f(r.delta_estimate, 3) if r.in_stored_figure else "-",
                     "status stored -> new": f"{r.stored_status} -> {r.new_status}" if r.in_stored_figure else r.new_status,
                     "n pos / bg (new)": f"{int(r.new_n_positive)} / {int(r.new_n_background)}"})
    print("### T1\n" + md(pd.DataFrame(rows)))


def t1b():
    d = pd.read_csv(f"{R}/t1b_seed_summary.csv")
    d = d[((d.arm == "GluC") & (d.rule == "gluc") & (d.background == "observed")) |
          ((d.arm == "Trypsin") & (d.background == "observed") & (d.band == "proximal_1_3")) |
          ((d.rule == "trypsin") & (d.arm != "Trypsin") & (d.background == "proteome") & (d.band == "distal_6_12"))]
    rows = [{"arm": ARMN[r.arm], "rule": r.rule, "background": r.background, "band": BAND[r.band],
             "default-seed interval": f"[{r.ci_low_default_seed:.3f}, {r.ci_high_default_seed:.3f}]",
             "lower bound range (50 seeds)": f"{r.ci_low_min:.3f} to {r.ci_low_max:.3f}",
             "upper bound range (50 seeds)": f"{r.ci_high_min:.3f} to {r.ci_high_max:.3f}",
             "seeds covering 0": f"{r.share_seeds_interval_covers_zero:.2f}",
             "status shares F/W/U/P": f"{r.share_seeds_status_FAIL:.2f}/{r.share_seeds_status_WARNING:.2f}/"
                                      f"{r.share_seeds_status_UNDECIDABLE:.2f}/{r.share_seeds_status_PASS:.2f}"}
            for _, r in d.iterrows()]
    print("### T1b\n" + md(pd.DataFrame(rows)))


def t2():
    d = pd.read_csv(f"{C.OUT}/t2_specific_counts.csv")
    d = d[d.match_key == "position"]
    rows = [{"arm": ARMN[r.arm], "HydN source": r.hydn_source, "HydP sites": r.hydp_cam_sites,
             "sites identified in HydN": r.hydp_sites_identified_in_hydn,
             "CAM in both (non-specific)": r.nonspecific_cam_in_both, "strict": r.strict, "lenient": r.lenient,
             "share of sites CAM in HydN": iv(r.frac_sites_also_cam_in_hydn, r.frac_sites_also_cam_in_hydn_ci_low,
                                               r.frac_sites_also_cam_in_hydn_ci_high),
             "share of testable sites CAM in HydN": iv(r.frac_testable_sites_cam_in_hydn, r.frac_testable_ci_low,
                                                       r.frac_testable_ci_high),
             "HydP cysteines testable": f(r.frac_hydp_identified_testable, 3)} for _, r in d.iterrows()]
    print("### T2\n" + md(pd.DataFrame(rows)))
    h = pd.read_csv(f"{C.OUT}/t2_hydn_summary.csv")
    rows = [{"arm": ARMN[r.arm], "HydP identified": r.hydp_identified, "HydP sites": r.hydp_cam,
             "HydP share": iv(r.hydp_share, r.hydp_share_ci_low, r.hydp_share_ci_high, 3),
             "HydN identified": r.hydn_identified, "HydN CAM": r.hydn_cam,
             "HydN share": iv(r.hydn_share, r.hydn_share_ci_low, r.hydn_share_ci_high, 3),
             "HydP/HydN identified": f(r.ratio_identified_hydp_over_hydn, 1),
             "HydP proteins never in any HydN arm": r.hydp_proteins_absent_from_all_hydn} for _, r in h.iterrows()]
    print("### T2b\n" + md(pd.DataFrame(rows)))


def t3():
    d = pd.read_csv(f"{R}/t3_cleavage_specific.csv")
    d = d[d.handling.isin(["none", "exclude", "restrict_proteins"])]
    sets = ["cam", "strict_matched", "lenient_matched", "nonspecific_matched", "strict_pooled", "lenient_pooled",
            "nonspecific_pooled", "cam_protein_absent_pooled", "in_trypsin_hydn", "not_in_trypsin_hydn",
            "in_pooled_hydn", "not_in_pooled_hydn"]
    for arm in C.ARMS:
        rows = []
        for s in sets:
            x = d[(d.arm == arm) & (d.positive_set == s)]
            if x.empty:
                continue
            r0 = x.iloc[0]
            row = {"positive set": s, "n pos": int(r0.n_positive)}
            for rt in (["own"] if arm == "Trypsin" else ["own", "trypsin"]):
                for bg in ("proteome", "observed"):
                    for b in ("proximal_1_3", "distal_6_12"):
                        y = x[(x.rule_type == rt) & (x.background == bg) & (x.band == b)].iloc[0]
                        row[f"{rt} {bg[:4]} {BAND[b]}"] = iv(y.estimate, y.ci_low, y.ci_high) if not y.below_minimum \
                            else f"<20 pos"
            rows.append(row)
        print(f"### T3 {arm}\n" + md(pd.DataFrame(rows)))
    h = pd.read_csv(f"{R}/t3_hydn_arm_audits.csv")
    rows = []
    for (arm, rule, bg), x in h.groupby(["arm", "rule", "background"], sort=False):
        row = {"HydN arm": ARMN[arm], "rule": rule, "background": bg, "n pos / bg": f"{int(x.n_positive.iloc[0])} / "
               f"{int(x.n_background.iloc[0])}" if np.isfinite(x.n_positive.iloc[0]) else "-"}
        for b in ("proximal_1_3", "distal_6_12"):
            y = x[x.band == b].iloc[0]
            row[BAND[b]] = iv(y.estimate, y.ci_low, y.ci_high) if np.isfinite(y.estimate) else "n/a (<20 pos)"
        rows.append(row)
    print("### T3 HydN arms\n" + md(pd.DataFrame(rows)))


def t4():
    d = pd.read_csv(f"{C.OUT}/t4_coincidence_specific.csv")
    rows = []
    for (dn, den), x in d.groupby(["definition", "denominator"], sort=False):
        row = {"definition": dn, "denominator": den}
        for arm in C.ARMS:
            y = x[x.arm == arm]
            row[ARMN[arm]] = (iv(y.coincidence.iloc[0], y.ci_low.iloc[0], y.ci_high.iloc[0], 3) +
                              f" ({int(y.n_sites.iloc[0])}/{int(y.n_identified.iloc[0])})") if len(y) else "-"
        rows.append(row)
    print("### T4\n" + md(pd.DataFrame(rows)))
    c = pd.read_csv(f"{C.OUT}/t4_arm_contrasts.csv")
    rows = []
    for (dn, den), x in c.groupby(["definition", "denominator"], sort=False):
        row = {"definition": dn, "denominator": den}
        for arm in C.ARMS[1:]:
            y = x[x.arm == arm]
            row[f"{ARMN[arm]} - trypsin"] = iv(y.difference.iloc[0], y.diff_ci_low.iloc[0], y.diff_ci_high.iloc[0], 3) \
                if len(y) else "-"
        rows.append(row)
    print("### T4b\n" + md(pd.DataFrame(rows)))


def t5():
    d = pd.read_csv(f"{C.OUT}/t5_deciles.csv")
    x = d[d.definition == "cam"]
    rows = []
    for dec, y in x.groupby("decile", sort=False):
        row = {"abundance decile": dec}
        for arm in C.ARMS:
            z = y[y.arm == arm]
            row[ARMN[arm]] = f"{z.coincidence.iloc[0]:.2f} ({int(z.n_sites.iloc[0])}/{int(z.n_identified.iloc[0])})" \
                if len(z) else "-"
        rows.append(row)
    print("### T5a\n" + md(pd.DataFrame(rows)))
    c = pd.read_csv(f"{C.OUT}/t5_depth_controls.csv")
    for dn in ("cam", "lenient_pooled", "lenient_matched"):
        rows = []
        for ctl, y in c[c.definition == dn].groupby("control", sort=False):
            row = {"control": ctl}
            for arm in C.ARMS[1:]:
                z = y[y.arm == arm].iloc[0]
                row[f"{ARMN[arm]}: arm vs trypsin(adj)"] = f"{z.coincidence_arm:.3f} vs {z.coincidence_trypsin_adjusted:.3f}"
                row[f"{ARMN[arm]}: difference"] = iv(z.difference, z.diff_ci_low, z.diff_ci_high, 3) + \
                    f" ({z.share_of_raw_difference_remaining:.0%})"
            rows.append(row)
        print(f"### T5b {dn}\n" + md(pd.DataFrame(rows)))


def t5_decile_check():
    """Deciles 1-10 and the missing-abundance class: is each non-tryptic arm's share above trypsin's?"""
    d = pd.read_csv(f"{C.OUT}/t5_deciles.csv", dtype={"decile": str})
    x = d[d.definition == "cam"]
    T = x[x.arm == "Trypsin"].set_index("decile").coincidence
    rows = []
    for arm in C.ARMS[1:]:
        y = x[x.arm == arm].set_index("decile").coincidence
        num = [k for k in y.index if k != "missing"]
        rows.append({"arm": ARMN[arm], "deciles 1-10 with arm > trypsin": f"{sum(y[k] > T[k] for k in num)}/{len(num)}",
                     "missing class: arm vs trypsin": f"{y['missing']:.3f} vs {T['missing']:.3f}"})
    print("### T5a check\n" + md(pd.DataFrame(rows)))


def t2_selection():
    d = pd.read_csv(f"{C.OUT}/t2_selection_summary.csv")
    print("### T2 selection\n" + md(d))


def t8():
    c = pd.read_csv(f"{C.OUT}/t8_hydn_vs_hydp_controls.csv")
    rows = []
    for _, r in c.iterrows():
        extra = ""
        if r.control.startswith("cysteine"):
            extra = (f"; both {int(r.cam_both)}, with only {int(r.cam_only_with_hydroxylamine)}, without only "
                     f"{int(r.cam_only_without_hydroxylamine)}, p = {r.exact_binomial_p_two_sided:.2g}")
        rows.append({"arm": ARMN[r.arm], "control": r.control,
                     "without (HydN)": f(r.cam_share_hydn, 3), "with (HydP, as compared)": f(r.cam_share_hydp_compared, 3),
                     "without - with": iv(r.difference_hydn_minus_hydp, r.diff_ci_low, r.diff_ci_high, 3) + extra,
                     "n with / without": f"{int(r.n_hydp)} / {int(r.n_hydn)}",
                     "share of raw": "-" if not np.isfinite(r.share_of_raw_difference_remaining)
                     else f"{r.share_of_raw_difference_remaining:.0%}"})
    print("### T8\n" + md(pd.DataFrame(rows)))
    dd = pd.read_csv(f"{C.OUT}/t8_hydn_vs_hydp_deciles.csv", dtype={"decile": str})
    x = dd[dd.arm == "Trypsin"]
    rows = []
    for dec, y in x.groupby("decile", sort=False):
        p = y[y.hydroxylamine == "HydP"].iloc[0]
        n = y[y.hydroxylamine == "HydN"].iloc[0]
        rows.append({"decile": dec, "with: share (k/n)": f"{p.cam_share:.2f} ({int(p.n_cam)}/{int(p.n_identified)})",
                     "without: share (k/n)": f"{n.cam_share:.2f} ({int(n.n_cam)}/{int(n.n_identified)})"})
    print("### T8 deciles (trypsin)\n" + md(pd.DataFrame(rows)))
    b = pd.read_csv(f"{C.OUT}/t8_bin_distribution.csv", dtype={"bin": str})
    b = b[b.arm == "Trypsin"].pivot(index="bin", columns="hydroxylamine", values="share_of_arm")
    print("### T8 bin distribution (trypsin)\n" + b.round(3).to_string())
    e = pd.read_csv(f"{C.OUT}/t8_evidence_depth.csv")
    print("### T8 evidence rows\n" + md(e[["arm", "hydroxylamine", "evidence_rows_in_arm_all_peptides",
                                           "identified_cysteines_stored", "evidence_rows_per_identified_cysteine"]]
                                        .round(1)))


def t9():
    s = pd.read_csv(f"{R}/t9_base_agreement.csv")
    print("### T9\n" + md(s[["table", "quantity", "class", "n_compared", "max_abs_difference", "row_at_max"]]
                         .assign(max_abs_difference=lambda z: z.max_abs_difference.round(5))))


def t11():
    """Revision round 2: stratum-matched split sets and candidate sets, flag shares, Mantel-Haenszel summary."""
    d = pd.read_csv(f"{R}/t11_stratum_matched_cleavage.csv")
    t3u = pd.read_csv(f"{R}/t3_cleavage_specific.csv")
    for design in ("split_stratum_matched", "candidate_in_stratum"):
        rows = []
        for (arm, s), x in d[d.design == design].groupby(["arm", "positive_set"], sort=False):
            r0 = x.iloc[0]
            row = {"arm": ARMN[arm], "set": s, "n pos / bg (proteome)": f"{int(r0.n_positive)} / {int(r0.n_background)}"}
            for rt in (["own"] if arm == "Trypsin" else ["own", "trypsin"]):
                for bg in ("proteome", "observed"):
                    for b in ("proximal_1_3", "distal_6_12"):
                        y = x[(x.rule_type == rt) & (x.background == bg) & (x.band == b)].iloc[0]
                        row[f"{rt} {bg[:4]} {BAND[b]}"] = iv(y.estimate, y.ci_low, y.ci_high) if not y.below_minimum \
                            else "<20"
            if design == "split_stratum_matched":
                u = t3u[(t3u.arm == arm) & (t3u.positive_set == s) & (t3u.handling == "exclude")
                        & (t3u.background == "proteome") & (t3u.band == "distal_6_12")]
                row["unmatched bg, distal own / trypsin"] = " / ".join(
                    f(u[u.rule_type == rt].estimate.iloc[0]) for rt in (["own"] if arm == "Trypsin" else ["own", "trypsin"]))
            rows.append(row)
        print(f"### T11 {design}\n" + md(pd.DataFrame(rows)))
    fl = pd.read_csv(f"{R}/t11_flag_shares_by_stratum.csv")
    x = fl[(fl.band == "distal_6_12") & (fl.stratum != "all")]
    x = x.assign(v=[f"{s:.2f} ({n})" for s, n in zip(x.share_flagged, x.n)])
    print("### T11 flag shares, distal band\n" + x.pivot_table(index=["arm", "rule"], columns=["stratum", "class"],
                                                               values="v", aggfunc="first").to_string())
    m = pd.read_csv(f"{R}/t11_mh_stratified.csv")
    m = m.assign(v=[iv(e, lo, hi) for e, lo, hi in zip(m.estimate, m.ci_low, m.ci_high)])
    print("### T11 Mantel-Haenszel (stored CAM sets)\n" +
          m.pivot_table(index=["arm", "rule_type", "background", "band"], columns="strata", values="v",
                        aggfunc="first").to_string())


def t12():
    c = pd.read_csv(f"{R}/t12_coverage_proxy.csv")
    c = c[c.base == "unfiltered"]
    rows = [{"pair": r.pair, "stratum": r.stratum, "cysteines": r.n_cysteines, "CAM A only / B only":
             f"{r.cam_A_only} / {r.cam_B_only}", "binomial p": f"{r.binomial_p_two_sided:.2g}",
             "median ratio": f(r.median_ratio), "B - A": iv(r.difference_B_minus_A, r.ci_low, r.ci_high, 3)
             if np.isfinite(r.difference_B_minus_A) else "-"} for _, r in c.iterrows()]
    print("### T12 coverage proxy (unfiltered)\n" + md(pd.DataFrame(rows)))
    t = pd.read_csv(f"{R}/t12_coverage_tertiles.csv")
    print("### T12 tertiles\n" + md(t[t.base == "unfiltered"].round(3)))


if __name__ == "__main__":
    t1(); t1b(); t2(); t2_selection(); t3(); t4(); t5(); t5_decile_check(); t8(); t9(); t11(); t12()
