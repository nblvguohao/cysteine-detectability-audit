"""H_artifact5_docs step 2 - POST HOC: what the ORIGINAL PXD015307 deposit reported for +32 (2026-09-30).

POST HOC revision analysis in response to a pre-submission review; not registered, not pre-specified.

Reads the repository's stored extraction of the deposit's own consensus files
(results/search_space_mass_accuracy.csv and _audit.json, produced by
scripts/audit_search_space_and_mass_accuracy.py from Fig-1D.msf and the four Fig4-L .pdResult files;
the consensus files themselves are not in the local copy) and re-reads the four Sulfide-assigned PSMs.
`observed_ppm` there is the deposit's own TargetPsms.DeltaMassInPPM (not recomputed); `neutral_mass` is a
recomputed approximation (residues + water + NEM x (marked Cys - 1) + sulfide) and only enters the
separation in ppm.

The stored audit called a Sulfide PSM 'resolvable' only if |observed ppm| < separation_ppm / 2, i.e. it
measured the error from 0 ppm although the same audit reports a systematic offset of about +5 ppm in
these files. The later re-search defined the rule relative to the file's own baseline
(|ppm - baseline median| < separation_ppm / 2) and states that this is the same rule; it is not. Here
both are applied.

Revision after verification (2026-09-30, round 2): the verifier showed that the four baseline definitions
of round 1 were all dominated by B-chain PSMs, and that the deposit's two non-sulfide A-chain PSMs sit
about 6 ppm below its B-chain PSMs. The sensitivity analysis therefore adds chain-matched and
chain-and-charge-matched baselines, reports the A- versus B-chain offset (with an m/z-matched comparison),
the tolerance-window selection effect under every baseline, and a per-PSM summary of whether the reading
survives every baseline. Nothing is refitted; no random numbers are drawn.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5_lib as L  # noqa: E402

CONS_CSV = os.path.join(L.REPO, "results", "search_space_mass_accuracy.csv")
CONS_AUDIT = os.path.join(L.REPO, "results", "search_space_mass_accuracy_audit.json")
CONS_SCRIPT = os.path.join(L.REPO, "scripts", "audit_search_space_and_mass_accuracy.py")
PREREG = os.path.join(L.W, "inputs", "repo_reports", "pxd015307_research_preregistration_2026-09-17.json")
VERIFIER_CSV = os.path.join(L.OUT, "verify_r1", "verify_r1_deposit_sulfide_psms_baseline_by_peptide.csv")
PROTON = 1.007276
WINDOW_EDGE_PPM = 10.0      # hypothetical: the deposit's tolerance is not recorded


def reading(zs, zo):
    if abs(zs) < 2 and abs(zo) >= 2:
        return "sulfide only"
    if abs(zs) < 2 and abs(zo) < 2:
        return "either"
    if abs(zo) < 2:
        return "dioxidation only"
    return "neither"


def main():
    df = pd.read_csv(L.register(CONS_CSV, "stored extraction of the deposit's targeted-arm PSMs"),
                     encoding="utf-8-sig")
    audit = L.read_json(CONS_AUDIT, "stored audit of the deposit's consensus files")
    script_txt = L.read_text(CONS_SCRIPT, "repository script that produced the consensus audit")
    prereg = L.read_json(PREREG, "re-search registration JSON (with amendments)")

    rule_in_script = "resolvable = (abs(obs) < ppm_sep / 2)" in script_txt
    rule_claim = prereg["primary_statistic_and_reading_rule"]["resolvable"]

    df["has_sulfide"] = df.has_sulfide.astype(str).str.lower() == "true"
    df["observed_ppm"] = df.observed_ppm.astype(float)
    df["chain"] = np.where(df.sequence.str.startswith("GIVEQ"), "A", "B")
    df["mz_approx"] = (df.neutral_mass + df.charge * PROTON) / df.charge
    ns = df[~df.has_sulfide]

    # ---------- baselines ----------
    def summ(v):
        v = np.asarray(v, dtype=float)
        return {"n": int(v.size), "median": float(np.median(v)) if v.size else np.nan,
                "mean": float(np.mean(v)) if v.size else np.nan,
                "sd": float(np.std(v, ddof=1)) if v.size > 1 else np.nan,
                "mad_sd": float(stats.median_abs_deviation(v, scale="normal")) if v.size > 1 else np.nan,
                "min": float(v.min()) if v.size else np.nan, "max": float(v.max()) if v.size else np.nan}

    base = []
    for arm, sub in df.groupby("arm"):
        s = summ(sub[~sub.has_sulfide].observed_ppm)
        base.append({"arm": arm, "subset": "all non-sulfide PSMs", "n_psms": int(len(sub)),
                     "n_sulfide_psms": int(sub.has_sulfide.sum()), **{f"{k}_ppm" if k != "n" else "n_baseline": v
                                                                       for k, v in s.items()}})
    s = summ(ns.observed_ppm)
    base.append({"arm": "ALL_targeted_arms", "subset": "all non-sulfide PSMs", "n_psms": int(len(df)),
                 "n_sulfide_psms": int(df.has_sulfide.sum()),
                 **{f"{k}_ppm" if k != "n" else "n_baseline": v for k, v in s.items()}})
    for chain in ("A", "B"):
        s = summ(ns[ns.chain == chain].observed_ppm)
        base.append({"arm": "ALL_targeted_arms", "subset": f"{chain}-chain non-sulfide PSMs",
                     "n_psms": int((df.chain == chain).sum()),
                     "n_sulfide_psms": int((df.has_sulfide & (df.chain == chain)).sum()),
                     **{f"{k}_ppm" if k != "n" else "n_baseline": v for k, v in s.items()}})
    base = pd.DataFrame(base)
    base["max_mass_da_separable_at_2sd"] = L.SEPARATION_DA * 1e6 / (4 * base.sd_ppm)
    base = base.rename(columns={"median_ppm": "median_ppm_non_sulfide", "mean_ppm": "mean_ppm_non_sulfide",
                                "sd_ppm": "sd_ppm_non_sulfide", "mad_sd_ppm": "mad_sd_ppm_non_sulfide",
                                "n_baseline": "n_non_sulfide_psms"})
    bmap = base[base.subset == "all non-sulfide PSMs"].set_index("arm")

    # A- versus B-chain offset among the deposit's own non-sulfide PSMs
    a_err = ns.loc[ns.chain == "A", "observed_ppm"].to_numpy()
    b_err = ns.loc[ns.chain == "B", "observed_ppm"].to_numpy()
    a_mz = ns.loc[ns.chain == "A", "mz_approx"].to_numpy()
    b_near = ns[(ns.chain == "B") & (ns.mz_approx >= a_mz.min() - 60) & (ns.mz_approx <= a_mz.max() + 60)]
    by_acc = ns.groupby(["arm", "chain", "charge"]).agg(
        n=("observed_ppm", "size"), median_ppm=("observed_ppm", "median"), mean_ppm=("observed_ppm", "mean"),
        sd_ppm=("observed_ppm", "std"), mz_approx=("mz_approx", "median")).reset_index()
    L.write_csv(by_acc, "deposit_nonsulfide_by_arm_chain_charge.csv")
    chain_offset = {
        "A_chain_non_sulfide_ppm": a_err.round(4).tolist(),
        "A_chain_arm_and_charge": ns.loc[ns.chain == "A", ["arm", "charge"]].astype(str).agg(" z=".join, axis=1).tolist(),
        "A_chain_mz_approx": a_mz.round(2).tolist(),
        "B_chain_n": int(b_err.size), "B_chain_median_ppm": float(np.median(b_err)),
        "B_chain_sd_ppm": float(np.std(b_err, ddof=1)),
        "A_minus_B_median_ppm": float(np.median(a_err) - np.median(b_err)),
        "mannwhitney_p_two_sided": float(stats.mannwhitneyu(a_err, b_err, alternative="two-sided").pvalue),
        "B_chain_within_60_mz_of_A_chain_ions": {
            "n": int(len(b_near)), "mz_range": [float(b_near.mz_approx.min()), float(b_near.mz_approx.max())],
            "median_ppm": float(b_near.observed_ppm.median()),
            "range_ppm": [float(b_near.observed_ppm.min()), float(b_near.observed_ppm.max())],
            "group_medians_by_arm_charge": {f"{a} z={z}": round(float(v), 3) for (a, z), v in
                                            b_near.groupby(["arm", "charge"]).observed_ppm.median().items()}},
        "note": ("the two A-chain PSMs are from the NaHS file (z = 3); the offset is not explained by m/z or charge "
                 "(B-chain ions at similar m/z sit near the B-chain median); its cause is not established, and "
                 "with n = 2 it cannot be characterised further"),
    }

    # ---------- the four sulfide PSMs ----------
    rows = []
    pooled = bmap.loc["ALL_targeted_arms"]
    for _, r in df[df.has_sulfide].iterrows():
        obs, sep = float(r.observed_ppm), float(r.separation_ppm)
        b = bmap.loc[r.arm]
        med, sd = b.median_ppm_non_sulfide, b.sd_ppm_non_sulfide
        mu_s, mu_o = med, med + sep
        z_s, z_o = (obs - mu_s) / sd, (obs - mu_o) / sd
        rows.append({
            "arm": r.arm, "chain": r.chain, "sequence": r.sequence, "charge": int(r.charge), "placed": r.placed,
            "neutral_mass_da_as_stored": float(r.neutral_mass), "mz_approx": float(r.mz_approx),
            "observed_ppm": obs, "separation_ppm": sep,
            "rule_stored_audit_abs_ppm_lt_half_sep": abs(obs) < sep / 2,
            "arm_baseline_median_ppm": med, "arm_baseline_sd_ppm": sd,
            "arm_baseline_chain_composition": "A {} / B {}".format(
                int(((ns.arm == r.arm) & (ns.chain == "A")).sum()), int(((ns.arm == r.arm) & (ns.chain == "B")).sum())),
            "deviation_from_arm_baseline_ppm": obs - med,
            "rule_baseline_corrected_arm": abs(obs - med) < sep / 2,
            "deviation_from_pooled_baseline_ppm": obs - pooled.median_ppm_non_sulfide,
            "rule_baseline_corrected_pooled": abs(obs - pooled.median_ppm_non_sulfide) < sep / 2,
            "expected_ppm_if_sulfide": mu_s, "expected_ppm_if_dioxidation": mu_o,
            "z_if_sulfide": z_s, "z_if_dioxidation": z_o,
            "gaussian_llr_sulfide_vs_dioxidation": float(stats.norm.logpdf(obs, mu_s, sd) - stats.norm.logpdf(obs, mu_o, sd)),
            "nearer_hypothesis": "sulfide" if abs(obs - mu_s) < abs(obs - mu_o) else "dioxidation",
            "reading_file_baseline": reading(z_s, z_o),
            "confidence_level_in_deposit": "not stored (consensus files not in the local copy)",
        })
    ps = pd.DataFrame(rows)

    # ---------- sensitivity to the baseline definition (now including chain-matched baselines) ----------
    def defs_for(r):
        f = ns[ns.arm == r.arm].observed_ppm
        fsd = float(f.std(ddof=1))
        same_file_chain = ns[(ns.arm == r.arm) & (ns.chain == r.chain)].observed_ppm
        same_chain = ns[ns.chain == r.chain].observed_ppm
        same_chain_z = ns[(ns.chain == r.chain) & (ns.charge == r.charge)].observed_ppm
        out = [
            ("D1 file, all non-sulfide PSMs (median, SD)", "B-dominated", f, float(f.median()), fsd, "SD of the same PSMs"),
            ("D2 file, all non-sulfide PSMs (median, MAD-scaled SD)", "B-dominated", f, float(f.median()),
             float(stats.median_abs_deviation(f, scale="normal")), "MAD-scaled SD of the same PSMs"),
            ("D3 four files pooled, all non-sulfide PSMs (median, SD)", "B-dominated", ns.observed_ppm,
             float(ns.observed_ppm.median()), float(ns.observed_ppm.std(ddof=1)), "SD of the same PSMs"),
            ("D4 file, all non-sulfide PSMs (mean, SD)", "B-dominated", f, float(f.mean()), fsd, "SD of the same PSMs"),
            ("D5 file, same chain", "chain-matched", same_file_chain,
             float(same_file_chain.median()) if len(same_file_chain) else np.nan,
             float(same_file_chain.std(ddof=1)) if len(same_file_chain) > 2 else np.nan, "SD of the same PSMs"),
            ("D6a four files, same chain (SD of those PSMs)", "chain-matched", same_chain,
             float(same_chain.median()) if len(same_chain) else np.nan,
             float(same_chain.std(ddof=1)) if len(same_chain) > 1 else np.nan,
             "SD of the same PSMs" + (" (n = 2)" if len(same_chain) == 2 else "")),
            ("D6b four files, same chain (file SD as stand-in)", "chain-matched", same_chain,
             float(same_chain.median()) if len(same_chain) else np.nan, fsd, "SD of all non-sulfide PSMs of the file"),
            ("D7 four files, same chain and charge", "chain-and-charge-matched", same_chain_z,
             float(same_chain_z.median()) if len(same_chain_z) else np.nan,
             float(same_chain_z.std(ddof=1)) if len(same_chain_z) > 2 else np.nan, "SD of the same PSMs"),
        ]
        return out

    sens = []
    for _, r in ps.iterrows():
        for label, family, subset, c, s_, s_src in defs_for(r):
            n = int(len(subset))
            rec = {"arm": r.arm, "chain": r.chain, "charge": r.charge, "placed": r.placed,
                   "observed_ppm": r.observed_ppm, "separation_ppm": r.separation_ppm,
                   "baseline_definition": label, "family": family, "n_baseline_psms": n,
                   "baseline_chains": "A {} / B {}".format(int((ns.loc[subset.index, "chain"] == "A").sum()),
                                                           int((ns.loc[subset.index, "chain"] == "B").sum())) if n else "",
                   "centre_ppm": c, "spread_ppm": s_, "spread_source": s_src}
            if n == 0 or not np.isfinite(c) or not np.isfinite(s_):
                rec.update({"reading": "no baseline (no PSM of this kind)" if n == 0 else "no spread (n < 3)"})
                sens.append(rec)
                continue
            dev = r.observed_ppm - c
            zs, zo = dev / s_, (dev - r.separation_ppm) / s_
            e_diox = c + r.separation_ppm
            rec.update({
                "deviation_ppm": dev, "z_if_sulfide": zs, "z_if_dioxidation": zo,
                "nearer_hypothesis": "sulfide" if abs(dev) < abs(dev - r.separation_ppm) else "dioxidation",
                "half_separation_rule_resolvable_as_sulfide": abs(dev) < r.separation_ppm / 2,
                "half_separation_rule_resolvable_as_dioxidation": abs(dev - r.separation_ppm) < r.separation_ppm / 2,
                "reading": reading(zs, zo),
                "expected_ppm_if_dioxidised_but_reported_as_sulfide": e_diox,
                "that_expectation_inside_plus10ppm_edge": bool(e_diox <= WINDOW_EDGE_PPM),
                "p_error_inside_plus10ppm_edge_if_dioxidised": float(stats.norm.cdf(WINDOW_EDGE_PPM, e_diox, s_)),
            })
            sens.append(rec)
    sens = pd.DataFrame(sens)
    L.write_csv(sens, "deposit_sulfide_psms_baseline_sensitivity.csv")

    robust = []
    for (arm, obs), g in sens.groupby(["arm", "observed_ppm"], sort=False):
        valid = g[g.reading.isin(["sulfide only", "dioxidation only", "either", "neither"])]
        readings = valid.reading.unique().tolist()
        bdom = valid[valid.family == "B-dominated"].reading.unique().tolist()
        cm = valid[valid.family != "B-dominated"].reading.unique().tolist()
        label = (f"{readings[0]} under every baseline" if len(readings) == 1 else
                 "baseline-dependent: " + "; ".join(
                     f"{rd} ({', '.join(x.split(' ')[0] for x in valid[valid.reading == rd].baseline_definition)})"
                     for rd in readings))
        first = g.iloc[0]
        robust.append({"arm": arm, "chain": first.chain, "charge": first.charge, "placed": first.placed,
                       "observed_ppm": obs, "n_baseline_definitions_with_a_value": int(len(valid)),
                       "readings_B_dominated_baselines": "; ".join(bdom),
                       "readings_chain_matched_baselines": "; ".join(cm) if cm else "no chain-matched PSMs",
                       "robust_summary": label,
                       "range_z_if_dioxidation": [float(valid.z_if_dioxidation.min()), float(valid.z_if_dioxidation.max())],
                       "expected_ppm_if_dioxidised_range": [float(valid.expected_ppm_if_dioxidised_but_reported_as_sulfide.min()),
                                                            float(valid.expected_ppm_if_dioxidised_but_reported_as_sulfide.max())]})
    robust = pd.DataFrame(robust)
    L.write_csv(robust, "deposit_sulfide_psms_robust_reading.csv")

    # ---------- cross-check against the verifier's independent table (if present) ----------
    xcheck = {"available": os.path.exists(VERIFIER_CSV)}
    if xcheck["available"]:
        v = pd.read_csv(L.register(VERIFIER_CSV, "verifier round-1 table (cross-check only)"))
        mapping = {"file, all non-sulfide (implementer)": "D1 file, all non-sulfide PSMs (median, SD)",
                   "file, same chain only": "D5 file, same chain",
                   "all files, same chain only": None,
                   "all files, same chain and charge": "D7 four files, same chain and charge"}
        diffs = []
        for _, vr in v.iterrows():
            lab = mapping.get(vr.baseline)
            if lab is None:
                # the verifier used the file SD as stand-in only when n < 3 (A chain), else the chain SD
                lab = ("D6b four files, same chain (file SD as stand-in)" if vr.chain == "A"
                       else "D6a four files, same chain (SD of those PSMs)")
            m = sens[(sens.arm == vr.arm) & (np.isclose(sens.observed_ppm, vr.observed_ppm)) & (sens.baseline_definition == lab)]
            if len(m) != 1:
                continue
            m = m.iloc[0]
            vr_read = {"both": "either"}.get(vr.reading_2sd, vr.reading_2sd)
            same = (vr_read == m.reading) or (vr.reading_2sd == "no baseline" and str(m.reading).startswith("no "))
            diffs.append({"arm": vr.arm, "observed_ppm": vr.observed_ppm, "verifier_baseline": vr.baseline,
                          "mine": lab, "verifier_reading": vr.reading_2sd, "my_reading": m.reading, "agree": bool(same),
                          "centre_diff": (float(vr.centre) - float(m.centre_ppm)) if pd.notna(vr.centre) and pd.notna(m.centre_ppm) else None})
        xcheck["rows_compared"] = len(diffs)
        xcheck["all_readings_agree"] = all(d["agree"] for d in diffs)
        xcheck["max_abs_centre_diff"] = max((abs(d["centre_diff"]) for d in diffs if d["centre_diff"] is not None), default=None)
        xcheck["rows"] = diffs

    summary = {
        "label": L.POSTHOC_LABEL,
        "proteome_scale_file_Fig1D": {
            "sulfide_ever_placed": audit["proteome_scale_search"]["sulfide_ever_placed"],
            "modifications_placed": [m["modification"] for m in audit["proteome_scale_search"]["modifications_actually_placed"]],
            "oxidation_residues": audit["proteome_scale_search"]["modifications_actually_placed"][0]["residues"],
            "table_counts": audit["proteome_scale_search"]["table_counts"],
            "high_confidence_peptides": "317 (report prose only; not in a stored machine-readable output)",
        },
        "targeted_arms": audit["targeted_arms"],
        "n_psms_by_arm": df.arm.value_counts().to_dict(),
        "control_file_psms": "0 (the Fig4-L_control consensus file carries no PSM)" if "control" not in set(df.arm) else int((df.arm == "control").sum()),
        "n_sulfide_psms": int(len(ps)),
        "sulfide_psms_by_arm": ps.arm.value_counts().to_dict(),
        "dioxidation_or_trioxidation_placed_anywhere": False,
        "n_targeted_psms_with_nem": int(df.placed.str.contains("Nethylmaleimide").sum()),
        "n_targeted_psms": int(len(df)),
        "insulin_sequences_assigned_by_deposit": sorted(df.sequence.unique().tolist()),
        "stored_rule_resolvable_count": int(ps.rule_stored_audit_abs_ppm_lt_half_sep.sum()),
        "baseline_corrected_arm_rule_resolvable_count": int(ps.rule_baseline_corrected_arm.sum()),
        "baseline_corrected_pooled_rule_resolvable_count": int(ps.rule_baseline_corrected_pooled.sum()),
        "readings_file_baseline": ps.reading_file_baseline.value_counts().to_dict(),
        "readings_by_baseline_definition": sens.groupby("baseline_definition").reading.value_counts().unstack(fill_value=0).to_dict(orient="index"),
        "robust_summary": robust[["arm", "chain", "observed_ppm", "robust_summary"]].to_dict(orient="records"),
        "chain_offset": chain_offset,
        "rule_in_consensus_audit_script_is_abs_ppm": rule_in_script,
        "rule_claimed_in_registration": rule_claim,
        "stored_verdict": audit["verdict"],
        "verifier_cross_check": {k: v for k, v in xcheck.items() if k != "rows"},
        "caveats": [
            "the deposit's precursor tolerance and searched modification list are not recorded in the consensus files",
            "if the window edge was +10 ppm, a dioxidised peptide reported as sulfide would sit at baseline + separation; "
            "depending on the baseline this is outside (+12.8 ppm, A chain against a B-dominated baseline) or inside "
            "(+6.1 ppm against the A-chain baseline; +9.6 to +10.5 ppm for the B chain) the window, so the window may "
            "truncate dioxidation-consistent errors and the errors of reported sulfide PSMs are not independent evidence for sulfide",
            "neutral masses in the stored extraction are recomputed approximations (residues + water + NEM/sulfide count), not the deposit's calculated masses",
            "the A-chain baseline rests on two PSMs from another file",
            "four PSMs; illustration, not a rate",
        ],
    }
    L.write_csv(ps, "deposit_sulfide_psms_reassessed.csv")
    L.write_csv(base, "deposit_targeted_arm_baselines.csv")
    L.write_json(summary, "deposit_plus32_summary.json")
    if xcheck["available"]:
        L.write_json(xcheck, "deposit_verifier_crosscheck.json")
    return summary, ps, base, sens, robust


if __name__ == "__main__":
    import json
    s, ps, base, sens, robust = main()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.max_colwidth", 80)
    print(base.to_string())
    print(sens[["arm", "chain", "observed_ppm", "baseline_definition", "n_baseline_psms", "centre_ppm", "spread_ppm",
                "z_if_sulfide", "z_if_dioxidation", "reading", "half_separation_rule_resolvable_as_sulfide",
                "half_separation_rule_resolvable_as_dioxidation", "expected_ppm_if_dioxidised_but_reported_as_sulfide",
                "p_error_inside_plus10ppm_edge_if_dioxidised"]].round(3).to_string())
    print(robust.to_string())
    print(json.dumps({k: s[k] for k in ("chain_offset", "verifier_cross_check", "readings_file_baseline")}, indent=1, default=str))
