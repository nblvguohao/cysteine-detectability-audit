"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific: revision after verification (round 1).

The CAM share of identified cysteines WITHOUT hydroxylamine (HydN arm) against WITH it (HydP arm of the same protease),
under the same depth controls that c05 applies to the between-digest contrast.

Why: the first version of this item (and the current manuscript, line 153) set the raw shares side by side (trypsin:
660/2,024 = 0.326 without against 1,359/5,328 = 0.255 with hydroxylamine) and read the higher share without
hydroxylamine as evidence that the readout is not specific. The verifier showed that this comparison carries the same
depth confound that c05 removes from the between-digest contrast: the HydN arm is 2.6-fold shallower and concentrated
on abundant proteins. This script applies c05's controls to it, with HydP (the deeper arm) in the reference role.

Quantity: share(HydN) - share(HydP), share = identified cysteines carrying CAM (localisation >= 0.75).
Controls (names as in c05):
  raw                                   all identified cysteines of each arm
  abundance_deciles_plus_missing        HydP directly standardised to HydN's distribution over the c05 decile bins of
                                        log10 Global protein intensity + a missing class (primary depth control)
  abundance_20bins_plus_missing         the same with 20 bins
  abundance_deciles_complete_case       deciles, cysteines of proteins with an abundance only
  reverse_hydn_std_to_hydp_deciles      HydN standardised to HydP's decile distribution (sensitivity; undefined when a
                                        HydP bin has no HydN rows, reported as NaN)
  protein_matched                       both arms restricted to proteins identified in both
  protein_matched_protein_equal         the same, every protein weighted equally
  cysteine_matched                      the same cysteines (identified in both arms), CAM compared within cysteine;
                                        discordant counts and an exact two-sided binomial (McNemar) test
  cysteine_matched_single_cys_peptides  the same, cysteines seen only on single-cysteine peptides in both arms
Decile edges: identical to c05 (deciles of log10 abundance over the union of proteins identified in the four HydP
arms), so the bins are the ones used for the between-digest contrast.
Intervals: paired protein bootstrap over the union of proteins (one draw applied to both arms; c_boot), 5000
replicates, seed 20260930 + 3000 + running offset (listed per row).

Context (not a control): evidence rows per arm from the stored Phase 4 conversion statistics
(repo/results/phase4_validation_2026-09-22_audit.json, V5.conversion). A cysteine counts as a site when >= 1 evidence
row carries CAM at localisation >= 0.75, so a cysteine covered by more evidence rows is more likely to be called.
Per-cysteine evidence counts are NOT in the converted tables, so within-cysteine depth cannot be controlled here.

Outputs (C.OUT): t8_hydn_vs_hydp_controls.csv, t8_hydn_vs_hydp_deciles.csv, t8_bin_distribution.csv,
t8_evidence_depth.csv
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402
import c_boot as B  # noqa: E402
import c05_depth_confound as D  # noqa: E402  (edges_from, assign_bins; main() is guarded)

PHASE4_AUDIT = f"{C.REPO}/results/phase4_validation_2026-09-22_audit.json"
SEED0 = C.NEW_SEED + 3000


def load_pair(a, e10, e20):
    """HydP (from the task-2 site table) and HydN (arm table) for one protease, with abundance bins."""
    P = pd.read_csv(f"{C.OUT}/t2_site_table_{a}.tsv", sep="\t", dtype={"protein": str})
    P = P.assign(cam=P.cam_hydp)
    N = C.read_arm(a, "HydN")
    N = N.assign(cam=N.label)
    for d in (P, N):
        d["log10_abundance"] = np.log10(d.abundance)
        d["bin10"] = D.assign_bins(d, e10)
        d["bin20"] = D.assign_bins(d, e20)
    # abundance is a protein property: it must agree between the two arms for every shared protein
    ab = pd.concat([P[["protein", "abundance"]], N[["protein", "abundance"]]]).dropna().drop_duplicates()
    assert ab.protein.is_unique, f"abundance differs between HydP and HydN for some protein ({a})"
    return P, N


def main():
    C.record_inputs(f"c08_hydn_vs_hydp_depth[{C.BASE}]",
                    [f"{C.OUT}/t2_site_table_{a}.tsv" for a in C.ARMS] + [C.arm_path(a, "HydN") for a in C.ARMS]
                    + [PHASE4_AUDIT])
    stP = D.load()                       # the four HydP site tables, as c05 reads them
    e10, _ = D.edges_from(stP, 10)
    e20, _ = D.edges_from(stP, 20)
    rows, dec_rows, dist_rows = [], [], []
    off = 0
    for a in C.ARMS:
        P, N = load_pair(a, e10, e20)
        # --- bin distribution (share of identified cysteines per abundance bin) and per-decile shares
        for b in range(11):
            for arm_name, d in (("HydP", P), ("HydN", N)):
                s = d[d.bin10 == b]
                dist_rows.append({"arm": a, "hydroxylamine": arm_name, "bin": b + 1 if b < 10 else "missing",
                                  "n_identified": len(s), "share_of_arm": len(s) / len(d),
                                  "n_proteins": s.protein.nunique()})
                if len(s):
                    pt, lo, hi = B.share_boot(s.protein, s.cam.values == 1, C.REPS, SEED0 + 500 + 20 * C.ARMS.index(a)
                                              + b + (0 if arm_name == "HydP" else 11))
                else:
                    pt = lo = hi = np.nan
                dec_rows.append({"arm": a, "hydroxylamine": arm_name, "decile": b + 1 if b < 10 else "missing",
                                 "n_identified": len(s), "n_cam": int(s.cam.sum()), "n_proteins": s.protein.nunique(),
                                 "cam_share": pt, "ci_low": lo, "ci_high": hi})
        # --- controls
        shared_p = set(P.protein) & set(N.protein)
        specs = [
            ("raw", P, N, None, 1, False, False),
            ("abundance_deciles_plus_missing", P, N, "bin10", 11, False, False),
            ("abundance_20bins_plus_missing", P, N, "bin20", 21, False, False),
            ("abundance_deciles_complete_case", P[P.bin10 < 10], N[N.bin10 < 10], "bin10", 11, False, False),
            ("reverse_hydn_std_to_hydp_deciles", P, N, "bin10", 11, False, True),
            ("protein_matched", P[P.protein.isin(shared_p)], N[N.protein.isin(shared_p)], None, 1, False, False),
            ("protein_matched_protein_equal", P[P.protein.isin(shared_p)], N[N.protein.isin(shared_p)], None, 1, True,
             False),
        ]
        for name, dP, dN, bc, nb, peq, reverse in specs:
            off += 1
            seed = SEED0 + off
            if not reverse:
                # reference T = HydP (standardised to HydN's bins when bc is given), X = HydN
                r = B.paired_contrast(dP, dN, "cam", C.REPS, seed, bin_col=bc, n_bins=nb, protein_equal=peq)
                share_n, share_p_raw, share_p_adj = r["cX"], r["cT"], r["cT_std"]
                adj_lo, adj_hi = r["cT_std_ci_low"], r["cT_std_ci_high"]
                diff, lo, hi, und = r["diff_std"], r["diff_std_ci_low"], r["diff_std_ci_high"], r["diff_std_n_undefined_reps"]
                n_p, n_n, np_p, np_n = r["n_T"], r["n_X"], r["n_proteins_T"], r["n_proteins_X"]
                adjusted_arm = "HydP"
            else:
                # reference T = HydN standardised to HydP's bins; X = HydP; report HydN_std - HydP (sign flipped)
                r = B.paired_contrast(dN, dP, "cam", C.REPS, seed, bin_col=bc, n_bins=nb, protein_equal=peq)
                share_n, share_p_raw, share_p_adj = r["cT_std"], r["cX"], r["cX"]
                adj_lo, adj_hi = r["cT_std_ci_low"], r["cT_std_ci_high"]
                diff, lo, hi, und = -r["diff_std"], -r["diff_std_ci_high"], -r["diff_std_ci_low"], r["diff_std_n_undefined_reps"]
                n_p, n_n, np_p, np_n = r["n_X"], r["n_T"], r["n_proteins_X"], r["n_proteins_T"]
                adjusted_arm = "HydN"
            rows.append({"arm": a, "control": name, "adjusted_arm": adjusted_arm if bc is not None else "",
                         "cam_share_hydn": share_n, "cam_share_hydp_raw": share_p_raw,
                         "cam_share_hydp_compared": share_p_adj, "adjusted_share_ci_low": adj_lo,
                         "adjusted_share_ci_high": adj_hi,
                         "difference_hydn_minus_hydp": diff, "diff_ci_low": lo, "diff_ci_high": hi,
                         "undefined_reps": und, "n_hydp": n_p, "n_hydn": n_n, "n_proteins_hydp": np_p,
                         "n_proteins_hydn": np_n, "seed": seed})
        j = P.merge(N, on=["protein", "position"], suffixes=("_P", "_N"))
        for name, jj in (("cysteine_matched", j),
                         ("cysteine_matched_single_cys_peptides",
                          j[(j.n_cys_peptide_P == 1) & (j.n_cys_peptide_N == 1)])):
            off += 1
            seed = SEED0 + off
            r = B.paired_matched(jj, "cam_P", "cam_N", C.REPS, seed)     # X = HydN, T = HydP: diff = HydN - HydP
            p_only, n_only = r["discordant_T1_X0"], r["discordant_T0_X1"]
            p_val = binomtest(p_only, p_only + n_only, 0.5).pvalue if (p_only + n_only) else np.nan
            # concordance of CAM status between the arms for the same cysteines (descriptive)
            n11, n00, nc = r["concordant_11"], r["concordant_00"], r["n_cys"]
            orr = (n11 + 0.5) * (n00 + 0.5) / ((p_only + 0.5) * (n_only + 0.5))       # Haldane-corrected odds ratio
            po = (n11 + n00) / nc
            pe = ((n11 + p_only) / nc) * ((n11 + n_only) / nc) + ((n00 + n_only) / nc) * ((n00 + p_only) / nc)
            kappa = (po - pe) / (1 - pe) if pe < 1 else np.nan
            rows.append({"arm": a, "control": name, "adjusted_arm": "",
                         "cam_share_hydn": r["cX"], "cam_share_hydp_raw": r["cT"], "cam_share_hydp_compared": r["cT"],
                         "adjusted_share_ci_low": r["cT_ci_low"], "adjusted_share_ci_high": r["cT_ci_high"],
                         "difference_hydn_minus_hydp": r["diff"], "diff_ci_low": r["diff_ci_low"],
                         "diff_ci_high": r["diff_ci_high"], "undefined_reps": 0, "n_hydp": r["n_cys"],
                         "n_hydn": r["n_cys"], "n_proteins_hydp": r["n_proteins"], "n_proteins_hydn": r["n_proteins"],
                         "seed": seed, "cam_both": r["concordant_11"], "cam_neither": r["concordant_00"],
                         "cam_only_with_hydroxylamine": p_only, "cam_only_without_hydroxylamine": n_only,
                         "exact_binomial_p_two_sided": p_val, "concordance_odds_ratio_haldane": orr,
                         "cohen_kappa": kappa,
                         "share_of_hydp_cam_also_cam_in_hydn": r["concordant_11"] / (r["concordant_11"] + p_only)
                         if (r["concordant_11"] + p_only) else np.nan})
    ctl = pd.DataFrame(rows)
    ctl["interval_excludes_zero"] = (ctl.diff_ci_low > 0) | (ctl.diff_ci_high < 0)
    raw = ctl[ctl.control == "raw"].set_index("arm").difference_hydn_minus_hydp
    ctl["share_of_raw_difference_remaining"] = [d / raw[a] for d, a in zip(ctl.difference_hydn_minus_hydp, ctl.arm)]
    ctl.to_csv(f"{C.OUT}/t8_hydn_vs_hydp_controls.csv", index=False)
    pd.DataFrame(dec_rows).to_csv(f"{C.OUT}/t8_hydn_vs_hydp_deciles.csv", index=False)
    dist = pd.DataFrame(dist_rows)
    dist.to_csv(f"{C.OUT}/t8_bin_distribution.csv", index=False)

    # evidence rows per arm (stored conversion statistics; context only)
    conv = json.load(open(PHASE4_AUDIT, encoding="utf-8"))["results"]["V5"]["conversion"]
    ev = []
    for a in C.ARMS:
        for h in ("HydP", "HydN"):
            c = conv[f"{a}_{h}"]
            ev.append({"arm": a, "hydroxylamine": h, "evidence_rows_in_arm_all_peptides": c["evidence_rows_in_arm"],
                       "identified_cysteines_stored": c["observed"], "cam_sites_stored": c["positive"],
                       "evidence_rows_per_identified_cysteine": c["evidence_rows_in_arm"] / c["observed"],
                       "note": "evidence rows of all peptides of the arm (not only cysteine peptides); per-cysteine "
                               "evidence counts are not in the converted tables"})
    evd = pd.DataFrame(ev)
    tr = evd[(evd.arm == "Trypsin") & (evd.hydroxylamine == "HydP")].evidence_rows_in_arm_all_peptides.iloc[0]
    evd["trypsin_hydp_rows_over_this_arm"] = tr / evd.evidence_rows_in_arm_all_peptides
    evd.to_csv(f"{C.OUT}/t8_evidence_depth.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.max_rows", 500)
    print(ctl[["arm", "control", "cam_share_hydn", "cam_share_hydp_raw", "cam_share_hydp_compared",
               "difference_hydn_minus_hydp", "diff_ci_low", "diff_ci_high", "share_of_raw_difference_remaining",
               "n_hydp", "n_hydn", "undefined_reps"]].round(4).to_string(index=False))
    print(ctl[ctl.control.str.startswith("cysteine")][["arm", "control", "cam_both", "cam_neither",
                                                       "cam_only_with_hydroxylamine", "cam_only_without_hydroxylamine",
                                                       "exact_binomial_p_two_sided", "concordance_odds_ratio_haldane",
                                                       "cohen_kappa",
                                                       "share_of_hydp_cam_also_cam_in_hydn"]].to_string(index=False))
    t = dist[dist.arm == "Trypsin"].pivot(index="bin", columns="hydroxylamine", values="share_of_arm")
    print(t.round(3))
    dd = pd.DataFrame(dec_rows)
    print(dd[dd.arm == "Trypsin"].pivot(index="decile", columns="hydroxylamine", values="cam_share").round(3))
    print(evd.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
