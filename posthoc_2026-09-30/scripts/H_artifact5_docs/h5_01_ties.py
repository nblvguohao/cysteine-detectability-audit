"""H_artifact5_docs step 1 - POST HOC verification of the PXD015307 score-tie numbers (2026-09-30).

POST HOC revision analysis in response to a pre-submission review; not registered, not pre-specified.
It re-reads the stored output results/pxd015307_posthoc_score_ties.csv (7,968 rows, one per spectrum
whose Comet top-5 output lines carry at least one +32-class cysteine modification) and checks every
number the manuscript quotes from it. It also asks what the tied spectra are:

  * decomposition of the 2,079 'dioxidation+sulfide' rows into (i) a single best interpretation that
    carries both chemistries on different cysteines (not a tie), (ii) ties between interpretations of
    identical composition (zero precursor-error span), (iii) ties between interpretations of different
    composition on one peptide backbone, (iv) ties between different peptide backbones;
  * target or decoy status of the best-scoring peptide(s), by matching the sequence to the mouse
    reference proteome (UniProt 2026_03, canonical; the search itself used an earlier Swiss-Prot
    canonical+isoform release, so isoform-only peptides stay 'unresolved') and, for decoys, to the
    reverse-all-but-C-terminal-residue construction used by Comet's internal decoys;
  * whether a non-zero span is simply k x 0.017758 Da expressed in ppm of a plausible peptide mass
    (arithmetic), with a shifted-mass control;
  * enzyme specificity and missed cleavages inferred from the target peptides in the output.
No random numbers are drawn.

Revision after verification (round 2): adds the count and median best E-value of the ties on one peptide
(categories B + C, 1,810 spectra), which is the set the proposed Results sentence describes.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5_lib as L  # noqa: E402

TIES = os.path.join(L.W, "inputs", "repo_results", "pxd015307_posthoc_score_ties.csv")
AUDIT = os.path.join(L.W, "inputs", "repo_results", "pxd015307_research_audit.json")
TIES_AUDIT = os.path.join(L.REPO, "results", "pxd015307_posthoc_score_ties_audit.json")
FASTA = os.path.join(L.W, "external", "UP000000589_10090.fasta.gz")
SRC_FIG3B = os.path.join(L.W, "inputs", "mcp_package", "Source_Data_Fig3b_tie_precursor_error_spans.csv")
SRC_FIG3 = os.path.join(L.W, "inputs", "mcp_package", "Source_Data_Fig3_search_space.csv")

# Bovine insulin chains as recorded internally (A chain: posthoc_pxd015307_score_ties.py docstring;
# B chain: bovine B30 Ala). Context residues are only those certain for proinsulins (KR before the
# A chain; RR after the B chain). P01317 itself is not available offline.
INSULIN_BOVINE = {"INS_BOVIN_B_context": "FVNQHLCGSHLVEALYLVCGERGFFYTPKARR",
                  "INS_BOVIN_A_context": "KRGIVEQCCASVCSLYQLENYCN"}
HUMAN_INSULIN_A = "GIVEQCCTSICSLYQLENYCN"
HUMAN_INSULIN_B = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"


def classify_sequences(peps, db):
    """Return dict peptide -> (is_target, is_decoy)."""
    out = {}
    for p in peps:
        t = db.find(p) >= 0
        d = db.find(L.comet_decoy(p)) >= 0
        out[p] = (t, d)
    return out


def td_label(t, d):
    if t and not d:
        return "target"
    if d and not t:
        return "decoy"
    if t and d:
        return "ambiguous"
    return "unresolved"


def row_td(peptide, cls):
    labs = {td_label(*cls[p]) for p in peptide.split("|")}
    if len(labs) == 1:
        return labs.pop()
    if "target" in labs and "decoy" in labs:
        return "mixed_target_decoy"
    return "mixed_other"


def category(r):
    if r.n_distinct_plus32_patterns_at_best_xcorr == 1:
        return "A_single_interpretation_not_a_tie"
    if r.n_backbones_at_best_xcorr > 1:
        return "D_tie_between_different_peptides"
    if r.ppm_span == 0:
        return "B_tie_same_composition"
    return "C_tie_different_composition_same_peptide"


def span_arithmetic(r, shift=0.0):
    """Does span == k * 0.017758 Da in ppm of M, for M = backbone + a plausible modification set?"""
    seq = r.peptide
    n_c, n_m = seq.count("C"), seq.count("M")
    m_bb = L.peptide_mass(seq) + shift
    span = r.ppm_span
    best = None
    for n32 in range(1, 4):
        for a in range(0, 4):          # iodoTMT6plex on C
            for b in range(0, 4):      # carbamidomethyl on C
                for c in range(0, 4):  # oxidation on M
                    if n32 + a + b + c > 3 or n32 + a + b > n_c or c > n_m:
                        continue
                    m = m_bb + n32 * 31.98 + a * L.IODOTMT + b * L.CAM + c * L.OX_M
                    for k in range(1, n32 + 1):
                        pred = k * L.SEPARATION_DA / m * 1e6
                        err = abs(pred - span)
                        if best is None or err < best[0]:
                            best = (err, k, n32, a, b, c, m, pred)
    err, k, n32, a, b, c, m, pred = best
    # stored ppm values are rounded to 3 decimals; allow 0.0015 ppm plus 0.05% of the span
    ok = err <= 0.0015 + 5e-4 * span
    return ok, k, n32, a, b, c, m, pred, err


def main():
    ties = pd.read_csv(L.register(TIES, "stored post hoc tie table (repo results)"),
                       dtype={"e_value_best": str})
    audit = L.read_json(AUDIT, "stored re-search audit (per-file thresholds, counts)")
    ties_audit = L.read_json(TIES_AUDIT, "stored tie-analysis audit (2,079; 18.829)")
    fig3b = pd.read_csv(L.register(SRC_FIG3B, "manuscript Source Data Fig 3b (MCP package copy)"))
    fig3 = pd.read_csv(L.register(SRC_FIG3, "manuscript Source Data Fig 3 (MCP package copy)"))
    ties["ev"] = ties.e_value_best.astype(float)
    thr = {arm: a["e_value_threshold_at_1pct_fdr"] for arm, a in audit["per_arm"].items()}
    ties["thr_1pct"] = ties.arm.map(thr)
    ties["targeted_arm"] = ties.arm != "Fig-1D"

    # ---------- target/decoy status of every peptide in the stored output ----------
    fasta = L.read_fasta_gz(FASTA, "UniProt mouse reference proteome 2026_03 (canonical)")
    db = "#".join(s for _, s in fasta) + "#" + "#".join(INSULIN_BOVINE.values()) + "#"
    peps = sorted({p for s in ties.peptide for p in s.split("|")})
    cls = classify_sequences(peps, db)
    ties["td_class"] = ties.peptide.map(lambda s: row_td(s, cls))
    uniq = pd.Series({p: td_label(*cls[p]) for p in peps}).value_counts()

    both = ties[ties.chemistries_at_best_xcorr == "dioxidation+sulfide"].copy()
    both["category"] = both.apply(category, axis=1)
    both["passes_1pct_fdr_evalue"] = both.ev <= both.thr_1pct
    both["is_bovine_insulin_A"] = both.peptide.str.contains("GIVEQCCASVCSLYQLENYCN")
    genuine = both[both.n_distinct_plus32_patterns_at_best_xcorr >= 2]
    # revision after verification: ties whose interpretations differ only in their sulfide/dioxidation
    # assignments on ONE peptide (categories B and C; excludes the 74 ties between different peptides)
    one_pep = both[both.category.isin(["B_tie_same_composition", "C_tie_different_composition_same_peptide"])]

    # ---------- claims ----------
    claims = []

    def claim(cid, text, claimed, recomputed, status, note=""):
        claims.append({"claim_id": cid, "manuscript_or_record_text": text, "claimed": claimed,
                       "recomputed": recomputed, "status": status, "note": note})

    claim("T1", "In 2,079 spectra the two chemistries tie at the best score", 2079, len(both),
          "count reproduced; wording not supported",
          f"{(both.n_distinct_plus32_patterns_at_best_xcorr == 1).sum()} of the 2,079 have a single best "
          f"interpretation carrying sulfide on one cysteine and dioxidation on another (no tie); spectra with "
          f">=2 distinct +32 patterns at the best xcorr and both chemistries present: {len(genuine)}")
    claim("T2", "none of these passes the 1% FDR threshold", 0, int(both.passes_1pct_fdr_evalue.sum()),
          "reproduced", f"lowest best E-value {both.ev.min():g} ({both.loc[both.ev.idxmin(), 'arm']}; "
          f"threshold {thr[both.loc[both.ev.idxmin(), 'arm']]:g})")
    claim("T3", "their median best expectation value is 17", 17, float(both.ev.median()), "reproduced",
          f"median over the {len(genuine)} spectra with >=2 patterns: {genuine.ev.median():g}; "
          f"Fig-1D {both[~both.targeted_arm].ev.median():g}, targeted arms {both[both.targeted_arm].ev.median():g}")
    hi = both[both.best_xcorr >= 1.5]
    claim("T4", "The 237 ties that reach xcorr >= 1.5 all lie in the four targeted arms", 237,
          f"{len(hi)} (Fig-1D: {(~hi.targeted_arm).sum()})", "count reproduced; interpretation qualified",
          f"{(hi.n_distinct_plus32_patterns_at_best_xcorr == 1).sum()} of the {len(hi)} are single "
          f"interpretations; {(hi.td_class == 'decoy').sum()} of the {len(hi)} best matches are decoy "
          f"peptides; xcorr from 1.0005-Da and 0.02-Da binning are not on one scale (max Fig-1D xcorr "
          f"{both[~both.targeted_arm].best_xcorr.max():g})")
    claim("T5", "which were acquired with ion-trap MS2 binned at 1.0005 Da", "prose", "not in a stored output",
          "consistent in prose records only",
          "scan counts stored only for Fig-1D (31,450 FTMS/3 ITMS) and Fig4-L_NaHS (1,080 FTMS/5,004 ITMS) "
          "in the registration amendment; fragment settings are in comet.params files not in the local copy")
    imax = both.ppm_span.idxmax()
    sb = both[both.n_backbones_at_best_xcorr == 1]
    claim("T6", "the precursor mass error spans up to 18.829 ppm", 18.829, float(both.ppm_span.max()),
          "number reproduced; interpretation not supported",
          f"the 18.829-ppm row is {both.loc[imax, 'arm']} scan {both.loc[imax, 'scan']}, a tie between two "
          f"different peptides ({both.loc[imax, 'n_backbones_at_best_xcorr']} backbones), best E-value "
          f"{both.loc[imax, 'ev']:g}; largest span on one backbone {sb.ppm_span.max():g} ppm; a span on one "
          f"backbone is k x 0.017758 Da in ppm of the peptide mass (see span check), i.e. arithmetic")
    z = both[both.ppm_span == 0]
    zc = z.apply(category, axis=1).value_counts()
    claim("T7", "in 398 of them it is exactly zero, meaning that the two interpretations share a composition "
                "and differ in position rather than in chemistry", 398, len(z), "count reproduced; wording wrong",
          f"{int(zc.get('A_single_interpretation_not_a_tie', 0))} have a single +32 pattern at the best xcorr "
          f"and are not ties ({(z.n_lines_at_best_xcorr == 1).sum()} of them are one output line); "
          f"{int(zc.get('B_tie_same_composition', 0))} are ties between interpretations of identical composition "
          f"that swap sulfide and dioxidation between cysteines; "
          f"{int(zc.get('D_tie_between_different_peptides', 0))} are ties between different peptides")
    claim("T8", "report: 877 tied spectra in the four targeted arms", 877, int(both.targeted_arm.sum()), "reproduced")
    bins = pd.cut(both.best_xcorr, [-np.inf, 0.5, 1.0, 1.5, 2.0, np.inf], right=False).value_counts().sort_index()
    claim("T9", "Fig. 3a bins <0.5, 0.5-1, 1-1.5, 1.5-2, >=2", "957/546/339/170/67",
          "/".join(str(int(x)) for x in bins.values), "reproduced", "bins pool files whose xcorr scales differ")
    ex = both[(both.arm == "Fig4-L_CTH_Cys") & (both.scan == 2124)].iloc[0]
    claim("T10", "report example: Fig4-L_CTH_Cys scan 2124, three chemistries at one score", "xcorr 1.9640",
          f"xcorr {ex.best_xcorr}, E {ex.ev:g}, peptide {ex.peptide}", "reproduced; example is not an identification",
          f"best E-value {ex.ev:g} vs 1%-FDR threshold {thr['Fig4-L_CTH_Cys']:g}; the peptide is the bovine "
          f"insulin A chain, whereas the deposit's own search of these files assigned the human A chain "
          f"({HUMAN_INSULIN_A})")
    claim("T11", "tie audit per-arm 'tied' counts (any two +32 patterns share the best xcorr)",
          "stored", int(sum(a['n_spectra_where_several_plus32_patterns_share_the_best_xcorr']
                            for a in ties_audit['per_arm'].values())),
          "different quantity", "the stored audit's per-arm 'tied' counts (1,287+273+254+16+337) include ties "
                                "between two positions of the same chemistry; they are not the 2,079")
    claim("T12", "Source Data Fig 3b rows", 2079, len(fig3b), "consistent with the 2,079 rows",
          "Fig 3b therefore includes the single-interpretation rows")
    claim("T13", "revision proposal (round 1): '1,884 spectra; median best E-value 16' for interpretations "
                 "differing only in their sulfide and dioxidation assignments", "1,884; 16",
          f"{len(one_pep)}; {one_pep.ev.median():g}", "corrected after verification",
          "74 of the 1,884 are ties between different peptides; ties on one peptide (categories B + C) number "
          f"{len(one_pep)} with median best E-value {one_pep.ev.median():g}")
    claims_df = pd.DataFrame(claims)

    # ---------- decomposition ----------
    dec = both.groupby(["arm", "category"]).size().unstack(fill_value=0)
    dec.loc["ALL"] = dec.sum()
    dec = dec.reset_index()
    zdec = pd.crosstab([z.n_lines_at_best_xcorr], [z.n_distinct_plus32_patterns_at_best_xcorr])
    zdec.columns = [f"n_patterns_{c}" for c in zdec.columns]
    zdec = zdec.reset_index()
    zcat = z.apply(category, axis=1).value_counts().rename_axis("category").reset_index(name="n")

    # ---------- target / decoy ----------
    td_arm = pd.crosstab(both.arm, both.td_class, margins=True, margins_name="ALL").reset_index()
    both["xcorr_bin"] = pd.cut(both.best_xcorr, [-np.inf, 0.5, 1.0, 1.5, 2.0, np.inf], right=False).astype(str)
    td_x = pd.crosstab([both.targeted_arm, both.xcorr_bin], both.td_class).reset_index()
    both["evalue_bin"] = pd.cut(both.ev, [0, 0.1, 1, 10, 100, 1000.1]).astype(str)
    td_e = pd.crosstab(both.evalue_bin, both.td_class).reset_index()
    td_hi = pd.crosstab(hi.arm, hi.td_class, margins=True, margins_name="ALL").reset_index()
    td_all_rows = pd.crosstab(ties.arm, ties.td_class, margins=True, margins_name="ALL").reset_index()
    td_cat = pd.crosstab(both.category, both.td_class, margins=True, margins_name="ALL").reset_index()

    def dfrac(df):
        t = (df.td_class == "target").sum()
        d = (df.td_class == "decoy").sum()
        return t, d, d / (t + d) if (t + d) else np.nan

    tdsum = []
    for name, sub in [("all 2,079", both), ("genuine ties (>=2 patterns)", genuine),
                      ("xcorr >= 1.5 (237)", hi), ("Fig-1D", both[~both.targeted_arm]),
                      ("targeted arms", both[both.targeted_arm]),
                      ("best E-value <= 1", both[both.ev <= 1]),
                      ("all 7,968 stored rows", ties)]:
        t, d, f = dfrac(sub)
        other = int(len(sub) - t - d)
        tdsum.append({"subset": name, "n_rows": len(sub), "n_target": int(t), "n_decoy": int(d),
                      "decoy_fraction_of_target_plus_decoy": f,
                      "decoy_fraction_lower_bound_other_counted_as_target": d / len(sub) if len(sub) else np.nan,
                      "decoy_fraction_upper_bound_other_counted_as_decoy": (d + other) / len(sub) if len(sub) else np.nan,
                      "n_other_classes": other})
    tdsum = pd.DataFrame(tdsum)

    # ---------- span arithmetic on single-backbone, different-composition ties ----------
    cset = both[both.category == "C_tie_different_composition_same_peptide"].copy()
    rows = []
    for _, r in cset.iterrows():
        ok, k, n32, a, b, c, m, pred, err = span_arithmetic(r)
        oks = [span_arithmetic(r, shift=s)[0] for s in (-7.7, -3.3, 3.3, 7.7)]
        rows.append({"arm": r.arm, "scan": r.scan, "peptide": r.peptide, "ppm_span": r.ppm_span,
                     "explained": ok, "k_substitutions": k, "n_plus32": n32, "n_iodoTMT": a, "n_CAM": b,
                     "n_oxM": c, "peptide_mass_used": m, "predicted_span_ppm": pred, "abs_error_ppm": err,
                     "explained_under_shifted_masses_fraction": float(np.mean(oks))})
    span_df = pd.DataFrame(rows)
    span_summary = {
        "n_single_backbone_different_composition_ties": int(len(span_df)),
        "n_explained_exactly_by_k_times_separation": int(span_df.explained.sum()),
        "fraction_explained": float(span_df.explained.mean()),
        "k_distribution_among_explained": span_df[span_df.explained].k_substitutions.value_counts().sort_index().to_dict(),
        "fraction_explained_under_shifted_masses_control_mean": float(span_df.explained_under_shifted_masses_fraction.mean()),
        "note": "a span on one peptide backbone equals k x 0.017758 Da / M x 1e6; it is set by peptide mass, not by the spectrum",
    }

    # ---------- precursor errors of single-composition tie rows vs the file baseline ----------
    pe_rows = []
    for arm, a in audit["per_arm"].items():
        sub = both[(both.arm == arm) & (both.ppm_span == 0)]
        vals = sub.ppm_min.values
        med, sd = a["baseline_median_ppm"], a["baseline_sd_ppm"]
        lo, hi2 = max(-10, med - 2 * sd), min(10, med + 2 * sd)
        within = float(np.mean((vals >= lo) & (vals <= hi2))) if len(vals) else np.nan
        ks = stats.kstest((vals + 10) / 20, "uniform") if len(vals) >= 5 else None
        pe_rows.append({"arm": arm, "n_zero_span_rows": int(len(vals)), "baseline_median_ppm": med,
                        "baseline_sd_ppm": sd, "window_lo": lo, "window_hi": hi2,
                        "fraction_within_baseline_2sd": within,
                        "expected_fraction_if_uniform_in_pm10ppm": (hi2 - lo) / 20,
                        "ks_uniform_statistic": None if ks is None else float(ks.statistic),
                        "ks_uniform_p": None if ks is None else float(ks.pvalue)})
    pe_df = pd.DataFrame(pe_rows)

    # ---------- enzyme / length inference from target peptides in the output ----------
    prot_seqs = [s for _, s in fasta]
    enz = []
    for p in peps:
        if td_label(*cls[p]) != "target" or "GIVEQCC" in p or "FVNQHLC" in p:
            continue
        full = False
        for s in prot_seqs:
            i = s.find(p)
            while i >= 0:
                n_ok = (i == 0) or (i == 1 and s[0] == "M") or (s[i - 1] in "KR" and p[0] != "P")
                j = i + len(p)
                c_ok = (j == len(s)) or (p[-1] in "KR" and s[j] != "P")
                if n_ok and c_ok:
                    full = True
                    break
                i = s.find(p, i + 1)
            if full:
                break
        mc = sum(1 for q in range(len(p) - 1) if p[q] in "KR" and p[q + 1] != "P")
        enz.append({"peptide": p, "length": len(p), "fully_tryptic_occurrence_found": full,
                    "missed_cleavages_P_rule": mc,
                    "internal_KP_or_RP": sum(1 for q in range(len(p) - 1) if p[q] in "KR" and p[q + 1] == "P")})
    enz = pd.DataFrame(enz)
    enz_summary = {
        "n_target_peptides_checked": int(len(enz)),
        "fraction_with_a_fully_tryptic_occurrence": float(enz.fully_tryptic_occurrence_found.mean()),
        "missed_cleavage_distribution": enz.missed_cleavages_P_rule.value_counts().sort_index().to_dict(),
        "max_missed_cleavages": int(enz.missed_cleavages_P_rule.max()),
        "n_with_internal_KP_or_RP": int((enz.internal_KP_or_RP > 0).sum()),
        "n_exceeding_2_missed_cleavages_if_trypsin_P": int(((enz.missed_cleavages_P_rule
                                                             + enz.internal_KP_or_RP) > 2).sum()),
        "length_min": int(enz.length.min()), "length_max": int(enz.length.max()),
        "length_min_all_output_peptides": int(min(len(p) for p in peps)),
        "length_max_all_output_peptides": int(max(len(p) for p in peps)),
        "note": "inferred from output sequences, POST HOC; the comet.params files are not in the local copy",
    }

    # ---------- insulin rows ----------
    ins = ties[ties.peptide.str.contains("GIVEQCC") | ties.peptide.str.contains("FVNQHLC")]
    ins_summary = {
        "n_rows_any_insulin_chain": int(len(ins)),
        "arms": ins.arm.value_counts().to_dict(),
        "best_evalue_min": float(ins.ev.min()) if len(ins) else None,
        "human_insulin_sequences_present_in_output": bool(ties.peptide.str.contains(HUMAN_INSULIN_A).any()
                                                          or ties.peptide.str.contains(HUMAN_INSULIN_B).any()),
        "bovine_A_chain_rows": int(ties.peptide.str.contains("GIVEQCCASVCSLYQLENYCN").sum()),
        "B_chain_rows": int(ties.peptide.str.contains("FVNQHLC").sum()),
    }

    # ---------- write ----------
    keep = ["arm", "scan", "peptide", "n_backbones_at_best_xcorr", "best_xcorr", "n_lines_at_best_xcorr",
            "n_distinct_plus32_patterns_at_best_xcorr", "chemistries_at_best_xcorr", "ppm_min", "ppm_max",
            "ppm_span", "e_value_best", "ev", "thr_1pct", "passes_1pct_fdr_evalue", "category", "td_class",
            "is_bovine_insulin_A"]
    L.write_csv(both[keep].rename(columns={"ev": "e_value_best_float", "thr_1pct": "evalue_threshold_1pct_fdr_file"}),
                "tie_rows_annotated.csv")
    L.write_csv(claims_df, "ties_claims_verification.csv")
    L.write_csv(dec, "ties_decomposition_by_arm.csv")
    L.write_csv(zdec, "ties_zero_span_lines_by_patterns.csv")
    L.write_csv(zcat, "ties_zero_span_categories.csv")
    L.write_csv(td_arm, "ties_target_decoy_by_arm.csv")
    L.write_csv(td_x, "ties_target_decoy_by_xcorr_bin.csv")
    L.write_csv(td_e, "ties_target_decoy_by_evalue_bin.csv")
    L.write_csv(td_hi, "ties_target_decoy_xcorr_ge_1p5.csv")
    L.write_csv(td_cat, "ties_target_decoy_by_category.csv")
    L.write_csv(td_all_rows, "all_stored_rows_target_decoy_by_arm.csv")
    L.write_csv(tdsum, "ties_target_decoy_summary.csv")
    L.write_csv(span_df, "ties_span_arithmetic_check.csv")
    L.write_csv(pe_df, "ties_precursor_error_vs_baseline.csv")
    L.write_csv(enz, "output_enzyme_inference_peptides.csv")
    summary = {
        "label": L.POSTHOC_LABEL,
        "n_rows_stored": int(len(ties)),
        "n_both_chemistries_rows": int(len(both)),
        "n_genuine_ties_ge2_patterns": int(len(genuine)),
        "categories": both.category.value_counts().to_dict(),
        "categories_targeted_arms": both[both.targeted_arm].category.value_counts().to_dict(),
        "categories_fig1d": both[~both.targeted_arm].category.value_counts().to_dict(),
        "zero_span_rows": int(len(z)), "zero_span_categories": zc.to_dict(),
        "zero_span_single_line": int((z.n_lines_at_best_xcorr == 1).sum()),
        "median_best_evalue_all": float(both.ev.median()),
        "median_best_evalue_genuine": float(genuine.ev.median()),
        "n_ties_one_peptide_B_plus_C": int(len(one_pep)),
        "median_best_evalue_ties_one_peptide_B_plus_C": float(one_pep.ev.median()),
        "median_best_evalue_ties_one_peptide_by_group": {
            "Fig-1D": float(one_pep[~one_pep.targeted_arm].ev.median()),
            "insulin files": float(one_pep[one_pep.targeted_arm].ev.median())},
        "decoy_share_ties_one_peptide": float((one_pep.td_class == "decoy").sum() /
                                              max(1, ((one_pep.td_class == "decoy") | (one_pep.td_class == "target")).sum())),
        "genuine_ties_by_group": {"Fig-1D": int((~genuine.targeted_arm).sum()),
                                  "insulin files": int(genuine.targeted_arm.sum())},
        "median_best_evalue_by_arm": both.groupby("arm").ev.median().to_dict(),
        "max_xcorr_by_arm": both.groupby("arm").best_xcorr.max().to_dict(),
        "n_xcorr_ge_1p5": int(len(hi)), "n_xcorr_ge_1p5_fig1d": int((~hi.targeted_arm).sum()),
        "n_xcorr_ge_1p5_genuine": int((hi.n_distinct_plus32_patterns_at_best_xcorr >= 2).sum()),
        "max_span_all": float(both.ppm_span.max()), "max_span_single_backbone": float(sb.ppm_span.max()),
        "max_span_row": both.loc[imax, ["arm", "scan", "peptide", "ev"]].to_dict(),
        "n_pass_1pct_fdr_evalue": int(both.passes_1pct_fdr_evalue.sum()),
        "unique_peptide_td_counts": uniq.to_dict(),
        "target_decoy_summary": tdsum.to_dict(orient="records"),
        "span_arithmetic": span_summary,
        "enzyme_inference": enz_summary,
        "insulin_rows": ins_summary,
        "ppm_range_all_rows": [float(ties.ppm_min.min()), float(ties.ppm_max.max())],
        "fig3_source_bins_match": bool(list(fig3[fig3.panel == "a"].value.astype(int)) == list(bins.values.astype(int))),
    }
    L.write_json(summary, "ties_summary.json")
    return summary


if __name__ == "__main__":
    import json
    s = main()
    print(json.dumps({k: v for k, v in s.items() if k not in ("target_decoy_summary",)}, indent=1, default=str)[:6000])
