"""H_artifact5_docs step 5 - POST HOC: can the fragment score separate sulfide from dioxidation? (2026-09-30).

POST HOC revision analysis in response to a pre-submission review; not registered, not pre-specified.
Added in the revision after verification: the verifier showed that the round-1 statement 'fragment bins
of 0.02 and 1.0005 Da cannot resolve 0.0178 Da' is false for 0.02-Da bins. Three checks:

(a) Arithmetic. A fragment ion that contains the modified cysteine moves by 0.017758/z in m/z between
    the two interpretations. With Comet's binning, BIN(m/z) = int(m/z / w + (1 - offset)), the shift
    changes the bin with probability (0.017758/z)/w for positions uniform within a bin.
(b) Theoretical fragments of the tied peptides. For every unique single-backbone peptide among the 2,079
    'dioxidation+sulfide' rows, b and y ions containing each cysteine (+31.972071 on that cysteine, no other
    modification, charges 1-3) are binned with w = 0.02 Da (offset 0.0 and 0.4; the Fig-1D offset was not
    recorded) and with w = 1.0005 Da, offset 0.4 (the insulin files); reported is the share whose bin
    changes when sulfide is replaced by dioxidation.
(c) What the search output shows. Among spectra whose best-scoring interpretation is a single-cysteine
    peptide with one +32 modification (one backbone at the best xcorr), the alternative chemistry was a
    scored candidate whenever its precursor error (best error -/+ separation) lies inside the +/-10-ppm
    window. If it was scored and does not share the best xcorr, the fragment score separated the two.
    Rows whose classification depends on unknown methionine oxidation (up to two) or that lie within
    0.3 ppm of the window edge are set aside as 'edge'.
No random numbers are drawn.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5_lib as L  # noqa: E402

TIES = os.path.join(L.W, "inputs", "repo_results", "pxd015307_posthoc_score_ties.csv")
PROTON = 1.007276
EDGE_MARGIN = 0.3
WINDOW = 10.0


def comet_bin(mz, w, offset):
    return np.floor(mz / w + (1.0 - offset)).astype(np.int64)


def fragments_with_cys(pep):
    """Neutral b (residue sum) and y (residue sum + H2O) fragment masses (Da) containing each cysteine,
    with +31.972071 on that cysteine and no other modification."""
    res = np.array([L.AA[c] for c in pep])
    n = len(pep)
    out = []
    for p in [i for i, c in enumerate(pep) if c == "C"]:
        mod = res.copy()
        mod[p] += L.SULFIDE
        csum = np.cumsum(mod)
        for k in range(1, n):
            if k - 1 >= p:            # b_k holds residues 0..k-1
                out.append(csum[k - 1])
            if n - k <= p:            # y_k holds residues n-k..n-1
                out.append(csum[-1] - csum[n - k - 1] + L.H2O)
    return np.array(out)


def main():
    ties = pd.read_csv(L.register(TIES, "stored post hoc tie table (repo results)"), dtype={"e_value_best": str})
    ties["ev"] = ties.e_value_best.astype(float)
    ties["group"] = np.where(ties.arm == "Fig-1D", "Fig-1D (FTMS MS2, 0.02-Da bins)",
                             "insulin files (ion-trap MS2, 1.0005-Da bins)")

    # ---------- (a) arithmetic ----------
    arith = []
    for w, lab in ((0.02, "0.02 Da (Fig-1D)"), (1.0005, "1.0005 Da (insulin files)")):
        for z in (1, 2, 3):
            shift = L.SEPARATION_DA / z
            arith.append({"bin_width": lab, "fragment_charge": z, "shift_mz": shift,
                          "p_bin_changes_uniform_positions": min(1.0, shift / w)})
    arith = pd.DataFrame(arith)

    # ---------- (b) theoretical fragments of the tied peptides ----------
    both = ties[(ties.chemistries_at_best_xcorr == "dioxidation+sulfide") & (ties.n_backbones_at_best_xcorr == 1)]
    theo = []
    for grp, sub in both.groupby("group"):
        peps = sorted({p for p in sub.peptide if all(c in L.AA for c in p)})
        frag = np.concatenate([fragments_with_cys(p) for p in peps])
        settings = ([(0.02, 0.0), (0.02, 0.4)] if grp.startswith("Fig-1D") else []) + [(1.0005, 0.4)]
        for w, off in settings:
            for z in (1, 2, 3):
                # b ions: residue sum; y ions: residue sum + H2O; add z protons and divide by z
                mz = (frag + z * PROTON) / z
                ch = comet_bin(mz + L.SEPARATION_DA / z, w, off) != comet_bin(mz, w, off)
                theo.append({"group": grp, "n_unique_peptides": len(peps), "n_cys_containing_fragments": int(frag.size),
                             "bin_width": w, "bin_offset": off, "fragment_charge": z,
                             "share_of_fragments_changing_bin": float(ch.mean())})
    theo = pd.DataFrame(theo)

    # ---------- (c) single-cysteine spectra: was the alternative scored, and did it tie? ----------
    t = ties[(ties.n_backbones_at_best_xcorr == 1)].copy()
    t = t[t.peptide.map(lambda p: p.count("C") == 1 and all(c in L.AA for c in p))].copy()
    t["n_met"] = t.peptide.str.count("M")
    m0 = t.peptide.map(L.peptide_mass) + 31.98
    s_hi = L.SEPARATION_DA / m0 * 1e6                                  # no Met oxidation
    s_lo = L.SEPARATION_DA / (m0 + L.OX_M * np.minimum(t.n_met, 2)) * 1e6   # up to two Met oxidations
    e = t.ppm_min.astype(float)
    sign = np.where(t.chemistries_at_best_xcorr == "sulfide", -1.0, 1.0)  # alt error = e - s (best S) or e + s (best 2O)
    alt_a, alt_b = e + sign * s_hi, e + sign * s_lo
    inside = (np.abs(alt_a) <= WINDOW - EDGE_MARGIN) & (np.abs(alt_b) <= WINDOW - EDGE_MARGIN)
    outside = (np.abs(alt_a) >= WINDOW + EDGE_MARGIN) & (np.abs(alt_b) >= WINDOW + EDGE_MARGIN)
    t["alt_error_ppm_no_metox"] = alt_a
    t["class"] = np.select(
        [t.chemistries_at_best_xcorr == "dioxidation+sulfide", inside, outside],
        ["both scored, tied at the best xcorr", "both scored, alternative scored lower",
         "alternative outside the precursor window"], default="edge (within 0.3 ppm of the window edge)")
    tab = pd.crosstab(t.group, t["class"], margins=True, margins_name="ALL").reset_index()
    scored = t[t["class"].str.startswith("both scored")]
    rates = []
    for (grp, eb), sub in scored.assign(evalue_bin=pd.cut(scored.ev, [0, 1, 10, 100, 1000.1]).astype(str)).groupby(["group", "evalue_bin"]):
        rates.append({"group": grp, "evalue_bin": eb, "n_both_scored": int(len(sub)),
                      "n_tied": int((sub["class"] == "both scored, tied at the best xcorr").sum()),
                      "share_tied": float((sub["class"] == "both scored, tied at the best xcorr").mean())})
    for grp, sub in scored.groupby("group"):
        rates.append({"group": grp, "evalue_bin": "all", "n_both_scored": int(len(sub)),
                      "n_tied": int((sub["class"] == "both scored, tied at the best xcorr").sum()),
                      "share_tied": float((sub["class"] == "both scored, tied at the best xcorr").mean())})
    rates = pd.DataFrame(rates)
    # sensitivity: no edge margin and no methionine allowance
    alt0 = e + sign * s_hi
    t["class_no_margin"] = np.select([t.chemistries_at_best_xcorr == "dioxidation+sulfide", np.abs(alt0) <= WINDOW],
                                     ["tied", "scored lower"], default="outside")
    sens = (t[t.class_no_margin != "outside"].groupby("group").class_no_margin
            .apply(lambda s: float((s == "tied").mean())).rename("share_tied_no_margin").reset_index())
    sens["n_both_scored_no_margin"] = t[t.class_no_margin != "outside"].groupby("group").size().values

    L.write_csv(arith, "fragment_bin_arithmetic.csv")
    L.write_csv(theo, "fragment_bin_theoretical_tie_peptides.csv")
    L.write_csv(tab, "single_cys_alternative_scored_classes.csv")
    L.write_csv(rates, "single_cys_tie_rate_when_both_scored.csv")
    L.write_csv(t[["arm", "scan", "peptide", "best_xcorr", "e_value_best", "chemistries_at_best_xcorr", "ppm_min",
                   "alt_error_ppm_no_metox", "n_met", "class", "class_no_margin"]], "single_cys_rows_classified.csv")
    summary = {
        "label": L.POSTHOC_LABEL,
        "arithmetic": arith.to_dict(orient="records"),
        "theoretical": theo.to_dict(orient="records"),
        "single_cys_classes": tab.to_dict(orient="records"),
        "single_cys_share_tied_when_both_scored": rates[rates.evalue_bin == "all"].to_dict(orient="records"),
        "single_cys_share_tied_by_evalue": rates[rates.evalue_bin != "all"].to_dict(orient="records"),
        "sensitivity_no_margin_no_met": sens.to_dict(orient="records"),
        "reading": ("with 0.02-Da bins a 0.017758-Da shift changes the bin of most singly charged cysteine-containing "
                    "fragments, and the fragment score separated the two candidates in most Fig-1D spectra in which "
                    "both were scored; with 1.0005-Da bins it rarely changes a bin and most such insulin-file "
                    "spectra tied. Ties in Fig-1D therefore reflect matches in which no discriminating fragment "
                    "was matched, not the bin width."),
    }
    L.write_json(summary, "fragment_bins_summary.json")
    return summary, arith, theo, tab, rates, sens


if __name__ == "__main__":
    import json
    s, arith, theo, tab, rates, sens = main()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(arith.to_string()); print(theo.to_string()); print(tab.to_string()); print(rates.to_string()); print(sens.to_string())
