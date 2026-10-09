"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific, task 5.

Is the higher coincidence of the non-tryptic arms (AspN, chymotrypsin, GluC) than of the trypsin arm explained by
their smaller depth (187, 724 and 880 identified cysteines against 5,293 after the FASTA filter)?

Abundance is the deposit's unenriched ('Global') protein intensity, carried on every row of the arm tables and
identical for a protein in every arm (checked below). It measures cellular abundance, not the amount of a protein in
the enriched sample, so it is only a partial proxy for depth.

(a) Coincidence by protein-abundance decile per arm. Decile edges are the deciles of log10 abundance over the union of
    proteins identified in the four HydP arms; proteins without abundance form an 11th bin ('missing').
(b) Trypsin standardised to each non-tryptic arm: (1) direct standardisation over the 10 decile bins + missing bin
    (primary); (2) 20 quantile bins + missing; (3) deciles, complete cases only; (4) protein-matched restriction (both
    arms restricted to proteins identified in both); (5) protein-matched with every protein weighted equally;
    (6) cysteine-matched (the same cysteines identified in both arms, CAM status compared within cysteine);
    (7) cysteine-matched, cysteines seen only on single-cysteine peptides in both arms (no localisation ambiguity).
    Intervals: paired protein bootstrap over the union of proteins, 5000 replicates, seed 20260930 + offset.
(c) Verdict per definition: does the difference survive each control (interval excludes 0)?
Definitions: cam (stored), lenient_pooled (hydroxylamine-specific, testability roughly equal across arms), and
lenient_matched (its matched-HydN filter is 2.6x to 13x weaker in the non-tryptic arms; reported for completeness).
Outputs: t5_deciles.csv, t5_bin_distribution.csv, t5_depth_controls.csv, t5_verdict.csv
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402
import c_boot as B  # noqa: E402

DEFS = ["cam", "lenient_pooled", "lenient_matched"]


def load():
    st = {}
    for a in C.ARMS:
        d = pd.read_csv(f"{C.OUT}/t2_site_table_{a}.tsv", sep="\t", dtype={"protein": str})
        d["cam"] = d.cam_hydp
        d["log10_abundance"] = np.log10(d.abundance)
        st[a] = d
    # abundance must be a protein property, identical across arms
    ab = pd.concat([st[a][["protein", "abundance"]] for a in C.ARMS]).dropna().drop_duplicates()
    assert ab.protein.is_unique, "abundance differs between arms for some protein"
    return st


def edges_from(st, q):
    prot = pd.concat([st[a][["protein", "log10_abundance"]] for a in C.ARMS]).dropna().drop_duplicates("protein")
    e = np.quantile(prot.log10_abundance.values, np.linspace(0, 1, q + 1))
    return e, len(prot)


def assign_bins(d, edges):
    q = len(edges) - 1
    b = np.full(len(d), q, dtype=int)          # missing bin = q
    m = d.log10_abundance.notna().values
    b[m] = np.clip(np.searchsorted(edges, d.log10_abundance.values[m], side="right") - 1, 0, q - 1)
    return b


def main():
    C.record_inputs(f"c05_depth_confound[{C.BASE}]", [f"{C.OUT}/t2_site_table_{a}.tsv" for a in C.ARMS])
    st = load()
    e10, nprot = edges_from(st, 10)
    e20, _ = edges_from(st, 20)
    for a in C.ARMS:
        st[a]["bin10"] = assign_bins(st[a], e10)
        st[a]["bin20"] = assign_bins(st[a], e20)

    # (a) deciles
    rows, dist = [], []
    for a in C.ARMS:
        d = st[a]
        for b in range(11):
            s = d[d.bin10 == b]
            dist.append({"arm": a, "bin": b if b < 10 else "missing", "n_identified": len(s),
                         "share_of_arm": len(s) / len(d), "n_proteins": s.protein.nunique(),
                         "log10_abundance_low": e10[b] if b < 10 else np.nan,
                         "log10_abundance_high": e10[b + 1] if b < 10 else np.nan})
            for dn in DEFS:
                if len(s) == 0:
                    continue
                pt, lo, hi = B.share_boot(s.protein, s[dn].values == 1, C.REPS, C.NEW_SEED + 100 + b)
                rows.append({"arm": a, "definition": dn, "decile": b + 1 if b < 10 else "missing",
                             "log10_abundance_low": e10[b] if b < 10 else np.nan,
                             "log10_abundance_high": e10[b + 1] if b < 10 else np.nan,
                             "n_identified": len(s), "n_sites": int(s[dn].sum()), "n_proteins": s.protein.nunique(),
                             "coincidence": pt, "ci_low": lo, "ci_high": hi})
    dec = pd.DataFrame(rows)
    dec.to_csv(f"{C.OUT}/t5_deciles.csv", index=False)
    pd.DataFrame(dist).to_csv(f"{C.OUT}/t5_bin_distribution.csv", index=False)

    # (b) controls
    ctl = []
    off = 0
    T = st["Trypsin"]
    for dn in DEFS:
        for x in C.ARMS[1:]:
            X = st[x]
            shared_p = set(T.protein) & set(X.protein)
            j = T.merge(X, on=["protein", "position"], suffixes=("_T", "_X"))
            specs = [
                ("raw", T, X, None, 1, False),
                ("abundance_deciles_plus_missing", T, X, "bin10", 11, False),
                ("abundance_20bins_plus_missing", T, X, "bin20", 21, False),
                ("abundance_deciles_complete_case", T[T.bin10 < 10], X[X.bin10 < 10], "bin10", 11, False),
                ("protein_matched", T[T.protein.isin(shared_p)], X[X.protein.isin(shared_p)], None, 1, False),
                ("protein_matched_protein_equal", T[T.protein.isin(shared_p)], X[X.protein.isin(shared_p)], None, 1, True),
            ]
            for name, dT, dX, bc, nb, peq in specs:
                off += 1
                r = B.paired_contrast(dT, dX, dn, C.REPS, C.NEW_SEED + 2000 + off, bin_col=bc, n_bins=nb, protein_equal=peq)
                ess = np.nan
                if bc is not None:
                    pX = np.bincount(dX[bc].values, minlength=nb) / len(dX)
                    pT = np.bincount(dT[bc].values, minlength=nb) / len(dT)
                    wts = (pX / np.where(pT > 0, pT, np.nan))[dT[bc].values]
                    ess = float(wts.sum() ** 2 / (wts ** 2).sum())
                ctl.append({"definition": dn, "arm": x, "control": name, "coincidence_arm": r["cX"],
                            "coincidence_arm_ci_low": r["cX_ci_low"], "coincidence_arm_ci_high": r["cX_ci_high"],
                            "coincidence_trypsin_raw": r["cT"], "coincidence_trypsin_adjusted": r["cT_std"],
                            "trypsin_adjusted_ci_low": r["cT_std_ci_low"], "trypsin_adjusted_ci_high": r["cT_std_ci_high"],
                            "difference": r["diff_std"], "diff_ci_low": r["diff_std_ci_low"],
                            "diff_ci_high": r["diff_std_ci_high"], "undefined_reps": r["diff_std_n_undefined_reps"],
                            "n_arm": r["n_X"], "n_trypsin": r["n_T"], "n_proteins_arm": r["n_proteins_X"],
                            "n_proteins_trypsin": r["n_proteins_T"], "trypsin_weight_ess": ess,
                            "seed": C.NEW_SEED + 2000 + off})
            for name, jj in (("cysteine_matched", j),
                             ("cysteine_matched_single_cys_peptides", j[(j.n_cys_peptide_T == 1) & (j.n_cys_peptide_X == 1)])):
                off += 1
                r = B.paired_matched(jj, f"{dn}_T", f"{dn}_X", C.REPS, C.NEW_SEED + 2000 + off)
                ctl.append({"definition": dn, "arm": x, "control": name, "coincidence_arm": r["cX"],
                            "coincidence_arm_ci_low": r["cX_ci_low"], "coincidence_arm_ci_high": r["cX_ci_high"],
                            "coincidence_trypsin_raw": r["cT"], "coincidence_trypsin_adjusted": r["cT"],
                            "trypsin_adjusted_ci_low": r["cT_ci_low"], "trypsin_adjusted_ci_high": r["cT_ci_high"],
                            "difference": r["diff"], "diff_ci_low": r["diff_ci_low"], "diff_ci_high": r["diff_ci_high"],
                            "undefined_reps": 0, "n_arm": r["n_cys"], "n_trypsin": r["n_cys"],
                            "n_proteins_arm": r["n_proteins"], "n_proteins_trypsin": r["n_proteins"],
                            "trypsin_weight_ess": np.nan, "seed": C.NEW_SEED + 2000 + off,
                            "discordant_trypsin1_arm0": r["discordant_T1_X0"],
                            "discordant_trypsin0_arm1": r["discordant_T0_X1"],
                            "concordant_11": r["concordant_11"], "concordant_00": r["concordant_00"]})
    cf = pd.DataFrame(ctl)
    cf["interval_excludes_zero"] = (cf.diff_ci_low > 0) | (cf.diff_ci_high < 0)
    raw = cf[cf.control == "raw"].set_index(["definition", "arm"]).difference
    cf["share_of_raw_difference_remaining"] = [d / raw[(dn, a)] for d, dn, a in zip(cf.difference, cf.definition, cf.arm)]
    cf.to_csv(f"{C.OUT}/t5_depth_controls.csv", index=False)
    ver = (cf.groupby(["definition", "control"])
             .apply(lambda g: pd.Series({"arms_difference_positive_and_excludes_zero": int(((g.diff_ci_low > 0)).sum()),
                                         "arms_tested": len(g),
                                         "min_difference": g.difference.min(), "max_difference": g.difference.max()}))
             .reset_index())
    ver.to_csv(f"{C.OUT}/t5_verdict.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    pd.set_option("display.max_columns", 40)
    print(pd.DataFrame(dist).pivot(index="bin", columns="arm", values="share_of_arm").round(3))
    print(dec[dec.definition == "cam"].pivot(index="decile", columns="arm", values="coincidence").round(3))
    print(dec[dec.definition == "cam"].pivot(index="decile", columns="arm", values="n_identified"))
    print(cf[["definition", "arm", "control", "coincidence_arm", "coincidence_trypsin_raw", "coincidence_trypsin_adjusted",
              "difference", "diff_ci_low", "diff_ci_high", "share_of_raw_difference_remaining", "n_arm", "n_trypsin",
              "trypsin_weight_ess"]].to_string(index=False))
    print(ver.to_string(index=False))
    print("decile edges (log10 Global intensity):", np.round(e10, 3).tolist(), "proteins used for edges:", nprot)


if __name__ == "__main__":
    main()
