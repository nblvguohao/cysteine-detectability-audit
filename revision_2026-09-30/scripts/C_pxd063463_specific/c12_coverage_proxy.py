"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific: revision after verification (round 2).

Within-cysteine depth cannot be measured in the converted tables (no per-cysteine evidence counts, no intensities).
The verifier (round 2) used a protein-level proxy: how many cysteines of the protein each arm identified. Where one
arm identifies more cysteines of a protein, it probably covers that protein's peptides with more evidence, and a
cysteine counts as a site when at least one evidence row carries CAM, so the proxy points to the direction of any
within-protein depth bias in a same-cysteine comparison. Reproduced here with this item's own code.

For each pair of arms (A, B) and the cysteines both identified (the same-cysteine controls of c05 and c08):
  coverage ratio r = (identified cysteines of the protein in A) / (identified cysteines of the protein in B)
  strata: all; A_more (r > 1); equal (r = 1); B_more (r < 1); A_not_more (r <= 1)
  per stratum: cysteines, CAM in A only, CAM in B only, exact two-sided binomial test on the discordant cysteines
  (scipy.stats.binomtest, p = 0.5), and for all / A_more / A_not_more the difference in CAM share B - A with a paired
  protein-clustered bootstrap 95% interval (c_boot.paired_matched, 5,000 replicates, seed 20260930 + 5000 + offset)
  tertiles of log r among discordant cysteines (pandas qcut; descriptive, as in the verifier's check)
Pairs: trypsin with vs without hydroxylamine (A = HydP, B = HydN; B - A is c08's 'without - with'); trypsin vs each
non-tryptic HydP arm (A = trypsin, B = the other arm; B - A is c05's 'arm - trypsin').
Both table bases are computed (unfiltered = base of the manuscript's counts; filtered = UniProt 2026_03).
Outputs: t12_coverage_proxy.csv, t12_coverage_tertiles.csv
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402
import c_boot as B  # noqa: E402

SEED0 = C.NEW_SEED + 5000
PAIRS = [("Trypsin", "HydP", "Trypsin", "HydN", "trypsin_with_vs_without_hydroxylamine")] + \
        [("Trypsin", "HydP", a, "HydP", f"trypsin_vs_{a}") for a in ("AspN", "CT", "GluC")]


def read(base, arm, hyd):
    p = C.raw_path(arm, hyd) if base == "unfiltered" else C.filtered_path(arm, hyd)
    d = pd.read_csv(p, sep="\t", dtype={"protein": str})
    d["position"] = d.position.astype(int)
    return d, p


def main():
    rows, trows, paths = [], [], []
    seed = SEED0
    for base in ("unfiltered", "filtered"):
        for aA, hA, aB, hB, name in PAIRS:
            A, pA = read(base, aA, hA)
            Bt, pB = read(base, aB, hB)
            paths += [pA, pB]
            cA, cB = A.groupby("protein").size(), Bt.groupby("protein").size()
            m = A[["protein", "position", "label"]].merge(Bt[["protein", "position", "label"]], on=["protein", "position"],
                                                          suffixes=("_A", "_B"))
            m["ratio"] = m.protein.map(cA).values / m.protein.map(cB).values
            shared = sorted(set(A.protein) & set(Bt.protein))
            strata = {"all": np.ones(len(m), bool), "A_more": m.ratio.values > 1, "equal": m.ratio.values == 1,
                      "B_more": m.ratio.values < 1, "A_not_more": m.ratio.values <= 1}
            for sn, sm in strata.items():
                g = m[sm]
                a1b0 = int(((g.label_A == 1) & (g.label_B == 0)).sum())
                a0b1 = int(((g.label_A == 0) & (g.label_B == 1)).sum())
                p = binomtest(a1b0, a1b0 + a0b1, 0.5).pvalue if a1b0 + a0b1 else np.nan
                rec = {"base": base, "pair": name, "arm_A": f"{aA}_{hA}", "arm_B": f"{aB}_{hB}", "stratum": sn,
                       "n_cysteines": len(g), "n_proteins": g.protein.nunique(), "cam_A_only": a1b0, "cam_B_only": a0b1,
                       "cam_both": int(((g.label_A == 1) & (g.label_B == 1)).sum()),
                       "share_A_only_of_discordant": a1b0 / (a1b0 + a0b1) if a1b0 + a0b1 else np.nan,
                       "binomial_p_two_sided": p,
                       "cam_share_A": float(g.label_A.mean()) if len(g) else np.nan,
                       "cam_share_B": float(g.label_B.mean()) if len(g) else np.nan,
                       "median_ratio": float(g.ratio.median()) if len(g) else np.nan,
                       "identified_cys_on_shared_proteins_A": int(cA.loc[shared].sum()),
                       "identified_cys_on_shared_proteins_B": int(cB.loc[shared].sum())}
                if sn in ("all", "A_more", "A_not_more") and len(g):
                    r = B.paired_matched(g.rename(columns={"label_A": "T", "label_B": "X"}), "T", "X", C.REPS, seed)
                    rec.update({"difference_B_minus_A": r["diff"], "ci_low": r["diff_ci_low"], "ci_high": r["diff_ci_high"],
                                "seed": seed})
                    seed += 1
                rows.append(rec)
            d = m[m.label_A != m.label_B].copy()
            if len(d) >= 20:
                d["tertile"] = pd.qcut(np.log(d.ratio), 3, labels=False, duplicates="drop")
                for t, g in d.groupby("tertile"):
                    trows.append({"base": base, "pair": name, "tertile": int(t) + 1, "n_discordant": len(g),
                                  "cam_A_only": int((g.label_A == 1).sum()),
                                  "share_A_only": float((g.label_A == 1).mean()),
                                  "ratio_min": float(g.ratio.min()), "ratio_max": float(g.ratio.max())})
    C.record_inputs("c12_coverage_proxy", sorted(set(paths)))
    out = pd.DataFrame(rows)
    out.to_csv(f"{C.RES}/t12_coverage_proxy.csv", index=False)
    tout = pd.DataFrame(trows)
    tout.to_csv(f"{C.RES}/t12_coverage_tertiles.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    pd.set_option("display.max_rows", 200)
    print(out[["base", "pair", "stratum", "n_cysteines", "cam_A_only", "cam_B_only", "binomial_p_two_sided",
               "median_ratio", "difference_B_minus_A", "ci_low", "ci_high"]].round(4).to_string(index=False))
    print(tout.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
