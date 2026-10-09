# -*- coding: utf-8 -*-
"""B_overadjust_sim step 1b (POST HOC): theoretical detectability per universe cysteine, with the
manuscript's own definition (Experimental Procedures): a residue is theoretically detectable if at
least one peptide of an in-silico tryptic digest (cleavage C-terminal to K/R, not before P; full
cleavage with up to 2 missed cleavages) contains it and has length 7-30 residues.

Used only for the Artifact-1 background-restriction control (within-protein z and Mantel-Haenszel
log2 OR against all cysteines vs against theoretically detectable cysteines), evaluated in the same
simulated tables. Writes results/B_overadjust_sim/universe_theor_detectable.csv.gz (same row order
as universe_cysteines.csv.gz).
"""
from __future__ import annotations

import bisect
import gzip
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_b as cb  # noqa: E402

import pandas as pd  # noqa: E402

import phase2b_claim_cohorts as p2b  # noqa: E402  (repo, read-only)
from run_cross_protease_detectability_probe import PROTEASES, boundaries_for  # noqa: E402

OUT = os.path.join(cb.RESULTS, "universe_theor_detectable.csv.gz")


def main():
    u = pd.read_csv(cb.UNIVERSE, dtype={"protein": str}, usecols=["protein", "position"])
    fasta = p2b.read_fasta_gz(cb.FASTA)
    rule = PROTEASES["Trypsin"]
    flags = []
    cache = {}
    for acc, pos in zip(u["protein"], u["position"]):
        if acc not in cache:
            cache[acc] = boundaries_for(fasta[acc], rule)
        b = cache[acc]
        i = pos - 1
        k = min(max(bisect.bisect_right(b, i) - 1, 0), len(b) - 2)  # fragment holding residue i
        ok = 0
        for lo in range(max(0, k - 2), k + 1):          # first fragment of the peptide
            for hi in range(k, min(len(b) - 2, lo + 2) + 1):  # last fragment, <= 2 missed cleavages
                L = b[hi + 1] - b[lo]
                if 7 <= L <= 30:
                    ok = 1
                    break
            if ok:
                break
        flags.append(ok)
    out = u.copy()
    out["theor_detectable_ms"] = flags
    out.to_csv(OUT, index=False, compression="gzip")
    print("theoretically detectable fraction", out["theor_detectable_ms"].mean(), "n", len(out))


if __name__ == "__main__":
    main()
