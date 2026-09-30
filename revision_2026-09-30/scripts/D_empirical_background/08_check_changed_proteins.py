"""Step 8 (POST HOC revision analysis, item D_empirical_background): check that the three proteins whose sequence
changed between the deposit-era FASTA and UniProt 2026_03 (P46718, Q3TX08, Q3UPF5; 9 identified cysteines) do not
move the restriction-design estimates. Point estimates only (Cys-Audit Haldane log2 OR).

Output: check_changed_proteins.csv
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

from common import RESULTS
from cys_audit import stats


def main():
    T = pd.read_csv(f"{RESULTS}/cys_scores_trypsin_arm.csv", dtype={"protein": str})
    rows = []
    for subset, U in (("all_1991_proteins", T), ("without_changed_sequences", T[T.seq_changed_since_2025 == 0])):
        lab, det = U.label.to_numpy(), U.detected.to_numpy()
        pos = lab == 1
        for d, bg in (("proteome", lab == 0), ("observed", (lab == 0) & (det == 1)),
                      ("theoretical", (lab == 0) & (U.theo_trypsin.to_numpy() == 1)),
                      ("empirical", (lab == 0) & (U.emp_obs.to_numpy() == 1))):
            r = {"subset": subset, "design": d, "n_positive": int(pos.sum()), "n_background": int(bg.sum())}
            for f in ("kr_dist", "kr_prox", "de_dist"):
                fl = U[f].to_numpy().astype(bool)
                r[f] = float(stats.log2_or((pos & fl).sum(), (pos & ~fl).sum(), (bg & fl).sum(), (bg & ~fl).sum()))
            rows.append(r)
    D = pd.DataFrame(rows)
    D.to_csv(f"{RESULTS}/check_changed_proteins.csv", index=False)
    print(D.round(4).to_string())


if __name__ == "__main__":
    main()
