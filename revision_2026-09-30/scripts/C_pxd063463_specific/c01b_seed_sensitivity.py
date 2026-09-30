"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific, task 1 supplement.

Monte Carlo (bootstrap-seed) variability of the Figure 2 interval bounds. The stored audits used one seed
(20260922); a 97.5% percentile bound from 5000 replicates carries Monte Carlo error, and a change of input rows
changes every draw. This re-runs Cys-Audit's cleavage check (unchanged) on the FASTA-filtered HydP tables with 50
seeds (20260922 + 7919*k, k = 0..49; k = 0 is the tool default) for all fourteen arm x rule x background cells and
reports, per band, the spread of the bounds and how often the interval covers zero or crosses the 0.5 margin.
Outputs: t1b_seed_bounds.csv (one row per cell x band x seed), t1b_seed_summary.csv
"""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402

SEEDS = [C.TOOL_SEED + 7919 * k for k in range(50)]
CELLS = []
for arm in C.ARMS:
    for rule in [C.OWN_RULE[arm]] + ([] if arm == "Trypsin" else ["trypsin"]):
        for bg in ("proteome", "observed"):
            CELLS.append((arm, rule, bg))


def run_cell(cell):
    import c_common as CC  # noqa: F401
    from cys_audit.checks import cleavage
    from cys_audit.io import read_sites
    arm, rule, bg = cell
    ds = read_sites(C.filtered_path(arm, "HydP"), fasta=C.FASTA, expand_background=True)
    rows = []
    for s in SEEDS:
        r = cleavage.run(ds, {"protease": rule, "background": bg, "seed": s, "reps": C.REPS})
        for b, t in r["details"]["bands"].items():
            rows.append({"arm": arm, "rule": rule, "background": bg, "band": b, "seed": s, "estimate": t["estimate"],
                         "ci_low": t["ci"][0], "ci_high": t["ci"][1], "status": t["status"]})
    return rows


def main():
    C.ensure_fasta()
    C.record_inputs("c01b_seed_sensitivity", [C.filtered_path(a, "HydP") for a in C.ARMS] + [C.FASTA],
                    extra={"seeds": SEEDS})
    rows = []
    with ProcessPoolExecutor(max_workers=4) as ex:
        for r in ex.map(run_cell, CELLS):
            rows.extend(r)
    df = pd.DataFrame(rows)
    df.to_csv(f"{C.RES}/t1b_seed_bounds.csv", index=False)
    g = df.groupby(["arm", "rule", "background", "band"], sort=False)
    summ = g.apply(lambda x: pd.Series({
        "estimate": x.estimate.iloc[0],
        "ci_low_default_seed": x.ci_low.iloc[0], "ci_high_default_seed": x.ci_high.iloc[0],
        "ci_low_min": x.ci_low.min(), "ci_low_max": x.ci_low.max(), "ci_low_sd": x.ci_low.std(ddof=1),
        "ci_high_min": x.ci_high.min(), "ci_high_max": x.ci_high.max(), "ci_high_sd": x.ci_high.std(ddof=1),
        "share_seeds_interval_covers_zero": ((x.ci_low <= 0) & (x.ci_high >= 0)).mean(),
        "share_seeds_status_FAIL": (x.status == "FAIL").mean(),
        "share_seeds_status_WARNING": (x.status == "WARNING").mean(),
        "share_seeds_status_UNDECIDABLE": (x.status == "UNDECIDABLE").mean(),
        "share_seeds_status_PASS": (x.status == "PASS").mean(),
        "n_seeds": len(x)}), include_groups=False).reset_index()
    summ.to_csv(f"{C.RES}/t1b_seed_summary.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(summ.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
