"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific: revision after verification (round 1).

How far do the counts/shares/differences computed on the FASTA-filtered tables (UniProt 2026_03; 43 HydP and 6 HydN
rows dropped) differ from those on the unfiltered deposit-converted tables? The first version of the report stated
'within 0.002 for every share and difference' without computing it over every table; the verifier found 0.0037. This
script computes the maximum absolute difference of every point estimate, per table and quantity, and names the row.

Inputs: the c02/c04/c05/c08 outputs of both bases (results/<ID>/ and results/<ID>/unfiltered_base/).
Output: t9_base_agreement.csv (one row per table x quantity), t9_base_agreement_rows.csv (every compared row).
Point estimates only; nothing is resampled.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402

F, U = C.RES, f"{C.RES}/unfiltered_base"

# table, key columns, value columns, class ('descriptive' marks quantities not proposed for the manuscript text)
SPECS = [
    ("t2_specific_counts.csv", ["arm", "hydn_source", "match_key"],
     ["frac_sites_also_cam_in_hydn", "frac_testable_sites_cam_in_hydn", "frac_hydp_identified_testable"]),
    ("t2_hydn_summary.csv", ["arm"], ["hydp_share", "hydn_share"]),
    ("t4_coincidence_specific.csv", ["arm", "definition", "denominator"], ["coincidence"]),
    ("t4_arm_contrasts.csv", ["definition", "denominator", "arm"], ["difference"]),
    ("t5_deciles.csv", ["arm", "definition", "decile"], ["coincidence"]),
    ("t5_depth_controls.csv", ["definition", "arm", "control"],
     ["difference", "coincidence_arm", "coincidence_trypsin_adjusted"]),
    ("t8_hydn_vs_hydp_controls.csv", ["arm", "control"],
     ["difference_hydn_minus_hydp", "cam_share_hydn", "cam_share_hydp_compared"]),
    ("t8_hydn_vs_hydp_deciles.csv", ["arm", "hydroxylamine", "decile"], ["cam_share"]),
]
DESCRIPTIVE_DENOMINATORS = {"proteins_absent_from_all_hydn", "proteins_present_in_some_hydn"}


def klass(table, row):
    if table in ("t5_deciles.csv", "t8_hydn_vs_hydp_deciles.csv"):
        return "per-decile (descriptive)"
    if "denominator" in row and row["denominator"] in DESCRIPTIVE_DENOMINATORS:
        return "protein-presence denominators (descriptive)"
    if "arm" in row and table.startswith("t8") and row["arm"] != "Trypsin":
        return "non-tryptic with vs without hydroxylamine (22-68 HydN cysteines)"
    return "main"


def main():
    C.record_inputs("c09_base_agreement", [f"{b}/{t}" for b in (F, U) for t, _, _ in SPECS])
    per_row, summ = [], []
    for t, keys, vals in SPECS:
        a = pd.read_csv(f"{F}/{t}", dtype={"decile": str})
        b = pd.read_csv(f"{U}/{t}", dtype={"decile": str})
        m = a.merge(b, on=keys, suffixes=("_filtered", "_unfiltered"), how="outer", indicator=True)
        assert (m["_merge"] == "both").all(), f"{t}: rows differ between bases"
        for v in vals:
            d = (m[f"{v}_filtered"] - m[f"{v}_unfiltered"]).abs()
            for i, r in m.iterrows():
                rec = {k: r[k] for k in keys}
                per_row.append({"table": t, "quantity": v, "key": "|".join(str(r[k]) for k in keys),
                                "filtered": r[f"{v}_filtered"], "unfiltered": r[f"{v}_unfiltered"],
                                "abs_difference": d.iloc[i], "class": klass(t, rec)})
    rows = pd.DataFrame(per_row)
    rows.to_csv(f"{C.RES}/t9_base_agreement_rows.csv", index=False)
    for (t, q, k), g in rows.groupby(["table", "quantity", "class"], sort=False):
        g2 = g.dropna(subset=["abs_difference"])
        if g2.empty:
            continue
        i = g2.abs_difference.idxmax()
        summ.append({"table": t, "quantity": q, "class": k, "n_compared": len(g2),
                     "max_abs_difference": g2.abs_difference.max(), "median_abs_difference": g2.abs_difference.median(),
                     "row_at_max": g2.loc[i, "key"], "filtered_at_max": g2.loc[i, "filtered"],
                     "unfiltered_at_max": g2.loc[i, "unfiltered"]})
    s = pd.DataFrame(summ)
    s.to_csv(f"{C.RES}/t9_base_agreement.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    pd.set_option("display.max_colwidth", 70)
    print(s.round(5).to_string(index=False))
    for k, g in rows.dropna(subset=["abs_difference"]).groupby("class"):
        i = g.abs_difference.idxmax()
        print(f"class {k}: max |filtered - unfiltered| = {g.abs_difference.max():.5f} at {g.loc[i, 'table']} "
              f"{g.loc[i, 'quantity']} {g.loc[i, 'key']}")


if __name__ == "__main__":
    main()
