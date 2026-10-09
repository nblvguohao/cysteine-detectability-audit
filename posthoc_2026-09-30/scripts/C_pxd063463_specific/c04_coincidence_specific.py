"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific, task 4.

Coincidence (share of identified HydP cysteines that carry a site) under the stored and the hydroxylamine-specific
site definitions of task 2, per arm, with protein-clustered 95% percentile intervals (5000 replicates; the tool's
seed 20260922 + 401, so the stored definition reproduces Cys-Audit's positive_detected_overlap exactly), and the
contrast of each non-tryptic arm with the trypsin arm (paired protein bootstrap, seed 20260930 + offset).

Denominators
  identified          every cysteine identified in the HydP arm (the stored denominator)
  exclude_ambiguous   identified minus HydP CAM sites that fail the definition (sensitivity)
  testable_<src>      identified cysteines that were also identified in the HydN source arm(s), i.e. the cysteines
                      for which the no-hydroxylamine arm could have shown CAM
  proteins_absent_from_all_hydn / proteins_present_in_some_hydn (stored definition only, descriptive) identified
                      cysteines of proteins never / at least once identified in a HydN arm
Outputs: t4_coincidence_specific.csv, t4_arm_contrasts.csv
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402
import c_boot as B  # noqa: E402

DEFS = ["cam", "lenient_matched", "strict_matched", "lenient_pooled", "strict_pooled"]


def frames(st, d):
    """Yield (denominator name, sub-frame, success column) for a definition."""
    cam = st.cam_hydp == 1
    succ = cam if d == "cam" else st[d] == 1
    yield "identified", st.assign(_s=succ.astype(int))
    if d != "cam":
        amb = cam & ~succ
        yield "exclude_ambiguous", st.loc[~amb].assign(_s=succ[~amb].astype(int))
    for src in ("matched", "pooled"):
        if d == "cam" or d.endswith(src):
            m = st[f"in_hydn_{src}"] == 1
            yield f"testable_{src}", st.loc[m].assign(_s=succ[m].astype(int))
    if d == "cam":
        # descriptive: proteins never identified in any HydN arm (presence/absence stand-in for a protein-level call)
        m = st.protein_in_hydn_pooled == 0
        yield "proteins_absent_from_all_hydn", st.loc[m].assign(_s=succ[m].astype(int))
        yield "proteins_present_in_some_hydn", st.loc[~m].assign(_s=succ[~m].astype(int))


def main():
    C.record_inputs(f"c04_coincidence_specific[{C.BASE}]", [f"{C.OUT}/t2_site_table_{a}.tsv" for a in C.ARMS])
    st = {a: pd.read_csv(f"{C.OUT}/t2_site_table_{a}.tsv", sep="\t", dtype={"protein": str}) for a in C.ARMS}
    rows, cons = [], []
    off = 0
    for d in DEFS:
        per = {}
        for a in C.ARMS:
            for den, f in frames(st[a], d):
                pt, lo, hi = B.share_boot(f.protein, f._s.values == 1, C.REPS, C.TOOL_SEED + 401)
                rows.append({"arm": a, "definition": d, "denominator": den, "n_identified": len(f),
                             "n_sites": int(f._s.sum()), "n_proteins": f.protein.nunique(), "coincidence": pt,
                             "ci_low": lo, "ci_high": hi})
                per[(a, den)] = f
        for den in sorted({k[1] for k in per}):
            if ("Trypsin", den) not in per:
                continue
            for x in C.ARMS[1:]:
                if (x, den) not in per:
                    continue
                off += 1
                r = B.paired_contrast(per[("Trypsin", den)], per[(x, den)], "_s", C.REPS, C.NEW_SEED + 1000 + off)
                cons.append({"definition": d, "denominator": den, "arm": x, "reference": "Trypsin",
                             "coincidence_arm": r["cX"], "coincidence_trypsin": r["cT"],
                             "difference": r["diff_raw"], "diff_ci_low": r["diff_raw_ci_low"],
                             "diff_ci_high": r["diff_raw_ci_high"], "ratio": r["cX"] / r["cT"] if r["cT"] else np.nan,
                             "n_arm": r["n_X"], "n_trypsin": r["n_T"], "n_proteins_union": r["n_proteins_union"],
                             "seed": C.NEW_SEED + 1000 + off})
    df = pd.DataFrame(rows)
    df.to_csv(f"{C.OUT}/t4_coincidence_specific.csv", index=False)
    cf = pd.DataFrame(cons)
    cf.to_csv(f"{C.OUT}/t4_arm_contrasts.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    print(df.to_string(index=False))
    print(cf.to_string(index=False))


if __name__ == "__main__":
    main()
