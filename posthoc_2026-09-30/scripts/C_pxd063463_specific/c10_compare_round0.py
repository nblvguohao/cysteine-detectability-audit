"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific: revision after verification (round 1).

Check that the re-run of the extended pipeline leaves every statistic of the first version unchanged: every table in
results/<ID>/_round0/ (snapshot taken before the re-run) is compared with the current table on its shared columns and
rows. New columns (c02 tryptic-HydN columns) and new rows (c03 selection-split sets) are allowed; any changed value in
a shared cell is reported. Output: t10_round0_comparison.csv (one row per table).
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402

KEYS = {"t3_cleavage_specific.csv": ["arm", "positive_set", "handling", "rule", "background", "band"]}


def main():
    rows = []
    if not os.path.isdir(f"{C.RES}/_round0"):
        print("no _round0 snapshot: nothing to compare (the snapshot exists only for the round-1 revision re-run)")
        return
    snap = sorted(glob.glob(f"{C.RES}/_round0/**/*", recursive=True))
    for p in snap:
        if not os.path.isfile(p) or not p.endswith((".csv", ".tsv")):
            continue
        rel = os.path.relpath(p, f"{C.RES}/_round0").replace("\\", "/")
        cur = f"{C.RES}/{rel}"
        sep = "\t" if p.endswith(".tsv") else ","
        a = pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False)
        b = pd.read_csv(cur, sep=sep, dtype=str, keep_default_na=False)
        name = os.path.basename(p)
        shared = [c for c in a.columns if c in b.columns]
        if name in KEYS:
            k = KEYS[name]
            m = a.merge(b[shared], on=k, how="left", suffixes=("_old", "_new"), indicator=True)
            missing = int((m["_merge"] != "both").sum())
            diff = sum(int((m[f"{c}_old"] != m[f"{c}_new"]).sum()) for c in shared if c not in k)
            rows.append({"table": rel, "rows_old": len(a), "rows_new": len(b), "old_rows_missing_now": missing,
                         "changed_shared_cells": diff, "columns_added": ";".join(c for c in b.columns if c not in a.columns)})
        else:
            same_len = len(a) == len(b)
            diff = int((a[shared].values != b[shared].values[:len(a)]).sum()) if same_len else -1
            rows.append({"table": rel, "rows_old": len(a), "rows_new": len(b), "old_rows_missing_now": 0 if same_len else -1,
                         "changed_shared_cells": diff, "columns_added": ";".join(c for c in b.columns if c not in a.columns)})
    out = pd.DataFrame(rows)
    out.to_csv(f"{C.RES}/t10_round0_comparison.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 80)
    print(out.to_string(index=False))
    print("ALL FIRST-VERSION VALUES UNCHANGED" if (out.changed_shared_cells == 0).all()
          and (out.old_rows_missing_now == 0).all() else "SOME VALUES CHANGED")


if __name__ == "__main__":
    main()
