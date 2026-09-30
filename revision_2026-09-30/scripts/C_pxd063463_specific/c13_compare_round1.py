"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific: revision after verification (round 2).

Check that adding c11/c12 and re-running the whole pipeline leaves every output of the round-1 revision unchanged.
The round-1 provenance (outputs_sha256) was copied to results/<ID>/_round1/provenance_round1.json before any
round-2 change; every output listed there (except logs, which record run order) is re-hashed and compared.
Output: t13_round1_comparison.csv (one row per round-1 output file) and a one-line verdict on stdout.
"""
from __future__ import annotations

import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402

SNAP = f"{C.RES}/_round1/provenance_round1.json"


def main():
    if not os.path.exists(SNAP):
        print("no _round1 snapshot: nothing to compare")
        return
    old = json.load(open(SNAP, encoding="utf-8"))["outputs_sha256"]
    rows = []
    for rel, h in sorted(old.items()):
        if rel.startswith("logs/"):
            continue
        p = f"{C.RES}/{rel}"
        cur = C.sha256_file(p) if os.path.exists(p) else None
        rows.append({"file": rel, "round1_sha256": h, "current_sha256": cur, "identical": cur == h})
    out = pd.DataFrame(rows)
    out.to_csv(f"{C.RES}/t13_round1_comparison.csv", index=False)
    changed = out[~out.identical]
    print(f"{len(out)} round-1 outputs compared (logs excluded); {int(out.identical.sum())} byte-identical")
    if len(changed):
        print(changed.to_string(index=False))
    print("ALL ROUND-1 OUTPUTS UNCHANGED" if changed.empty else "SOME ROUND-1 OUTPUTS CHANGED")


if __name__ == "__main__":
    main()
