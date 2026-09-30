"""Phase 6c: source-data files renamed to the R27 figure numbers (copies, hash-verified; originals untouched).

Mapping is the one in scripts/revise_manuscript_submission_r27_2026-09-22.py: old 2->3, 3->5, 4->6, 5->7, 6->4, 7->2.
Figure 1 is a schematic and has no source data.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "source_data_r27")
AUD = os.path.join(ROOT, "results", "source_data_r27_2026-09-22_audit.json")
MOVES = [("source_data_r24", "Source_Data_Fig2_search_space.csv", "Source_Data_Fig3_search_space.csv"),
         ("source_data_r24", "Source_Data_Fig3_three_axes.csv", "Source_Data_Fig5_three_axes.csv"),
         ("source_data_r24", "Source_Data_Fig4_claim_retests.csv", "Source_Data_Fig6_claim_retests.csv"),
         ("source_data_r24", "Source_Data_Fig5_self_audit_public.csv", "Source_Data_Fig7_self_audit_public.csv"),
         ("source_data_r25", "Source_Data_Fig6_synthetic_ground_truth.csv", "Source_Data_Fig4_synthetic_ground_truth.csv"),
         ("source_data_r25", "Source_Data_Fig7_public_four_protease.csv", "Source_Data_Fig2_public_four_protease.csv")]


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    if os.path.exists(OUT) or os.path.exists(AUD):
        sys.exit("REFUSE: source_data_r27 or its audit exists")
    os.makedirs(OUT)
    rows = {}
    for d, old, new in MOVES:
        s = os.path.join(ROOT, d, old)
        if not os.path.exists(s):
            sys.exit(f"REFUSE: missing {s}")
        t = os.path.join(OUT, new)
        shutil.copy2(s, t)
        if sha(s) != sha(t):
            sys.exit(f"REFUSE: copy differs for {new}")
        rows[new] = {"from": os.path.join(d, old), "sha256": sha(t)}
    nums = sorted(int(n.split("_Fig")[1][0]) for n in rows)
    if nums != [2, 3, 4, 5, 6, 7]:
        sys.exit(f"REFUSE: figure numbers covered are {nums}, expected 2..7")
    with open(AUD, "w", encoding="utf-8") as fh:
        json.dump({"script": "scripts/build_source_data_r27_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
                   "files": rows, "figure_numbers_covered": nums,
                   "note": "Figure 1 is a schematic and has no source data"}, fh, indent=1)
    print(json.dumps({"files": len(rows), "numbers": nums}, indent=1))


if __name__ == "__main__":
    main()
