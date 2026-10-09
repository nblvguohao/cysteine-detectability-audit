# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC side check (2026-09-30). Not registered.

The released tool's claim_verdict (cys_audit.verdict, documented as the rule pre-declared for the
published-claim re-tests) requires the BASELINE interval to exclude 0 before a positive claim can
survive; the pipeline rule that produced Fig. 6 (classify + phase-2b guards) does not read the
baseline interval for positive claims. This applies the tool's rule to the stored intervals
(controlled = matched interval, the pipeline's primary caliber) and lists the differences.
Writes results/E_bootstrap_coverage/rule_difference.csv
"""
import os, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import simlib as S
from cys_audit.verdict import claim_verdict

OUT = r"/path/to/local/_cys_repo_work/public/revision_2026-09-30/results/E_bootstrap_coverage"
ci = pd.read_csv(os.path.join(OUT, "claims_intervals.csv"))
rows = []
for _, r in ci.iterrows():
    v = claim_verdict(r.baseline_point, [r.baseline_lo, r.baseline_hi], [r.matched_lo, r.matched_hi],
                      [r.random_lo, r.random_hi], r.matched_point, r.claim_direction)
    tool = {"reverses": ("null_broken_by_control" if r.claim_direction == "null_no_preference" else "reverses"),
            "survives": ("null_survives" if r.claim_direction == "null_no_preference" else "survives")}.get(v["verdict"], v["verdict"])
    rows.append({"claim_id": r.claim_id, "stored_verdict": r.stored_verdict, "tool_rule_verdict": tool,
                 "tool_reason": v["reason"], "baseline_excludes0": r.baseline_excludes0})
d = pd.DataFrame(rows)
d["differs"] = d.stored_verdict != d.tool_rule_verdict
d.to_csv(os.path.join(OUT, "rule_difference.csv"), index=False)
pd.set_option("display.width", 200)
print(d[d.differs].to_string())
print("differences:", int(d.differs.sum()))
