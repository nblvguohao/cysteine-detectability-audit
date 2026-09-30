"""POST-HOC (not pre-registered; written after reading the X1 cross-check): which stored verdicts are decided by
the baseline alone and which by the control step, and how that lines up with Cys-Audit agreement.

Classification (from the shared verdict rule's branch order, applied to stored values; fixed before running):
  baseline_decided  claim_direction positive: baseline CI covers 0 (undecidable at baseline) or baseline CI excludes
                    0 with the opposite sign (baseline_contradicts_claim);
                    claim_direction null: baseline CI excludes 0 (the author's null does not reproduce)
  control_decided   everything else: the verdict depends on the controlled / random-control intervals
Cross-tabulated against results/phase5_cys_audit_crosscheck_2026-09-22.csv agreement (site-level claims only).
Output: results/phase5_claim_control_dependence_2026-09-22.csv and _audit.json. Standard library only.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M = os.path.join(ROOT, "results", "claim_audit_2026-09-22.csv")
X = os.path.join(ROOT, "results", "phase5_cys_audit_crosscheck_2026-09-22.csv")
OUT = os.path.join(ROOT, "results", "phase5_claim_control_dependence_2026-09-22.csv")
AUD = os.path.join(ROOT, "results", "phase5_claim_control_dependence_2026-09-22_audit.json")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    for p in (OUT, AUD):
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    xm = {r["claim_id"]: r for r in csv.DictReader(open(X, encoding="utf-8"))}
    rows = []
    for r in csv.DictReader(open(M, encoding="utf-8")):
        if r["verdict_status"] == "NO_BASELINE":
            continue
        b = float(r["effect_original"])
        lo, hi = json.loads(r["effect_original_ci"])
        covers = lo <= 0 <= hi
        if r["claim_direction"] == "null_no_preference":
            decided = "baseline_decided" if not covers else "control_decided"
        else:
            decided = "baseline_decided" if (covers or b < 0) else "control_decided"
        x = xm.get(r["claim_id"], {})
        rows.append({"claim_id": r["claim_id"], "verdict": r["verdict"], "verdict_status": r["verdict_status"],
                     "claim_direction": r["claim_direction"], "decided_by": decided,
                     "cys_audit_status": x.get("tool_status", ""), "agreement": x.get("agreement", "")})
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    site = [r for r in rows if r["agreement"]]
    tab = Counter((r["decided_by"], r["agreement"]) for r in site)
    ctrl_non_pass = [r["claim_id"] for r in site if r["decided_by"] == "control_decided" and r["verdict_status"] != "PASS"]
    audit = {"script": "scripts/posthoc_claim_control_dependence_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
             "posthoc": True, "inputs": {os.path.relpath(p, ROOT): sha(p) for p in (M, X)},
             "decided_by_all28": dict(Counter(r["decided_by"] for r in rows)),
             "site_level_crosstab": {f"{a}|{b}": n for (a, b), n in sorted(tab.items())},
             "control_decided_non_pass_site_claims": ctrl_non_pass,
             "outputs": {os.path.relpath(OUT, ROOT): sha(OUT)}}
    json.dump(audit, open(AUD, "w"), indent=1)
    print(json.dumps(audit, indent=1))


if __name__ == "__main__":
    main()
