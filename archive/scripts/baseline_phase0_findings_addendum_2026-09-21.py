"""Phase 0 baseline, part 5b: one text finding found after part 5 had written its outputs (appended).

Same rules as scripts/baseline_phase0_summarise_2026-09-21.py (checks imported, not rewritten):
the finding is written only if the manuscript literal is on its line and every data check passes.

F20 (should) Methods 'Per-residue structural features' says secondary structure was assigned by
DSSP and the Software paragraph lists pydssp. The structural claim re-tests the Results report
(round phase2e: SFE-001, SNO-001, SNO-002, SNO-009, SNO-012) read features built by
scripts/build_structural_site_features.py, whose audit records biotite annotate_sse (P-SEA) and
biotite Shrake-Rupley because the DSSP binary was unavailable. pydssp is used by the collaborator
tree's large-scale structural scripts (15, 21), whose results R24 no longer reports.

Output (new name): baseline/manuscript_text_findings_2026-09-21_addendum.csv
"""
from __future__ import annotations

import csv
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("s", os.path.join(ROOT, "scripts", "baseline_phase0_summarise_2026-09-21.py"))
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
OUT = os.path.join(ROOT, "baseline", "manuscript_text_findings_2026-09-21_addendum.csv")

F20 = ("F20", "should", "Methods, Per-residue structural features", 366, "Secondary structure was assigned by DSSP",
       "Secondary structure was assigned by DSSP; Software lists pydssp 0.9.1.",
       "The structural claim re-tests reported in Results (phase2e) used biotite annotate_sse (P-SEA) and biotite "
       "Shrake-Rupley, because the DSSP binary was unavailable; pydssp appears only in the collaborator's large-scale "
       "structural scripts (15, 21), whose results R24 does not report.",
       [("json", "this", "results/structural_site_features_audit.json", ["sse_method"],
         "biotite annotate_sse (P-SEA); DSSP binary unavailable in this sandbox"),
        ("json", "this", "results/phase2e_claim_retest_audit.json", ["structural_features", "script"],
         "scripts/build_structural_site_features.py"),
        ("text", "collab", "scripts/15_structural_features_large.py", "import pydssp", True)],
       "results/structural_site_features_audit.json; results/phase2e_claim_retest_audit.json")


def main():
    if os.path.exists(OUT):
        sys.exit(f"REFUSE: {OUT} exists")
    with open(S.MS, encoding="utf-8") as fh:
        L = fh.read().split("\n")
    fid, sev, sec, ln, lit, says, shows, checks, ev = F20
    ms_ok = any(lit in L[i] for i in range(max(0, ln - 2), min(len(L), ln + 1)))
    res = [bool(S.run_check(c)) for c in checks]
    if not ms_ok or not all(res):
        sys.exit(f"F20 not written: manuscript literal {ms_ok}, checks {res}")
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "severity", "section", "md_line", "manuscript_literal",
                                           "what_the_text_says", "what_the_products_show",
                                           "n_mechanical_checks_passed", "evidence"])
        w.writeheader()
        w.writerow({"id": fid, "severity": sev, "section": sec, "md_line": ln, "manuscript_literal": lit,
                    "what_the_text_says": says, "what_the_products_show": shows,
                    "n_mechanical_checks_passed": len(res), "evidence": ev})
    print("F20 written; checks", res, "script sha256", S.sha(os.path.abspath(__file__)))


if __name__ == "__main__":
    main()
