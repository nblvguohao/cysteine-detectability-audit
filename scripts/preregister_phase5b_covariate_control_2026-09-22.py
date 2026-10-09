"""Phase 5b pre-registration: rerun the 18 site-level claims with Cys-Audit v0.2.0's covariate-based claim control.

Registered before scripts/run_phase5b_covariate_control_2026-09-22.py exists. NOT BLIND in one respect, stated up
front: the tree's stored stratified estimates (MH over propensity-score quintiles on VIS10, the same construction as
the new tool control) were read before this registration, for the five claims in question:
  PERS-009 0.3991 [0.2649, 0.5308]; SFE-006 0.1824 [-0.0333, 0.4017]; SFE-008 -0.0379 [-0.3442, 0.2845];
  SNO-004 0.5408 [0.0224, 1.1261]; SNO-012 -0.5799 [-1.5473, 0.3160].
The predictions below therefore test whether an independent implementation reproduces that construction and what
verdicts it yields; they are not blind forecasts. The stored verdicts use the 1:1 caliper-MATCHED estimate, which
the tool does not implement.

Frozen: every file under cys-audit/src and cys-audit/tests at v0.2.0 (sha256 below); v0.1.0 is archived in
cys-audit/archive/v0.1.0/ (hash-identical to the Phase 4 frozen list).

B0 backward compatibility (a failure stops): v0.2.0 WITHOUT covariates on the 18 Phase 5 exports
   (external/intake/phase5_claim_inputs_20260922/) must reproduce results/phase5_cys_audit_crosscheck_2026-09-22.csv:
   same claim status and verdict, and baseline and controlled estimates equal to 6 decimals, for all 18.
G1 cohort rebuild as in Phase 5 (a failure on > 3 claims stops); covariates exported = the cohort's VIS10 matrix
   (10 columns, names from FEATURE_NAMES filtered by VISIBILITY_ONLY, in FEATURE_NAMES order).
G3 implementation check (reported; failure = implementation differs, not a stop): |tool propensity AUC - stored
   propensity_auc| <= 0.001 for every claim (same objective: ridge logistic, C = 1, standardised covariates).
G4 implementation check (reported): |tool controlled estimate - stored stratified_log2_or| <= 0.05 for every claim.
Run: --protease trypsin --background proteome --claim-feature claim_attr --claim-covariates <10 VIS10 columns>,
   direction as in Phase 5, 5000 replicates, default seed; canonical interpreter 3.11.16.
Predictions (status = the four-level claim status; comparison target = the stored verdict_status):
   P1 of the five control-decided non-PASS claims, at least 4 now match the stored status; expected matches
      PERS-009 ATTENUATED, SFE-006 FAIL, SFE-008 UNDECIDABLE, SNO-004 FAIL; expected mismatch SNO-012 (its stratified
      interval covers 0, so a stratified control reads the null as surviving -> PASS, while the matched estimate
      breaks it).
   P2 of the 13 claims that agreed under v0.1.0, at least 12 still agree.
   P3 no claim where the tool gives FAIL or ATTENUATED and the stored status is PASS.
Reading fixed now: if P1 holds with SNO-012 as the mismatch, the stored SNO-012 verdict rests on matching versus
stratification and is reported as specification-dependent; nothing in the stored tables is changed.
Standard library only.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, "cys-audit")
OUT = os.path.join(ROOT, "protocols", "phase5b_covariate_control_preregistration_2026-09-22.json")
AUD = os.path.join(ROOT, "results", "phase5b_covariate_control_prereg_2026-09-22_audit.json")
RUNNER = os.path.join(ROOT, "scripts", "run_phase5b_covariate_control_2026-09-22.py")
FUTURE = [os.path.join(ROOT, "results", f) for f in ("phase5b_covariate_control_2026-09-22.csv",
                                                     "phase5b_covariate_control_2026-09-22_audit.json")]


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    for p in [OUT, AUD, RUNNER] + FUTURE:
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    frozen = {}
    for sub in ("src", "tests"):
        for dp, _, fs in os.walk(os.path.join(PKG, sub)):
            for f in sorted(fs):
                if f.endswith(".py"):
                    p = os.path.join(dp, f)
                    frozen[os.path.relpath(p, ROOT)] = sha(p)
    sys.path.insert(0, os.path.join(PKG, "src"))
    from cys_audit import constants as C
    protocol = {
        "phase": "5b", "date": "2026-09-22", "tool_version": C.VERSION, "frozen_files": frozen,
        "decision_constants": json.loads(json.dumps({k: getattr(C, k) for k in dir(C) if k.isupper()}, default=list)),
        "not_blind": {"stored_stratified_read_before_registration": {
            "PERS-009": [0.3991, 0.2649, 0.5308], "SFE-006": [0.1824, -0.0333, 0.4017],
            "SFE-008": [-0.0379, -0.3442, 0.2845], "SNO-004": [0.5408, 0.0224, 1.1261],
            "SNO-012": [-0.5799, -1.5473, 0.3160]}},
        "B0": {"reference": "results/phase5_cys_audit_crosscheck_2026-09-22.csv", "decimals": 6},
        "G1": {"stop_if_fail_more_than": 3}, "G3": {"tolerance": 0.001}, "G4": {"tolerance": 0.05},
        "run": {"protease": "trypsin", "background": "proteome", "bootstrap_reps": 5000, "interpreter": "3.11.16"},
        "predictions": {
            "P1": {"claims": ["PERS-009", "SFE-006", "SFE-008", "SNO-004", "SNO-012"], "min_matches": 4,
                   "expected_mismatch": "SNO-012"},
            "P2": {"min_still_agree": 12, "of": 13}, "P3": {"stricter_max": 0}},
        "archive_v010": "cys-audit/archive/v0.1.0/",
    }
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(protocol, fh, indent=1)
    audit = {"script": "scripts/preregister_phase5b_covariate_control_2026-09-22.py",
             "script_sha256": sha(os.path.abspath(__file__)), "protocol_sha256": sha(OUT), "n_frozen_files": len(frozen),
             "runner_absent_when_written": True}
    with open(AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=1)
    print(json.dumps(audit, indent=1))


if __name__ == "__main__":
    main()
