"""Phase 6a pre-registration: the text corrections R24 -> R25 (figures and packaging are Phase 6b).

Registered before scripts/revise_manuscript_submission_r25_2026-09-22.py exists. Nine edits, each with the finding
it closes and the product its numbers come from. Exact old/new strings live in the generator; this file fixes what
may change, what may not, and how the result is checked.

E1  F01 (must) Artefact 2: the headline pair +2.54 / -2.61 is presented as "multi-cysteine peptides invert the
    feature". Products: -2.61 is the ASYMMETRIC comparison; filtering both classes gives +3.34 (689 sites), the sign
    of the unfiltered estimate, and in the human dataset -8.78 asymmetric against +2.67 symmetric. Rewritten as
    decision D1 states: asymmetric single-Cys filtering reverses the sign.
    Numbers from baseline/manuscript_text_findings_2026-09-21.csv (F01), whose evidence column names the
    collaborator products 35_symmetric_ncys_test.json / 36_ncys_decompose.json / Q1c正式版_实测对照集.md.
E2  F02 (must) claim PERS-002 provenance: 1,923 negatives was corrected to 1,289 on 2026-09-15, so the stated
    "629-negative shortfall" does not exist; Round B differs from the settled cohort by five in each class.
    Product: results/qpers_sid_label_criteria_audit.json (corrections).
E3  F03 (must) Data availability carries an empty citation "refs **[ ]**". Rewritten to name the source without
    inventing reference numbers.
E4  F06 (should) Abstract attributes the 2,079 ties to ion-trap binning; 1,202 of them come from a high-resolution
    file (0.02 Da bins) and 877 from ion-trap arms (1.0005 Da). Both bin widths exceed 0.0178 Da, so the mechanism
    stands and only the attribution changes. Product: results/pxd015307_posthoc_score_ties.csv.
E5  F09 (should) Figure 1 legend says all five paired controls carry a same-size random control; only the claim
    re-tests do. Product: the legend's own source, COLLAB 12_matched_*.json as recorded in the Phase 0 finding.
E6  P1 + Phase 4: the manuscript says a non-tryptic digest "was not available" and lists its absence as a limit.
    Phase 4 ran the audit on PXD063463, where one enrichment was digested with four proteases: the cleavage
    signature follows the protease used (own-rule against trypsin-rule log2 OR: AspN 1.539092 / 0.118877, CT
    1.194606 / 0.401074, GluC 0.944881 / 0.383989; trypsin arm 0.972971 [0.711516, 1.275616]). The screen found no
    proteome-scale non-tryptic deposit for persulfidation, S-nitrosylation or sulfenylation ("not found" is not
    "absent"). Product: results/phase4_public_demo_2026-09-22.csv; PHASE_1_REPORT.md for the screen.
E7  P4: Artefact 4 reads as a general property. In the same four-protease deposit the trypsin arm's coincidence is
    0.255068 [0.238042, 0.272426] and the other three arms 0.56456-0.614973, so the coincidence is design-dependent.
    Product: results/phase4_public_demo_2026-09-22.csv.
E8+E9 P5/P5b + peer intel: `undecidable` is stated as an audit output carrying a reason per claim, and two verdicts
    (SNO-012, SFE-008) are shown to depend on the control specification rather than on power.
    Product: results/phase5b_covariate_control_2026-09-22.csv and its audit.

May not change: any other sentence; the abstract's four deliverables; Limitation 10's self-audit ceiling; the
verdict counts 12/11/2/1/1/1; the survey figures; SFE-006's author-cohort numbers (1.3840 / 0.5950 / 1.3918 /
100.6%) which E9 only cites earlier, never restates differently.

Gates in the generator (all must pass or it refuses to write):
  G1 R24 sha256 equals the value recorded here.
  G2 every number token that appears in R25 but not in R24 is in the generator's declared ADDED list, and each
     entry names the product file it comes from; every token that disappears is in the declared REMOVED list.
  G3 reverse substitution of the nine edits reconstructs R24 byte for byte.
  G4 required sentinels present in R25 (the corrected clauses).
  G5 forbidden strings absent from R25: "1,923 negatives", "629-negative", "below ion-trap fragment binning",
     "refs **[ ]**", "but none was available", "each with a same-size random control".
  G6 abstract <= 200 words; body and Methods word counts reported, not gated.
Standard library only.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(ROOT, "reports", "MANUSCRIPT_SUBMISSION_R24_2026-09-21.md")
OUT = os.path.join(ROOT, "protocols", "phase6a_manuscript_preregistration_2026-09-22.json")
AUD = os.path.join(ROOT, "results", "phase6a_manuscript_prereg_2026-09-22_audit.json")
RUNNER = os.path.join(ROOT, "scripts", "revise_manuscript_submission_r25_2026-09-22.py")
FUTURE = [os.path.join(ROOT, "reports", "MANUSCRIPT_SUBMISSION_R25_2026-09-22.md"),
          os.path.join(ROOT, "results", "manuscript_revise_r25_2026-09-22_audit.json")]
PRODUCTS = ["baseline/manuscript_text_findings_2026-09-21.csv", "results/qpers_sid_label_criteria_audit.json",
            "results/pxd015307_posthoc_score_ties.csv", "results/phase4_public_demo_2026-09-22.csv",
            "results/phase5b_covariate_control_2026-09-22.csv", "results/claim_audit_2026-09-22.csv"]


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    for p in [OUT, AUD, RUNNER] + FUTURE:
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    protocol = {
        "phase": "6a", "date": "2026-09-22", "source": os.path.relpath(MS, ROOT), "source_sha256": sha(MS),
        "target": "reports/MANUSCRIPT_SUBMISSION_R25_2026-09-22.md",
        "edits": ["E1_F01_artefact2_symmetry", "E2_F02_pers002_provenance", "E3_F03_empty_citation",
                  "E4_F06_abstract_binning", "E5_F09_figure1_legend", "E6_P1_nontryptic_control",
                  "E7_P4_artefact4_design_dependent", "E8_E9_undecidable_and_specification_dependence"],
        "products_sha256": {p: sha(os.path.join(ROOT, p)) for p in PRODUCTS},
        "gates": ["G1_source_unchanged", "G2_number_tokens_declared", "G3_reverse_reconstructs_source",
                  "G4_sentinels_present", "G5_forbidden_absent", "G6_abstract_le_200_words"],
        "may_not_change": ["four deliverables in the abstract", "Limitation 10 self-audit ceiling",
                           "verdict counts 12/11/2/1/1/1", "survey figures 6.8% and 47.3%",
                           "SFE-006 author-cohort values 1.3840 / 0.5950 / 1.3918 / 100.6%"],
        "figures_and_packaging": "Phase 6b, not this round",
    }
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(protocol, fh, indent=1, ensure_ascii=False)
    audit = {"script": "scripts/preregister_phase6a_manuscript_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
             "protocol_sha256": sha(OUT), "runner_absent_when_written": True}
    with open(AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=1)
    print(json.dumps(audit, indent=1))


if __name__ == "__main__":
    main()
