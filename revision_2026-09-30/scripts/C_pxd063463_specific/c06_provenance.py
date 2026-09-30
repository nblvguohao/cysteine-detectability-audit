"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific: write results/<ID>/provenance.json.

Merges the per-step input records (sha256 of every input read), and adds library versions, seeds, the sha256 of
every script of this item, of the Cys-Audit source files that computed the statistics, and of every output file.
"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402


def main():
    parts = [json.loads(l) for l in open(C.PROV_FILE, encoding="utf-8")]
    inputs, steps = {}, []
    for p in parts:
        for path, h in p["inputs"].items():
            if path in inputs and inputs[path] != h:
                raise SystemExit(f"input changed between steps: {path}")
            inputs[path] = h
        steps.append({k: v for k, v in p.items() if k != "inputs"} | {"n_inputs": len(p["inputs"])})
    tool_files = sorted(glob.glob(f"{C.TOOL_SRC}/cys_audit/**/*.py", recursive=True))
    v010_files = sorted(glob.glob(f"{C.TOOL_V010}/cys_audit/**/*.py", recursive=True))
    # this item's own pipeline scripts; the adversarial verifier's scripts (verify_r1_*.py, verify_r2_*.py) and
    # outputs (results/<ID>/verify_r1/, verify_r2/) live in the same folders but are not part of this pipeline and
    # carry their own provenance (results/<ID>/verify_r*/provenance*.json); they are listed separately, not as outputs
    all_py = sorted(glob.glob(f"{C.SCRIPTS}/*.py"))
    scripts = [s for s in all_py if not os.path.basename(s).startswith("verify_")]
    verifier_scripts = [s for s in all_py if os.path.basename(s).startswith("verify_")]

    def excluded(p):
        q = p.replace("\\", "/")
        return ("/verify_r" in q or "/_round0/" in q or "/_round1/" in q
                or q.endswith(("provenance.json", "_provenance_parts.jsonl")))

    outputs = sorted(p for p in glob.glob(f"{C.RES}/**/*", recursive=True) if os.path.isfile(p) and not excluded(p))
    # snapshot of the first version's tables, taken before the revision after verification (round 1) re-run
    round0 = sorted(p for p in glob.glob(f"{C.RES}/_round0/**/*", recursive=True) if os.path.isfile(p))
    # round-1 provenance (sha256 of every round-1 output), copied before the round-2 changes
    round1 = f"{C.RES}/_round1/provenance_round1.json"
    import cys_audit
    rec = {
        "item": C.ITEM,
        "label": "POST HOC revision analysis in response to pre-submission review criticism (2026-09-30); "
                 "not registered, not pre-specified",
        "date": "2026-09-30",
        "environment": C.env_record(),
        "cys_audit": {"version_used": cys_audit.__version__, "source": C.TOOL_SRC,
                      "files_sha256": {os.path.relpath(f, C.TOOL_SRC).replace("\\", "/"): C.sha256_file(f)
                                       for f in tool_files},
                      "archived_v0_1_0_used_for_equivalence_check": C.TOOL_V010,
                      "v0_1_0_files_sha256": {os.path.relpath(f, C.TOOL_V010).replace("\\", "/"): C.sha256_file(f)
                                              for f in v010_files}},
        "seeds": {"cys_audit_default_seed_all_tool_statistics": C.TOOL_SEED,
                  "cleavage_band_seed_rule": "seed + 101 + band index (tool)",
                  "coincidence_seed_rule": "seed + 401 (tool); reused for every single-arm share in t4",
                  "new_analyses_base_seed": C.NEW_SEED,
                  "t2_share_intervals": "20260930 + running offset",
                  "t4_arm_contrasts": "20260930 + 1000 + running offset (listed per row)",
                  "t5_deciles": "20260930 + 100 + decile index",
                  "t5_controls": "20260930 + 2000 + running offset (listed per row)",
                  "t8_hydn_vs_hydp_controls": "20260930 + 3000 + running offset (listed per row)",
                  "t8_hydn_vs_hydp_deciles": "20260930 + 3500 + 20 * arm index + bin index (+11 for HydN)",
                  "t11_stratum_matched_cleavage": "tool seed 20260922 (cleavage.run; bands seed + 101 + j), as c03",
                  "t11_mh_stratified": "20260930 + 4000 + running offset (listed per row)",
                  "t12_coverage_proxy": "20260930 + 5000 + running offset (listed per row)",
                  "t1b_seed_sensitivity": [C.TOOL_SEED + 7919 * k for k in range(50)],
                  "bootstrap_replicates": C.REPS},
        "revision_after_verification_round_1": {
            "date": "2026-09-30",
            "added_scripts": ["c08_hydn_vs_hydp_depth.py", "c09_base_agreement.py", "c10_compare_round0.py"],
            "changed_scripts": ["c02_specific_sites.py (tryptic-HydN columns, t2_selection_summary.csv)",
                                "c03_cleavage_specific.py (selection split sets)",
                                "c07_print_tables.py (new tables)", "c06_provenance.py", "run_all.py"],
            "unchanged_statistics": "every statistic of the first version is recomputed unchanged (same seeds)"},
        "revision_after_verification_round_2": {
            "date": "2026-09-30",
            "problem": "selection split (c03 in_/not_in_trypsin_hydn) compared selected site subsets with an unselected "
                       "background; re-run with backgrounds restricted to the same identification stratum",
            "added_scripts": ["c11_stratum_matched.py", "c12_coverage_proxy.py", "c13_compare_round1.py"],
            "changed_scripts": ["run_all.py (steps c11-c13)", "c06_provenance.py (round-2 record, verify_r2 and _round1 "
                                "excluded from outputs)", "c07_print_tables.py (new tables)"],
            "unchanged_statistics": "every round-1 output is re-created byte-identical (t13_round1_comparison.csv); "
                                    "the c03 split rows with the unmatched background are kept, not overwritten",
            "round1_provenance_snapshot_sha256": C.sha256_file(round1) if os.path.exists(round1) else None},
        "verifier_scripts_sha256_not_part_of_pipeline": {os.path.basename(s): C.sha256_file(s)
                                                         for s in verifier_scripts},
        "round0_snapshot_sha256_first_version_outputs": {os.path.relpath(p, C.RES).replace("\\", "/"):
                                                         C.sha256_file(p) for p in round0},
        "steps": steps,
        "inputs_sha256": inputs,
        "scripts_sha256": {os.path.basename(s): C.sha256_file(s) for s in scripts},
        "outputs_sha256": {os.path.relpath(p, C.RES).replace("\\", "/"): C.sha256_file(p) for p in outputs},
        "notes": ["The UniProt 2026_03 FASTA is decompressed to results/<ID>/_work/ because Cys-Audit reads plain "
                  "FASTA; its sha256 is recorded both compressed (inputs) and decompressed (step c01 extra).",
                  "Stored originals used for comparison: repo/results/phase4_reports/*/audit.json and "
                  "inputs/repo_results/phase4_public_demo_2026-09-22.csv (read-only)."],
    }
    with open(f"{C.RES}/provenance.json", "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=1, ensure_ascii=False)
    print("provenance written:", len(inputs), "inputs,", len(outputs), "outputs")


if __name__ == "__main__":
    main()
