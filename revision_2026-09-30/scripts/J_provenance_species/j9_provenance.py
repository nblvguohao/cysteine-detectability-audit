# -*- coding: utf-8 -*-
"""J9. Merge the per-script provenance fragments of item J into results/J_provenance_species/provenance.json.

POST HOC revision analysis (2026-09-30); not registered, not pre-specified.
Run order: j1_species_evidence.py, j2_artifact12_registry.py, j3_ablation_statistic.py,
j4_deposit_manifest.py, j5_singlecys_construction.py and j6_threshold_free_branch_test.py (both added in
revision round 2), j7_round3_checks.py (added in revision round 3), j9_provenance.py. j5 reads outputs of j1 and
j3, so it runs after them; j6 tests the deposit-time branch of j3 on synthetic data held in memory; j7 reads outputs
of j1 and j5. The verify_r*.py scripts in the same folder belong to the adversarial verifier; they are hashed here
but not run by this item.

Revision round 3: the merge now asserts that every fragment records the same sha256 for a shared input and that
each recorded hash equals the file's hash at merge time. In round 2 the j2 fragment still carried the hash of an
earlier manuscript version (j2 had not been re-run after the 05:37 edit), and the merge silently kept it.
"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jcommon import OUT, POSTHOC_LABEL, SCRIPTS, sha256_file, versions, write_json  # noqa: E402

RUN_ORDER = ["j1_species_evidence.py", "j2_artifact12_registry.py", "j3_ablation_statistic.py",
             "j4_deposit_manifest.py", "j5_singlecys_construction.py", "j6_threshold_free_branch_test.py",
             "j7_round3_checks.py", "j9_provenance.py"]


def main():
    frags = []
    for p in sorted(glob.glob(os.path.join(OUT, "_provenance_fragment_*.json"))):
        with open(p, encoding="utf-8") as fh:
            frags.append(json.load(fh))
    ran = {f["script"] for f in frags}
    missing = [s for s in RUN_ORDER[:-1] if s not in ran]
    assert not missing, f"fragments missing for {missing}"
    inputs = {}
    conflicts = []
    for f in frags:
        for k, v in f["inputs"].items():
            if k in inputs and inputs[k]["sha256"] != v["sha256"]:
                conflicts.append((k, inputs[k]["read_by"], f["script"]))
            inputs.setdefault(k, dict(v))
            inputs[k].setdefault("read_by", [])
            if f["script"] not in inputs[k]["read_by"]:
                inputs[k]["read_by"].append(f["script"])
    assert not conflicts, f"fragments disagree on input hashes (re-run the stale script): {conflicts}"
    live_prefix = os.path.normpath(os.path.join(os.path.dirname(OUT), "..", "reports")).replace("\\", "/")
    stale, live_changed = [], []
    for k, v in inputs.items():
        changed = sha256_file(k) != v["sha256"]
        if k.startswith(live_prefix + "/"):
            # other items' reports, read by j7 part A4 only: live documents of concurrent items
            v["live_document_of_another_item"] = True
            v["changed_since_read"] = changed
            if changed:
                live_changed.append(k)
        elif changed:
            stale.append(k)
    assert not stale, f"inputs changed after their fragment was written (re-run): {stale}"
    if live_changed:
        print("note: other items' reports changed after j7 read them (expected; recorded):", live_changed)
    all_py = sorted(glob.glob(os.path.join(SCRIPTS, "*.py")))
    prov = {
        "item": "J_provenance_species",
        "date": "2026-09-30",
        "label": POSTHOC_LABEL,
        "network_access": "none",
        "run_order": RUN_ORDER,
        "n_inputs": len(inputs),
        "scripts": {os.path.basename(p): sha256_file(p) for p in all_py if not os.path.basename(p).startswith("verify_")},
        "verifier_scripts_not_run_by_this_item": {os.path.basename(p): sha256_file(p) for p in all_py
                                                  if os.path.basename(p).startswith("verify_")},
        "inputs": dict(sorted(inputs.items())),
        "outputs": {k: v for f in frags for k, v in f["outputs"].items()},
        "seeds": {f["script"]: f["seeds"] for f in frags},
        "notes": {f["script"]: f["notes"] for f in frags},
        "versions": versions(),
        "worker_processes": "j1 uses 4 worker processes for the substring search; all other steps are single-process",
        "other_items_read_not_modified": ["results/F_rice_artifact3/fetched/fetched_sequences.tsv (rice sequences of active "
                                          "UniProtKB 2026_03 entries outside UP000059680, fetched by item F)",
                                          "results/F_rice_artifact3/s05_filter_chain.csv (cited in round 2; superseded by item F's "
                                          "round-2 chain)",
                                          "results/F_rice_artifact3/s11_rice_chain.csv (cited in round 3: 1,769 -> 1,753 -> 1,377 -> "
                                          "1,377 -> 1,340 -> 817 sites)",
                                          "reports/<other items>.md (round 3, j7 part A4: line pointers compared read-only; these "
                                          "are live documents of concurrent items and may change after this run)"],
        "hash_checks_at_merge": "every input recorded by more than one fragment has one sha256, and every recorded sha256 "
                                "equals the file's sha256 when this merge ran",
        "determinism": "fixed seeds (20260930 for the j1 bootstraps, the j3 self-test and the j3 rice protein-draw reference; "
                       "20260921 for the j3 threshold-free bootstrap, the seed of the registered refit analysis; 20260930, +1 and +2 "
                       "for the j5 protein-clustered bootstraps; 20260930 for the synthetic j6 test; 20260930 + 3 for the j7 "
                       "protein-clustered bootstrap); "
                       "multiprocessing in j1 only parallelises a pure membership test, so outputs do not depend on scheduling",
        "revision": "fourth version of the report (revision round 3, after adversarial verification round 3), 2026-09-30; "
                    "round 2 added j5_singlecys_construction.py and j6_threshold_free_branch_test.py and a "
                    "coverable-only (symmetric-filter) variant to the deposit-time branch of j3; round 3 made j2 locate "
                    "manuscript lines by anchor phrases, re-ran j2 and j3 on the manuscript version saved at 05:37 "
                    "(sha256 79ce7159...), added j7_round3_checks.py and the hash checks of this merge",
    }
    write_json(os.path.join(OUT, "provenance.json"), prov)
    print("inputs:", len(inputs), "outputs:", len(prov["outputs"]))


if __name__ == "__main__":
    main()
