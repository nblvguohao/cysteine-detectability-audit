"""H_artifact5_docs - run every step and write provenance.json (POST HOC revision analysis, 2026-09-30).

Order: step 1 (ties) -> step 2 (deposit +32) -> step 3 (mass accuracy, phospho/sulfo; reads step 2 output)
-> step 4 (parameter table; reads step 1 and 2 outputs) -> step 5 (fragment bins) -> step 6 (records consulted);
steps 5 and 6 were added in the revision after verification. No random numbers are drawn anywhere, so no seed
is needed; the run is deterministic given the inputs listed in provenance.json.

    python h5_run_all.py
"""
from __future__ import annotations

import datetime as dt
import glob
import importlib
import json
import os
import platform
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import h5_lib as L  # noqa: E402

STEPS = ["h5_01_ties", "h5_02_deposit_plus32", "h5_03_mass_accuracy_phospho", "h5_04_parameters",
         "h5_05_fragment_bins", "h5_06_records"]

# recorded script hashes in the repository's stored audits; the Windows checkout converts line endings,
# so the LF-normalised hash is the one that must match
RECORDED = {
    "scripts/run_pxd015307_research_mass_accuracy.py": ("results/pxd015307_research_audit.json", "script_sha256"),
    "scripts/posthoc_pxd015307_score_ties.py": ("results/pxd015307_posthoc_score_ties_audit.json", "script_sha256"),
    "scripts/audit_search_space_and_mass_accuracy.py": ("results/search_space_mass_accuracy_audit.json", "script_sha256"),
}


def main():
    os.makedirs(L.OUT, exist_ok=True)
    started = dt.datetime.now(dt.timezone.utc)
    for name in STEPS:
        mod = importlib.import_module(name)
        mod.main()
        print("done:", name)

    checks = []
    for script, (audit_rel, key) in RECORDED.items():
        spath = os.path.join(L.REPO, script)
        apath = os.path.join(L.REPO, audit_rel)
        L.register(spath, "repository script (hash check)", text=True)
        rec = L.read_json(apath, "repository audit (recorded script hash)")[key]
        checks.append({"script": script, "recorded_sha256": rec, "sha256_on_disk": L.sha256_file(spath),
                       "sha256_lf_normalised": L.sha256_lf(spath), "match_after_lf_normalisation": rec == L.sha256_lf(spath)})
    report_md = os.path.join(L.REPO, "reports", "PXD015307_RESEARCH_SEARCH_SPACE.md")
    L.register(report_md, "repository copy of the round report (hash check)", text=True)
    checks.append({"script": "reports/PXD015307_RESEARCH_SEARCH_SPACE.md",
                   "recorded_sha256": "371448867c42e14bbb7e77994dd1f000649be06ce487b54d0069fb9f5d30b376 (archived Supplementary Note 8)",
                   "sha256_on_disk": L.sha256_file(report_md), "sha256_lf_normalised": L.sha256_lf(report_md),
                   "match_after_lf_normalisation": L.sha256_lf(report_md) == "371448867c42e14bbb7e77994dd1f000649be06ce487b54d0069fb9f5d30b376"})

    import numpy
    import pandas
    import scipy
    sys.path.insert(0, L.CYS_AUDIT_SRC)
    import cys_audit
    outputs = []
    for p in sorted(glob.glob(os.path.join(L.OUT, "*"))):
        if (not os.path.isfile(p)) or os.path.basename(p) == "provenance.json" or p.endswith(".tmp"):
            continue
        outputs.append({"file": os.path.basename(p), "bytes": os.path.getsize(p), "sha256": L.sha256_file(p)})
    scripts = [{"file": os.path.basename(p), "sha256": L.sha256_file(p)}
               for p in sorted(glob.glob(os.path.join(L.SCRIPTS, "h5_*.py")))]
    verifier = {"scripts": [{"file": os.path.basename(p), "sha256": L.sha256_file(p)}
                            for p in sorted(glob.glob(os.path.join(L.SCRIPTS, "verify_r1_*.py")))],
                "outputs": [{"file": "verify_r1/" + os.path.basename(p), "sha256": L.sha256_file(p)}
                            for p in sorted(glob.glob(os.path.join(L.OUT, "verify_r1", "*"))) if os.path.isfile(p)],
                "note": "written by the adversarial verifier (round 1), not by this item's pipeline; "
                        "only verify_r1_deposit_sulfide_psms_baseline_by_peptide.csv is read, as a cross-check in step 2"}
    prov = {
        "item": L.ITEM,
        "status": L.POSTHOC_LABEL,
        "run_started_utc": started.isoformat(timespec="seconds"),
        "run_finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "host": platform.node(), "platform": platform.platform(),
        "python": sys.version,
        "libraries": {"numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__,
                      "cys_audit": getattr(cys_audit, "__version__", "unknown")},
        "seeds": "none: no step draws random numbers (all computations are deterministic counts, medians and closed-form masses)",
        "network": "not used",
        "steps": STEPS,
        "scripts": scripts,
        "inputs": L.inputs_registry(),
        "repository_hash_checks": checks,
        "outputs": outputs,
        "verifier_round1_files_present": verifier,
        "revision": "round 2 (revision after verification, 2026-09-30): steps 1-4 amended, steps 5 and 6 added",
        "notes": [
            "The raw files and the Comet output tables of the PXD015307 re-search are not in the local copy; nothing was re-searched.",
            "Target/decoy status of tie peptides uses UniProt mouse 2026_03 canonical sequences (the search used an earlier Swiss-Prot canonical+isoform release); isoform-only peptides remain 'unresolved'.",
            "Bovine insulin chain sequences used for classification are the A chain recorded in the repository and the canonical bovine B chain; P01317 itself was not available offline.",
        ],
    }
    L.write_json(prov, "provenance.json")
    print(json.dumps(checks, indent=1))


if __name__ == "__main__":
    main()
