"""Step 5 - provenance.json for item I_fig5a_estimand (POST HOC revision analysis, 2026-09-30).

Merges the per-step input records, hashes every output and every script of this item, records
library versions and seeds, and checks that the original repository scripts behind the stored
outputs are unchanged apart from line endings.

Run:  python -B s05_provenance.py   (after s01-s04)
"""
from __future__ import annotations

import datetime
import json
import platform
import sys

import matplotlib
import numpy
import pandas
import scipy
import sklearn

import common as C


def main():
    steps = {}
    for p in sorted(C.OUT.glob("_inputs_*.json")):
        rec = json.loads(p.read_text(encoding="utf-8"))
        steps[rec["step"]] = rec
    inputs = {}
    for rec in steps.values():
        for path, meta in rec["inputs"].items():
            inputs[path] = meta
    census_audit = json.loads(C.CENSUS_AUDIT.read_text(encoding="utf-8"))
    sites_audit = json.loads(C.CENSUS_SITES_AUDIT.read_text(encoding="utf-8"))
    recorded = {
        C.CENSUS_SCRIPT: census_audit["script_sha256"],
        C.PROBE_SCRIPT: census_audit["imported_feature_definitions"]["sha256"],
        C.INGEST_SCRIPT: census_audit["ingest_script_sha256"],
    }
    repo_scripts = {}
    for path in (C.CENSUS_SCRIPT, C.PROBE_SCRIPT, C.INGEST_SCRIPT, C.DIAGNOSTICS_SCRIPT):
        d = {"sha256_raw": C.sha256(path), "sha256_lf_normalised": C.sha256_lf(path)}
        if path in recorded:
            d["sha256_recorded_in_census_audit"] = recorded[path]
            d["unchanged_except_line_endings"] = d["sha256_lf_normalised"] == recorded[path]
        repo_scripts[str(path)] = d
    diag = json.loads(C.DIAGNOSTICS_AUDIT.read_text(encoding="utf-8"))
    repo_scripts[str(C.DIAGNOSTICS_SCRIPT)]["sha256_recorded_in_diagnostics_audit"] = diag["script_sha256"]
    repo_scripts[str(C.DIAGNOSTICS_SCRIPT)]["unchanged_except_line_endings"] = (
        repo_scripts[str(C.DIAGNOSTICS_SCRIPT)]["sha256_lf_normalised"] == diag["script_sha256"])
    for extra in (C.CENSUS_AUDIT, C.CENSUS_SITES_AUDIT, C.DIAGNOSTICS_AUDIT, C.CENSUS_VERIFICATION):
        inputs[str(extra)] = {"exists": extra.exists(), "sha256": C.sha256(extra), "bytes": extra.stat().st_size}
    outputs = {p.name: {"sha256": C.sha256(p), "bytes": p.stat().st_size}
               for p in sorted(C.OUT.iterdir()) if p.is_file() and not p.name.startswith("_inputs_")
               and p.name != "provenance.json"}
    scripts_dir = C.W / "scripts" / C.ITEM
    scripts = {p.name: C.sha256(p) for p in sorted(scripts_dir.glob("*.py"))}
    prov = {
        "item": C.ITEM,
        "label": ("POST HOC revision analysis in response to pre-submission review criticism of Figure 5a; "
                  "not registered and not pre-specified"),
        "generated": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "analysis_date": "2026-09-30",
        "what_was_computed": [
            "s01: tables of stored absolute AUCs (global and within-protein) for VIS10, DIG25, CPL15 under NEG_A/NEG_B; "
            "share recomputation and rounding bounds; Source Data and manuscript-number checks; derived point quantities "
            "(within-protein share, complement ratio, unit composition); task-3 table; proposed Source Data rows",
            "s02: feature-definition identities (CPL15 encodes peptide length) using the original feature code on "
            "synthetic sequences and the first 3,000 mouse reference-proteome entries (not a Fig. 5a cohort)",
            "s03: availability of the inputs a chemistry-aware denominator would need; no model fitted",
            "s04: draft figure of the stored absolute AUCs",
        ],
        "models_fitted": "none",
        "network_access": "none",
        "seeds": {"synthetic_sequences_s02": C.SEED,
                  "original_census_bootstrap_seed_not_rerun": census_audit["seed"]},
        "python": sys.version, "platform": platform.platform(),
        "libraries": {"numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__,
                      "scikit-learn": sklearn.__version__, "matplotlib": matplotlib.__version__},
        "inputs_read": inputs,
        "original_repository_scripts": repo_scripts,
        "original_census_inputs_recorded_but_absent_locally": {
            "results/ptm_census_sites.csv": census_audit["sites_table_sha256"],
            "external/proteomes/hsa.fasta.gz": sites_audit["inputs"]["hsa_proteome"]["sha256"],
            "external/proteomes/ath.fasta.gz": sites_audit["inputs"]["ath_proteome"]["sha256"],
        },
        "item_scripts_sha256": scripts,
        "outputs_sha256": outputs,
        "steps": {k: {kk: vv for kk, vv in v.items() if kk != "inputs"} for k, v in steps.items()},
        "reproduce": "cd scripts/I_fig5a_estimand && python -B run_all.py",
    }
    (C.OUT / "provenance.json").write_text(json.dumps(prov, indent=1, ensure_ascii=False), encoding="utf-8")
    print("provenance written;", len(inputs), "inputs,", len(outputs), "outputs")
    print(json.dumps({k: v.get("unchanged_except_line_endings") for k, v in repo_scripts.items()}, indent=1,
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
