"""G_plmsnosite_rescore, final step (system python): write results/G_plmsnosite_rescore/provenance.json.

POST HOC revision analysis (2026-09-30); not registered. Collects the SHA-256 of every input read
(git blobs of the pLMSNOSite clone at e9158af, the repository script imported for the features and
classifier configuration, this item's intermediate files), every script and every output, the
library versions of both Python environments and all seeds.

Run order: 10_detectability_baseline.py, 11_repro_diagnostics.py (system python);
20_run_plmsnosite.py (W/.venv_tf/Scripts/python.exe); 21_numpy_forward_pass.py, 30_compare.py,
31_detectability_tracking.py, 32_kr_stratified.py, 33_within_protein_mean.py, 35_proposed_source_data.py,
40_provenance.py (system python). Steps 21 and 32 and the long-format 35 were added after adversarial
verification round 1 (2026-09-30); step 33 (within-protein AUROC under the manuscript's estimand) and the
revised 35 after verification round 2. Files named verify_r1_* / verify_r2_* (scripts) and
results/G_plmsnosite_rescore/verify_r1/ and verify_r2/ were written by the verifier, not by this item's
implementer; they are hashed here but not used.
"""
from __future__ import annotations

import json
import os
import platform
import sys

import numpy as np
import pandas as pd
import scipy
import sklearn
import statsmodels

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402


def main():
    s1 = json.loads((C.OUT / "detectability_reproduction.json").read_text(encoding="utf-8"))
    s2 = json.loads((C.OUT / "plmsnosite_run.json").read_text(encoding="utf-8"))
    s3 = json.loads((C.OUT / "comparison_inputs_sha256.json").read_text(encoding="utf-8"))
    s4 = json.loads((C.OUT / "numpy_forward_pass_check.json").read_text(encoding="utf-8"))
    s5 = json.loads((C.OUT / "kr_stratified.json").read_text(encoding="utf-8"))
    s6 = json.loads((C.OUT / "withinprot_mean.json").read_text(encoding="utf-8"))
    inputs = {}
    for d in (s1.get("inputs_sha256", {}), s2.get("inputs_sha256_git_blobs", {}), s3, s4.get("inputs_sha256", {}),
              s5.get("inputs_sha256", {}), s6.get("inputs_sha256", {})):
        inputs.update(d)
    verifier_files = {}
    for vdir in ("verify_r1", "verify_r2"):
        if (C.OUT / vdir).exists():
            verifier_files.update({f"{vdir}/{p.name}": C.sha256_file(p) for p in sorted((C.OUT / vdir).glob("*"))
                                   if p.is_file()})
    outputs = {p.name: C.sha256_file(p) for p in sorted(C.OUT.iterdir())
               if p.is_file() and p.name != "provenance.json"}
    scripts = {p.name: C.sha256_file(p) for p in sorted(C.SCRIPTS.glob("*.py"))}
    prov = {
        "item": "G_plmsnosite_rescore",
        "label": ("POST HOC revision analysis (2026-09-30) in response to pre-submission review criticism; "
                  "not registered, not pre-specified"),
        "external_source": {"repository": "https://github.com/KCLabMTU/pLMSNOSite", "commit": C.EXPECTED_COMMIT,
                            "git_head_verified": C.git_head() == C.EXPECTED_COMMIT, "local_clone": str(C.EXT),
                            "read_via": "git show HEAD:<path> (LF blobs; the Windows working tree is CRLF)",
                            "registered_sha256_match": {k: inputs.get(f"pLMSNOSite@{C.EXPECTED_COMMIT[:7]}:{k}") == v
                                                        for k, v in C.REGISTERED_SHA256.items()}},
        "inputs_sha256": inputs,
        "scripts_sha256": scripts,
        "outputs_sha256": outputs,
        "verifier_outputs_sha256_not_used": verifier_files,
        "seeds": {"bootstrap_main": 20260930, "bootstrap_stability_check": 20260924,
                  "reproduction_bootstrap": "20260924 (one generator shared DIG25 -> VIS10 -> label-permuted)",
                  "label_permutation": "numpy default_rng(20260924).permutation(y_train)",
                  "hgb_random_state": "0 as in HGB_KWARGS (verified inert: identical predictions with 20260924)",
                  "tensorflow": "tf.random.set_seed(20260930), np.random.seed(20260930), TF_DETERMINISTIC_OPS=1; "
                                "inference only, repeated prediction identical",
                  "kr_stratified_bootstrap": "20260930 (all 278 test proteins; replicates identical to 30_compare.py)",
                  "within_protein_mean_bootstrap": ("20260930 (all 278 test proteins; replicates identical to "
                                                    "30_compare.py, asserted against the stored replicate file); "
                                                    "sensitivities: 20260930 over the 220 proteins with both labels "
                                                    "(Supplemental Note 12 scheme) and 20260924 (Monte Carlo check)"),
                  "sanity_replicate_seeds": "1 (30_compare.py and 33_within_protein_mean.py) and 7 (32_kr_stratified.py): "
                                            "one replicate each, used only to check weighted statistics against "
                                            "explicitly expanded rows"},
        "environment_system_python": {"python": platform.python_version(), "numpy": np.__version__,
                                      "pandas": pd.__version__, "scipy": scipy.__version__,
                                      "scikit-learn": sklearn.__version__, "statsmodels": statsmodels.__version__,
                                      "platform": platform.platform(), "OMP_NUM_THREADS": "4"},
        "environment_tf_venv": s2.get("environment"),
        "environment_numpy_forward_pass": s4.get("environment"),
        "numpy_forward_pass_max_abs_diff_vs_tensorflow": s4.get("max_abs_diff_vs_tensorflow"),
        "notes": [
            "TensorFlow cannot open non-ASCII paths on Windows; the verbatim evaluate_model.py run used a staging copy "
            "of the git blobs in an ASCII-only temporary directory (session scratchpad), deleted after the run.",
            "evaluate_model.py imports tqdm, matplotlib, seaborn and mpl_toolkits without using them; inert stub "
            "modules stood in for these packages, which are absent from the venv.",
            "The original detectability-baseline script and registration JSON are not in the repository; features and "
            "classifier were imported from repo/scripts/run_cross_protease_detectability_probe.py.",
            "Added after verification round 1: 21_numpy_forward_pass.py recomputes the released models' probabilities "
            "without TensorFlow and without any staging copy on disk (h5py reads the git blobs from memory); the "
            "TensorFlow run of 20_run_plmsnosite.py was not repeated.",
            "Added after verification round 2: 33_within_protein_mean.py computes the within-protein AUROC under the "
            "manuscript's estimand (unweighted mean of per-protein AUROCs over test proteins carrying both labels, "
            "computed by repo/scripts/common.py protein_rows imported by path with bytecode writing disabled); the "
            "pooled-pair within-protein AUROC of 30_compare.py is kept as a labelled sensitivity. No earlier step was "
            "re-run except 35_proposed_source_data.py and this script.",
        ],
    }
    (C.OUT / "provenance.json").write_text(json.dumps(prov, indent=1), encoding="utf-8")
    print(json.dumps({k: prov[k] for k in ("external_source", "seeds")}, indent=1))
    print(len(inputs), "inputs;", len(scripts), "scripts;", len(outputs), "outputs")


if __name__ == "__main__":
    main()
