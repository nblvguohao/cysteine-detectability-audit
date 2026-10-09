"""Run every test on one dataset and assemble the machine-readable audit record."""
from __future__ import annotations

import hashlib
import math
import platform

import numpy as np

from . import constants as C
from .checks import ORDER, abundance, cleavage, multicys, overlap, searchspace
from .stats import bh_qvalues
from .verdict import claim_retest, dataset_summary

RUNNERS = {"cleavage_geometry": cleavage.run, "multi_cysteine": multicys.run, "protein_abundance": abundance.run,
           "positive_detected_overlap": overlap.run, "search_space": searchspace.run}

INTERPRETATION = {
    C.PASS: "an artefact signature of meaningful size is excluded on this axis",
    C.WARNING: "an artefact signature is detectable on this axis but its size is not shown to be meaningful",
    C.FAIL: "an artefact signature of meaningful size is established on this axis; site-preference claims that "
            "touch it need an explicit control",
    C.UNDECIDABLE: "this axis cannot be assessed with the inputs given (see reason)",
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def clean(obj):
    """Recursively replace NaN/inf with None and numpy scalars with Python scalars (strict JSON)."""
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def config_record():
    return {"bootstrap_reps": C.BOOTSTRAP_REPS, "ci_level": C.CI_LEVEL, "min_positives": C.MIN_POSITIVES,
            "min_background": C.MIN_BACKGROUND, "mmd_log2_or": C.MMD_LOG2_OR, "mmd_auc": C.MMD_AUC,
            "cleavage_bands": {k: list(v) for k, v in C.CLEAVAGE_BANDS.items()},
            "cleavage_ci_level": C.CLEAVAGE_CI_LEVEL, "overlap_fail_lower": C.OVERLAP_FAIL_LOWER,
            "overlap_warning_lower": C.OVERLAP_WARNING_LOWER, "overlap_pass_upper": C.OVERLAP_PASS_UPPER,
            "abundance_min_coverage": C.ABUNDANCE_MIN_COVERAGE, "searchspace_fail_share": C.SEARCHSPACE_FAIL_SHARE,
            "attenuation_floor": C.ATTENUATION_FLOOR, "orphaned_limit": C.ORPHANED_LIMIT}


def run_audit(ds, cfg, inputs=None):
    tests = {}
    for name in ORDER:
        if name in cfg.get("skip", ()):
            continue
        tests[name] = RUNNERS[name](ds, cfg)
        tests[name]["interpretation"] = INTERPRETATION[tests[name]["status"]]
    names = list(tests)
    q = bh_qvalues([tests[n].get("p_value") for n in names])
    for n, qv in zip(names, q):
        tests[n]["q_value_bh"] = qv
    claim = claim_retest(ds, cfg)
    rec = {
        "tool": "cys-audit", "version": C.VERSION,
        "dataset": cfg.get("dataset", "unnamed"), "modification": cfg.get("modification"),
        "protease": cfg["protease"], "background": cfg["background"], "seed": cfg["seed"],
        "input": dict(inputs or {}, n_rows=len(ds), n_positive=int((ds.label == 1).sum()),
                      n_detected=None if ds.detected is None else int((ds.detected == 1).sum()),
                      n_proteins=len(set(ds.protein)), n_clusters=len(set(ds.cluster)), notes=ds.notes),
        "decision_constants": config_record(),
        "tests": tests,
        "multiple_testing": "Benjamini-Hochberg across the tests of this audit that report a bootstrap p "
                            "(q_value_bh); statuses are set by intervals and margins, never by p or q",
        "dataset_summary": dataset_summary(tests),
        "claim": claim,
        "overall_claim_status": claim["status"],
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
    }
    return clean(rec)
