"""Phase 0 baseline, part 2b: correction run for the second-tree rows of the first rerun.

WHY THIS FILE EXISTS (correction C1, appended; run-1 outputs are kept unchanged)
The first rerun (scripts/baseline_phase0_rerun_2026-09-21.py, outputs baseline/rerun_*_2026-09-21.*)
declared results/24_extra_sequences.json as an OUTPUT of the second-tree script
24_build_persulfidation_table.py, because the filename appears in that script as a literal, and so
deleted it from the mirror before the run. It is an INPUT: line 29 of that script reads it
("sequences retrieved one by one by accession for proteins outside the reference proteome"), and
no script on disk writes it. Deleting it removed ~87 sequences, which is why run 1 reported
24_mapped_osa_persulfidation.json with 581 proteins instead of 640 and why every downstream row of
36_ncys_decompose differed. Those seven run-1 rows (tree2_24_persulf_table x4,
tree2_36_ncys x3) are superseded by this run; they are not edited.

WHAT THIS RUN DOES (declared before it runs)
C1  Re-run second-tree scripts 24 and 36 in a fresh, guarded mirror of the second tree with
    24_extra_sequences.json left in place; declared outputs are only the files the scripts WRITE:
      24 -> 24_mapped_osa_persulfidation.json, 24_mapped_osa_allpeptides.json, 24_peptide_stats.json
      36 -> 36_ncys_decompose.csv, 36_ncys_decompose.json, 36_ncys_split.json
    Comparison rules and verdict words are imported unchanged from the run-1 harness.
C2  Characterise the run-1 DIFFERENT verdicts for second-tree scripts 10 and 12, which looked like
    last-digit float noise: over every numeric leaf of the stored and regenerated JSON, report the
    maximum absolute and relative difference, and whether any leaf that the manuscript prints
    (z at 2 decimals, p at 2 decimals) would change at printed precision. The run-1 mirror files
    are read from the run-1 scratch directory (PHASE0_RUN1_SCRATCH). This characterises; it does
    not change the run-1 verdicts.

OUTPUTS (new names)
  baseline/rerun_comparison_2026-09-21_correction.csv
  baseline/rerun_float_noise_tree2_10_12_2026-09-21.csv
  baseline/rerun_logs/correction_<label>.log
  baseline/rerun_correction_2026-09-21_audit.json   (with a 'corrections' key)
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRATCH = os.environ.get("PHASE0_SCRATCH") or sys.exit("set PHASE0_SCRATCH (fresh directory)")
RUN1 = os.environ.get("PHASE0_RUN1_SCRATCH") or sys.exit("set PHASE0_RUN1_SCRATCH (run-1 scratch)")
HARNESS = os.path.join(ROOT, "scripts", "baseline_phase0_rerun_2026-09-21.py")
spec = importlib.util.spec_from_file_location("h", HARNESS)
H = importlib.util.module_from_spec(spec)
spec.loader.exec_module(H)          # the harness has a __main__ guard; importing runs nothing

TREE2 = H.TREE2
M = os.path.join(SCRATCH, "mirror_tree2")
GUARD = os.path.join(SCRATCH, "guard")
GLOG = os.path.join(SCRATCH, "guard_log.tsv")
OUT = os.path.join(ROOT, "baseline")
PLAN = [
    ("tree2_24_persulf_table", "scripts/24_build_persulfidation_table.py",
     ["results/24_mapped_osa_persulfidation.json", "results/24_mapped_osa_allpeptides.json",
      "results/24_peptide_stats.json"], 900),
    ("tree2_36_ncys", "scripts/36_ncys_decompose.py",
     ["results/36_ncys_decompose.csv", "results/36_ncys_decompose.json", "results/36_ncys_split.json"], 2400),
]


def leaves(o, p=""):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaves(v, f"{p}/{k}")
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from leaves(v, f"{p}[{i}]")
    else:
        yield p, o


def main():
    outs = [os.path.join(OUT, f) for f in ("rerun_comparison_2026-09-21_correction.csv",
                                           "rerun_float_noise_tree2_10_12_2026-09-21.csv",
                                           "rerun_correction_2026-09-21_audit.json")]
    for p in outs:
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    if os.path.exists(M):
        sys.exit("REFUSE: fresh PHASE0_SCRATCH required")
    t0 = time.time()
    os.makedirs(M)
    for d in ("scripts", "results"):
        shutil.copytree(os.path.join(TREE2, d), os.path.join(M, d), symlinks=True,
                        ignore=shutil.ignore_patterns("__pycache__"))
    os.symlink(os.path.join(TREE2, "data"), os.path.join(M, "data"))
    os.makedirs(GUARD)
    with open(os.path.join(GUARD, "sitecustomize.py"), "w") as fh:
        fh.write(H.GUARD_SRC)
    env = dict(os.environ)
    env.update({"PYTHONHASHSEED": "17", "PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "4",
                "PHASE0_PROTECT": os.pathsep.join([os.path.realpath(ROOT), os.path.realpath(TREE2)]),
                "PHASE0_GUARD_LOG": GLOG,
                "PYTHONPATH": os.pathsep.join([GUARD, os.path.join(M, "scripts")])})
    comps, runs = [], []
    for label, script, outputs, timeout in PLAN:
        for o in outputs:
            if os.path.lexists(os.path.join(M, o)):
                os.remove(os.path.join(M, o))
        st = time.time()
        proc = subprocess.run([H.PY_STD, os.path.join(M, script)], cwd=M, env=env,
                              capture_output=True, text=True, timeout=timeout)
        wall = round(time.time() - st, 1)
        with open(os.path.join(OUT, "rerun_logs", f"correction_{label}.log"), "w", encoding="utf-8") as fh:
            fh.write(f"# {label} rc {proc.returncode} wall_s {wall}\n{proc.stdout[-20000:]}\n# stderr\n{proc.stderr[-20000:]}\n")
        runs.append({"label": label, "returncode": proc.returncode, "wall_s": wall})
        for o in outputs:
            if proc.returncode != 0:
                v, d = "RUN_FAILED", f"returncode {proc.returncode}"
            else:
                v, d = H.compare(os.path.join(TREE2, o), os.path.join(M, o), M, TREE2)
            comps.append({"label": label, "tree": "tree2", "output": o, "verdict": v, "detail": d[:600]})
    with open(outs[0], "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(comps[0]))
        w.writeheader()
        w.writerows(comps)
    # C2: float-noise characterisation of the run-1 second-tree 10/12 JSONs
    noise = []
    for f in ("10_large_scale_motif.json", "12_matched_sno.json", "12_matched_so.json", "12_matched_ox.json"):
        with open(os.path.join(TREE2, "results", f), encoding="utf-8") as fh:
            a = dict(leaves(json.load(fh)))
        with open(os.path.join(RUN1, "mirror_tree2", "results", f), encoding="utf-8") as fh:
            b = dict(leaves(json.load(fh)))
        mabs = mrel = 0.0
        n_num = n_diff = n_print_change = 0
        keys_mismatch = sorted(set(a) ^ set(b))
        for k in set(a) & set(b):
            x, y = a[k], b[k]
            if isinstance(x, bool) or isinstance(y, bool):
                n_diff += x != y
                continue
            if isinstance(x, (int, float)) and isinstance(y, (int, float)):
                n_num += 1
                if x != y:
                    n_diff += 1
                    d = abs(x - y)
                    mabs = max(mabs, d)
                    mrel = max(mrel, d / max(abs(x), 1e-300))
                    if round(x, 2) != round(y, 2):
                        n_print_change += 1
            elif x != y:
                n_diff += 1
                n_print_change += 1
        noise.append({"file": f, "numeric_leaves": n_num, "differing_leaves": n_diff,
                      "max_abs_diff": f"{mabs:.3g}", "max_rel_diff": f"{mrel:.3g}",
                      "leaves_changing_at_2_decimals": n_print_change,
                      "keys_present_in_only_one": len(keys_mismatch)})
    with open(outs[1], "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(noise[0]))
        w.writeheader()
        w.writerows(noise)
    glog = open(GLOG).read().splitlines() if os.path.exists(GLOG) else []
    audit = {
        "script": "scripts/baseline_phase0_rerun_correction_2026-09-21.py",
        "script_sha256": H.sha(os.path.abspath(__file__)),
        "imports_harness": "scripts/baseline_phase0_rerun_2026-09-21.py",
        "imported_harness_sha256": H.sha(HARNESS),
        "python": sys.version,
        "runs": runs,
        "comparison": comps,
        "float_noise_10_12": noise,
        "guard_refusals": len(glog),
        "corrections": [{
            "id": "C1",
            "what": "run 1 declared results/24_extra_sequences.json as an output of second-tree script 24 and deleted it before running; it is an input (read at line 29, written by no script on disk)",
            "effect_on_run1": "rows tree2_24_persulf_table (4) and tree2_36_ncys (3) in baseline/rerun_comparison_2026-09-21.csv are superseded by this file; they are left as written",
            "what_this_correction_does_NOT_change": ["any verdict for this tree's 18 producer scripts",
                                                    "the run-1 verdicts for second-tree scripts 10 and 12",
                                                    "the comparison rules or verdict vocabulary (imported unchanged)"],
        }],
        "outputs": {os.path.relpath(p, ROOT): H.sha(p) for p in outs[:2]},
        "total_wall_s": round(time.time() - t0, 1),
    }
    with open(outs[2], "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1, allow_nan=False)
    for c in comps:
        print(c["verdict"].ljust(18), c["output"], c["detail"][:200])
    for n in noise:
        print(n)


if __name__ == "__main__":
    main()
