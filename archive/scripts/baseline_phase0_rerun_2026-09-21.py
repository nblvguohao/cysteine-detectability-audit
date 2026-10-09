"""Phase 0 baseline, part 2: re-execute producer scripts in an isolated mirror and compare outputs.

(Declared before the run. Nothing below is tuned after looking at a comparison.)

QUESTION
Can the stored products behind the current manuscript (R24) be regenerated, byte for byte or value
for value, by the code on disk today? This is a minimal-scale reproduction: only producer scripts
that are cheap, deterministic in design, and need no network or remote host are re-run.

ISOLATION (nothing in either real tree can be written)
* A mirror is built under a scratch folder (PHASE0_SCRATCH). For this tree, the directories a
  producer may write into are COPIED (scripts, results, protocols, reports, shared,
  source_data*, supplementary*); large read-only inputs are SYMLINKED (inputs, external, features,
  plm, models_v2, interaction, structures, tabular). For the second tree
  (/path/to/project) scripts/ and results/ are copied and data/ is
  symlinked.
* Every child interpreter starts with a sitecustomize.py that installs a PEP 578 audit hook. Any
  open-for-write, rename, replace, remove, rmdir, new mkdir, truncate, chmod, utime, link, symlink,
  copyfile/copytree destination or rmtree whose real path lies inside either real tree raises
  PermissionError and is logged. Reads are allowed. (Tectonic is a subprocess, not Python; it is
  given mirror paths only.)
* Before a script runs, every output it is declared to write is deleted from the mirror, so a
  script that silently fails to write cannot be mistaken for a reproduction.

COMPARISON, per declared output
  IDENTICAL          sha256 equal to the stored file in the real tree
  EQUIVALENT_VALUES  sha256 differs but: CSV -> header and every cell equal as strings;
                     JSON -> equal after (a) replacing the mirror root string with the real root
                     string and (b) dropping VOLATILE keys (below); PNG -> decoded pixel arrays equal
  DIFFERENT          anything else (first differences are listed; for CSVs the maximum absolute
                     numeric difference among differing numeric cells is also reported, so float
                     noise can be told apart from a changed value; it does not change the verdict)
  NOT_REGENERATED    the output does not exist after the run
  NO_STORED_COPY     the real tree has no stored copy to compare against
  RUN_FAILED         non-zero exit or timeout (outputs are then not compared)
VOLATILE keys (JSON only): any key whose lower-cased name contains one of
  time, date, elapsed, started, finished, runtime, wall, timestamp, generated, utc, host, created
  -- plus 'pdf_sha256' and 'log_sha256' for the LaTeX build (Tectonic embeds a creation date).
  No CSV column is ever dropped. Vector/PDF/SVG files are compared by sha256 only and reported as
  informational, because matplotlib and Tectonic embed creation dates.

RUN SETTINGS
  PYTHONHASHSEED=17 (the original runs did not fix it; a hash-seed-sensitive script shows up as
  DIFFERENT, which is itself a finding), PYTHONDONTWRITEBYTECODE=1, OMP/MKL threads 4.
  Interpreter per script = the class recorded in its audit: stdlib -> /usr/bin/python3 (3.9.6),
  numpy/sklearn/lightgbm -> the project venv (3.9.6), matplotlib -> the lnrna venv (3.9.6).

NOT RE-RUN, AND WHY (reported in the output, not hidden)
  - refit of the ablation / detectability-only arms: ran on amax for ~29 min on 80 threads; the
    stored out-of-fold scores are re-analysed instead (analyse_public_refit, posthoc cost).
  - Comet re-search of PXD015307: Comet is not installed here and the raw files live on amax; the
    stored search tables are re-analysed instead.
  - literature-survey ingest: reads the external artifact store, not a file in this tree.
  - claim re-test rounds phase2/2B/2C/2D/2E: canonical interpreter differs (3.11.16) and some
    rounds read paywalled intake; their stored verdict tables are re-tallied instead.
  - AlphaFold structural features, PLM embeddings, the v2 model training: hours of compute.

OUTPUTS (new names; nothing in the real tree is overwritten)
  baseline/rerun_comparison_2026-09-21.csv
  baseline/rerun_runs_2026-09-21.csv
  baseline/rerun_logs/<script>.log
  baseline/rerun_guard_log_2026-09-21.tsv
  baseline/baseline_figures/  (regenerated PNG/PDF copied out of the mirror)
  baseline/rerun_2026-09-21_audit.json   (this script's sha256, settings, counts)
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREE2 = "/path/to/project"
SCRATCH = os.environ.get("PHASE0_SCRATCH") or sys.exit("set PHASE0_SCRATCH to a scratch directory")
M_THIS = os.path.join(SCRATCH, "mirror_this")
M_TREE2 = os.path.join(SCRATCH, "mirror_tree2")
GUARD_DIR = os.path.join(SCRATCH, "guard")
GUARD_LOG = os.path.join(SCRATCH, "guard_log.tsv")
OUT = os.path.join(ROOT, "baseline")
LOGS = os.path.join(OUT, "rerun_logs")
FIGS = os.path.join(OUT, "baseline_figures")

PY_STD = "/usr/bin/python3"
PY_VENV = "/path/to/venv/bin/python"
PY_MPL = "/path/to/venv2/bin/python"

COPY_THIS = ["scripts", "results", "protocols", "reports", "shared", "source_data",
             "source_data_r22", "source_data_r23", "source_data_r24", "supplementary",
             "supplementary_r21"]
LINK_THIS = ["inputs", "external", "features", "plm", "models_v2", "interaction", "structures",
             "tabular", "deliverables", "handoff", "logs"]
VOLATILE = ("time", "date", "elapsed", "started", "finished", "runtime", "wall", "timestamp",
            "generated", "utc", "host", "created")
VOLATILE_EXACT = {"pdf_sha256", "log_sha256"}


def figs_in(tree_dir, prefix=None, exclude_prefix=None):
    d = os.path.join(ROOT, tree_dir)
    return sorted(os.path.join(tree_dir, f) for f in os.listdir(d)
                  if (prefix is None or f.startswith(prefix))
                  and (exclude_prefix is None or not f.startswith(exclude_prefix))
                  and not f.startswith("."))


def plan():
    """(label, tree, interpreter, script relpath, [declared outputs relpaths], timeout_s)"""
    r = "results/"
    P = [
        ("claim_verdict_tally", "this", PY_STD, "scripts/recount_claim_verdicts_2026-09-17.py",
         [r + "claim_verdict_tally_2026-09-17.csv", r + "claim_verdict_tally_2026-09-17_audit.json"], 300),
        ("diagnostics_regression", "this", PY_VENV, "scripts/ptm_detectability_diagnostics.py",
         [r + "label_semantics_tiers.csv", r + "ptm_detectability_diagnostics_audit.json"], 600),
        ("artefact3_abundance", "this", PY_VENV, "scripts/recompute_artefact3_osa_abundance_2026-09-19.py",
         [r + "artefact3_osa_abundance_recomputed_2026-09-19.csv", r + "artefact3_osa_abundance_recomputed_2026-09-19_audit.json"], 900),
        ("artefact3_effect_sizes", "this", PY_VENV, "scripts/artefact3_effect_sizes_and_matching_2026-09-19.py",
         [r + "artefact3_effect_sizes_2026-09-19.csv", r + "artefact3_effect_sizes_2026-09-19_audit.json"], 1800),
        ("artefact3_positive_set", "this", PY_VENV, "scripts/artefact3_positive_set_sensitivity_2026-09-19.py",
         [r + "artefact3_positive_set_sensitivity_2026-09-19.csv", r + "artefact3_positive_set_sensitivity_2026-09-19_audit.json"], 1800),
        ("artefact3_original_arms", "this", PY_VENV, "scripts/artefact3_reconstructed_original_arms_2026-09-19.py",
         [r + "artefact3_reconstructed_original_arms_2026-09-19.csv", r + "artefact3_reconstructed_original_arms_2026-09-19_audit.json"], 1800),
        ("artefact4_pxd048216", "this", PY_VENV, "scripts/analyse_public_artefact4_2026-09-21.py",
         [r + "public_artefact4_coincidence_2026-09-21.csv", r + "public_artefact4_calibration_2026-09-21.csv",
          r + "public_artefact4_multicys_2026-09-21.csv", r + "public_artefact4_2026-09-21_audit.json"], 1800),
        ("artefact5_score_ties", "this", PY_STD, "scripts/posthoc_pxd015307_score_ties.py",
         [r + "pxd015307_posthoc_score_ties.csv", r + "pxd015307_posthoc_score_ties_audit.json"], 900),
        ("artefact5_mass_accuracy", "this", PY_STD, "scripts/run_pxd015307_research_mass_accuracy.py",
         [r + "pxd015307_research_summary.csv", r + "pxd015307_research_psms.csv", r + "pxd015307_research_audit.json"], 900),
        ("protease_reach", "this", PY_STD, "scripts/recompute_protease_reach_public_2026-09-17.py",
         [r + "protease_reach_public_2026-09-17.csv", r + "protease_reach_public_audit.json"], 900),
        ("brg3_ring", "this", PY_STD, "scripts/check_brg3_ring_all_cysteines_2026-09-17.py",
         [r + "brg3_ring_all_cysteines_2026-09-17.csv", r + "brg3_ring_all_cysteines_audit.json"], 300),
        ("positional_kr_public", "this", PY_VENV, "scripts/positional_kr_profile_public_2026-09-17.py",
         [r + "positional_kr_profile_public_2026-09-17.csv", r + "positional_kr_profile_public_summary_2026-09-17.csv",
          r + "positional_kr_profile_public_audit.json"], 3600),
        ("self_audit_public", "this", PY_VENV, "scripts/run_self_audit_public_cohorts_2026-09-20.py",
         [r + "self_audit_public_cohorts_2026-09-20.csv", r + "self_audit_public_five_proteases_2026-09-20.csv",
          r + "self_audit_public_cohorts_2026-09-20_audit.json"], 3600),
        ("public_refit_analysis", "this", PY_VENV, "scripts/analyse_public_refit_2026-09-21.py",
         [r + "public_refit_ablation_2026-09-21.csv", r + "public_refit_detectability_only_2026-09-21.csv",
          r + "public_refit_ablation_2026-09-21_audit.json"], 3600),
        ("public_refit_posthoc_cost", "this", PY_VENV, "scripts/posthoc_public_refit_performance_cost_2026-09-21.py",
         [r + "public_refit_performance_cost_posthoc_2026-09-21.csv", r + "public_refit_performance_cost_posthoc_2026-09-21_audit.json"], 1800),
        ("figs_2to5", "this", PY_MPL, "scripts/plot_figs_r23_2026-09-21.py",
         figs_in("figures_r23", exclude_prefix="Fig1_") + figs_in("source_data_r23")
         + [r + "figs_r23_2026-09-21_audit.json"], 900),
        ("fig1_overview", "this", PY_MPL, "scripts/plot_fig1_overview_r24_2026-09-21.py",
         figs_in("figures_r24", "Fig1_overview") + [r + "fig1_overview_r24_2026-09-21_audit.json"], 600),
        ("latex_r24", "this", PY_STD, "scripts/md_to_latex_r24_2026-09-21.py",
         ["latex_r24/manuscript_R24.tex", r + "latex_r24_build_audit_2026-09-21.json"], 900),
        ("tree2_10_motif", "tree2", PY_STD, "scripts/10_large_scale_motif.py",
         ["results/10_large_scale_motif.json"], 900),
        ("tree2_12_matched", "tree2", PY_STD, "scripts/12_detectability_matched.py",
         ["results/12_matched_sno.json", "results/12_matched_so.json", "results/12_matched_ox.json"], 900),
        ("tree2_24_persulf_table", "tree2", PY_STD, "scripts/24_build_persulfidation_table.py",
         ["results/24_mapped_osa_persulfidation.json", "results/24_mapped_osa_allpeptides.json",
          "results/24_peptide_stats.json", "results/24_extra_sequences.json"], 900),
        ("tree2_36_ncys", "tree2", PY_STD, "scripts/36_ncys_decompose.py",
         ["results/36_ncys_decompose.csv", "results/36_ncys_decompose.json", "results/36_ncys_split.json"], 2400),
    ]
    return P


GUARD_SRC = r'''
import os, sys
_PROTECT = [os.path.realpath(p) for p in os.environ.get("PHASE0_PROTECT", "").split(os.pathsep) if p]
_LOG = os.environ.get("PHASE0_GUARD_LOG")
_WF = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND | os.O_TRUNC
def _real(p):
    try:
        return os.path.realpath(os.fsdecode(p))
    except Exception:
        return None
def _under(p):
    if p is None or isinstance(p, int):
        return False
    rp = _real(p)
    return rp is not None and any(rp == r or rp.startswith(r + os.sep) for r in _PROTECT)
def _deny(ev, p):
    if _LOG:
        with open(_LOG, "a") as fh:
            fh.write("%s\t%s\t%s\n" % (ev, os.fsdecode(p), sys.argv[0] if sys.argv else ""))
    raise PermissionError("PHASE0 guard refused %s on protected tree: %s" % (ev, os.fsdecode(p)))
def _hook(event, args):
    if event == "open":
        path, mode, flags = (list(args) + [None, None, None])[:3]
        if isinstance(mode, str):
            w = any(c in mode for c in "wax+")
        else:
            w = isinstance(flags, int) and bool(flags & _WF)
        if w and _under(path):
            _deny(event, path)
    elif event in ("os.rename", "os.replace", "shutil.move"):
        for p in args[:2]:
            if _under(p):
                _deny(event, p)
    elif event in ("os.link", "os.symlink", "shutil.copyfile", "shutil.copytree"):
        if len(args) > 1 and _under(args[1]):
            _deny(event, args[1])
    elif event == "os.mkdir":
        p = args[0] if args else None
        if _under(p) and not os.path.exists(os.fsdecode(p)):
            _deny(event, p)
    elif event in ("os.remove", "os.rmdir", "os.truncate", "os.chmod", "os.utime", "os.chown",
                   "shutil.rmtree", "os.unlink"):
        p = args[0] if args else None
        if _under(p):
            _deny(event, p)
if _PROTECT:
    sys.addaudithook(_hook)
'''


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def build_mirror():
    if os.path.exists(M_THIS) or os.path.exists(M_TREE2):
        sys.exit(f"REFUSE: mirror already exists under {SCRATCH}; use a fresh PHASE0_SCRATCH")
    os.makedirs(M_THIS)
    for d in COPY_THIS:
        src = os.path.join(ROOT, d)
        if os.path.isdir(src):
            shutil.copytree(src, os.path.join(M_THIS, d), symlinks=True,
                            ignore=shutil.ignore_patterns("__pycache__", "*.tiff", "*.tif"))
    for d in LINK_THIS:
        if os.path.exists(os.path.join(ROOT, d)):
            os.symlink(os.path.join(ROOT, d), os.path.join(M_THIS, d))
    # output folders of the figure and LaTeX producers: real, empty, except the figure PDFs the
    # LaTeX converter reads as inputs (copied; Fig1 is regenerated before the converter runs)
    for d in ("figures_r23", "figures_r24"):
        os.makedirs(os.path.join(M_THIS, d))
    for f in os.listdir(os.path.join(ROOT, "figures_r24")):
        if f.endswith(".pdf") and not f.startswith("Fig1_"):
            shutil.copy2(os.path.join(ROOT, "figures_r24", f), os.path.join(M_THIS, "figures_r24", f))
    os.makedirs(os.path.join(M_THIS, "latex_r24"))
    os.makedirs(M_TREE2)
    for d in ("scripts", "results"):
        shutil.copytree(os.path.join(TREE2, d), os.path.join(M_TREE2, d), symlinks=True,
                        ignore=shutil.ignore_patterns("__pycache__"))
    os.symlink(os.path.join(TREE2, "data"), os.path.join(M_TREE2, "data"))
    os.makedirs(GUARD_DIR)
    with open(os.path.join(GUARD_DIR, "sitecustomize.py"), "w") as fh:
        fh.write(GUARD_SRC)


def strip_volatile(o):
    if isinstance(o, dict):
        return {k: strip_volatile(v) for k, v in o.items()
                if not (any(t in k.lower() for t in VOLATILE) or k in VOLATILE_EXACT)}
    if isinstance(o, list):
        return [strip_volatile(v) for v in o]
    return o


def json_diff(a, b, path="", out=None, limit=20):
    if out is None:
        out = []
    if len(out) >= limit:
        return out
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{path}/{k}: only in {'stored' if k in a else 'rerun'}")
            else:
                json_diff(a[k], b[k], f"{path}/{k}", out, limit)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: list length {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            json_diff(x, y, f"{path}[{i}]", out, limit)
    elif a != b:
        out.append(f"{path}: {str(a)[:60]!r} vs {str(b)[:60]!r}")
    return out


def compare(stored, rerun, mirror_root, real_root):
    if not os.path.exists(stored):
        return "NO_STORED_COPY", ""
    if not os.path.exists(rerun):
        return "NOT_REGENERATED", ""
    if sha(stored) == sha(rerun):
        return "IDENTICAL", ""
    ext = os.path.splitext(stored)[1].lower()
    if ext == ".csv":
        with open(stored, newline="", encoding="utf-8") as f1, open(rerun, newline="", encoding="utf-8") as f2:
            a, b = list(csv.reader(f1)), list(csv.reader(f2))
        b = [[c.replace(mirror_root, real_root) for c in row] for row in b]
        if a == b:
            return "EQUIVALENT_VALUES", "differs only in the mirror path string"
        diffs, maxabs = [], 0.0
        if len(a) != len(b):
            diffs.append(f"rows {len(a)} vs {len(b)}")
        for i, (ra, rb) in enumerate(zip(a, b)):
            for j, (ca, cb) in enumerate(zip(ra, rb)):
                if ca != cb:
                    try:
                        maxabs = max(maxabs, abs(float(ca) - float(cb)))
                    except ValueError:
                        pass
                    if len(diffs) < 8:
                        head = a[0][j] if a and j < len(a[0]) else j
                        diffs.append(f"row {i} col {head}: {ca[:40]!r} vs {cb[:40]!r}")
        n = sum(1 for ra, rb in zip(a, b) for ca, cb in zip(ra, rb) if ca != cb)
        return "DIFFERENT", f"{n} differing cells; max_abs_numeric_diff={maxabs:.3g}; " + "; ".join(diffs)
    if ext == ".json":
        with open(stored, encoding="utf-8") as fh:
            a = json.load(fh)
        with open(rerun, encoding="utf-8") as fh:
            b = json.loads(fh.read().replace(mirror_root, real_root))
        sa, sb = strip_volatile(a), strip_volatile(b)
        if sa == sb:
            return "EQUIVALENT_VALUES", "equal after dropping volatile keys"
        return "DIFFERENT", "; ".join(json_diff(sa, sb))
    if ext == ".png":
        try:
            import matplotlib.image as mpimg
            import numpy as np
            x, y = mpimg.imread(stored), mpimg.imread(rerun)
            if x.shape == y.shape and np.array_equal(x, y):
                return "EQUIVALENT_VALUES", "pixel arrays equal"
            if x.shape != y.shape:
                return "DIFFERENT", f"shape {x.shape} vs {y.shape}"
            return "DIFFERENT", f"pixels differ: max abs {float(np.abs(x - y).max()):.3g}, frac {float((x != y).mean()):.3g}"
        except Exception as e:  # noqa: BLE001
            return "DIFFERENT", f"png compare failed: {e}"
    if ext == ".tex":
        with open(stored, encoding="utf-8") as fh:
            a = fh.read()
        with open(rerun, encoding="utf-8") as fh:
            b = fh.read().replace(mirror_root, real_root)
        return ("EQUIVALENT_VALUES", "equal after path normalisation") if a == b else ("DIFFERENT", "tex text differs")
    return "DIFFERENT", "binary differs (informational for pdf/svg: creation dates are embedded)"


def main():
    for p in (os.path.join(OUT, "rerun_comparison_2026-09-21.csv"), os.path.join(OUT, "rerun_2026-09-21_audit.json")):
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    os.makedirs(LOGS, exist_ok=True)
    os.makedirs(FIGS, exist_ok=True)
    t0 = time.time()
    build_mirror()
    env_base = dict(os.environ)
    env_base.update({
        "PYTHONHASHSEED": "17", "PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "4",
        "MKL_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "4", "MPLBACKEND": "Agg",
        "PHASE0_PROTECT": os.pathsep.join([os.path.realpath(ROOT), os.path.realpath(TREE2)]),
        "PHASE0_GUARD_LOG": GUARD_LOG,
    })
    runs, comps = [], []
    for label, tree, py, script, outputs, timeout in plan():
        mroot = M_THIS if tree == "this" else M_TREE2
        rroot = ROOT if tree == "this" else TREE2
        for o in outputs:
            mp = os.path.join(mroot, o)
            if os.path.lexists(mp):
                os.remove(mp)
        env = dict(env_base)
        env["PYTHONPATH"] = os.pathsep.join([GUARD_DIR, os.path.join(mroot, "scripts")])
        log_path = os.path.join(LOGS, f"{label}.log")
        started = time.time()
        try:
            proc = subprocess.run([py, os.path.join(mroot, script)], cwd=mroot, env=env,
                                  capture_output=True, text=True, timeout=timeout)
            rc, so, se = proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired as e:
            rc = "TIMEOUT"
            so = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
            se = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or "")
        wall = round(time.time() - started, 1)
        with open(log_path, "w", encoding="utf-8") as fh:
            fh.write(f"# {label}\n# interpreter {py}\n# script {script}\n# rc {rc}\n# wall_s {wall}\n")
            fh.write("# ---- stdout (last 200 lines)\n" + "\n".join(so.splitlines()[-200:]) + "\n")
            fh.write("# ---- stderr (last 200 lines)\n" + "\n".join(se.splitlines()[-200:]) + "\n")
        runs.append({"label": label, "tree": tree, "interpreter": py, "script": script,
                     "returncode": rc, "wall_s": wall, "log": os.path.relpath(log_path, ROOT)})
        print(f"[{label}] rc={rc} wall={wall}s", flush=True)
        for o in outputs:
            if rc != 0:
                v, d = "RUN_FAILED", f"returncode {rc}"
            else:
                v, d = compare(os.path.join(rroot, o), os.path.join(mroot, o), mroot, rroot)
            comps.append({"label": label, "tree": tree, "output": o, "verdict": v, "detail": d[:600]})
            ext = os.path.splitext(o)[1].lower()
            if ext in (".png", ".pdf") and os.path.exists(os.path.join(mroot, o)):
                shutil.copy2(os.path.join(mroot, o), os.path.join(FIGS, os.path.basename(o)))
    with open(os.path.join(OUT, "rerun_runs_2026-09-21.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(runs[0]))
        w.writeheader()
        w.writerows(runs)
    comp_csv = os.path.join(OUT, "rerun_comparison_2026-09-21.csv")
    with open(comp_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(comps[0]))
        w.writeheader()
        w.writerows(comps)
    guard_rows = open(GUARD_LOG).read().splitlines() if os.path.exists(GUARD_LOG) else []
    with open(os.path.join(OUT, "rerun_guard_log_2026-09-21.tsv"), "w", encoding="utf-8") as fh:
        fh.write("event\tpath\tscript\n" + "\n".join(guard_rows) + ("\n" if guard_rows else ""))
    from collections import Counter
    audit = {
        "script": "scripts/baseline_phase0_rerun_2026-09-21.py",
        "script_sha256": sha(os.path.abspath(__file__)),
        "python": sys.version,
        "scratch_mirror": SCRATCH,
        "settings": {"PYTHONHASHSEED": "17", "threads": 4, "volatile_json_keys": list(VOLATILE) + sorted(VOLATILE_EXACT)},
        "interpreters": {"stdlib": PY_STD, "numpy_sklearn_lightgbm": PY_VENV, "matplotlib": PY_MPL},
        "n_scripts": len(runs),
        "run_returncodes": {r["label"]: r["returncode"] for r in runs},
        "verdict_counts": dict(Counter(c["verdict"] for c in comps)),
        "guard_refusals": len(guard_rows),
        "comparison_csv_sha256": sha(comp_csv),
        "total_wall_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(OUT, "rerun_2026-09-21_audit.json"), "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1, allow_nan=False)
    print(json.dumps(audit["verdict_counts"], indent=1))


if __name__ == "__main__":
    main()
