"""Run every step of item C_pxd063463_specific in order (POST HOC revision analysis, 2026-09-30).

  python run_all.py
Each step is a separate process; logs go to results/C_pxd063463_specific/logs/.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c_common as C  # noqa: E402

STEPS = [("c01_reproduce_fig2.py", "filtered"), ("c01b_seed_sensitivity.py", "filtered"),
         ("c02_specific_sites.py", "filtered"), ("c03_cleavage_specific.py", "filtered"),
         ("c04_coincidence_specific.py", "filtered"), ("c05_depth_confound.py", "filtered"),
         ("c08_hydn_vs_hydp_depth.py", "filtered"),
         ("c11_stratum_matched.py", "filtered"),       # c11-c13 added in revision after verification (round 2)
         ("c12_coverage_proxy.py", "filtered"),        # computes both table bases itself
         ("c02_specific_sites.py", "unfiltered"), ("c04_coincidence_specific.py", "unfiltered"),
         ("c05_depth_confound.py", "unfiltered"), ("c08_hydn_vs_hydp_depth.py", "unfiltered"),
         ("c09_base_agreement.py", "filtered"),        # c08-c10 added in revision after verification (round 1)
         ("c10_compare_round0.py", "filtered"),
         ("c13_compare_round1.py", "filtered"),
         ("c06_provenance.py", "filtered")]


def main():
    os.makedirs(f"{C.RES}/logs", exist_ok=True)
    if os.path.exists(C.PROV_FILE):
        os.remove(C.PROV_FILE)          # regenerated below; holds only this item's own provenance lines
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1")
    for s, base in STEPS:
        log = f"{C.RES}/logs/{s.replace('.py', '')}_{base}.log"
        with open(log, "w", encoding="utf-8") as fh:
            r = subprocess.run([sys.executable, os.path.join(HERE, s)], stdout=fh, stderr=subprocess.STDOUT,
                               env=dict(env, C_BASE=base), cwd=HERE)
        print(s, base, "exit", r.returncode, flush=True)
        if r.returncode != 0:
            sys.exit(f"step {s} failed; see {log}")


if __name__ == "__main__":
    main()
