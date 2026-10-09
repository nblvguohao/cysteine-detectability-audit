"""Run item I_fig5a_estimand end to end (POST HOC revision analysis, 2026-09-30).

    cd scripts/I_fig5a_estimand
    python -B run_all.py

-B keeps Python from writing bytecode caches, in particular into the read-only repository
whose feature module s02 imports.
"""
import runpy
import sys

sys.dont_write_bytecode = True
for step in ("s01_tables", "s02_feature_identity", "s03_input_availability", "s04_plot_draft", "s05_provenance"):
    print(f"==== {step}", flush=True)
    runpy.run_path(f"{step}.py", run_name="__main__")
