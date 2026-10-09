"""POST HOC CORRECTION (2026-09-24): the public refit read with a single LightGBM member instead of the blend.

WHY
The pre-submission review of R55 (reports/PRESUBMISSION_REVIEW_R55_SUMMARY_2026-09-24.md, item A3) found that
every arm of the public refit -- including `detect_only`, described in the manuscript as "a 25-feature set that
contains no chemistry at all" -- was scored by the BLEND of three members (run_v2_chem_basic_ablation.run_arm):
lgb_rank, lgb_bin and fusion_mlp, and that fusion_mlp receives, besides the arm's columns, the 1,024-dimensional
ProtBERT embeddings (raw and protein-centred). The registered analysis (scripts/analyse_public_refit_2026-09-21.py,
protocols/public_refit_3_3_and_3_4_preregistration_2026-09-21.json) read the `blend` column. A blend with an
embedding-fed member is not a detectability-only model, and not the "single LightGBM binary classifier" the
manuscript describes, so the reported recovery ratios and Table 2 do not measure what the text says.

WHAT THIS SCRIPT DOES
Nothing is refitted. The stored out-of-fold files already hold each member's score. The registered analysis
module and the post hoc performance-cost module are imported UNCHANGED and run with their score column switched:
  lgb_bin   PRIMARY: a single LightGBM binary member, the same model class as the audited ranker; it sees only
            the arm's columns, so detect_only is genuinely the 25 detectability columns
  lgb_rank  sensitivity: the LightGBM ranking member
Outputs go to new files (results/public_refit_*_{member}_2026-09-24*). The 2026-09-21 outputs are kept as the
registered record; the manuscript reports the lgb_bin reading and states the change.

GATES
  G0 outputs do not exist before the run
  G1 positive control: the unchanged modules run through this wrapper with score column `blend` reproduce the
     stored registered outputs byte for byte (results/public_refit_ablation_2026-09-21.csv,
     results/public_refit_detectability_only_2026-09-21.csv, results/public_refit_performance_cost_posthoc_2026-09-21.csv)
  G2 the registered NC1 (label-permuted) gate passes for each member reading, or nothing from that member is used
INTERPRETER: project venv (Python 3.9.6 / numpy 2.0.2), as the registered analysis.
"""
import csv
import hashlib
import importlib.util
import json
import pathlib
import sys
import tempfile

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
RES = ROOT / "results"
A_PATH = ROOT / "scripts/analyse_public_refit_2026-09-21.py"
B_PATH = ROOT / "scripts/posthoc_public_refit_performance_cost_2026-09-21.py"
SELF = pathlib.Path(__file__).read_bytes()
MEMBERS = ("lgb_bin", "lgb_rank")
AUD = RES / "public_refit_lgb_member_reanalysis_2026-09-24_audit.json"
STORED = {"ablation": RES / "public_refit_ablation_2026-09-21.csv",
          "detect": RES / "public_refit_detectability_only_2026-09-21.csv",
          "cost": RES / "public_refit_performance_cost_posthoc_2026-09-21.csv"}


def sha(b):
    return hashlib.sha256(b).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run(column, out_ablation, out_detect, out_cost, audit_a, audit_b):
    A = load(A_PATH, "refit_a_" + column)
    orig_read = A.read_oof

    def read_oof(cohort, arm):
        d = orig_read(cohort, arm)
        f = A.OUTDIR / cohort / f"{cohort}_{arm}_oof.csv"
        rows = list(csv.DictReader(open(f, encoding="utf-8-sig")))
        d["score"] = np.asarray([float(r[column]) for r in rows])
        return d
    A.read_oof = read_oof
    A.OUT_CSV, A.OUT_DET, A.OUT_AUDIT, A.SELF = out_ablation, out_detect, audit_a, SELF
    A.main()

    B = load(B_PATH, "refit_b_" + column)

    def per_protein(cohort, arm):
        f = B.OUTDIR / cohort / f"{cohort}_{arm}_oof.csv"
        rows = list(csv.DictReader(open(f, encoding="utf-8-sig")))
        acc = np.asarray([r["accession"] for r in rows])
        y = np.asarray([int(r["label"]) for r in rows])
        s = np.asarray([float(r[column]) for r in rows])
        return {r["accession"]: r["auc"] for r in B.protein_rows(y, s, acc, acc)}
    B.per_protein = per_protein
    B.OUT, B.AUDIT, B.SELF = out_cost, audit_b, SELF
    B.main()


def main():
    outs = {m: dict(ablation=RES / f"public_refit_ablation_{m}_2026-09-24.csv",
                    detect=RES / f"public_refit_detectability_only_{m}_2026-09-24.csv",
                    cost=RES / f"public_refit_performance_cost_posthoc_{m}_2026-09-24.csv",
                    audit_a=RES / f"public_refit_ablation_{m}_2026-09-24_audit.json",
                    audit_b=RES / f"public_refit_performance_cost_posthoc_{m}_2026-09-24_audit.json")
            for m in MEMBERS}
    if AUD.exists() or any(p.exists() for o in outs.values() for p in o.values()):
        sys.exit("REFUSE: G0 outputs exist")

    # G1: the wrapper with the registered column reproduces the registered outputs
    with tempfile.TemporaryDirectory() as tmp:
        t = pathlib.Path(tmp)
        run("blend", t / "abl.csv", t / "det.csv", t / "cost.csv", t / "a.json", t / "b.json")
        g1 = {k: (t / f).read_bytes() == STORED[k].read_bytes()
              for k, f in (("ablation", "abl.csv"), ("detect", "det.csv"), ("cost", "cost.csv"))}
    if not all(g1.values()):
        sys.exit("REFUSE: G1 wrapper does not reproduce the registered outputs %s" % g1)

    summary = {}
    for m in MEMBERS:
        o = outs[m]
        run(m, o["ablation"], o["detect"], o["cost"], o["audit_a"], o["audit_b"])
        a = json.loads(o["audit_a"].read_text(encoding="utf-8"))
        a.update(wrapper_script=pathlib.Path(__file__).name, wrapper_sha256=sha(SELF), score_column=m,
                 status="POST HOC correction of the registered reading (which used `blend`); see wrapper docstring")
        o["audit_a"].write_text(json.dumps(a, ensure_ascii=False, indent=2), encoding="utf-8")
        nc1_ok = all(v.get("crosses_zero") and v.get("wp_covers_half")
                     for k, v in a["gates"].items() if k.startswith("NC1_"))
        summary[m] = dict(nc1_pass=bool(nc1_ok), gates_failed=a.get("gates_failed"),
                          ablation=list(csv.DictReader(open(o["ablation"], encoding="utf-8"))),
                          detect=list(csv.DictReader(open(o["detect"], encoding="utf-8"))),
                          cost=list(csv.DictReader(open(o["cost"], encoding="utf-8"))))
    AUD.write_text(json.dumps(dict(
        script=pathlib.Path(__file__).name, script_sha256=sha(SELF),
        imported_unchanged={A_PATH.name: sha(A_PATH.read_bytes()), B_PATH.name: sha(B_PATH.read_bytes())},
        G1_blend_reproduces_registered_outputs=g1, primary_member="lgb_bin", sensitivity_member="lgb_rank",
        outputs={m: {k: str(p.relative_to(ROOT)) for k, p in o.items()} for m, o in outs.items()},
        summary=summary, all_pass=all(s["nc1_pass"] and not s["gates_failed"] for s in summary.values())),
        ensure_ascii=False, indent=2), encoding="utf-8")
    for m in MEMBERS:
        print("==", m, "NC1 pass:", summary[m]["nc1_pass"], "failed:", summary[m]["gates_failed"])
        for r in summary[m]["detect"]:
            print("  detect", r)
        for r in summary[m]["ablation"]:
            print("  abl", {k: r[k] for k in ("cohort", "arm", "log2_or", "ci_low", "ci_high", "within_protein_auc")})
        for r in summary[m]["cost"]:
            print("  cost", r)


if __name__ == "__main__":
    main()
