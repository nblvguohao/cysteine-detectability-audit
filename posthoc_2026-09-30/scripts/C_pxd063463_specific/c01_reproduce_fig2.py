"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific, task 1.

Reproduce the stored Figure 2 values of PXD063463 (four HydP arms x own/trypsin rule x proteome/observed background x
both bands, plus coincidence) with the UniProt 2026_03 mouse FASTA and the current Cys-Audit source (v0.2.2).

Steps
  1. Check every (protein, position) of the eight arm tables against the new FASTA; drop rows whose protein is absent
     or whose position is not a cysteine; write the filtered tables and a per-row list of what was dropped.
  2. Re-run the eleven stored audits (same command line as Phase 4: --expand-background, --modification
     S-palmitoylation, --identity-readout chemistry_inferred, default seed 20260922, 5000 replicates) on the filtered
     tables with v0.2.2; three further audits (trypsin rule, observed background, non-tryptic arms) are added and
     flagged as not part of the stored figure.
  3. Re-run the eleven audits on the unfiltered tables (the tool's native handling of absent proteins) and with the
     archived v0.1.0 source (subprocess) to separate FASTA, filtering and version effects.
  4. Compare every cell with the stored audit.json of repo/results/phase4_reports.
Outputs (results/C_pxd063463_specific/): fasta_check.csv, dropped_rows.csv, inputs_filtered/*.tsv,
  t1_reports*/<name>/, t1_fig2_reproduction.csv, t1_summary.json
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402

from cys_audit import cli  # noqa: E402
from cys_audit.io import read_fasta  # noqa: E402

PLAN = []
for arm in C.ARMS:
    PLAN.append((arm, C.OWN_RULE[arm], "proteome", True))
    PLAN.append((arm, C.OWN_RULE[arm], "observed", True))
    if arm != "Trypsin":
        PLAN.append((arm, "trypsin", "proteome", True))
        PLAN.append((arm, "trypsin", "observed", False))   # not in the stored figure


def name_of(arm, rule, bg):
    return f"pxd063463_{arm}_HydP_{rule}_{bg}"


def step1_filter(seqs):
    os.makedirs(f"{C.RES}/inputs_filtered", exist_ok=True)
    rows, dropped = [], []
    for arm in C.ARMS:
        for hyd in ("HydP", "HydN"):
            d = pd.read_csv(C.raw_path(arm, hyd), sep="\t", dtype={"protein": str}, keep_default_na=False)
            absent = ~d.protein.isin(set(seqs))
            beyond = np.array([p in seqs and pos > len(seqs[p]) for p, pos in zip(d.protein, d.position)])
            notc = np.array([p in seqs and pos <= len(seqs[p]) and seqs[p][pos - 1] != "C"
                             for p, pos in zip(d.protein, d.position)])
            bad = absent.values | beyond | notc
            for i in np.flatnonzero(bad):
                r = d.iloc[i]
                dropped.append({"arm": arm, "hyd": hyd, "protein": r.protein, "position": int(r.position),
                                "label": int(r.label), "reason": "protein absent from 2026_03 FASTA" if absent.iloc[i]
                                else ("position beyond sequence" if beyond[i] else "residue is not C")})
            keep = d.loc[~bad]
            # write exactly the original columns and cell text (keep_default_na=False keeps blanks as blanks)
            keep.to_csv(C.filtered_path(arm, hyd), sep="\t", index=False, lineterminator="\n")
            rows.append({"arm": arm, "hyd": hyd, "rows": len(d), "cam_rows": int((d.label == 1).sum()),
                         "proteins": d.protein.nunique(), "protein_absent_rows": int(absent.sum()),
                         "protein_absent_proteins": d.protein[absent].nunique(),
                         "position_beyond_rows": int(beyond.sum()), "not_cys_rows": int(notc.sum()),
                         "dropped_rows": int(bad.sum()), "dropped_cam_rows": int((d.label[bad] == 1).sum()),
                         "rows_kept": len(keep), "cam_rows_kept": int((keep.label == 1).sum())})
    pd.DataFrame(rows).to_csv(f"{C.RES}/fasta_check.csv", index=False)
    pd.DataFrame(dropped).to_csv(f"{C.RES}/dropped_rows.csv", index=False)
    return pd.DataFrame(rows)


def run_cli(inp, arm, rule, bg, out_root):
    name = name_of(arm, rule, bg)
    out = f"{out_root}/{name}"
    rc = cli.main(["audit", "--input", inp, "--protease", rule, "--background", bg, "--fasta", C.FASTA,
                   "--expand-background", "--modification", "S-palmitoylation", "--identity-readout",
                   "chemistry_inferred", "--output", out, "--dataset", name])
    if rc != 0:
        raise RuntimeError(f"cys-audit returned {rc} for {name}")
    return json.load(open(f"{out}/audit.json", encoding="utf-8"))


V010_RUNNER = r'''
import json, sys
sys.path.insert(0, sys.argv[1])
from cys_audit import cli
from cys_audit import constants
assert constants.VERSION == "0.1.0", constants.VERSION
args = json.loads(sys.argv[2])
sys.exit(cli.main(args))
'''


def run_v010(inp, arm, rule, bg, out_root):
    name = name_of(arm, rule, bg)
    out = f"{out_root}/{name}"
    args = ["audit", "--input", inp, "--protease", rule, "--background", bg, "--fasta", C.FASTA,
            "--expand-background", "--modification", "S-palmitoylation", "--identity-readout",
            "chemistry_inferred", "--output", out, "--dataset", name]
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, "-c", V010_RUNNER, C.TOOL_V010, json.dumps(args)], capture_output=True,
                       text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"v0.1.0 failed for {name}: {r.stderr[-2000:]}")
    return json.load(open(f"{out}/audit.json", encoding="utf-8"))


def cells(rec):
    """Flatten an audit record into the Figure 2 quantities."""
    out = {}
    cg = rec["tests"]["cleavage_geometry"]
    bands = (cg.get("details") or {}).get("bands") or {}
    for b in ("proximal_1_3", "distal_6_12"):
        t = bands.get(b)
        out[b] = None if t is None else {"estimate": t["estimate"], "ci_low": t["ci"][0], "ci_high": t["ci"][1],
                                         "status": t["status"], "n_positive": cg.get("n_positive"),
                                         "n_background": cg.get("n_background"), "n_clusters": cg.get("n_clusters"),
                                         "share_flagged_positive": t.get("share_flagged_positive"),
                                         "share_flagged_background": t.get("share_flagged_background")}
    ov = rec["tests"]["positive_detected_overlap"]
    out["coincidence"] = {"estimate": ov["estimate"], "ci_low": ov["ci"][0], "ci_high": ov["ci"][1],
                          "status": ov["status"], "n_positive": ov["details"]["positive_and_detected"],
                          "n_background": ov["details"]["detected"], "n_clusters": ov.get("n_clusters"),
                          "share_flagged_positive": None, "share_flagged_background": None}
    out["_meta"] = {"n_rows": rec["input"]["n_rows"], "notes": rec["input"]["notes"],
                    "headline_band": (cg.get("details") or {}).get("headline_band"), "cg_status": cg["status"],
                    "version": rec["version"], "fasta_sha256": rec["input"].get("fasta_sha256")}
    return out


def input_identity():
    """The arm tables in inputs/ carry CRLF line endings (Windows checkout); after LF normalisation their sha256 must
    equal the sites_sha256 recorded in each stored audit.json, i.e. the audited content."""
    out = {}
    for arm in C.ARMS:
        rec = json.load(open(f"{C.STORED_REPORTS}/{name_of(arm, C.OWN_RULE[arm], 'proteome')}/audit.json",
                             encoding="utf-8"))
        raw = open(C.raw_path(arm, "HydP"), "rb").read()
        lf = hashlib.sha256(raw.replace(b"\r\n", b"\n")).hexdigest()
        out[arm] = {"stored_sites_sha256": rec["input"]["sites_sha256"], "input_lf_normalised_sha256": lf,
                    "identical_content": lf == rec["input"]["sites_sha256"]}
    return out


def source_data_check():
    """Every value of Source Data Fig. 2 must equal the stored audit.json / phase4 demo value it cites."""
    rows = list(csv.DictReader(open(C.SOURCE_FIG2, encoding="utf-8")))
    demo = list(csv.DictReader(open(C.STORED_DEMO, encoding="utf-8")))
    bad = []
    for r in rows:
        if r["panel"] in ("a", "b"):
            arm, rule, bg, band = r["row"].split("|")
            a = json.load(open(f"{C.STORED_REPORTS}/{name_of(arm, rule, bg)}/audit.json", encoding="utf-8"))
            t = a["tests"]["cleavage_geometry"]["details"]["bands"][band]
            v = {"estimate": t["estimate"], "ci_low": t["ci"][0], "ci_high": t["ci"][1]}[r["field"]]
        else:
            d = [x for x in demo if x["arm"] == r["row"] and x["test"] == "positive_detected_overlap"][0]
            v = float(d[r["field"]])
        if abs(float(r["value"]) - v) > 1e-9:
            bad.append(r)
    return {"fields_checked": len(rows), "mismatches": bad}


def main():
    fasta = C.ensure_fasta()
    C.record_inputs("c01_reproduce_fig2", [C.FASTA_GZ] + [C.raw_path(a, h) for a in C.ARMS for h in ("HydP", "HydN")]
                    + [C.STORED_DEMO, C.SOURCE_FIG2]
                    + [f"{C.STORED_REPORTS}/{name_of(a, r, b)}/audit.json" for a, r, b, s in PLAN if s],
                    extra={"fasta_decompressed_sha256": C.sha256_file(fasta)})
    seqs = read_fasta(fasta)
    fc = step1_filter(seqs)
    print(fc.to_string(index=False), flush=True)

    rows = []
    summary = {"plan": [list(p) for p in PLAN], "input_identity": input_identity(),
               "source_data_fig2_check": source_data_check()}
    print(json.dumps(summary["input_identity"], indent=1), json.dumps(summary["source_data_fig2_check"]), flush=True)
    for arm, rule, bg, stored in PLAN:
        name = name_of(arm, rule, bg)
        new = cells(run_cli(C.filtered_path(arm, "HydP"), arm, rule, bg, f"{C.RES}/t1_reports"))
        nat = cells(run_cli(C.raw_path(arm, "HydP"), arm, rule, bg, f"{C.RES}/t1_reports_native")) if stored else None
        v010 = cells(run_v010(C.filtered_path(arm, "HydP"), arm, rule, bg, f"{C.RES}/t1_reports_v010")) if stored else None
        old = cells(json.load(open(f"{C.STORED_REPORTS}/{name}/audit.json", encoding="utf-8"))) if stored else None
        for q in ("proximal_1_3", "distal_6_12", "coincidence"):
            if q == "coincidence" and (rule != C.OWN_RULE[arm] or bg != "observed"):
                continue          # coincidence does not depend on rule or background; report it once per arm
            r = {"arm": arm, "rule": rule, "background": bg, "quantity": q, "in_stored_figure": stored}
            for tag, src in (("stored", old), ("new", new), ("native", nat), ("v010", v010)):
                c = None if src is None else src[q]
                for f in ("estimate", "ci_low", "ci_high", "status", "n_positive", "n_background", "n_clusters"):
                    r[f"{tag}_{f}"] = None if c is None else c[f]
            if old is not None:
                for f in ("estimate", "ci_low", "ci_high"):
                    r[f"delta_{f}"] = r[f"new_{f}"] - r[f"stored_{f}"]
                    r[f"delta_native_{f}"] = r[f"native_{f}"] - r[f"stored_{f}"]
                    r[f"v010_minus_new_{f}"] = r[f"v010_{f}"] - r[f"new_{f}"]
                if q != "coincidence":
                    ex_old = r["stored_ci_low"] > 0 or r["stored_ci_high"] < 0
                    ex_new = r["new_ci_low"] > 0 or r["new_ci_high"] < 0
                    r["stored_excludes_zero"], r["new_excludes_zero"] = ex_old, ex_new
                    r["same_zero_conclusion"] = ex_old == ex_new
                    r["same_status"] = r["stored_status"] == r["new_status"]
            rows.append(r)
        summary[name] = {"new_meta": new["_meta"], "stored_meta": None if old is None else old["_meta"],
                         "native_meta": None if nat is None else nat["_meta"]}
        print("done", name, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(f"{C.RES}/t1_fig2_reproduction.csv", index=False)
    st = df[df.in_stored_figure]
    band = st[st.quantity != "coincidence"]
    summary["agreement"] = {
        "cells_compared_bands": int(len(band)), "cells_compared_coincidence": int((st.quantity == "coincidence").sum()),
        "max_abs_delta_estimate_bands": float(band.delta_estimate.abs().max()),
        "max_abs_delta_ci_bound_bands": float(np.nanmax(np.abs(band[["delta_ci_low", "delta_ci_high"]].values))),
        "median_abs_delta_estimate_bands": float(band.delta_estimate.abs().median()),
        "max_abs_delta_estimate_coincidence": float(st[st.quantity == "coincidence"].delta_estimate.abs().max()),
        "status_changes": band.loc[~band.same_status.astype(bool), ["arm", "rule", "background", "quantity", "stored_status",
                                                        "new_status"]].to_dict("records"),
        "zero_conclusion_changes": band.loc[~band.same_zero_conclusion.astype(bool), ["arm", "rule", "background", "quantity"]]
        .to_dict("records"),
        "max_abs_v010_minus_new_all": float(np.nanmax(np.abs(st[["v010_minus_new_estimate", "v010_minus_new_ci_low",
                                                                 "v010_minus_new_ci_high"]].values))),
        "max_abs_native_minus_stored_estimate_bands": float(band.delta_native_estimate.abs().max()),
    }
    with open(f"{C.RES}/t1_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, default=str)
    print(json.dumps(summary["agreement"], indent=1, default=str))


if __name__ == "__main__":
    main()
