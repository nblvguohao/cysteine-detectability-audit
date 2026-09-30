"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific, task 3.

Cleavage-geometry audit (Cys-Audit v0.2.2 checks.cleavage.run, unchanged: Haldane log2 OR of a cleavage-competent
residue within 1-3 and 6-12 residues, protein-clustered percentile bootstrap, 5000 replicates, 97.5% intervals, seed
20260922 as in the stored audits) re-run with hydroxylamine-specific positive sets from task 2.

Positive sets per arm (HydP arm; filtered to the 2026_03 FASTA):
  cam            the stored definition (CAM, localisation >= 0.75), re-run here as the reference
  strict_<src>   CAM in HydP and identified without CAM in HydN
  lenient_<src>  CAM in HydP and not CAM in HydN
  nonspecific_<src> (control) CAM in HydP and CAM in HydN
  <src> = matched (same-protease HydN arm, primary) or pooled (union of the four HydN arms, sensitivity)
  cam_protein_absent_pooled (descriptive) the stored definition restricted to proteins never identified in any HydN
                 arm (handling 'restrict_proteins': both classes restricted to those proteins)
  Selection split (added in revision after verification, round 1; 'exclude' handling only):
  in_trypsin_hydn / not_in_trypsin_hydn   CAM sites whose cysteine was / was not identified (with or without CAM) in
                 the TRYPTIC HydN arm. The pooled HydN identifications are 95.5% tryptic, and identification in a tryptic
                 digest itself selects cysteines with tryptic cleavage geometry, so every pooled definition and the
                 strict trypsin set are conditioned on it; this split shows how much.
  in_pooled_hydn / not_in_pooled_hydn     the same split by identification in any HydN arm.
Handling of HydP CAM sites that fail a definition:
  exclude (primary)  removed from both classes, so the background is exactly the stored background
  demote (sensitivity) kept as label-0 rows (they join the observed and proteome backgrounds)
Backgrounds: proteome (every other cysteine of the identified proteins, --expand-background) and observed
(identified cysteines of the HydP arm that are not positives). Rules: the arm's own rule and the trypsin rule.
A set with fewer than 20 positives is UNDECIDABLE by the tool's fixed minimum; its flagged shares are still listed.

Also audits the HydN arms themselves (CAM in the arm without hydroxylamine as 'positives') with the tool CLI.
Outputs: t3_cleavage_specific.csv, t3_hydn_arm_audits.csv, t3_hydn_reports/<name>/
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C  # noqa: E402

from cys_audit import cli  # noqa: E402
from cys_audit.checks import cleavage  # noqa: E402
from cys_audit.constants import CLEAVAGE_BANDS  # noqa: E402
from cys_audit.io import background_mask, read_sites  # noqa: E402

SETS = ["cam", "strict_matched", "lenient_matched", "nonspecific_matched", "strict_pooled", "lenient_pooled",
        "nonspecific_pooled"]
SPLIT_SETS = ["in_trypsin_hydn", "not_in_trypsin_hydn", "in_pooled_hydn", "not_in_pooled_hydn"]


def flagged_shares(ds, rule, bg):
    """Share of positives / background with a cleavage-competent residue in each band (descriptive)."""
    flags, has_seq = cleavage.band_flags(ds, rule)
    bgm = background_mask(ds, bg) & has_seq
    pos = (ds.label == 1) & has_seq
    out = {}
    for b in CLEAVAGE_BANDS:
        out[b] = (float(flags[b][pos].mean()) if pos.any() else np.nan,
                  float(flags[b][bgm].mean()) if bgm.any() else np.nan)
    return out, int(pos.sum()), int(bgm.sum())


def run_set(ds, rule, bg):
    res = cleavage.run(ds, {"protease": rule, "background": bg, "seed": C.TOOL_SEED, "reps": C.REPS})
    shares, npos, nbg = flagged_shares(ds, rule, bg)
    rows = []
    bands = (res.get("details") or {}).get("bands") or {}
    for b in CLEAVAGE_BANDS:
        t = bands.get(b)
        rows.append({"band": b, "n_positive": npos, "n_background": nbg, "n_clusters": res.get("n_clusters"),
                     "estimate": None if t is None else t["estimate"],
                     "ci_low": None if t is None else t["ci"][0], "ci_high": None if t is None else t["ci"][1],
                     "status": res["status"] if t is None else t["status"],
                     "reason": res["reason"] if t is None else t["reason"],
                     "share_flagged_positive": shares[b][0], "share_flagged_background": shares[b][1],
                     "below_minimum": t is None})
    return rows


def main():
    C.ensure_fasta()
    C.record_inputs("c03_cleavage_specific", [C.filtered_path(a, h) for a in C.ARMS for h in ("HydP", "HydN")]
                    + [f"{C.RES}/t2_site_table_{a}.tsv" for a in C.ARMS] + [C.FASTA])
    out = []
    for arm in C.ARMS:
        ds = read_sites(C.filtered_path(arm, "HydP"), fasta=C.FASTA, expand_background=True)
        st = pd.read_csv(f"{C.RES}/t2_site_table_{arm}.tsv", sep="\t", dtype={"protein": str})
        idx = {(p, int(q)): i for i, (p, q) in enumerate(zip(st.protein, st.position))}
        rowmap = np.array([idx.get((p, int(q)), -1) for p, q in zip(ds.protein, ds.position)])
        listed = rowmap >= 0
        assert listed.sum() == len(st), "every HydP-identified cysteine must be a listed row"
        cam = ds.label == 1
        assert np.array_equal(cam[listed], st.cam_hydp.values[rowmap[listed]] == 1)
        prot_in_hydn = set(st.protein[st.protein_in_hydn_pooled == 1])
        # selection split (revision after verification, round 1)
        camst = st.cam_hydp == 1
        st["in_trypsin_hydn"] = (camst & (st.in_hydn_trypsin_arm == 1)).astype(int)
        st["not_in_trypsin_hydn"] = (camst & (st.in_hydn_trypsin_arm == 0)).astype(int)
        st["in_pooled_hydn"] = (camst & (st.in_hydn_pooled == 1)).astype(int)
        st["not_in_pooled_hydn"] = (camst & (st.in_hydn_pooled == 0)).astype(int)
        for sname in SETS + SPLIT_SETS + ["cam_protein_absent_pooled"]:
            if sname in ("cam", "cam_protein_absent_pooled"):
                spec = cam.copy()
            else:
                spec = np.zeros(len(ds), dtype=bool)
                spec[listed] = st[sname].values[rowmap[listed]] == 1
            ambiguous = cam & ~spec
            if sname == "cam":
                handlings = ["none"]
            elif sname == "cam_protein_absent_pooled":
                handlings = ["restrict_proteins"]
            elif sname in SPLIT_SETS:
                handlings = ["exclude"]
            else:
                handlings = ["exclude", "demote"]
            for h in handlings:
                if h == "restrict_proteins":
                    # keep only proteins never identified in any HydN arm (both classes, within-protein pairing kept)
                    keep = np.array([p not in prot_in_hydn for p in ds.protein])
                    d2 = ds.subset(keep)
                elif h == "exclude":
                    d2 = ds.subset(~ambiguous)
                    d2.label = spec[~ambiguous].astype(ds.label.dtype)
                elif h == "demote":
                    d2 = ds.subset(np.ones(len(ds), dtype=bool))
                    d2.label = spec.astype(ds.label.dtype)
                else:
                    d2 = ds
                rules = [C.OWN_RULE[arm]] + ([] if arm == "Trypsin" else ["trypsin"])
                for rule in rules:
                    for bg in ("proteome", "observed"):
                        for r in run_set(d2, rule, bg):
                            r.update({"arm": arm, "positive_set": sname, "handling": h, "rule": rule,
                                      "rule_type": "own" if rule == C.OWN_RULE[arm] else "trypsin", "background": bg,
                                      "n_sites_in_set": int((d2.label == 1).sum()),
                                      "n_ambiguous_cam": int(ambiguous.sum()) if h in ("exclude", "demote") else 0})
                            out.append(r)
        print("done", arm, flush=True)
    cols = ["arm", "positive_set", "handling", "rule", "rule_type", "background", "band", "n_sites_in_set",
            "n_ambiguous_cam", "n_positive", "n_background", "n_clusters", "estimate", "ci_low", "ci_high", "status",
            "below_minimum", "share_flagged_positive", "share_flagged_background", "reason"]
    df = pd.DataFrame(out)[cols]
    df.to_csv(f"{C.RES}/t3_cleavage_specific.csv", index=False)

    # HydN arms audited on their own (CLI, full report)
    hrows = []
    for arm in C.ARMS:
        rules = [C.OWN_RULE[arm]] + ([] if arm == "Trypsin" else ["trypsin"])
        for rule in rules:
            for bg in ("proteome", "observed"):
                name = f"pxd063463_{arm}_HydN_{rule}_{bg}"
                o = f"{C.RES}/t3_hydn_reports/{name}"
                rc = cli.main(["audit", "--input", C.filtered_path(arm, "HydN"), "--protease", rule, "--background", bg,
                               "--fasta", C.FASTA, "--expand-background", "--modification", "S-palmitoylation",
                               "--identity-readout", "chemistry_inferred", "--output", o, "--dataset", name])
                assert rc == 0
                rec = json.load(open(f"{o}/audit.json", encoding="utf-8"))
                cg = rec["tests"]["cleavage_geometry"]
                bands = (cg.get("details") or {}).get("bands") or {}
                for b in CLEAVAGE_BANDS:
                    t = bands.get(b)
                    hrows.append({"arm": arm, "arm_type": "HydN", "rule": rule,
                                  "rule_type": "own" if rule == C.OWN_RULE[arm] else "trypsin", "background": bg,
                                  "band": b, "n_positive": cg.get("n_positive"), "n_background": cg.get("n_background"),
                                  "estimate": None if t is None else t["estimate"],
                                  "ci_low": None if t is None else t["ci"][0], "ci_high": None if t is None else t["ci"][1],
                                  "status": cg["status"] if t is None else t["status"],
                                  "reason": cg["reason"] if t is None else t["reason"],
                                  "coincidence": rec["tests"]["positive_detected_overlap"]["estimate"]})
    pd.DataFrame(hrows).to_csv(f"{C.RES}/t3_hydn_arm_audits.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    pd.set_option("display.max_rows", 500)
    show = df[(df.handling != "demote")]
    print(show[["arm", "positive_set", "rule_type", "background", "band", "n_positive", "n_background", "estimate",
                "ci_low", "ci_high", "status"]].to_string(index=False))
    print(pd.DataFrame(hrows)[["arm", "rule", "background", "band", "n_positive", "n_background", "estimate", "ci_low",
                               "ci_high", "status"]].to_string(index=False))


if __name__ == "__main__":
    main()
