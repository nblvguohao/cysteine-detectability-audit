# -*- coding: utf-8 -*-
"""E_bootstrap_coverage -- POST HOC revision analysis (2026-09-30). Not registered.

Assembles the 28 claims with a measured baseline (Source_Data_Fig6 rows) with every interval the
verdict reads (baseline, matched, stratified = secondary caliber, same-size random control), from
Supplemental Data 4/5/6 in W/inputs/mcp_package, plus the number of bootstrap clusters (from the
stored audit records of the re-test rounds in the read-only repository, because Supplemental Data 4
does not carry it).

Then:
  * re-derives every verdict from the stored intervals with the pipeline's own rule and checks it
    equals the stored verdict for all 28. The rule's functions (classify, DOWNGRADE, classify_null,
    crosses_zero) are compiled VERBATIM from the repository by AST extraction (simlib.pipeline()); the
    FLOW that combines them is RE-TYPED in verdict_from_intervals below: the primary/secondary caliber
    disagreement downgrade of run_contingency_claim (run_phase2_claims_under_detectability_control.py),
    and the null-claim branch and baseline-direction guard of run_phase2b_claims_backfill.run_one /
    run_phase2d_claims.run_one (run_phase2e_claims.py imports the latter). [Docstring corrected in the
    revision after verification, 2026-09-30; the code is unchanged and its outputs are identical.]
  * for every interval computes w* = |point| / (half-width toward 0): the factor by which that
    half-width would have to be multiplied for the bound to reach 0 (w* > 1: excludes 0);
  * marks an interval as DECISIVE when toggling its zero-exclusion alone changes the verdict.
Writes results/E_bootstrap_coverage/claims_intervals.csv and claims_table_check.json.
"""
from __future__ import annotations

import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd

import simlib as S

W = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30"
PKG = os.path.join(W, "inputs", "mcp_package")
REPO_RES = r"C:/Users/admin/Desktop/小论文/_cys_repo_work/repo/results"
OUT = os.path.join(W, "results", "E_bootstrap_coverage")
AUDITS = ["phase2_claim_retest_audit.json", "phase2b_claim_retest_audit.json",
          "phase2d_claim_retest_audit.json", "phase2e_claim_retest_audit.json"]

# Which stored row carries each Fig. 6 claim (checked against the Fig. 6 baseline value below)
ROW_SELECT = {"SD5": {"PERS-008": "zf_background"}, "SD6": {c: "primary" for c in
                                                            ("SFE-001", "SNO-001", "SNO-002", "SNO-009", "SNO-012")}}
AUDIT_KEY = {"PERS-008": "PERS-008|zf_background"}


def load_groups():
    groups = {}
    for f in AUDITS:
        a = json.load(open(os.path.join(REPO_RES, f), encoding="utf-8"))
        recs = a.get("records") or a.get("claims")
        it = recs.items() if isinstance(recs, dict) else [(r["claim_id"], r) for r in recs]
        for k, r in it:
            pdg = r.get("precheck_permutation_degeneracy") or {}
            if pdg.get("groups") is not None:
                groups[k] = {"groups": int(pdg["groups"]), "mean_group_size": float(pdg["mean_group_size"]),
                             "singleton_groups": int(pdg["singleton_groups"]), "source": f}
                if "baseline_table" in r:
                    groups[k]["baseline_table"] = r["baseline_table"]
    return groups


def verdict_from_intervals(r, P):
    """The stored rule, applied to stored intervals (no recomputation of any estimate)."""
    b = float(r["baseline_log2_or"])
    bi = [float(r["baseline_ci_low"]), float(r["baseline_ci_high"])]
    c = float(r["matched_log2_or"])
    ci = [float(r["matched_ci_low"]), float(r["matched_ci_high"])]
    mh = float(r["stratified_log2_or"])
    mi = [float(r["stratified_ci_low"]), float(r["stratified_ci_high"])]
    ri = [float(r["random_control_ci_low"]), float(r["random_control_ci_high"])]
    blocked = bool(r.get("precheck_blocked", False))
    primary, _ = P["classify"](b, c, ci, ri, blocked)
    secondary, _ = P["classify"](b, mh, mi, ri, blocked)
    final = P["DOWNGRADE"].get(primary, primary) if primary != secondary else primary
    if P["crosses_zero"](bi):
        agrees = "no_baseline_effect"
    elif r["claim_direction"] == "null_no_preference":
        agrees = "no"
    else:
        agrees = "yes" if b > 0 else "no"
    if r["claim_direction"] == "null_no_preference":
        final, _ = P["classify_null"](bi, ci)
    elif agrees == "no":
        final = "baseline_contradicts_claim"
    return final, primary, secondary


def ratio(point, lo, hi):
    return float(S.ratio_to_zero(point, lo, hi))


INTERVALS = {"baseline": ("baseline_log2_or", "baseline_ci_low", "baseline_ci_high"),
             "matched": ("matched_log2_or", "matched_ci_low", "matched_ci_high"),
             "stratified": ("stratified_log2_or", "stratified_ci_low", "stratified_ci_high"),
             "random": ("random_control_log2_or", "random_control_ci_low", "random_control_ci_high")}


def toggle(r, which):
    """Copy of the row with interval `which` moved to the other side of 0 (point unchanged)."""
    r2 = dict(r)
    p, lo, hi = (float(r[k]) for k in INTERVALS[which])
    excludes = not (lo <= 0 <= hi)
    if excludes:  # make it cover 0 by pulling the near bound just across 0
        if p >= 0:
            r2[INTERVALS[which][1]] = -1e-9
        else:
            r2[INTERVALS[which][2]] = 1e-9
    else:  # make it exclude 0 by pulling the near bound just across 0 toward the point
        if p >= 0:
            r2[INTERVALS[which][1]] = 1e-9
        else:
            r2[INTERVALS[which][2]] = -1e-9
    return r2


def main():
    P = S.pipeline()
    fig6 = pd.read_csv(os.path.join(PKG, "Source_Data_Fig6_claim_retests.csv"))
    wide = fig6.pivot_table(index="row", columns="field", values="value", aggfunc="first")
    src = fig6.groupby("row").source_table.first()
    sd = {"SD4": pd.read_csv(os.path.join(PKG, "Supplemental_Data_4_retest_round_b.csv")),
          "SD5": pd.read_csv(os.path.join(PKG, "Supplemental_Data_5_retest_round_d.csv")),
          "SD6": pd.read_csv(os.path.join(PKG, "Supplemental_Data_6_retest_round_e.csv"))}
    tally = pd.read_csv(os.path.join(PKG, "Supplemental_Data_2_claim_verdict_tally.csv"))
    groups = load_groups()
    rows, mismatches = [], []
    for claim in wide.index:
        tab = {"Supplemental_Data_4_retest_round_b.csv": "SD4", "Supplemental_Data_5_retest_round_d.csv": "SD5",
               "Supplemental_Data_6_retest_round_e.csv": "SD6"}[src[claim]]
        t = sd[tab]
        cand = t[t.claim_id == claim]
        if tab in ROW_SELECT and claim in ROW_SELECT[tab]:
            cand = cand[cand.specification == ROW_SELECT[tab][claim]]
        cand = cand[np.isclose(cand.baseline_log2_or.astype(float), float(wide.loc[claim, "baseline_log2_or"]),
                               atol=5e-4)]
        assert len(cand) == 1, (claim, len(cand))
        r = cand.iloc[0].to_dict()
        for f in ("matched_log2_or", "matched_ci_low", "random_control_ci_high"):
            assert abs(float(r[f]) - float(wide.loc[claim, f])) < 5e-4, (claim, f)
        final, primary, secondary = verdict_from_intervals(r, P)
        stored = str(r["verdict"])
        tally_v = tally.loc[tally.claim_id == claim, "verdict"].iloc[0]
        if final != stored or stored != tally_v:
            mismatches.append({"claim": claim, "recomputed": final, "stored": stored, "tally": tally_v})
        key = AUDIT_KEY.get(claim, claim)
        g = groups.get(key) or groups.get(f"{claim}|primary") or {}
        rec = {"claim_id": claim, "source_table": tab, "specification": r.get("specification", "primary"),
               "claim_direction": r["claim_direction"], "unit": r["unit"],
               "n_observations": int(float(r["n_observations"])), "n_positive": int(float(r["n_positive"])),
               "n_attribute": (int(float(r["n_attribute"])) if pd.notna(r.get("n_attribute")) else None),
               "n_matched_pairs": int(float(r["n_matched_pairs"])),
               "bootstrap_clusters": g.get("groups"), "mean_cluster_size": g.get("mean_group_size"),
               "singleton_clusters": g.get("singleton_groups"), "clusters_source": g.get("source"),
               "is_transfer": bool(str(r.get("is_transfer", "False")) == "True"),
               "baseline_reproduction": r.get("baseline_reproduction"),
               "stored_verdict": stored, "recomputed_verdict": final,
               "primary_verdict_recomputed": primary, "secondary_verdict_recomputed": secondary}
        bt = g.get("baseline_table")
        if bt:
            rec.update({f"baseline_{k}": v for k, v in bt.items()})
        for name, (pk, lk, hk) in INTERVALS.items():
            p_, lo, hi = float(r[pk]), float(r[lk]), float(r[hk])
            rec[f"{name}_point"], rec[f"{name}_lo"], rec[f"{name}_hi"] = p_, lo, hi
            rec[f"{name}_excludes0"] = not (lo <= 0 <= hi)
            rec[f"{name}_wstar"] = ratio(p_, lo, hi)
            v_t, _, _ = verdict_from_intervals(toggle(r, name), P)
            rec[f"{name}_decisive"] = v_t != final
            rec[f"{name}_verdict_if_toggled"] = v_t
        rows.append(rec)
    df = pd.DataFrame(rows).sort_values("claim_id")
    df.to_csv(os.path.join(OUT, "claims_intervals.csv"), index=False)
    chk = {"label": "POST HOC revision analysis E_bootstrap_coverage (not registered)",
           "n_claims": int(len(df)), "verdicts_reproduced_from_stored_intervals": int((df.stored_verdict == df.recomputed_verdict).sum()),
           "mismatches": mismatches,
           "stored_verdict_counts": df.stored_verdict.value_counts().to_dict(),
           "inputs_sha256": {f: S.sha256(os.path.join(PKG, f)) for f in (
               "Source_Data_Fig6_claim_retests.csv", "Supplemental_Data_4_retest_round_b.csv",
               "Supplemental_Data_5_retest_round_d.csv", "Supplemental_Data_6_retest_round_e.csv",
               "Supplemental_Data_2_claim_verdict_tally.csv")},
           "repo_audit_sha256": {f: S.sha256(os.path.join(REPO_RES, f)) for f in AUDITS}}
    with open(os.path.join(OUT, "claims_table_check.json"), "w", encoding="utf-8") as fh:
        json.dump(chk, fh, indent=1)
    print(json.dumps({k: v for k, v in chk.items() if "sha" not in k}, indent=1))
    cols = ["claim_id", "unit", "bootstrap_clusters", "n_positive", "n_matched_pairs", "stored_verdict"]
    for n in INTERVALS:
        cols += [f"{n}_wstar", f"{n}_decisive"]
    pd.set_option("display.width", 250)
    print(df[cols].round(3).to_string())


if __name__ == "__main__":
    main()
