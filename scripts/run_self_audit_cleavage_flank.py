"""Does OUR deployed model rank cysteines by the protease's own cleavage residues?

`reports/PHASE2E_STRUCTURAL_CLAIMS.md` and `reports/PHASE2D_CLAIMS_AFTER_MATERIAL_DROP.md`
established, inside one paper's own claim, that how much of a published site-preference survives
detectability control tracks whether the named residue participates in the protease's cleavage
rule: on the Yang 2014 cohort the glutamate half retains 0.829 of its baseline while the lysine
half retains 0.435, and SFE-006's K/R version retained 0.430. That finding indicts a mechanism,
and the same mechanism is available to OUR model, whose training label is a capture-workflow
annotation (tier T3 by `results/label_semantics_tiers.csv`). This script therefore turns the
instrument on ourselves before the manuscript claims anything about others.

**What is audited.** The two v2 blends as actually deployed, from stored out-of-fold scores - no
refitting:
  * `results/v2_stack_chem_oof.csv` blend - the WET-LAB model (deployment weights 0.64/0.20/0.16;
    this is what `collaboration/locked_predictions_2026-09-14/` was produced from),
  * `results/v2_stack_full_oof.csv` blend - the BENCHMARK model, expected to be worse here because
    it absorbs peptide visibility by construction,
  * the training LABEL itself, as the reference: a model can only be charged with amplifying an
    association the label does not already carry.
The 90 locked predictions handed to the collaborators are audited separately from their own
`context_31aa` windows.

**Attribute, predeclared.** Two flags, both computed from `inputs/sly_proteome.tsv`:
  * primary `kr_sfe006_offsets` - K or R at any of [-10,-8,-7,-6,-4,-2,4,5,6,7,8], the exact offset
    set of the claim that generated this hypothesis, so the two are numerically comparable;
  * secondary `kr_within_5` - K or R anywhere in +-5, a simpler offset-free version.

**Statistics, predeclared.** Cluster bootstrap over homology components, 5000 replicates, seed
20260915, all reusing `run_phase2_claims_under_detectability_control` unchanged:
  * `score -> attribute` weighted AUC: how strongly the model ranks K/R-flanked cysteines;
  * the same AUC computed WITHIN each label stratum (positives only, negatives only). This is the
    decomposition that matters: an association present overall but absent within strata is
    inherited from the label, not added by the model;
  * top-k attribute rate at k = 50, 100, 500, 1000 and 90 (the locked-set size) against the rate
    among all remaining sites, as a log2 odds ratio.

**Reading rule, fixed before the run.** Let A be the score->attribute AUC within label strata and
L the label's own attribute association.
  1. If A's interval covers 0.5 in BOTH strata - the deployed model does not rank by cleavage-rule
     residues beyond what the label carries. The self-audit clears it and the manuscript may say so.
  2. If A > 0.5 within strata AND the model's overall association exceeds L - the model AMPLIFIES
     the cleavage-rule signal. This is the same defect the manuscript documents in others and must
     be disclosed as such, with the locked-prediction top-k reported alongside its K/R rate.
  3. If the overall association exceeds 0.5 but A does not within strata - the association is
     INHERITED from the T3 label rather than added by the model. Disclose as inherited.
The one thing this script may not do is pick the branch after seeing the numbers.

Writes results/self_audit_cleavage_flank.csv, _topk.csv and _audit.json.
"""
from __future__ import annotations

import csv
import json
import os
import platform
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from phase2_claim_cohorts import SFE006_KR_OFFSETS
from run_phase2_claims_under_detectability_control import (
    REPLICATES, SEED, bootstrap_interval, log_odds_ratio, make_weighted_auc,
    make_weighted_table, sha256_of,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
INPUTS = os.path.join(ROOT, "inputs")
LOCKED = os.path.join(ROOT, "collaboration", "locked_predictions_2026-09-14",
                      "locked_site_predictions.csv")
CHEM = os.path.join(RESULTS, "v2_stack_chem_oof.csv")
FULL = os.path.join(RESULTS, "v2_stack_full_oof.csv")
PROTEOME = os.path.join(INPUTS, "sly_proteome.tsv")
SCRIPT = os.path.abspath(__file__)
TOPK = (50, 90, 100, 500, 1000)
WITHIN_5 = list(range(-5, 0)) + list(range(1, 6))


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def sequences():
    out = {}
    with open(PROTEOME, encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("Entry") or "").strip()
            seq = (row.get("Sequence") or "").strip().upper()
            if acc and seq:
                out[acc] = seq
    return out


def flank_flag(seq, position, offsets, residues=("K", "R")):
    """1 when any offset around the 1-based position carries K or R, else 0. None off-protein."""
    i = position - 1
    if not (0 <= i < len(seq)) or seq[i] != "C":
        return None
    for d in offsets:
        j = i + d
        if 0 <= j < len(seq) and seq[j] in residues:
            return 1
    return 0


def audit_scores(name, y, attribute, score, groups, rows_out):
    auc_fn = make_weighted_auc(attribute, score)
    point = auc_fn(np.ones(len(score)))
    interval, _ = bootstrap_interval(auc_fn, groups)
    row = {"model": name, "stratum": "all", "n": len(score),
           "attribute_rate": round(float(attribute.mean()), 4),
           "auc_score_to_attribute": round(float(point), 4),
           "ci_low": round(interval[0], 4), "ci_high": round(interval[1], 4)}
    rows_out.append(row)
    for label_value, stratum in ((1, "label_positive"), (0, "label_negative")):
        mask = y == label_value
        if mask.sum() < 20 or len(np.unique(attribute[mask])) < 2:
            rows_out.append({"model": name, "stratum": stratum, "n": int(mask.sum()),
                             "attribute_rate": round(float(attribute[mask].mean()), 4)
                             if mask.any() else None,
                             "auc_score_to_attribute": None, "ci_low": None, "ci_high": None})
            continue
        fn = make_weighted_auc(attribute[mask], score[mask])
        p = fn(np.ones(int(mask.sum())))
        iv, _ = bootstrap_interval(fn, groups[mask])
        rows_out.append({"model": name, "stratum": stratum, "n": int(mask.sum()),
                         "attribute_rate": round(float(attribute[mask].mean()), 4),
                         "auc_score_to_attribute": round(float(p), 4),
                         "ci_low": round(iv[0], 4), "ci_high": round(iv[1], 4)})
    return row


def audit_topk(name, attribute, score, groups, rows_out):
    order = np.argsort(-score, kind="stable")
    for k in TOPK:
        if k >= len(score):
            continue
        in_top = np.zeros(len(score), dtype=int)
        in_top[order[:k]] = 1
        table = make_weighted_table(in_top, attribute)
        point = log_odds_ratio(*table(np.ones(len(score))))
        interval, _ = bootstrap_interval(lambda w: log_odds_ratio(*table(w)), groups)
        rows_out.append({
            "model": name, "k": k,
            "attribute_rate_in_top_k": round(float(attribute[order[:k]].mean()), 4),
            "attribute_rate_rest": round(float(attribute[order[k:]].mean()), 4),
            "log2_or": round(float(point), 4),
            "ci_low": round(interval[0], 4), "ci_high": round(interval[1], 4)})


def main():
    started = time.time()
    seqs = sequences()
    chem = read_csv(CHEM)
    full = read_csv(FULL)
    assert len(chem) == len(full), (len(chem), len(full))
    for a, b in zip(chem[:200], full[:200]):
        assert (a["accession"], a["position"]) == (b["accession"], b["position"])

    keep, y, groups, chem_score, full_score, primary, secondary = [], [], [], [], [], [], []
    missing_protein = off_site = 0
    for a, b in zip(chem, full):
        seq = seqs.get(a["accession"])
        if seq is None:
            missing_protein += 1
            continue
        pos = int(a["position"])
        p = flank_flag(seq, pos, SFE006_KR_OFFSETS)
        s = flank_flag(seq, pos, WITHIN_5)
        if p is None or s is None:
            off_site += 1
            continue
        keep.append((a["accession"], pos))
        y.append(int(a["label"]))
        groups.append(a["component"])
        chem_score.append(float(a["blend"]))
        full_score.append(float(b["blend"]))
        primary.append(p)
        secondary.append(s)
    y = np.asarray(y, dtype=int)
    groups = np.asarray(groups)
    chem_score = np.asarray(chem_score, dtype=float)
    full_score = np.asarray(full_score, dtype=float)
    primary = np.asarray(primary, dtype=int)
    secondary = np.asarray(secondary, dtype=int)
    print(f"sites audited {len(y)} | positives {int(y.sum())} | components {len(set(groups))} "
          f"| dropped: protein absent {missing_protein}, site off-protein or not Cys {off_site}",
          flush=True)

    rows, topk_rows = [], []
    label_reference = {}
    for attr_name, attr in (("kr_sfe006_offsets", primary), ("kr_within_5", secondary)):
        # the label's own association, as the reference the model is judged against
        fn = make_weighted_auc(attr, y.astype(float))
        point = fn(np.ones(len(y)))
        iv, _ = bootstrap_interval(fn, groups)
        label_reference[attr_name] = {"auc_label_to_attribute": round(float(point), 4),
                                      "ci": [round(iv[0], 4), round(iv[1], 4)],
                                      "rate_in_positives": round(float(attr[y == 1].mean()), 4),
                                      "rate_in_negatives": round(float(attr[y == 0].mean()), 4)}
        rows.append({"model": f"LABEL[{attr_name}]", "stratum": "all", "n": len(y),
                     "attribute_rate": round(float(attr.mean()), 4),
                     "auc_score_to_attribute": round(float(point), 4),
                     "ci_low": round(iv[0], 4), "ci_high": round(iv[1], 4)})
        for model_name, score in (("v2_chem_blend_wetlab", chem_score),
                                  ("v2_full_blend_benchmark", full_score)):
            audit_scores(f"{model_name}[{attr_name}]", y, attr, score, groups, rows)
            audit_topk(f"{model_name}[{attr_name}]", attr, score, groups, topk_rows)
            print(f"  done {model_name}[{attr_name}] ({round(time.time()-started, 1)}s)", flush=True)

    # the 90 locked predictions, from their own windows
    locked = read_csv(LOCKED)
    locked_rows = []
    for r in locked:
        w = (r.get("context_31aa") or "").strip().upper()
        if len(w) != 31 or w[15] != "C":
            locked_rows.append({"ok": False}); continue
        centre = 16
        prim = 1 if any(w[centre - 1 + d] in ("K", "R") for d in SFE006_KR_OFFSETS
                        if 0 <= centre - 1 + d < 31) else 0
        sec = 1 if any(w[centre - 1 + d] in ("K", "R") for d in WITHIN_5
                       if 0 <= centre - 1 + d < 31) else 0
        locked_rows.append({"ok": True, "primary": prim, "secondary": sec,
                            "rank": int(r["chem_blend_rank"])})
    usable = [r for r in locked_rows if r["ok"]]
    locked_summary = {
        "n_locked": len(locked), "n_with_usable_window": len(usable),
        "kr_sfe006_rate": round(sum(r["primary"] for r in usable) / max(1, len(usable)), 4),
        "kr_within_5_rate": round(sum(r["secondary"] for r in usable) / max(1, len(usable)), 4),
        "cohort_rate_kr_sfe006": round(float(primary.mean()), 4),
        "cohort_rate_kr_within_5": round(float(secondary.mean()), 4),
    }

    with open(os.path.join(RESULTS, "self_audit_cleavage_flank.csv"), "w",
              encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(RESULTS, "self_audit_cleavage_flank_topk.csv"), "w",
              encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(topk_rows[0]))
        w.writeheader()
        w.writerows(topk_rows)

    audit = {
        "script": "scripts/run_self_audit_cleavage_flank.py", "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.node(),
        "round": "self-audit: does our deployed model rank by the protease's cleavage residues",
        "no_refit": ("stored out-of-fold scores only; results/v2_stack_chem_oof.csv and "
                     "results/v2_stack_full_oof.csv are read as produced on 2026-09-13"),
        "hypothesis_source": ("reports/PHASE2D_CLAIMS_AFTER_MATERIAL_DROP.md section 1: Glu half "
                              "retains 0.829, Lys half 0.435, K/R version 0.430 on one cohort"),
        "instrument": {"replicates": REPLICATES, "seed": SEED,
                       "grouping": "homology component from the oof table's component column",
                       "source": "run_phase2_claims_under_detectability_control (unchanged)"},
        "attributes": {"kr_sfe006_offsets": SFE006_KR_OFFSETS, "kr_within_5": WITHIN_5},
        "cohort": {"sites": len(y), "positives": int(y.sum()),
                   "components": len(set(groups)),
                   "dropped_protein_absent": missing_protein,
                   "dropped_site_off_protein_or_not_cys": off_site},
        "label_reference": label_reference,
        "locked_predictions": locked_summary,
        "inputs": {"results/v2_stack_chem_oof.csv": sha256_of(CHEM),
                   "results/v2_stack_full_oof.csv": sha256_of(FULL),
                   "inputs/sly_proteome.tsv": sha256_of(PROTEOME),
                   "collaboration/locked_predictions_2026-09-14/locked_site_predictions.csv":
                       sha256_of(LOCKED)},
        "versions": {"python": sys.version, "numpy": np.__version__},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    json.dump(audit, open(os.path.join(RESULTS, "self_audit_cleavage_flank_audit.json"), "w"),
              ensure_ascii=False, indent=2)
    for r in rows:
        print(f"  {r['model']:44s} {r['stratum']:15s} n={r['n']:>6d} "
              f"rate={r['attribute_rate']} auc={r['auc_score_to_attribute']} "
              f"[{r['ci_low']}, {r['ci_high']}]")
    print("label reference:", json.dumps(label_reference, ensure_ascii=False))
    print("locked:", json.dumps(locked_summary, ensure_ascii=False))
    print(f"done in {audit['elapsed_minutes']} min")


if __name__ == "__main__":
    main()
