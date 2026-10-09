"""Does the peptide-length signature hold for all five proteases, or only trypsin?

`scripts/run_self_audit_cleavage_flank.py` asks whether a deployed model ranks cysteines by
trypsin's cleavage structure beyond what the label carries. A post-hoc offset profile suggested the pattern
is not K/R enrichment but PEPTIDE-LENGTH logic: cleavage residues DEPLETED immediately beside
the cysteine and enriched 6 to 12 residues away - which is what a model learns if it is predicting "this cysteine lands in an
identifiable tryptic peptide".

That profile was post-hoc in the trypsin round. **Here it is the predeclared hypothesis, tested on
the other four proteases**, which is the right way round: one cohort generates, the replication
tests.

**Predeclared, fixed before the run.**

* The five proteases and their cleavage rules are taken UNCHANGED from
  `run_cross_protease_detectability_probe.PROTEASES` - the same definitions behind the project's
  standing conclusion that the detectability bias holds across five proteases. Cut positions come
  from that module's `boundaries_for`, so cleavage side (`after` / `before`) and proline blocking
  are respected rather than approximated by counting residues at offsets.
* Two attributes per protease, from the distance between the cysteine and the nearest cleavage
  boundary:
  - primary `cut_at_6_to_12` - a boundary exists 6 to 12 residues away (the workable-length band);
  - control `cut_within_5` - a boundary exists within 5 residues (the too-short band).
* **Why the proximal band is a control rather than a second test.** Three of these proteases cleave
  at residues with genuine chemical meaning next to a cysteine: GluC at glutamate and AspN at
  aspartate, whose proximity lowers cysteine pKa, and chymotrypsin at F/W/Y/L/M, which report
  hydrophobic packing. A proximal association could therefore be real chemistry. Through-space
  chemistry acts at SHORT range; peptide-length logic acts at 6-12 residues and is DEPLETED at
  short range. The contrast between the two bands is what separates them, and neither band alone
  would.
* Statistics are the registered instrument, unchanged: weighted AUC of score against attribute,
  cluster bootstrap over homology components, 5000 replicates, seed 20260915, plus top-100
  enrichment as a log2 odds ratio. The AUC is also computed WITHIN each label stratum, where the
  label is constant and an association therefore cannot be inherited from it - that is the clause
  the trypsin round's verdict rested on.
* The label's own association is reported per attribute, WITH the analytic cap for a binary
  predictor (0.5 + (rate_positive - rate_negative)/2). The trypsin round recorded that comparing a
  continuous score against an uncapped label AUC is unfair; the cap is printed here so the same
  mistake is not repeated.
* Descriptive, no verdict attached: the mean length of the containing peptide in the top 100 versus
  the cohort, per protease.

**Reading rule per protease, fixed before the run.** Using the within-stratum AUCs:
  1. distal interval entirely above 0.5 AND distal point estimate above proximal - the
     peptide-length signature REPLICATES for that protease.
  2. both intervals cover 0.5 - no cleavage-structure association; report as such.
  3. proximal above distal with its interval above 0.5 - OPPOSITE to trypsin. That is the pattern
     genuine short-range chemistry would produce, and it refutes the peptide-length reading for
     that protease.
**Aggregate rule.** Signature in 4 or 5 proteases: the finding upgrades from a trypsin-specific
artefact to a structural bias of the capture workflow, and belongs in the manuscript's main text.
In trypsin only: it stays trypsin-specific. In 2 or 3: partial, and the report must name which
proteases and not average over them.

Writes results/self_audit_five_proteases.csv, _topk.csv and _audit.json.
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

from run_cross_protease_detectability_probe import PROTEASES, boundaries_for
from run_phase2_claims_under_detectability_control import (
    REPLICATES, SEED, bootstrap_interval, log_odds_ratio, make_weighted_auc,
    make_weighted_table, sha256_of,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
CHEM = os.path.join(RESULTS, "v2_stack_chem_oof.csv")
FULL = os.path.join(RESULTS, "v2_stack_full_oof.csv")
PROTEOME = os.path.join(ROOT, "inputs", "cohort_proteome.tsv")
SCRIPT = os.path.abspath(__file__)
DISTAL = (6, 12)
PROXIMAL = 5
TOP_K = 100


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


def cut_distances(sequence, rule):
    """Interior cleavage boundaries as 0-based sequence offsets."""
    cuts = boundaries_for(sequence, rule)
    return np.asarray([c for c in cuts if 0 < c < len(sequence)], dtype=int)


def peptide_length_at(cuts_full, position_index):
    """Length of the peptide containing the 0-based residue index."""
    idx = np.searchsorted(cuts_full, position_index, side="right")
    return int(cuts_full[idx] - cuts_full[idx - 1])


def band_flags(sequence, position, rule):
    """(distal 6-12, proximal <=5, containing peptide length) for a 1-based position."""
    i = position - 1
    interior = cut_distances(sequence, rule)
    cuts_full = np.concatenate([[0], interior, [len(sequence)]])
    length = peptide_length_at(cuts_full, i)
    if interior.size == 0:
        return 0, 0, length
    d = np.abs(interior - i)
    return int(((d >= DISTAL[0]) & (d <= DISTAL[1])).any()), int((d <= PROXIMAL).any()), length


def auc_rows(model, protease, band, y, attribute, score, groups, out):
    for stratum, mask in (("all", np.ones(len(y), dtype=bool)),
                          ("label_positive", y == 1), ("label_negative", y == 0)):
        sub_attr, sub_score, sub_groups = attribute[mask], score[mask], groups[mask]
        if mask.sum() < 20 or len(np.unique(sub_attr)) < 2:
            out.append({"model": model, "protease": protease, "band": band, "stratum": stratum,
                        "n": int(mask.sum()), "attribute_rate": None,
                        "auc": None, "ci_low": None, "ci_high": None})
            continue
        fn = make_weighted_auc(sub_attr, sub_score)
        point = fn(np.ones(int(mask.sum())))
        interval, _ = bootstrap_interval(fn, sub_groups)
        out.append({"model": model, "protease": protease, "band": band, "stratum": stratum,
                    "n": int(mask.sum()),
                    "attribute_rate": round(float(sub_attr.mean()), 4),
                    "auc": round(float(point), 4),
                    "ci_low": round(interval[0], 4), "ci_high": round(interval[1], 4)})


def main():
    started = time.time()
    seqs = sequences()
    chem = read_csv(CHEM)
    full = read_csv(FULL)
    assert len(chem) == len(full)
    y = np.asarray([int(r["label"]) for r in chem], dtype=int)
    groups = np.asarray([r["component"] for r in chem])
    scores = {"v2_chem_blend": np.asarray([float(r["blend"]) for r in chem], dtype=float),
              "v2_full_blend_benchmark": np.asarray([float(r["blend"]) for r in full], dtype=float)}

    rows, topk_rows, label_ref, descriptive = [], [], {}, {}
    for protease, rule in PROTEASES.items():
        distal, proximal, peplen = [], [], []
        for r in chem:
            seq = seqs[r["accession"]]
            a, b, L = band_flags(seq, int(r["position"]), rule)
            distal.append(a); proximal.append(b); peplen.append(L)
        distal = np.asarray(distal, dtype=int)
        proximal = np.asarray(proximal, dtype=int)
        peplen = np.asarray(peplen, dtype=float)

        for band, attr in (("distal_6_12", distal), ("proximal_le_5", proximal)):
            fn = make_weighted_auc(attr, y.astype(float))
            point = fn(np.ones(len(y)))
            iv, _ = bootstrap_interval(fn, groups)
            cap = 0.5 + (attr[y == 1].mean() - attr[y == 0].mean()) / 2
            label_ref[f"{protease}|{band}"] = {
                "auc_label_to_attribute": round(float(point), 4),
                "ci": [round(iv[0], 4), round(iv[1], 4)],
                "binary_predictor_analytic_cap": round(float(cap), 4),
                "rate_in_positives": round(float(attr[y == 1].mean()), 4),
                "rate_in_negatives": round(float(attr[y == 0].mean()), 4)}
            for model, score in scores.items():
                auc_rows(model, protease, band, y, attr, score, groups, rows)
                order = np.argsort(-score, kind="stable")
                in_top = np.zeros(len(score), dtype=int)
                in_top[order[:TOP_K]] = 1
                table = make_weighted_table(in_top, attr)
                lor = log_odds_ratio(*table(np.ones(len(score))))
                liv, _ = bootstrap_interval(lambda w: log_odds_ratio(*table(w)), groups)
                topk_rows.append({
                    "model": model, "protease": protease, "band": band, "k": TOP_K,
                    "rate_in_top_k": round(float(attr[order[:TOP_K]].mean()), 4),
                    "rate_rest": round(float(attr[order[TOP_K:]].mean()), 4),
                    "log2_or": round(float(lor), 4),
                    "ci_low": round(liv[0], 4), "ci_high": round(liv[1], 4)})
        descriptive[protease] = {
            "mean_peptide_length_cohort": round(float(peplen.mean()), 2),
            "mean_peptide_length_top100_chem": round(
                float(peplen[np.argsort(-scores["v2_chem_blend"], kind="stable")[:TOP_K]].mean()), 2),
            "mean_peptide_length_top100_full": round(
                float(peplen[np.argsort(-scores["v2_full_blend_benchmark"], kind="stable")[:TOP_K]].mean()), 2),
            "distal_rate_cohort": round(float(distal.mean()), 4),
            "proximal_rate_cohort": round(float(proximal.mean()), 4)}
        print(f"  {protease} done ({round(time.time()-started, 1)}s)", flush=True)

    # per-protease verdict by the predeclared rule, on the within-stratum AUCs
    R = {(r["model"], r["protease"], r["band"], r["stratum"]): r for r in rows}
    verdicts = {}
    for protease in PROTEASES:
        per_model = {}
        for model in scores:
            calls = []
            for stratum in ("label_positive", "label_negative"):
                d = R[(model, protease, "distal_6_12", stratum)]
                p = R[(model, protease, "proximal_le_5", stratum)]
                if d["auc"] is None or p["auc"] is None:
                    calls.append("not_measurable"); continue
                distal_above = d["ci_low"] > 0.5
                proximal_above = p["ci_low"] > 0.5
                if distal_above and d["auc"] > p["auc"]:
                    calls.append("signature_replicates")
                elif proximal_above and p["auc"] > d["auc"]:
                    calls.append("opposite_short_range")
                elif not distal_above and not proximal_above:
                    calls.append("no_association")
                else:
                    calls.append("mixed")
            per_model[model] = calls
        chem_calls = per_model["v2_chem_blend"]
        verdicts[protease] = {
            "per_model": per_model,
            "chem_model": ("signature_replicates" if all(c == "signature_replicates" for c in chem_calls)
                             else ("opposite_short_range" if all(c == "opposite_short_range" for c in chem_calls)
                                   else ("no_association" if all(c == "no_association" for c in chem_calls)
                                         else "mixed_across_label_strata")))}
    replicating = [p for p, v in verdicts.items() if v["chem_model"] == "signature_replicates"]
    aggregate = ("structural_bias_of_the_capture_workflow" if len(replicating) >= 4
                 else ("trypsin_specific" if replicating == ["Trypsin"]
                       else f"partial_in_{len(replicating)}_of_5"))

    with open(os.path.join(RESULTS, "self_audit_five_proteases.csv"), "w",
              encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(os.path.join(RESULTS, "self_audit_five_proteases_topk.csv"), "w",
              encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(topk_rows[0])); w.writeheader(); w.writerows(topk_rows)

    audit = {
        "script": "scripts/run_self_audit_five_proteases.py", "script_sha256": sha256_of(SCRIPT),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.node(),
        "round": "self-audit part 2: the peptide-length signature across all five proteases",
        "hypothesis_status": ("post-hoc in the trypsin round "
                              "(results/self_audit_offset_profile_posthoc.json), PREDECLARED here"),
        "no_refit": "stored out-of-fold scores from 2026-09-13 only",
        "protease_rules": {k: dict(v) for k, v in PROTEASES.items()},
        "protease_rules_source": "run_cross_protease_detectability_probe.PROTEASES (unchanged)",
        "bands": {"distal_6_12": list(DISTAL), "proximal_le_5": PROXIMAL,
                  "why_proximal_is_a_control": (
                      "GluC cleaves at E and AspN at D, whose proximity genuinely lowers cysteine "
                      "pKa, and chymotrypsin at F/W/Y/L/M reports hydrophobic packing; through-space "
                      "chemistry acts at short range while peptide-length logic acts at 6-12 and is "
                      "depleted at short range, so only the contrast separates them")},
        "instrument": {"replicates": REPLICATES, "seed": SEED, "top_k": TOP_K,
                       "grouping": "homology component",
                       "source": "run_phase2_claims_under_detectability_control (unchanged)"},
        "cohort": {"sites": len(y), "positives": int(y.sum()), "components": len(set(groups))},
        "label_reference": label_ref,
        "descriptive_peptide_length": descriptive,
        "verdicts_per_protease": verdicts,
        "proteases_where_signature_replicates": replicating,
        "aggregate_verdict": aggregate,
        "inputs": {"results/v2_stack_chem_oof.csv": sha256_of(CHEM),
                   "results/v2_stack_full_oof.csv": sha256_of(FULL),
                   "inputs/cohort_proteome.tsv": sha256_of(PROTEOME)},
        "versions": {"python": sys.version, "numpy": np.__version__},
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }
    json.dump(audit, open(os.path.join(RESULTS, "self_audit_five_proteases_audit.json"), "w"),
              ensure_ascii=False, indent=2)
    for protease, v in verdicts.items():
        print(f"  {protease:14s} chem={v['chem_model']:26s} calls={v['per_model']['v2_chem_blend']}")
    print("replicates in:", replicating)
    print("AGGREGATE:", aggregate)
    print("peptide length (cohort / top100 chem):",
          {p: (d["mean_peptide_length_cohort"], d["mean_peptide_length_top100_chem"])
           for p, d in descriptive.items()})
    print(f"done in {audit['elapsed_minutes']} min")


if __name__ == "__main__":
    main()
