"""Does the QTRP negative-class contamination move the struct53 result?

`scripts/ingest_frbm_compilation.py` found that "detected as free thiol only"
negatives are not clean: for the two arms with the most positives, 45-46% of
the negatives are called persulfidated somewhere else - most often in another
QTRP arm (different pH or replicate of the same paper), a smaller remainder
only in the external Free Radic Biol Med 2024 compilation. That is label noise
in the negative class, and it should bias every AUC in
`reports/QTRP_PERSULFIDATION_TEST.md` downward, not upward - a positive judged
against contaminated negatives is a harder problem than one judged against
clean negatives, in expectation.

Three negative definitions are compared, all evaluated with the exact struct53
design that produced the reported 0.6164 / 0.6892:

* `raw`         - as reported: negative = free thiol only, this arm;
* `arm_clean`   - negative = free thiol only in this arm AND never called SSH
                  in any other human QTRP arm (removes the within-paper
                  inconsistency, which is most of the flagged mass);
* `full_clean`  - arm_clean AND not in the external compilation's SSH list
                  either (removes the remaining outside-literature overlap).

Positives are unchanged throughout - only sites removed from the negative
class differ. If struct53's discrimination survives or improves as the
negative class gets cleaner, the 0.6164 / 0.6892 floor was conservative.
"""
from __future__ import annotations

import collections
import json
import platform
import sys
import time

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json
from run_qtrp_persulfidation_test import (
    ARMS, SITES, grouped_scores, load_proteome, probe_matrix, structural_matrix,
)
from run_cross_protease_detectability_probe import FEATURE_NAMES, VISIBILITY_ONLY

COMPILED = RESULTS / "frbm_compiled_sites.csv"


def all_arm_labels(rows):
    """Every (accession, site) called SSH in any human QTRP arm, pooled."""
    pooled_ssh = set()
    per_key_arms = collections.defaultdict(set)
    for row in rows:
        if row.get("species") != "human":
            continue
        label = row["class_from_adduct"] or row["class_from_sheet"]
        if label not in ("SH", "SSH"):
            continue
        key = (row["accession"], int(row["site"]))
        arm = f"{row['source_table']}_pH{row['ph']}"
        if label == "SSH":
            pooled_ssh.add(key)
        per_key_arms[key].add((arm, label))
    return pooled_ssh, per_key_arms


def arm_sites_raw(rows, table, ph):
    classes = collections.defaultdict(set)
    for row in rows:
        if row["source_table"] != table or row["ph"] != ph:
            continue
        label = row["class_from_adduct"] or row["class_from_sheet"]
        if label not in ("SH", "SSH"):
            continue
        classes[label].add((row["accession"], int(row["site"])))
    positives = classes["SSH"]
    negatives = classes["SH"] - classes["SSH"]
    return positives, negatives


def main():
    started = time.time()
    rows = read_csv(SITES)
    sequences = load_proteome()
    visibility_columns = [i for i, n in enumerate(FEATURE_NAMES) if n in VISIBILITY_ONLY]

    pooled_ssh, per_key_arms = all_arm_labels(rows)
    compiled_ssh = {
        (r["accession"], int(r["site"]))
        for r in read_csv(COMPILED) if r["compiled_ssh"] == "1"
    }
    print(f"pooled human QTRP SSH (any arm): {len(pooled_ssh)}", flush=True)
    print(f"external compiled SSH: {len(compiled_ssh)}", flush=True)

    summary = []
    audit_arms = {}
    for table, ph, description in ARMS:
        arm_name = f"{table}_pH{ph}"
        positives, negatives_raw = arm_sites_raw(rows, table, ph)
        if not negatives_raw:
            continue

        contaminating_other_arm = negatives_raw & pooled_ssh
        negatives_arm_clean = negatives_raw - pooled_ssh
        contaminating_compilation_only = (negatives_arm_clean & compiled_ssh)
        negatives_full_clean = negatives_arm_clean - compiled_ssh

        variants = {
            "raw": negatives_raw,
            "arm_clean": negatives_arm_clean,
            "full_clean": negatives_full_clean,
        }
        counts = {k: len(v) for k, v in variants.items()}
        print(json.dumps({"arm": arm_name, "positives": len(positives), **{
            f"negatives_{k}": v for k, v in counts.items()}}), flush=True)

        arm_result = {
            "description": description,
            "positives": len(positives),
            "negatives_raw": len(negatives_raw),
            "negatives_flagged_by_other_qtrp_arm": len(contaminating_other_arm),
            "negatives_additionally_flagged_by_external_compilation_only": len(contaminating_compilation_only),
            "negatives_arm_clean": len(negatives_arm_clean),
            "negatives_full_clean": len(negatives_full_clean),
            "designs": {},
        }

        for variant_name, negatives in variants.items():
            keys = sorted(positives) + sorted(negatives)
            y = np.asarray([1] * len(positives) + [0] * len(negatives))
            groups = np.asarray([a for a, _ in keys])

            def is_cysteine(accession, position):
                sequence = sequences.get(accession)
                if not sequence or not (1 <= position <= len(sequence)):
                    return False
                return sequence[position - 1] == "C"

            residue_ok = np.asarray([is_cysteine(a, p) for a, p in keys], dtype=bool)
            S, S_names, stats = structural_matrix(keys, sequences)
            P, probe_ok = probe_matrix(keys, sequences)
            covered = np.isfinite(S).any(axis=1) & residue_ok & probe_ok
            if covered.sum() < 100 or y[covered].sum() < 20 or (y[covered] == 0).sum() < 20:
                arm_result["designs"][variant_name] = {"skipped": "too few usable sites after cleaning"}
                continue

            designs = {
                "struct53": S[covered],
                "probe_visibility10": P[np.ix_(covered, visibility_columns)],
            }
            yy, gg = y[covered], groups[covered]
            variant_result = {"n_sites": int(covered.sum()), "n_positive": int(yy.sum())}
            for label, x in designs.items():
                scores = grouped_scores(x, yy, gg)
                finite = np.isfinite(scores)
                auc = float(roc_auc_score(yy[finite], scores[finite]))
                ap = float(average_precision_score(yy[finite], scores[finite]))
                variant_result[label] = {"roc_auc": round(auc, 4), "average_precision": round(ap, 4)}
                summary.append({
                    "arm": arm_name, "negative_definition": variant_name, "design": label,
                    "n_sites": int(finite.sum()), "n_positive": int(yy[finite].sum()),
                    "roc_auc": round(auc, 4), "average_precision": round(ap, 4),
                })
                print(f"  {variant_name:12s} {label:20s} AUC {auc:.4f} AP {ap:.4f}", flush=True)
            arm_result["designs"][variant_name] = variant_result

        audit_arms[arm_name] = arm_result

    write_csv(RESULTS / "qtrp_negative_cleaning_test.csv", summary)
    write_json(RESULTS / "qtrp_negative_cleaning_test_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": (
            "does struct53's discrimination between persulfidated and "
            "free-thiol-only cysteines survive removing negatives that are "
            "known SSH sites elsewhere (another QTRP arm, or the external "
            "Free Radic Biol Med 2024 compilation)"
        ),
        "negative_definitions": {
            "raw": "as originally reported: detected as free thiol only, this arm",
            "arm_clean": "raw, minus sites called SSH in any other human QTRP arm",
            "full_clean": "arm_clean, minus sites in the external compilation's SSH list",
        },
        "pooled_human_qtrp_ssh_any_arm": len(pooled_ssh),
        "external_compiled_ssh": len(compiled_ssh),
        "arms": audit_arms,
        "results": summary,
        "input_hashes": {
            "results/qtrp_sites_normalised.csv": sha256(SITES),
            "results/frbm_compiled_sites.csv": sha256(COMPILED),
        },
        "versions": {"python": sys.version, "numpy": np.__version__},
    })
    print(json.dumps(audit_arms, indent=2, default=str)[:2000])


if __name__ == "__main__":
    main()
