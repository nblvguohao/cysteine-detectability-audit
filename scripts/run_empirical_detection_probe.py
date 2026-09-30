"""Can the visibility features predict whether a cysteine was actually detected?

Everything in this project's detectability line so far shows that visibility
features predict an *annotation* label.  The mechanism behind that claim is one
step further back: those features should predict whether mass spectrometry sees
the cysteine at all.  That step has never been measured directly here, because
no dataset in the delivery records which cysteines were detected and which were
not.

PXD015307 supplies it.  Its Fig-1D search result is a mouse iodoTMT6plex thiol
labelling run: 6,746 iodoTMT assignments, every one on a cysteine, no
persulfide channel involved.  A cysteine carrying iodoTMT was seen; the other
cysteines of the same protein were not.  The comparison is therefore
"within a protein that was observed, which of its cysteines was observed",
the same within-protein framing as the project's own metric, with the label
being detection itself rather than a modification.

Trypsin is confirmed two ways: the deposited protocol states MS-grade trypsin,
and 98.04% of the 25,014 identified peptides end in K or R.
"""
from __future__ import annotations

import collections
import json
import platform
import re
import sqlite3
import sys
import time

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from common import RESULTS, ROOT, sha256, write_csv, write_json
from run_cross_protease_detectability_probe import (
    FEATURE_NAMES, HGB_KWARGS, PROTEASES, VISIBILITY_ONLY, boundaries_for, features,
)

MSF = ROOT / "external/sulfhydrome_bta_tmt/Fig-1D.msf"
RULE = "Trypsin"
MODIFICATION = "iodoTMT6plex"
N_SPLITS = 5
ACCESSION_PATTERN = re.compile(r"(?:sp|tr)\|([A-Z0-9\-]+)\|")


def accession_of(description):
    if not description:
        return None
    match = ACCESSION_PATTERN.search(description)
    if match:
        return match.group(1)
    first = description.lstrip(">").split()[0]
    return first.split("|")[1] if "|" in first else first


def load_cohort():
    connection = sqlite3.connect(MSF)
    proteins = {}
    for protein_id, sequence, sequence_hash in connection.execute(
        "select ProteinID, Sequence, SequenceHashCode from Proteins"
    ):
        proteins[protein_id] = {
            "sequence": (sequence or "").upper(),
            "hash": sequence_hash,
            "accession": None,
        }
    for protein_id, description in connection.execute(
        "select ProteinID, Description from ProteinAnnotations"
    ):
        if protein_id in proteins:
            proteins[protein_id]["accession"] = accession_of(description)

    detected = collections.defaultdict(set)
    ambiguous = 0
    considered = 0
    rows = connection.execute(
        """
        select p.Sequence, pm.Position, pp.ProteinID
        from PeptidesAminoAcidModifications pm
        join AminoAcidModifications m
             on m.AminoAcidModificationID = pm.AminoAcidModificationID
        join Peptides p on p.PeptideID = pm.PeptideID
        join PeptidesProteins pp on pp.PeptideID = p.PeptideID
        where m.ModificationName = ?
        """,
        (MODIFICATION,),
    ).fetchall()
    for peptide, position, protein_id in rows:
        considered += 1
        peptide = (peptide or "").upper()
        record = proteins.get(protein_id)
        if record is None or not peptide:
            continue
        if not (0 <= position < len(peptide)) or peptide[position] != "C":
            continue
        sequence = record["sequence"]
        starts = [
            index for index in range(len(sequence))
            if sequence.startswith(peptide, index)
        ]
        if len(starts) != 1:
            # a repeated peptide cannot be localised to one residue
            ambiguous += 1
            continue
        detected[protein_id].add(starts[0] + position + 1)

    # collapse identical protein sequences so shared entries are counted once
    by_hash = collections.defaultdict(list)
    for protein_id, record in proteins.items():
        by_hash[record["hash"]].append(protein_id)
    cohort = []
    for group in by_hash.values():
        union = set()
        for protein_id in group:
            union |= detected.get(protein_id, set())
        if not union:
            continue
        representative = min(
            group, key=lambda pid: (proteins[pid]["accession"] or "zzz", pid)
        )
        sequence = proteins[representative]["sequence"]
        cysteines = [i + 1 for i, residue in enumerate(sequence) if residue == "C"]
        if len(cysteines) < 2:
            continue
        cohort.append({
            "accession": proteins[representative]["accession"] or f"PID{representative}",
            "sequence": sequence,
            "cysteines": cysteines,
            "detected": sorted(union & set(cysteines)),
            "n_duplicate_entries": len(group),
        })
    assignments = connection.execute(
        """
        select count(*) from PeptidesAminoAcidModifications pm
        join AminoAcidModifications m
             on m.AminoAcidModificationID = pm.AminoAcidModificationID
        where m.ModificationName = ?
        """,
        (MODIFICATION,),
    ).fetchone()[0]
    return cohort, {
        "iodotmt_modification_assignments": assignments,
        "iodotmt_peptide_protein_join_rows": considered,
        "ambiguous_repeated_peptide_rows": ambiguous,
        "proteins_in_database": len(proteins),
        "unique_sequences_with_a_detected_cysteine": sum(
            1 for group in by_hash.values()
            if any(detected.get(pid) for pid in group)
        ),
    }


def build(cohort, subset):
    keep = [
        i for i, name in enumerate(FEATURE_NAMES)
        if subset == "all" or name in VISIBILITY_ONLY
    ]
    rule = PROTEASES[RULE]
    rows, labels, groups = [], [], []
    for record in cohort:
        bounds = boundaries_for(record["sequence"], rule)
        detected = set(record["detected"])
        for position in record["cysteines"]:
            row = features(record["sequence"], position - 1, bounds, rule)
            rows.append([row[i] for i in keep])
            labels.append(1 if position in detected else 0)
            groups.append(record["accession"])
    return (
        np.asarray(rows, dtype=np.float64),
        np.asarray(labels, dtype=np.int8),
        np.asarray(groups),
    )


def grouped_oof(x, y, groups):
    scores = np.full(len(y), np.nan)
    for train, test in GroupKFold(n_splits=N_SPLITS).split(x, y, groups):
        model = HistGradientBoostingClassifier(**HGB_KWARGS)
        model.fit(x[train], y[train])
        scores[test] = model.predict_proba(x[test])[:, 1]
    if not np.isfinite(scores).all():
        raise RuntimeError("Incomplete out-of-fold predictions")
    return scores


def within_protein(y, scores, groups):
    order = np.argsort(groups, kind="stable")
    values, proteins, sites = [], 0, 0
    start = 0
    while start < len(order):
        stop = start
        while stop < len(order) and groups[order[stop]] == groups[order[start]]:
            stop += 1
        index = order[start:stop]
        block = y[index]
        if 0 < block.sum() < len(block):
            values.append(roc_auc_score(block, scores[index]))
            proteins += 1
            sites += len(index)
        start = stop
    return float(np.mean(values)), proteins, sites


def main():
    started = time.time()
    if not MSF.exists():
        raise SystemExit(f"Missing {MSF}")
    cohort, counts = load_cohort()
    rows = []
    detail = {}
    for subset in ("all", "visibility_only"):
        x, y, groups = build(cohort, subset)
        scores = grouped_oof(x, y, groups)
        mean_within, n_proteins, n_sites = within_protein(y, scores, groups)
        summary = {
            "label": "cysteine detected with iodoTMT in this run",
            "feature_subset": subset,
            "n_features": int(x.shape[1]),
            "n_sites": int(len(y)),
            "n_detected": int(y.sum()),
            "detection_rate": round(float(y.mean()), 4),
            "n_proteins": int(len(set(groups))),
            "roc_auc": round(float(roc_auc_score(y, scores)), 4),
            "average_precision": round(float(average_precision_score(y, scores)), 4),
            "within_protein_auc": round(mean_within, 4),
            "within_protein_n_proteins": n_proteins,
            "within_protein_n_sites": n_sites,
        }
        rows.append(summary)
        detail[subset] = summary
        print(json.dumps({**summary, "seconds": round(time.time() - started, 1)}), flush=True)

    def share(visibility, full):
        return round((visibility - 0.5) / (full - 0.5), 4) if full > 0.5 else None

    write_csv(RESULTS / "empirical_detection_probe.csv", rows)
    write_json(RESULTS / "empirical_detection_probe_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "question": "do the visibility features predict whether a cysteine was detected at all",
        "deposition": "PXD015307, Fig-1D search result, Quantitative profiling of protein sulfhydrome by the BTA-TMT",
        "species": "Mus musculus",
        "modification_used_as_detection_evidence": MODIFICATION,
        "no_persulfide_channel_in_this_run": True,
        "protease": "Trypsin, stated in the deposited protocol and consistent with 98.04% of peptides ending in K or R",
        "cohort_construction": (
            "identical protein sequences collapsed by hash; a protein enters only if at least one "
            "of its cysteines was detected and it carries at least two cysteines; positives are "
            "iodoTMT-labelled cysteines, negatives are the remaining cysteines of the same protein"
        ),
        "counts": counts,
        "splits": f"GroupKFold by protein, {N_SPLITS} folds",
        "classifier": HGB_KWARGS,
        "results": detail,
        "visibility_share": share(
            detail["visibility_only"]["roc_auc"], detail["all"]["roc_auc"]
        ),
        "visibility_share_within_protein": share(
            detail["visibility_only"]["within_protein_auc"],
            detail["all"]["within_protein_auc"],
        ),
        "limits": [
            "detection is instrument, method and abundance dependent; a single run understates how many cysteines are detectable in principle",
            "negatives are undetected in this run only, not undetectable",
            "repeated peptides that cannot be localised to one residue were dropped",
            "mouse; no persulfidation label is involved anywhere in this test",
        ],
        "source_sha256": sha256(MSF),
        "hardware": {"platform": platform.platform()},
        "versions": {"python": sys.version, "numpy": np.__version__},
    })
    print(json.dumps({
        "visibility_share_global": share(detail["visibility_only"]["roc_auc"], detail["all"]["roc_auc"]),
        "visibility_share_within_protein": share(
            detail["visibility_only"]["within_protein_auc"], detail["all"]["within_protein_auc"]
        ),
        "counts": counts,
    }, indent=2))


if __name__ == "__main__":
    main()
