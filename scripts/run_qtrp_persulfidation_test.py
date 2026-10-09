"""The test the project has been blocked on: real persulfidation, matched negatives.

Low-pH QTRP labels free thiols and persulfides in the same run, and the adduct on
the cysteine separates them (`C(IPM)` against `C(IPMS)`). The ingestion
(`results/qtrp_class_pairing.csv`) shows every persulfide site is also seen as a
free thiol and no site is persulfide-only, which is what partial modification of
a detected cysteine should look like. So the design is:

* positive - a cysteine whose persulfide form was detected in this run;
* negative - a cysteine detected in the same run as a free thiol only.

Both classes were seen by the same probe in the same experiment, so neither is
defined by being invisible. This is the first time in this project that a
persulfidation label comes with a defensible negative set.

Four feature sets are compared, all evaluated with protein-grouped folds:

* struct53      - the AlphaFold descriptors that reached Spearman 0.4422 on the
                  ABPP reactivity label and did nothing on the annotation label;
* probe_vis10   - the ten pure-visibility features;
* probe_digest25- the full in-silico digest block;
* combined      - struct53 plus the digest block.

If the structural block discriminates here while the visibility block does not,
the project finally has evidence of persulfidation chemistry that is not an
artefact of what mass spectrometry can see.
"""
from __future__ import annotations

import collections
import gzip
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from common import RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import build_structural_chemistry_features as builder
from run_cross_protease_detectability_probe import (
    FEATURE_NAMES, HGB_KWARGS, PROTEASES, VISIBILITY_ONLY, boundaries_for, features,
)

SITES = RESULTS / "qtrp_sites_normalised.csv"
PDB_DIR = ROOT / "external/alphafold_human/pdb"
PROTEOME = ROOT / "external/proteomes/hsa.fasta.gz"
N_SPLITS = 5
RULE = "Trypsin"
ARMS = (
    ("S1", "5.0", "A549 lysate, low-pH QTRP, primary arm"),
    ("S2", "5.0", "A549 lysate, reactivity profiling, replication arm"),
    ("S1", "7.6", "A549 lysate at pH 7.6, contrast arm where the method loses sensitivity"),
)


def load_proteome():
    sequences, name, parts = {}, None, []
    with gzip.open(PROTEOME, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(parts)
                header = line[1:]
                name = header.split("|")[1] if "|" in header else header.split()[0]
                parts = []
            else:
                parts.append(line)
    if name:
        sequences[name] = "".join(parts)
    return {k: v.upper() for k, v in sequences.items()}


def arm_sites(rows, table, ph):
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
    return sorted(positives), sorted(negatives)


def structural_matrix(keys, sequences):
    by_protein = collections.defaultdict(list)
    for index, (accession, position) in enumerate(keys):
        by_protein[accession].append((index, position))
    names, matrix = None, None
    stats = {"parsed": 0, "no_model": 0, "no_sulfur": 0}
    for accession, entries in sorted(by_protein.items()):
        path = PDB_DIR / f"AF-{accession}-F1-model_v6.pdb"
        structure = builder.read_pdb(path) if path.exists() else None
        if structure is None:
            stats["no_model"] += 1
            continue
        stats["parsed"] += 1
        for index, position in entries:
            values = builder.site_features(structure, position)
            if values is None:
                stats["no_sulfur"] += 1
                continue
            if names is None:
                names = sorted(values)
                matrix = np.full((len(keys), len(names)), np.nan, dtype=np.float32)
            matrix[index] = [values[name] for name in names]
    return matrix, names, stats


def probe_matrix(keys, sequences):
    rule = PROTEASES[RULE]
    cache, rows, ok = {}, [], []
    for accession, position in keys:
        sequence = sequences.get(accession)
        if not sequence or not (1 <= position <= len(sequence)) or sequence[position - 1] != "C":
            rows.append([np.nan] * len(FEATURE_NAMES))
            ok.append(False)
            continue
        if accession not in cache:
            cache[accession] = boundaries_for(sequence, rule)
        rows.append(features(sequence, position - 1, cache[accession], rule))
        ok.append(True)
    return np.asarray(rows, dtype=np.float64), np.asarray(ok)


def grouped_scores(x, y, groups):
    scores = np.full(len(y), np.nan)
    for train, test in GroupKFold(n_splits=N_SPLITS).split(x, y, groups):
        model = HistGradientBoostingClassifier(**HGB_KWARGS)
        model.fit(x[train], y[train])
        scores[test] = model.predict_proba(x[test])[:, 1]
    return scores


def main():
    started = time.time()
    rows = read_csv(SITES)
    sequences = load_proteome()
    visibility_columns = [i for i, n in enumerate(FEATURE_NAMES) if n in VISIBILITY_ONLY]

    summary = []
    audit_arms = {}
    for table, ph, description in ARMS:
        positives, negatives = arm_sites(rows, table, ph)
        keys = positives + negatives
        y = np.asarray([1] * len(positives) + [0] * len(negatives))
        groups = np.asarray([a for a, _ in keys])
        def is_cysteine(accession, position):
            sequence = sequences.get(accession)
            if not sequence or not (1 <= position <= len(sequence)):
                return False
            return sequence[position - 1] == "C"

        residue_ok = np.asarray(
            [is_cysteine(a, p) for a, p in keys], dtype=bool
        )
        S, S_names, stats = structural_matrix(keys, sequences)
        P, probe_ok = probe_matrix(keys, sequences)
        covered = np.isfinite(S).any(axis=1) & residue_ok & probe_ok
        print(json.dumps({
            "arm": f"{table}_pH{ph}", "positives": int(len(positives)),
            "negatives": int(len(negatives)), "usable_sites": int(covered.sum()),
            **stats, "minutes": round((time.time() - started) / 60, 1),
        }), flush=True)
        if covered.sum() < 100 or y[covered].sum() < 20:
            audit_arms[f"{table}_pH{ph}"] = {"skipped": "too few usable sites", **stats}
            continue

        designs = {
            "struct53": S[covered],
            "probe_visibility10": P[np.ix_(covered, visibility_columns)],
            "probe_digest25": P[covered],
            "struct53_plus_digest25": np.hstack([S[covered], P[covered]]),
        }
        yy, gg = y[covered], groups[covered]
        for label, x in designs.items():
            scores = grouped_scores(x, yy, gg)
            finite = np.isfinite(scores)
            row = {
                "arm": f"{table}_pH{ph}",
                "arm_description": description,
                "design": label,
                "n_features": int(x.shape[1]),
                "n_sites": int(finite.sum()),
                "n_positive": int(yy[finite].sum()),
                "positive_rate": round(float(yy[finite].mean()), 4),
                "n_proteins": int(len(set(gg[finite]))),
                "roc_auc": round(float(roc_auc_score(yy[finite], scores[finite])), 4),
                "average_precision": round(
                    float(average_precision_score(yy[finite], scores[finite])), 4),
            }
            summary.append(row)
            print(f"  {label:24s} AUC {row['roc_auc']:.4f} AP {row['average_precision']:.4f}", flush=True)
        audit_arms[f"{table}_pH{ph}"] = {
            "description": description,
            "positives": int(len(positives)),
            "negatives": int(len(negatives)),
            "usable_sites": int(covered.sum()),
            "structure_stats": stats,
        }

    write_csv(RESULTS / "qtrp_persulfidation_test.csv", summary)
    write_json(RESULTS / "qtrp_persulfidation_test_audit.json", {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "question": "among cysteines detected in the same run, does anything predict which ones carry a persulfide",
        "label": "positive = persulfide form detected; negative = detected as free thiol only, same run",
        "why_this_negative_set_is_defensible": "both classes were detected by the same probe in the same experiment, so neither class is defined by invisibility",
        "source": "Fu et al., Antioxid Redox Signal 2020, doi 10.1089/ars.2019.7777, Supplementary Tables via institutional subscription",
        "splits": f"GroupKFold by protein, {N_SPLITS} folds",
        "classifier": HGB_KWARGS,
        "arms": audit_arms,
        "results": summary,
        "limits": [
            "positives are modest in number, a few hundred at best",
            "a persulfide form can go undetected for abundance reasons, so negatives are 'no persulfide seen here' rather than chemically unmodified",
            "AlphaFold single conformer; the experiment is a lysate at pH 5.0",
            "NaHS-treated lysate is a strong exogenous sulfide load, not physiological persulfidation",
        ],
        "input_hashes": {"results/qtrp_sites_normalised.csv": sha256(SITES)},
        "versions": {"python": sys.version, "numpy": np.__version__},
    })
    print(json.dumps(summary, indent=2, default=str)[:1200])


if __name__ == "__main__":
    main()
