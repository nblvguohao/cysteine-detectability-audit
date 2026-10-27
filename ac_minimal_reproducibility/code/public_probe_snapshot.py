"""Does a peptide-detectability probe alone predict PTM site labels?

This generalises the v2 finding beyond persulfidation and beyond one dataset.
On an annotation benchmark, 25 in-silico tryptic digest features alone ranked
the annotation label well.  If that is a
property of proteomics-derived PTM labels rather than of one cohort, the same
probe should work on somebody else's benchmark, for a different PTM, and it
should stop working once negatives are restricted to peptides the experiment
actually observed.

Substrate: the phosphorylation benchmarks released by Zuallaert et al. with
PhosphoLingo (doi 10.1101/2022.11.28.518163).  They are the only published
PTM-site datasets that ship the same modification split five ways by protease,
each in two variants:

* NPF, negatives not filtered: every S/T in the protein is a candidate.
* PF, negatives protease-filtered: candidates are restricted to peptides
  matched in the experiment, which is the fix those authors proposed.

Three predeclared probes per dataset:

* matched probe - all digest features computed with that dataset's own protease.
* visibility-only probe - the same protease, but restricted to features that
  describe whether the peptide can be seen at all: length, mass, GRAVY and the
  detectability windows.  Position inside the peptide and the local count of
  cleavage residues are dropped, because a K/R count is also a real feature of
  basophilic kinase motifs and would let motif biology masquerade as
  detectability.
* trypsin probe - all features computed with the trypsin rule regardless of the
  dataset's protease.

Two predeclared readings:

1. NPF minus PF for the matched probe.  This is how much of the benchmark's
   apparent signal is peptide visibility rather than modification chemistry.
2. matched minus trypsin on the non-trypsin datasets.  If the mechanism really
   is digestion, a probe built on the wrong enzyme should lose most of its
   edge.  This is a specificity control, not a performance claim.

No tuning: one fixed HistGradientBoosting configuration, copied from the v1
baseline in inputs/source_configuration.json, fitted on the released train
split and evaluated on the released test split.  Cleavage rules are stated
explicitly in PROTEASES below and are assumptions of this analysis.

Usage:
    python scripts/run_cross_protease_detectability_probe.py
Inputs are staged in external/phospholingo/, downloaded from
https://raw.githubusercontent.com/jasperzuallaert/PhosphoLingo/master/data/
"""
from __future__ import annotations

import bisect
import csv
import hashlib
import json
import pathlib

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "external/phospholingo"
OUT_SUMMARY = ROOT / "external/cross_protease_probe_summary.csv"
OUT_AUDIT = ROOT / "external/cross_protease_probe_audit.json"

# Cleavage rules.  "after" cuts C-terminal to the listed residues, "before"
# cuts N-terminal to them.  block_proline suppresses a cut before proline.
PROTEASES = {
    "Trypsin": {"residues": "KR", "side": "after", "block_proline": True},
    "LysC": {"residues": "K", "side": "after", "block_proline": True},
    "GluC": {"residues": "E", "side": "after", "block_proline": False},
    "AspN": {"residues": "D", "side": "before", "block_proline": False},
    "Chymotrypsin": {"residues": "FWYLM", "side": "after", "block_proline": True},
}
DATASETS = list(PROTEASES) + ["multi_protease"]
VARIANTS = ("NPF", "PF")

RESIDUE_MASS = {
    "A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
    "F": 147.06841, "G": 57.02146, "H": 137.05891, "I": 113.08406,
    "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
    "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203,
    "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y": 163.06333,
    "X": 110.0,
}
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2, "X": 0.0,
}
WATER_MASS = 18.010565
DETECTABLE_LENGTH = (7, 30)
DETECTABLE_MASS = (700.0, 3500.0)

HGB_KWARGS = dict(
    max_iter=300, learning_rate=0.06, max_leaf_nodes=31, min_samples_leaf=20,
    l2_regularization=0.0, max_bins=255, early_stopping=False, random_state=0,
)

VISIBILITY_ONLY = {
    "pep_len", "pep_mass", "pep_gravy", "pep_log_len",
    "pep_detectable_length", "pep_detectable_mass", "pep_detectable_both",
    "pep_mc1_detectable", "pep_mc2_detectable",
    "pep_detectable_any_missed_cleavage",
}

FEATURE_NAMES = [
    "pep_len", "pep_mass", "pep_gravy", "pep_cys_count",
    "pep_offset_in_peptide", "pep_relative_offset",
    "pep_is_protein_n_terminal", "pep_is_protein_c_terminal",
    "pep_detectable_length", "pep_detectable_mass", "pep_detectable_both",
    "pep_log_len",
    "pep_mc1_left_len", "pep_mc1_right_len", "pep_mc1_min_len", "pep_mc1_detectable",
    "pep_mc2_left_len", "pep_mc2_right_len", "pep_mc2_min_len", "pep_mc2_detectable",
    "pep_detectable_any_missed_cleavage",
    "pep_n_peptides_in_protein", "pep_cleavage_sites_within_20",
    "pep_dist_to_peptide_n_term", "pep_dist_to_peptide_c_term",
]


def boundaries_for(sequence: str, rule: dict) -> list:
    n = len(sequence)
    cuts = [0]
    residues, side, block = rule["residues"], rule["side"], rule["block_proline"]
    for i, residue in enumerate(sequence):
        if residue not in residues:
            continue
        cut = i + 1 if side == "after" else i
        if cut <= 0 or cut >= n:
            continue
        if block and sequence[cut] == "P":
            continue
        cuts.append(cut)
    if cuts[-1] != n:
        cuts.append(n)
    return sorted(set(cuts))


def descriptors(sequence: str, start: int, stop: int):
    peptide = sequence[start:stop]
    if not peptide:
        return 0.0, 0.0, 0.0, 0.0
    mass = WATER_MASS + sum(RESIDUE_MASS.get(r, RESIDUE_MASS["X"]) for r in peptide)
    gravy = sum(HYDROPATHY.get(r, 0.0) for r in peptide) / len(peptide)
    return float(len(peptide)), mass, gravy, float(peptide.count("C"))


def features(sequence: str, index: int, bounds: list, rule: dict) -> list:
    n = len(sequence)
    k = min(max(bisect.bisect_right(bounds, index) - 1, 0), len(bounds) - 2)
    start, stop = bounds[k], bounds[k + 1]
    length, mass, gravy, cys = descriptors(sequence, start, stop)

    det_len = 1.0 if DETECTABLE_LENGTH[0] <= length <= DETECTABLE_LENGTH[1] else 0.0
    det_mass = 1.0 if DETECTABLE_MASS[0] <= mass <= DETECTABLE_MASS[1] else 0.0
    row = [
        length, mass, gravy, cys,
        float(index - start), float(index - start) / max(1.0, length - 1),
        1.0 if start == 0 else 0.0, 1.0 if stop == n else 0.0,
        det_len, det_mass, 1.0 if (det_len and det_mass) else 0.0,
        float(np.log1p(length)),
    ]
    detectable_any = det_len
    for missed in (1, 2):
        left = bounds[max(0, k - missed)]
        right = bounds[min(len(bounds) - 1, k + 1 + missed)]
        left_len, left_mass, _, _ = descriptors(sequence, left, stop)
        right_len, right_mass, _, _ = descriptors(sequence, start, right)
        detectable = 0.0
        for cand_len, cand_mass in ((left_len, left_mass), (right_len, right_mass)):
            if (DETECTABLE_LENGTH[0] <= cand_len <= DETECTABLE_LENGTH[1]
                    and DETECTABLE_MASS[0] <= cand_mass <= DETECTABLE_MASS[1]):
                detectable = 1.0
        row += [left_len, right_len, min(left_len, right_len), detectable]
        detectable_any = max(detectable_any, detectable)
    window = sequence[max(0, index - 20):min(n, index + 21)]
    row += [
        detectable_any,
        float(len(bounds) - 1),
        float(sum(1 for r in window if r in rule["residues"])),
        float(index - start), float(stop - 1 - index),
    ]
    return row


def read_fasta(path: pathlib.Path):
    entries = []
    name, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                entries.append((name, "".join(chunks)))
            name, chunks = line[1:].strip(), []
        elif line:
            chunks.append(line.strip())
    if name is not None:
        entries.append((name, "".join(chunks)))

    records = []
    for name, annotated in entries:
        sequence, sites = [], []
        for char in annotated:
            if char == "#":
                sites.append((len(sequence) - 1, 1))
            elif char == "@":
                sites.append((len(sequence) - 1, 0))
            else:
                sequence.append(char)
        records.append((name, "".join(sequence), sites))
    return records


def build(records, rule_name: str, subset: str = "all"):
    rule = PROTEASES[rule_name]
    keep = [i for i, name in enumerate(FEATURE_NAMES)
            if subset == "all" or name in VISIBILITY_ONLY]
    rows, labels = [], []
    for _, sequence, sites in records:
        bounds = boundaries_for(sequence, rule)
        for index, label in sites:
            row = features(sequence, index, bounds, rule)
            rows.append([row[i] for i in keep])
            labels.append(label)
    return np.asarray(rows, dtype=np.float64), np.asarray(labels, dtype=np.int8)


def auc_interval(y_true, scores, replicates=400, seed=0):
    rng = np.random.default_rng(seed)
    n = len(y_true)
    draws = []
    for _ in range(replicates):
        idx = rng.integers(0, n, n)
        if y_true[idx].min() == y_true[idx].max():
            continue
        draws.append(roc_auc_score(y_true[idx], scores[idx]))
    draws.sort()
    lo = draws[int(0.025 * (len(draws) - 1))]
    hi = draws[int(np.ceil(0.975 * (len(draws) - 1)))]
    return float(lo), float(hi)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    summary, manifest = [], {}
    cache = {}
    for dataset in DATASETS:
        folder = "single_protease" if dataset in PROTEASES else "."
        for variant in VARIANTS:
            base = DATA / folder / dataset / variant if folder != "." else DATA / dataset / variant
            train_path, test_path = base / "train.fasta", base / "test.fasta"
            if not train_path.exists():
                print("missing", train_path)
                continue
            manifest[str(train_path.relative_to(ROOT))] = sha256(train_path)
            manifest[str(test_path.relative_to(ROOT))] = sha256(test_path)
            train_records = read_fasta(train_path)
            test_records = read_fasta(test_path)
            cache[(dataset, variant)] = (train_records, test_records)

            own = "Trypsin" if dataset == "multi_protease" else dataset
            probes = [(own, "all"), (own, "visibility_only")]
            if dataset != "Trypsin":
                probes.append(("Trypsin", "all"))
            if dataset == "multi_protease":
                probes = [("Trypsin", "all"), ("Trypsin", "visibility_only"), ("GluC", "all")]

            for probe, subset in probes:
                x_train, y_train = build(train_records, probe, subset)
                x_test, y_test = build(test_records, probe, subset)
                model = HistGradientBoostingClassifier(**HGB_KWARGS)
                model.fit(x_train, y_train)
                scores = model.predict_proba(x_test)[:, 1]
                lo, hi = auc_interval(y_test, scores)
                summary.append({
                    "dataset": dataset,
                    "negatives": variant,
                    "probe_protease": probe,
                    "feature_subset": subset,
                    "probe_is_matched": "yes" if probe == dataset or (dataset == "multi_protease" and probe == "Trypsin") else "no",
                    "n_train_sites": int(len(y_train)),
                    "n_train_pos": int(y_train.sum()),
                    "n_test_sites": int(len(y_test)),
                    "n_test_pos": int(y_test.sum()),
                    "test_positive_rate": round(float(y_test.mean()), 4),
                    "roc_auc": round(float(roc_auc_score(y_test, scores)), 4),
                    "roc_auc_ci": f"[{lo:.4f}, {hi:.4f}]",
                    "average_precision": round(float(average_precision_score(y_test, scores)), 4),
                    "n_features": int(x_train.shape[1]),
                })
                print(f"{dataset:<14}{variant:<4}{probe:<13}{subset:<17}"
                      f"n={len(y_test):>6} pos={y_test.mean():.3f} "
                      f"AUC={summary[-1]['roc_auc']:.4f} [{lo:.4f}, {hi:.4f}]")

    with OUT_SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    OUT_AUDIT.write_text(json.dumps({
        "substrate": "PhosphoLingo phosphorylation benchmarks, Zuallaert et al. doi 10.1101/2022.11.28.518163",
        "source_url": "https://raw.githubusercontent.com/jasperzuallaert/PhosphoLingo/master/data/",
        "label_convention": "'#' after a residue is a positive site, '@' a negative candidate",
        "variants": {"NPF": "negatives not filtered", "PF": "negatives restricted to observed peptides"},
        "feature_block": FEATURE_NAMES,
        "n_features": len(FEATURE_NAMES),
        "visibility_only_subset": sorted(VISIBILITY_ONLY),
        "visibility_only_rationale": "drops position-in-peptide and local cleavage-residue counts, which also encode basophilic kinase motifs",
        "auc_interval": {"method": "bootstrap over test sites", "replicates": 400, "seed": 0},
        "cleavage_rules": PROTEASES,
        "cleavage_rule_caveat": "GluC is treated as cleaving after E only; LysC and chymotrypsin block cleavage before proline; these are assumptions of this analysis",
        "classifier": {"name": "HistGradientBoostingClassifier", **HGB_KWARGS},
        "tuning": "none; single fixed configuration copied from the v1 baseline",
        "split": "the released train.fasta and test.fasta of each dataset, unchanged",
        "input_sha256": manifest,
        "interpretation_limit": "probe-only performance measures how much of a benchmark's signal is peptide visibility; it is not a phosphorylation predictor",
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
