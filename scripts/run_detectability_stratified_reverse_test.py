"""Retrospective version of the pre-registered reverse test.

The v2 blind-test protocol (protocols/BLIND_TEST_PROTOCOL_V2.md, section 4)
declares one check as the only evidence that could support a chemistry claim:
inside the subset where the Cys-containing peptide IS detectable by a standard
trypsin workflow, does the chemistry model still rank observed sites above
chance?  That check does not need new wet-lab data.  Both external cohorts
already carry full protein sequences, so peptide detectability can be computed
and conditioned on.

Pre-declared design, fixed before looking at any result:

* Candidates are every Cys of every external protein, exactly the set already
  scored in external/v2_retrospective_site_scores.csv.
* Strata.  `all` reproduces the published unstratified numbers.  `detectable`
  keeps Cys whose peptide satisfies the frozen project criterion
  T:pep_detectable_any_missed_cleavage.  `detectable_strict` keeps Cys whose
  fully cleaved peptide satisfies BOTH the length and the mass window, which
  removes a quirk in the frozen criterion (its fully-cleaved term tests length
  only, while its missed-cleavage terms test length and mass).  Both are
  reported so that the quirk cannot drive the conclusion.
  `not_detectable` is the complement of `detectable`.
* A protein enters a stratum when the stratum holds at least two candidates and
  at least one observed site.  Proteins whose ranking problem is trivial inside
  the stratum are therefore excluded, not counted as successes.
* Ranking metrics use the project's tie rule, mean annotation among tied
  top-ranked sites, extended exactly to Top-2 and to MRR by averaging over
  orderings consistent with the ties.
* The random baseline is computed per protein from that protein's own stratum
  size and observed count, so it adapts to the stratum instead of being a
  global constant.
* Uncertainty is a paired protein bootstrap, 400 replicates, seed 0, matching
  the bootstrap unit and replicate count used throughout the project.

Limits carried over unchanged: an unobserved Cys is not a confirmed chemical
negative, both cohorts are retrospective rather than blinded, and the external
deployment is the single LightGBM binary member because no ProtBERT embeddings
exist for these sequences.
"""
from __future__ import annotations

import bisect
import csv
import json
import math
import pathlib
import random

ROOT = pathlib.Path(__file__).resolve().parent.parent
SCORES = ROOT / "external/v2_retrospective_site_scores.csv"
COHORTS = (
    ("human_PXD044043", ROOT / "inputs/PXD044043_mapped_human.json", "hsa_persulfidation"),
    ("rice_PXD072089", ROOT / "inputs/PXD072089_mapped_rice.json", "osa_persulfidation"),
)
OUT_SUMMARY = ROOT / "external/v2_detectability_stratified_summary.csv"
OUT_INTERVALS = ROOT / "external/v2_detectability_stratified_intervals.json"
OUT_SITES = ROOT / "external/v2_detectability_stratified_sites.csv"
OUT_AUDIT = ROOT / "external/v2_detectability_stratified_audit.json"

MODELS = {
    "v2_chem": "v2_chem_score",
    "v2_full": "v2_full_score",
    "v1_corrected_mlp": "mlp_corrected_score",
}
STRATA = ("all", "detectable", "detectable_strict", "not_detectable")
BOOTSTRAP = 400
SEED = 0

RESIDUE_MASS = {
    "A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
    "F": 147.06841, "G": 57.02146, "H": 137.05891, "I": 113.08406,
    "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
    "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203,
    "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y": 163.06333,
    "X": 110.0,
}
WATER_MASS = 18.010565
DETECTABLE_LENGTH = (7, 30)
DETECTABLE_MASS = (700.0, 3500.0)


def cleavage_boundaries(sequence: str) -> list[int]:
    n = len(sequence)
    boundaries = [0]
    for i, residue in enumerate(sequence):
        if residue in ("K", "R") and (i + 1 >= n or sequence[i + 1] != "P"):
            boundaries.append(i + 1)
    if boundaries[-1] != n:
        boundaries.append(n)
    return boundaries


def peptide_length_mass(sequence: str, start: int, stop: int) -> tuple[float, float]:
    peptide = sequence[start:stop]
    if not peptide:
        return 0.0, 0.0
    mass = WATER_MASS + sum(RESIDUE_MASS.get(r, RESIDUE_MASS["X"]) for r in peptide)
    return float(len(peptide)), mass


def detectability(sequence: str, position: int, boundaries: list[int]) -> dict:
    """Replicates the frozen project flags, including their asymmetry."""
    index = position - 1
    peptide_index = bisect.bisect_right(boundaries, index) - 1
    peptide_index = min(max(peptide_index, 0), len(boundaries) - 2)
    start, stop = boundaries[peptide_index], boundaries[peptide_index + 1]
    length, mass = peptide_length_mass(sequence, start, stop)

    detectable_length = DETECTABLE_LENGTH[0] <= length <= DETECTABLE_LENGTH[1]
    detectable_mass = DETECTABLE_MASS[0] <= mass <= DETECTABLE_MASS[1]
    # The frozen feature seeds the running flag with the length test only.
    detectable_any = 1.0 if detectable_length else 0.0
    for missed in (1, 2):
        left = boundaries[max(0, peptide_index - missed)]
        right = boundaries[min(len(boundaries) - 1, peptide_index + 1 + missed)]
        left_length, left_mass = peptide_length_mass(sequence, left, stop)
        right_length, right_mass = peptide_length_mass(sequence, start, right)
        for cand_length, cand_mass in ((left_length, left_mass), (right_length, right_mass)):
            if (DETECTABLE_LENGTH[0] <= cand_length <= DETECTABLE_LENGTH[1]
                    and DETECTABLE_MASS[0] <= cand_mass <= DETECTABLE_MASS[1]):
                detectable_any = 1.0
    return {
        "pep_len": length,
        "pep_mass": mass,
        "pep_cys_count": float(sequence[start:stop].count("C")),
        "detectable_any_missed_cleavage": detectable_any,
        "detectable_strict_fully_cleaved": 1.0 if (detectable_length and detectable_mass) else 0.0,
    }


def groups_by_score(items: list[tuple[float, int]]) -> list[tuple[int, int]]:
    """Descending score groups as (size, n_observed), ties kept together."""
    ordered = sorted(items, key=lambda t: -t[0])
    out = []
    i = 0
    while i < len(ordered):
        j = i
        while j < len(ordered) and ordered[j][0] == ordered[i][0]:
            j += 1
        block = ordered[i:j]
        out.append((len(block), sum(o for _, o in block)))
        i = j
    return out


def top1_tie_averaged(groups) -> float:
    size, observed = groups[0]
    return observed / size


def top2_tie_averaged(groups) -> float:
    size, observed = groups[0]
    if size >= 2:
        misses = size - observed
        return 1.0 - (math.comb(misses, 2) / math.comb(size, 2) if misses >= 2 else 0.0)
    if observed == 1:
        return 1.0
    size2, observed2 = groups[1]
    return observed2 / size2


def expected_reciprocal_rank(groups) -> float:
    """E[1/rank of first observed], averaging over within-tie orderings."""
    offset = 0
    for size, observed in groups:
        if observed == 0:
            offset += size
            continue
        total = math.comb(size, observed)
        value = 0.0
        for r in range(1, size - observed + 2):
            higher = math.comb(size - r + 1, observed)
            lower = math.comb(size - r, observed) if size - r >= observed else 0
            value += ((higher - lower) / total) / (offset + r)
        return value
    return 0.0


def random_metrics(size: int, observed: int) -> tuple[float, float, float]:
    single = [(size, observed)]
    return top1_tie_averaged(single), top2_tie_averaged(single), expected_reciprocal_rank(single)


def bootstrap_interval(values: list[float], seed: int) -> tuple[float, float, float]:
    rng = random.Random(seed)
    n = len(values)
    mean = sum(values) / n
    draws = []
    for _ in range(BOOTSTRAP):
        total = 0.0
        for _ in range(n):
            total += values[rng.randrange(n)]
        draws.append(total / n)
    draws.sort()
    lo = draws[int(0.025 * (BOOTSTRAP - 1))]
    hi = draws[int(math.ceil(0.975 * (BOOTSTRAP - 1)))]
    return mean, lo, hi


def main() -> None:
    with SCORES.open(encoding="utf-8-sig") as handle:
        score_rows = list(csv.DictReader(handle))

    site_records = []
    for dataset, path, key in COHORTS:
        payload = json.loads(path.read_text())[key]
        proteins = {p["acc"]: p for p in payload["proteins"]}
        boundaries = {acc: cleavage_boundaries(p["seq"]) for acc, p in proteins.items()}
        for row in score_rows:
            if row["dataset"] != dataset:
                continue
            protein = proteins[row["accession"]]
            position = int(row["position"])
            flags = detectability(protein["seq"], position, boundaries[row["accession"]])
            record = {
                "dataset": dataset,
                "accession": row["accession"],
                "position": position,
                "observed": int(row["observed_in_retrospective_dataset"]),
            }
            record.update(flags)
            for label, column in MODELS.items():
                record[label] = float(row[column])
            site_records.append(record)

    with OUT_SITES.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(site_records[0]))
        writer.writeheader()
        writer.writerows(site_records)

    by_protein = {}
    for record in site_records:
        by_protein.setdefault((record["dataset"], record["accession"]), []).append(record)

    summary_rows = []
    intervals = {}
    for dataset, _, _ in COHORTS:
        for stratum in STRATA:
            def in_stratum(record):
                if stratum == "all":
                    return True
                if stratum == "detectable":
                    return record["detectable_any_missed_cleavage"] == 1.0
                if stratum == "detectable_strict":
                    return record["detectable_strict_fully_cleaved"] == 1.0
                return record["detectable_any_missed_cleavage"] == 0.0

            kept = []
            for (ds, acc), records in by_protein.items():
                if ds != dataset:
                    continue
                subset = [r for r in records if in_stratum(r)]
                observed = sum(r["observed"] for r in subset)
                if len(subset) >= 2 and observed >= 1:
                    kept.append((acc, subset))
            if not kept:
                continue

            per_model = {}
            for label in MODELS:
                rows = []
                for acc, subset in kept:
                    groups = groups_by_score([(r[label], r["observed"]) for r in subset])
                    size = len(subset)
                    observed = sum(r["observed"] for r in subset)
                    r1, r2, rm = random_metrics(size, observed)
                    rows.append({
                        "accession": acc,
                        "n": size,
                        "k": observed,
                        "top1": top1_tie_averaged(groups),
                        "top2": top2_tie_averaged(groups),
                        "mrr": expected_reciprocal_rank(groups),
                        "rand_top1": r1,
                        "rand_top2": r2,
                        "rand_mrr": rm,
                    })
                per_model[label] = rows

            reference = per_model["v2_chem"]
            n_proteins = len(reference)
            n_sites = sum(r["n"] for r in reference)
            n_observed = sum(r["k"] for r in reference)
            for label, rows in per_model.items():
                row = {
                    "dataset": dataset,
                    "stratum": stratum,
                    "model": label,
                    "n_proteins": n_proteins,
                    "n_candidate_sites": n_sites,
                    "n_observed_sites": n_observed,
                }
                for metric in ("top1", "top2", "mrr"):
                    row[metric] = round(sum(r[metric] for r in rows) / n_proteins, 4)
                    row["random_" + metric] = round(sum(r["rand_" + metric] for r in rows) / n_proteins, 4)
                    diffs = [r[metric] - r["rand_" + metric] for r in rows]
                    mean, lo, hi = bootstrap_interval(diffs, SEED)
                    row[metric + "_minus_random"] = round(mean, 5)
                    row[metric + "_minus_random_ci"] = f"[{lo:.5f}, {hi:.5f}]"
                    intervals[f"{dataset}|{stratum}|{label}|{metric}_vs_random"] = {
                        "difference": mean, "ci_95": [lo, hi], "n_proteins": n_proteins,
                    }
                summary_rows.append(row)

            for metric in ("top1", "top2", "mrr"):
                chem = {r["accession"]: r[metric] for r in per_model["v2_chem"]}
                for other in ("v2_full", "v1_corrected_mlp"):
                    paired = [chem[r["accession"]] - r[metric] for r in per_model[other]]
                    mean, lo, hi = bootstrap_interval(paired, SEED)
                    intervals[f"{dataset}|{stratum}|v2_chem_minus_{other}|{metric}"] = {
                        "difference": mean, "ci_95": [lo, hi], "n_proteins": n_proteins,
                    }

    with OUT_SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    OUT_INTERVALS.write_text(json.dumps(intervals, indent=2), encoding="utf-8")
    OUT_AUDIT.write_text(json.dumps({
        "purpose": "retrospective form of BLIND_TEST_PROTOCOL_V2 section 4 reverse test",
        "scores_read_from": str(SCORES.relative_to(ROOT)),
        "scores_refitted": False,
        "external_deployment": "single LightGBM binary member, no ProtBERT embeddings for external sequences",
        "strata": list(STRATA),
        "inclusion_rule": "at least two candidates and at least one observed site inside the stratum",
        "tie_rule": "mean annotation among tied top-ranked sites, extended to Top-2 and MRR",
        "random_baseline": "per protein, from that protein's own stratum size and observed count",
        "bootstrap": {"unit": "external protein", "replicates": BOOTSTRAP, "seed": SEED},
        "n_site_records": len(site_records),
        "interpretation_limit": "an unobserved Cys is not a confirmed chemical negative; both cohorts are retrospective",
    }, indent=2), encoding="utf-8")

    header = f"{'dataset':<18}{'stratum':<20}{'model':<18}{'prot':>6}{'top1':>8}{'rand':>8}{'diff':>9}  ci"
    print(header)
    for row in summary_rows:
        print(f"{row['dataset']:<18}{row['stratum']:<20}{row['model']:<18}"
              f"{row['n_proteins']:>6}{row['top1']:>8.4f}{row['random_top1']:>8.4f}"
              f"{row['top1_minus_random']:>9.4f}  {row['top1_minus_random_ci']}")


if __name__ == "__main__":
    main()
