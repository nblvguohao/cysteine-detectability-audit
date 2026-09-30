"""Per-protein independence screen and v2 prediction for the user's 15 proteins.

The question this answers is not only "which Cys ranks first" but first of all
"is this protein independent of the training cohort".  Three checks run per
protein, reusing the rule of the 2026-09-12 cross-conversation screen:

1. exact sequence identity against every training-cohort protein;
2. strong homology: top-five shared-5mer training proteins, then local
   BLOSUM62 alignment, strong means at least 30 percent identity and at least
   80 percent coverage on both sides;
3. local-window duplication: the 31-aa window centred on each candidate Cys
   compared with the windows of annotated training positives.

For a protein that is itself in the training cohort the annotated positions are
reported as data, not as prediction, because the model was trained on them.
"""
from __future__ import annotations

import json
import platform
import sys
import time
from collections import Counter, defaultdict

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

from common import INPUTS, RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import score_experimental_candidates as v1
import v2_apply
from v2_stack import FEATURES, load_v2_matrix

SEEDS = [11, 29, 47, 71, 97]
FLANK = 15
ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")


def centred_window(sequence, position, flank=FLANK):
    return sequence[max(0, position - 1 - flank): min(len(sequence), position + flank)]


def build_screen(candidates, train_sequences, train_rows):
    accessions = sorted({row["accession"] for row in train_rows})
    sequences = [train_sequences[a] for a in accessions]
    by_sequence = defaultdict(list)
    for accession, sequence in zip(accessions, sequences):
        by_sequence[sequence].append(accession)

    positive_windows = defaultdict(list)
    cohort_sites = defaultdict(list)
    for row in train_rows:
        accession, position, label = row["accession"], int(row["position"]), int(row["label"])
        cohort_sites[accession].append((position, label))
        if label == 1:
            positive_windows[
                centred_window(train_sequences[accession], position)
            ].append((accession, position))

    kmers = defaultdict(list)
    for index, sequence in enumerate(sequences):
        for word in {sequence[i:i + 5] for i in range(max(0, len(sequence) - 4))}:
            if "X" not in word:
                kmers[word].append(index)

    aligner = PairwiseAligner(
        mode="local",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-10,
        extend_gap_score=-0.5,
    )
    alphabet = set(aligner.substitution_matrix.alphabet)
    clean = lambda value: "".join(c if c in alphabet else "X" for c in value)

    results = {}
    for name, sequence in candidates.items():
        counts = Counter()
        for word in {sequence[i:i + 5] for i in range(max(0, len(sequence) - 4))}:
            if "X" in word:
                continue
            for index in kmers.get(word, ()):
                counts[index] += 1
        shortlist = [
            index for index in sorted(counts, key=lambda i: (-counts[i], accessions[i]))[:5]
            if counts[index] >= 4
        ]
        best = None
        strong = []
        for index in shortlist:
            target = sequences[index]
            alignment = aligner.align(clean(sequence), clean(target))[0]
            query_blocks, target_blocks = alignment.aligned
            paired = sum(int(end - start) for start, end in query_blocks)
            identical = 0
            for (qs, qe), (ts, te) in zip(query_blocks, target_blocks):
                identical += sum(
                    1 for a, b in zip(sequence[qs:qe], target[ts:te]) if a == b
                )
            identity = identical / paired if paired else 0.0
            query_coverage = paired / len(sequence)
            target_coverage = paired / len(target)
            record = {
                "training_accession": accessions[index],
                "shared_5mers": int(counts[index]),
                "identity": float(identity),
                "query_coverage": float(query_coverage),
                "target_coverage": float(target_coverage),
                "alignment_score": float(alignment.score),
            }
            if best is None or record["alignment_score"] > best["alignment_score"]:
                best = record
            if identity >= 0.30 and query_coverage >= 0.80 and target_coverage >= 0.80:
                strong.append(record)

        exact = sorted(by_sequence.get(sequence, []))
        window_hits = []
        for index, residue in enumerate(sequence):
            if residue != "C":
                continue
            position = index + 1
            hits = positive_windows.get(centred_window(sequence, position), [])
            if hits:
                window_hits.append({
                    "candidate_position": position,
                    "training_matches": [
                        {"accession": a, "position": p} for a, p in hits
                    ],
                })

        annotated = []
        for accession in exact + [r["training_accession"] for r in strong]:
            for position, label in sorted(cohort_sites.get(accession, [])):
                if label == 1:
                    annotated.append({"accession": accession, "position": position})

        if exact:
            status = "in_training_cohort_exact_sequence"
        elif strong:
            status = "strong_homolog_of_training_protein"
        elif window_hits:
            status = "local_window_duplicates_training_positive"
        else:
            status = "independent_by_this_screen"

        results[name] = {
            "status": status,
            "exact_training_accessions": exact,
            "strong_homologs": strong,
            "best_screened_hit": best,
            "positive_window_duplications": window_hits,
            "training_annotated_positions_on_matched_proteins": annotated,
        }
    return results


def score_distribution_context(scores_by_set):
    """Empirical percentile of a site score inside the cohort OOF distributions."""
    X, names, y, folds, proteins, components, positions, meta = load_v2_matrix()
    context = {}
    for feature_set, path in (("chem", "v2_stack_chem_oof.npz"), ("full", "v2_stack_full_oof.npz")):
        archive = np.load(RESULTS / path, allow_pickle=False)
        oof = archive["lgb_bin"]
        context[feature_set] = {
            "annotated": np.sort(oof[y == 1]),
            "unlabelled": np.sort(oof[y == 0]),
        }
    return context


def percentile_and_ratio(context, feature_set, score):
    annotated = context[feature_set]["annotated"]
    unlabelled = context[feature_set]["unlabelled"]
    above_annotated = float(np.searchsorted(annotated, score) / len(annotated))
    above_unlabelled = float(np.searchsorted(unlabelled, score) / len(unlabelled))
    tail_annotated = max(1.0 - above_annotated, 1.0 / len(annotated))
    tail_unlabelled = max(1.0 - above_unlabelled, 1.0 / len(unlabelled))
    return {
        "percentile_among_annotated": above_annotated,
        "percentile_among_unlabelled": above_unlabelled,
        "tail_likelihood_ratio": tail_annotated / tail_unlabelled,
    }


def main():
    started = time.time()
    encoder = v2_apply.load_source_encoder()
    reference_names = np.load(FEATURES / "v2_features.npz", allow_pickle=False)["names"].astype(str)
    source_sequences, structures, _, _ = encoder.read_inputs()
    train_rows = read_csv(INPUTS / "primary_site_folds.csv")
    train_accessions = {row["accession"] for row in train_rows}
    train_sequences = {a: source_sequences[a] for a in train_accessions}

    targets = {
        name: sequence.strip().upper().rstrip("*")
        for name, sequence in v1.read_fasta(INPUTS / "candidate_proteins.fasta").items()
    }
    screen = build_screen(targets, train_sequences, train_rows)
    print(json.dumps({k: v["status"] for k, v in screen.items()}, indent=1), flush=True)

    reverse = defaultdict(list)
    for accession, sequence in source_sequences.items():
        reverse[sequence].append(accession)
    structure_lookup = {}
    scorable = {}
    withheld = {}
    for name, sequence in targets.items():
        unknown = sorted({c for c in sequence if c not in ALPHABET})
        if unknown:
            withheld[name] = unknown
        scorable[name] = sequence
        matches = sorted(reverse.get(sequence, []))
        usable = [a for a in matches if structures.get(a, {}).get("status") == "ok"]
        if usable:
            structure_lookup[name] = {
                int(entry["pos"]): entry for entry in structures[usable[0]]["cys"]
            }

    records = [{"accession": n, "sequence": s} for n, s in sorted(scorable.items())]
    matrix, proteins, positions = v2_apply.build_matrix(
        records, encoder, reference_names, structure_lookup=structure_lookup
    )

    predictions = {}
    deployment_audit = {}
    for feature_set in ("chem", "full"):
        boosters, names, columns, deployment = v2_apply.fit_deployed_model(
            feature_set, seeds=SEEDS
        )
        stacked = np.stack([
            booster.predict(np.ascontiguousarray(matrix[:, columns])) for booster in boosters
        ])
        predictions[feature_set] = {"mean": stacked.mean(axis=0), "sd": stacked.std(axis=0, ddof=1)}
        deployment_audit[feature_set] = deployment
        print(json.dumps({"feature_set": feature_set, "rounds": deployment["rounds"]}), flush=True)

    context = score_distribution_context(predictions)

    rows = []
    for name in sorted(scorable):
        mask = np.flatnonzero(proteins == name)
        site_positions = positions[mask]
        chem = predictions["chem"]["mean"][mask]
        chem_sd = predictions["chem"]["sd"][mask]
        full = predictions["full"]["mean"][mask]
        full_sd = predictions["full"]["sd"][mask]
        _, chem_rank = v1.ranks(chem, site_positions)
        _, full_rank = v1.ranks(full, site_positions)
        unknown_positions = {
            i + 1 for i, residue in enumerate(scorable[name]) if residue not in ALPHABET
        }
        for index, position in enumerate(site_positions):
            window = centred_window(scorable[name], int(position))
            chem_context = percentile_and_ratio(context, "chem", float(chem[index]))
            full_context = percentile_and_ratio(context, "full", float(full[index]))
            nearest_unknown = (
                min(abs(int(position) - u) for u in unknown_positions)
                if unknown_positions else None
            )
            rows.append({
                "protein": name,
                "position": int(position),
                "independence_status": screen[name]["status"],
                "v2_chem_score": float(chem[index]),
                "v2_chem_rank": chem_rank[index],
                "v2_chem_seed_sd": float(chem_sd[index]),
                "v2_chem_percentile_among_annotated": chem_context["percentile_among_annotated"],
                "v2_chem_tail_likelihood_ratio": chem_context["tail_likelihood_ratio"],
                "v2_full_score": float(full[index]),
                "v2_full_rank": full_rank[index],
                "v2_full_seed_sd": float(full_sd[index]),
                "v2_full_percentile_among_annotated": full_context["percentile_among_annotated"],
                "v2_full_tail_likelihood_ratio": full_context["tail_likelihood_ratio"],
                "rank_gap_full_minus_chem": full_rank[index] - chem_rank[index],
                "context_31aa": window,
                "window_contains_unknown_residue": bool(
                    any(c not in ALPHABET for c in window)
                ),
                "distance_to_nearest_unknown_residue": nearest_unknown,
                "primary_mutant": f"C{int(position)}S",
                "orthogonal_mutant": f"C{int(position)}A",
                "scope": "annotation-ranking priority; not a chemical occupancy probability",
            })

    write_csv(RESULTS / "v2_user_protein_sites.csv", rows)
    write_json(RESULTS / "v2_user_protein_independence.json", screen)
    write_json(RESULTS / "v2_user_protein_run_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "n_proteins": len(scorable),
        "n_sites": len(rows),
        "withheld_unknown_residues": withheld,
        "unknown_residue_policy": "scored but flagged; every site keeps a window-contains-unknown flag and the distance to the nearest unknown residue",
        "screen_rule": "exact sequence identity; top-five shared-5mer plus local BLOSUM62 with identity >= 0.30 and bilateral coverage >= 0.80; 31-aa positive window duplication",
        "deployment": deployment_audit,
        "seeds": SEEDS,
        "likelihood_ratio_definition": "ratio of the upper-tail mass of the annotated OOF score distribution to that of the unlabelled OOF score distribution at this site score; an annotation-enrichment factor, not a chemical probability",
        "hardware": {"platform": platform.platform()},
        "versions": {
            "python": sys.version, "numpy": np.__version__,
            "lightgbm": v2_apply.members.lgb.__version__,
        },
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "candidate_proteins.fasta",
                INPUTS / "primary_site_folds.csv",
                FEATURES / "v2_features.npz",
            ]
        },
    })
    print("sites written:", len(rows))


if __name__ == "__main__":
    main()
