"""Rescore the wet-lab candidate proteins with the v2 deployed models.

The v1 candidate table stays in place; this program adds the v2 chemistry and
full models, their seed spread, the cross-model rank span and an updated assay
shortlist.  Scores remain full-fit prioritisation numbers: they are neither
chemical occupancy probabilities nor independent validation.
"""
from __future__ import annotations

import json
import platform
import sys
import time
from collections import defaultdict

import numpy as np

from common import INPUTS, RESULTS, ROOT, read_csv, sha256, write_csv, write_json
import score_experimental_candidates as v1
import v2_apply
from v2_stack import FEATURES

SEEDS = [11, 29, 47, 71, 97]


def stability(chem_rank, full_rank, v1_rank, chem_sd, full_sd):
    """Predeclared three-model agreement rule, mirroring the v1 definition."""
    span = max(chem_rank, full_rank, v1_rank) - min(chem_rank, full_rank, v1_rank)
    if max(chem_rank, full_rank) <= 2 and span <= 2 and max(chem_sd, full_sd) <= 0.08:
        return "high"
    if max(chem_rank, full_rank) <= 4 and span <= 4:
        return "moderate"
    return "low"


def main():
    started = time.time()
    encoder = v2_apply.load_source_encoder()
    reference_names = np.load(FEATURES / "v2_features.npz", allow_pickle=False)[
        "names"
    ].astype(str)
    source_sequences, structures, _, _ = encoder.read_inputs()
    reverse = defaultdict(list)
    for accession, sequence in source_sequences.items():
        reverse[sequence].append(accession)

    targets = v1.read_fasta(INPUTS / "candidate_proteins.fasta")
    previous = {
        (row["protein"], int(row["position"])): row
        for row in read_csv(RESULTS / "experimental_candidate_all_cys.csv")
    }

    scorable = {}
    structure_lookup = {}
    withheld = {}
    for protein, sequence in targets.items():
        sequence = sequence.strip().upper()
        if any(residue not in "ACDEFGHIKLMNPQRSTVWY" for residue in sequence):
            withheld[protein] = "sequence contains unknown residues; no silent imputation"
            continue
        scorable[protein] = sequence
        matches = sorted(reverse.get(sequence, []))
        usable = [
            accession for accession in matches
            if structures.get(accession, {}).get("status") == "ok"
        ]
        if usable:
            structure_lookup[protein] = {
                int(entry["pos"]): entry for entry in structures[usable[0]]["cys"]
            }

    records = [
        {"accession": protein, "sequence": sequence}
        for protein, sequence in sorted(scorable.items())
    ]
    matrix, proteins, positions = v2_apply.build_matrix(
        records, encoder, reference_names, structure_lookup=structure_lookup
    )

    scores = {}
    deployment_audit = {}
    for feature_set in ("chem", "full"):
        boosters, names, columns, deployment = v2_apply.fit_deployed_model(
            feature_set, seeds=SEEDS
        )
        stacked = np.stack([
            booster.predict(np.ascontiguousarray(matrix[:, columns]))
            for booster in boosters
        ])
        scores[feature_set] = {
            "mean": stacked.mean(axis=0),
            "sd": stacked.std(axis=0, ddof=1),
        }
        deployment_audit[feature_set] = deployment
        print(json.dumps({
            "feature_set": feature_set,
            "rounds": deployment["rounds"],
            "seeds": SEEDS,
            "n_sites": int(stacked.shape[1]),
            "seconds": round(time.time() - started, 1),
        }), flush=True)

    rows = []
    summaries = []
    for protein in sorted(scorable):
        mask = np.flatnonzero(proteins == protein)
        site_positions = positions[mask]
        chem_mean = scores["chem"]["mean"][mask]
        chem_sd = scores["chem"]["sd"][mask]
        full_mean = scores["full"]["mean"][mask]
        full_sd = scores["full"]["sd"][mask]
        chem_order, chem_rank = v1.ranks(chem_mean, site_positions)
        full_order, full_rank = v1.ranks(full_mean, site_positions)
        protein_rows = []
        for index, position in enumerate(site_positions):
            old = previous[(protein, int(position))]
            v1_rank = int(old["corrected_mlp_rank"])
            row = {
                "protein": protein,
                "position": int(position),
                "v2_chem_score": float(chem_mean[index]),
                "v2_chem_seed_sd": float(chem_sd[index]),
                "v2_chem_rank": chem_rank[index],
                "v2_full_score": float(full_mean[index]),
                "v2_full_seed_sd": float(full_sd[index]),
                "v2_full_rank": full_rank[index],
                "v1_corrected_mlp_score": float(old["corrected_mlp_score"]),
                "v1_corrected_mlp_rank": v1_rank,
                "v1_histgb_rank": int(old["histgb_rank"]),
                "rank_span_v2chem_v2full_v1mlp": max(chem_rank[index], full_rank[index], v1_rank)
                - min(chem_rank[index], full_rank[index], v1_rank),
                "rank_stability": stability(
                    chem_rank[index], full_rank[index], v1_rank,
                    float(chem_sd[index]), float(full_sd[index]),
                ),
                "structure_present": old["structure_present"],
                "relative_sasa": old["relative_sasa"],
                "plddt": old["plddt"],
                "predicted_disulfide": old["predicted_disulfide"],
                "in_theoretical_digest_window": old["in_theoretical_digest_window"],
                "peptide": old["peptide"],
                "peptide_cys_count": old["peptide_cys_count"],
                "context": old["context"],
                "primary_mutant": f"C{int(position)}S",
                "orthogonal_mutant": f"C{int(position)}A",
                "mutation_warning": "interpret with expression, folding and activity controls",
                "scope": "full-fit experimental prioritisation; not a chemical probability and not validation",
            }
            rows.append(row)
            protein_rows.append(row)
        top = [protein_rows[index] for index in chem_order[:2]]
        assay = sorted(
            [
                row for row in protein_rows
                if v1.bool_value(row["in_theoretical_digest_window"])
                and int(row["peptide_cys_count"]) == 1
            ],
            key=lambda row: (row["v2_chem_rank"], row["position"]),
        )[:2]
        summaries.append({
            "protein": protein,
            "status": "scored",
            "n_cys": len(site_positions),
            "v2_chem_top1_position": top[0]["position"],
            "v2_chem_top1_peptide": top[0]["peptide"],
            "v2_chem_top1_stability": top[0]["rank_stability"],
            "v2_chem_top2_position": top[1]["position"] if len(top) > 1 else "",
            "v2_chem_top2_peptide": top[1]["peptide"] if len(top) > 1 else "",
            "v2_chem_top2_stability": top[1]["rank_stability"] if len(top) > 1 else "",
            "v2_full_top1_position": protein_rows[full_order[0]]["position"],
            "v1_mlp_top1_position": min(
                protein_rows, key=lambda row: row["v1_corrected_mlp_rank"]
            )["position"],
            "assay_priority1_position": assay[0]["position"] if assay else "",
            "assay_priority1_peptide": assay[0]["peptide"] if assay else "",
            "assay_priority1_mutant": assay[0]["primary_mutant"] if assay else "",
            "assay_priority2_position": assay[1]["position"] if len(assay) > 1 else "",
            "assay_priority2_peptide": assay[1]["peptide"] if len(assay) > 1 else "",
            "assay_priority2_mutant": assay[1]["primary_mutant"] if len(assay) > 1 else "",
            "assay_rule": "v2 chemistry-model rank restricted to theoretical digest peptides carrying one Cys; model ranks are not modified",
        })
    for protein, reason in sorted(withheld.items()):
        summaries.append({
            "protein": protein,
            "status": "withheld_unknown_residues",
            "n_cys": targets[protein].count("C"),
            "v2_chem_top1_position": "", "v2_chem_top1_peptide": "", "v2_chem_top1_stability": "",
            "v2_chem_top2_position": "", "v2_chem_top2_peptide": "", "v2_chem_top2_stability": "",
            "v2_full_top1_position": "", "v1_mlp_top1_position": "",
            "assay_priority1_position": "", "assay_priority1_peptide": "", "assay_priority1_mutant": "",
            "assay_priority2_position": "", "assay_priority2_peptide": "", "assay_priority2_mutant": "",
            "assay_rule": reason,
        })

    write_csv(RESULTS / "v2_candidate_all_cys.csv", rows)
    write_csv(RESULTS / "v2_candidate_top2_and_assay.csv", summaries)
    write_json(RESULTS / "v2_candidate_run_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "primary_model": "v2 chemistry LightGBM binary member, five-seed full-fit mean",
        "secondary_model": "v2 full LightGBM binary member, five-seed full-fit mean",
        "secondary_model_meaning": "includes in-silico tryptic detectability features; a high rank partly reflects how visible the site is to a trypsin-based workflow",
        "reference_model": "v1 corrected MLP and frozen 23-feature HGB ranks from results/experimental_candidate_all_cys.csv",
        "deployment": deployment_audit,
        "seeds": SEEDS,
        "unknown_residue_policy": "withhold entire protein",
        "structure_features": "used where the candidate sequence matches a source protein with a usable structure; missing otherwise",
        "candidate_scores_used_for_validation": False,
        "hardware": {"platform": platform.platform()},
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "lightgbm": v2_apply.members.lgb.__version__,
        },
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "candidate_proteins.fasta",
                INPUTS / "frozen_benchmark_data.npz",
                FEATURES / "v2_features.npz",
                RESULTS / "experimental_candidate_all_cys.csv",
            ]
        },
        "source_encoder": {
            "path": str(v2_apply.SOURCE_PIPELINE),
            "sha256": sha256(v2_apply.SOURCE_PIPELINE),
        },
    })
    print(json.dumps(summaries, indent=1, default=str)[:800])


if __name__ == "__main__":
    main()
