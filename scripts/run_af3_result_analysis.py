"""Analyse AlphaFold Server (AF3) results for the three-protein interaction question.

Reads the five diffusion samples per job, ranks them by AF3's own ranking_score
and pairwise ipTM, and for the best-supported pair extracts the confident
sub-interface using per-residue cross-chain PAE (not just static contact
geometry from PLIP, which cannot tell a confident contact from a spurious one).

A coupled folding-and-binding check compares complex-model pLDDT against the
isolated AlphaFold Server monomer at the same residues: a large pLDDT gain
localised to the interface, together with low cross-chain PAE, is the
strongest available signal for a real interaction without wet-lab data.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser

from common import INPUTS, RESULTS, ROOT, sha256, write_csv, write_json

AF3_DIR = ROOT / "interaction" / "af3_results"
MONOMER_DIR = Path(
    "/Users/lyuguohao/.codex/.chatgpt-projects/"
    "g-p-6aa27e760cf88191b6a228fe3c180b34/analysis_outputs/af_server_results_20260912"
)
JOBS = {
    "mapk_plus_wrky": ("Solyc06g068990.4.1", "XP_004233015.1", "solyc06g068990_4_1_wt", "xp_004233015_1_wt"),
    "mapk_plus_pp2a_b56": ("Solyc06g068990.4.1", "Solyc12g006920.2.1", "solyc06g068990_4_1_wt", "solyc12g006920_2_1_wt"),
    "wrky_plus_pp2a_b56": ("XP_004233015.1", "Solyc12g006920.2.1", "xp_004233015_1_wt", "solyc12g006920_2_1_wt"),
    "mapk_plus_wrky_domain_171_234": ("Solyc06g068990.4.1", "XP_004233015.1 (171-234)", "solyc06g068990_4_1_wt", None),
}
RAW_RESULTS_DIR = Path("/tmp/af3_results")  # the five-sample archive extracted from the user's zip


def monomer_plddt(folder, positions):
    if folder is None:
        return {}
    path = MONOMER_DIR / folder / f"fold_{folder}_model_0.cif"
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("m", str(path))
    chain = next(iter(next(iter(structure)).get_chains()))
    return {
        residue.id[1]: float(residue["CA"].get_bfactor())
        for residue in chain
        if residue.id[1] in positions and "CA" in residue
    }


def complex_plddt(job, chain_id, positions):
    path = AF3_DIR / job / "best_model_0.cif"
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("m", str(path))
    chain = next(iter(structure))[chain_id]
    return {
        residue.id[1]: float(residue["CA"].get_bfactor())
        for residue in chain
        if residue.id[1] in positions and "CA" in residue
    }


def main():
    started = time.time()
    ranking_rows = []
    for job in JOBS:
        for sample in range(5):
            path = RAW_RESULTS_DIR / job / f"fold_{job}_summary_confidences_{sample}.json"
            data = json.loads(path.read_text())
            ranking_rows.append({
                "job": job,
                "sample": sample,
                "ranking_score": data["ranking_score"],
                "iptm": data["iptm"],
                "ptm": data["ptm"],
                "chain_pair_iptm_off_diagonal": data["chain_pair_iptm"][0][1],
                "chain_pair_pae_min_off_diagonal": min(
                    data["chain_pair_pae_min"][0][1], data["chain_pair_pae_min"][1][0]
                ),
                "fraction_disordered": data["fraction_disordered"],
                "has_clash": data["has_clash"],
            })
    ranking_rows.sort(key=lambda row: -row["ranking_score"])
    write_csv(RESULTS / "v2_af3_model_ranking.csv", ranking_rows)

    best_per_job = {}
    for job in JOBS:
        rows = [row for row in ranking_rows if row["job"] == job]
        best_per_job[job] = max(rows, key=lambda row: row["ranking_score"])

    # confident-interface analysis for the best-supported pair
    job = "wrky_plus_pp2a_b56"
    receptor, ligand, receptor_monomer_folder, ligand_monomer_folder = JOBS[job]
    full_data = json.loads(
        (RAW_RESULTS_DIR / job / f"fold_{job}_full_data_0.json").read_text()
    )
    pae = np.asarray(full_data["pae"])
    n_a = 317  # chain A (WRKY) length in this job

    interface_rows = []
    seen = set()
    plip_path = AF3_DIR / job / "plip_report.txt"
    import re

    # PLIP's fixed-width markdown tables: "| RESNR | RESTYPE | RESCHAIN | ..."
    # column 1 is the receptor-chain (A, i.e. WRKY here) residue number whenever
    # RESCHAIN (column 3) is A; only the first two columns are needed here.
    row_pattern = re.compile(r"^\|\s*(\d+)\s*\|\s*([A-Z]{3})\s*\|\s*([AB])\s*\|")
    for line in plip_path.read_text().splitlines():
        match = row_pattern.match(line)
        if not match:
            continue
        residue_number, _, chain = match.groups()
        if chain != "A":
            continue
        residue_number = int(residue_number)
        if residue_number not in seen:
            seen.add(residue_number)
            idx = residue_number - 1
            min_pae = float(min(pae[idx, n_a:].min(), pae[n_a:, idx].min()))
            interface_rows.append({"residue": residue_number, "min_cross_chain_pae": min_pae})

    interface_rows.sort(key=lambda row: row["min_cross_chain_pae"])
    confident = [row["residue"] for row in interface_rows if row["min_cross_chain_pae"] <= 10.0]

    monomer_values = monomer_plddt(receptor_monomer_folder, set(range(1, 40)))
    complex_values = complex_plddt(job, "A", set(range(1, 40)))
    coupled_folding_rows = []
    for position in sorted(set(monomer_values) & set(complex_values)):
        coupled_folding_rows.append({
            "residue": position,
            "monomer_plddt": monomer_values[position],
            "complex_plddt": complex_values[position],
            "gain": complex_values[position] - monomer_values[position],
        })
    write_csv(RESULTS / "v2_af3_wrky_pp2a_coupled_folding.csv", coupled_folding_rows)
    write_json(RESULTS / "v2_af3_wrky_pp2a_confident_interface.json", {
        "job": job,
        "confident_residues_wrky_side": confident,
        "confident_threshold_angstrom": 10.0,
        "all_plip_interface_residues_with_pae": interface_rows,
        "interpretation": "PLIP reports geometric contacts from one static model; a contact is only trustworthy where AF3's own cross-chain PAE is also low. Contacts with high cross-chain PAE are excluded here even though PLIP lists them.",
    })

    write_json(RESULTS / "v2_af3_result_analysis_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "n_jobs": len(JOBS),
        "n_samples_per_job": 5,
        "best_ranking_score_by_job": {job: row["ranking_score"] for job, row in best_per_job.items()},
        "best_iptm_by_job": {job: row["iptm"] for job, row in best_per_job.items()},
        "best_cross_chain_pae_min_by_job": {
            job: row["chain_pair_pae_min_off_diagonal"] for job, row in best_per_job.items()
        },
        "confident_interface_method": "per-residue cross-chain PAE threshold plus coupled-folding pLDDT comparison against the isolated AlphaFold Server monomer",
        "conventional_confidence_thresholds_note": "ipTM >= 0.8 confident, 0.6-0.8 moderate, < 0.6 low; these are AlphaFold-Multimer community heuristics, not an AF3-validated cutoff, applied here descriptively",
        "versions": {"python": sys.version, "numpy": np.__version__},
        "input_hashes": {
            str((INPUTS / "candidate_proteins.fasta").relative_to(ROOT)): sha256(
                INPUTS / "candidate_proteins.fasta"
            )
        },
    })
    print(json.dumps(best_per_job, indent=1, default=str))
    print("confident WRKY-side interface residues:", confident)


if __name__ == "__main__":
    main()
