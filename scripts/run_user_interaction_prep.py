"""Free-route preparation for the three-protein interaction question.

Reads the AlphaFold Server monomer models that already exist locally, derives
the ordered regions from per-residue pLDDT, writes trimmed single-chain PDB
files for rigid-body docking, and emits ready-to-submit AlphaFold Server job
requests for the three pairwise complexes.  Nothing here is paid and nothing is
uploaded; the job files are written to disk for the user to submit.
"""
from __future__ import annotations

import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select

from common import RESULTS, ROOT, sha256, write_csv, write_json
import score_experimental_candidates as v1
from common import INPUTS

AF_DIR = Path(
    "/Users/lyuguohao/.codex/.chatgpt-projects/"
    "g-p-6aa27e760cf88191b6a228fe3c180b34/analysis_outputs/af_server_results_20260912"
)
STRUCTURES = ROOT / "structures"
TARGETS = {
    "XP_004233015.1": "xp_004233015_1_wt",
    "Solyc06g068990.4.1": "solyc06g068990_4_1_wt",
    "Solyc12g006920.2.1": "solyc12g006920_2_1_wt",
}
PLDDT_CUTOFF = 70.0
MIN_SEGMENT = 30

# Regions the motif scan flagged, checked against the ordered segments.
FEATURES_OF_INTEREST = {
    "XP_004233015.1": {
        "T98_T101_MAPK_phosphosite_cluster": (98, 101),
        "S184_MAPK_phosphosite": (184, 184),
        "C193_C198_candidate_persulfidation": (193, 198),
        "WRKY_domain_core": (171, 234),
    },
    "Solyc06g068990.4.1": {
        "activation_loop_TDY": (270, 272),
        "C246_evidence_site": (246, 246),
        "C284_evidence_site": (284, 284),
        "C599_model_top1": (599, 599),
    },
    "Solyc12g006920.2.1": {
        "D_motif_RRRLTWERLEI": (470, 480),
        "D_motif_KRQVLLEL": (116, 123),
        "C379_evidence_site": (379, 379),
    },
}
PAIRS = [
    ("Solyc06g068990.4.1", "XP_004233015.1", "MAPK_plus_WRKY"),
    ("Solyc06g068990.4.1", "Solyc12g006920.2.1", "MAPK_plus_PP2A_B56"),
    ("XP_004233015.1", "Solyc12g006920.2.1", "WRKY_plus_PP2A_B56"),
]


class Ordered(Select):
    def __init__(self, keep):
        self.keep = keep

    def accept_residue(self, residue):
        return residue.id[1] in self.keep

    def accept_atom(self, atom):
        return atom.element != "H"


def plddt_profile(path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", str(path))
    chain = next(iter(next(iter(structure)).get_chains()))
    profile = {}
    for residue in chain:
        if "CA" not in residue:
            continue
        profile[residue.id[1]] = float(residue["CA"].get_bfactor())
    return structure, chain, profile


def segments(profile, cutoff=PLDDT_CUTOFF, minimum=MIN_SEGMENT):
    ordered = sorted(position for position, value in profile.items() if value >= cutoff)
    blocks = []
    for position in ordered:
        if blocks and position == blocks[-1][-1] + 1:
            blocks[-1].append(position)
        else:
            blocks.append([position])
    return [(block[0], block[-1]) for block in blocks if len(block) >= minimum]


def main():
    started = time.time()
    STRUCTURES.mkdir(parents=True, exist_ok=True)
    sequences = {
        name: value.strip().upper().rstrip("*")
        for name, value in v1.read_fasta(INPUTS / "candidate_proteins.fasta").items()
    }

    rows = []
    prepared = {}
    for name, folder in TARGETS.items():
        source = AF_DIR / folder / f"fold_{folder}_model_0.cif"
        if not source.exists():
            raise SystemExit(f"Missing AlphaFold model: {source}")
        structure, chain, profile = plddt_profile(source)
        blocks = segments(profile)
        keep = {p for start, stop in blocks for p in range(start, stop + 1)}
        target = STRUCTURES / f"{folder}_ordered.pdb"
        writer = PDBIO()
        writer.set_structure(structure)
        writer.save(str(target), select=Ordered(keep))
        full = STRUCTURES / f"{folder}_full.pdb"
        writer.save(str(full), select=Ordered(set(profile)))
        prepared[name] = {
            "ordered_pdb": str(target.relative_to(ROOT)),
            "full_pdb": str(full.relative_to(ROOT)),
            "ordered_segments": blocks,
            "n_residues_total": len(profile),
            "n_residues_ordered": len(keep),
            "mean_plddt": float(np.mean(list(profile.values()))),
            "fraction_ordered": len(keep) / len(profile),
            "source_model": str(source),
        }
        for label, (start, stop) in FEATURES_OF_INTEREST[name].items():
            values = [profile[p] for p in range(start, stop + 1) if p in profile]
            inside = any(
                block_start <= start and stop <= block_stop
                for block_start, block_stop in blocks
            )
            rows.append({
                "protein": name,
                "feature": label,
                "residue_range": f"{start}-{stop}" if start != stop else str(start),
                "mean_plddt": round(float(np.mean(values)), 2) if values else "",
                "min_plddt": round(float(np.min(values)), 2) if values else "",
                "inside_an_ordered_segment": inside,
                "usable_for_rigid_docking": inside,
                "note": (
                    "ordered, can be used as a docking restraint"
                    if inside else
                    "low-confidence or disordered; a rigid-body docking restraint here is not meaningful"
                ),
            })
        print(json.dumps({
            "protein": name,
            "residues": len(profile),
            "fraction_ordered": round(len(keep) / len(profile), 3),
            "mean_plddt": round(float(np.mean(list(profile.values()))), 1),
            "ordered_segments": blocks,
        }), flush=True)

    write_csv(RESULTS / "v2_user_interaction_region_check.csv", rows)

    jobs_dir = STRUCTURES / "alphafold_server_jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for first, second, label in PAIRS:
        payload = [{
            "name": label,
            "modelSeeds": [],
            "sequences": [
                {"proteinChain": {"sequence": sequences[first], "count": 1, "useStructureTemplate": True}},
                {"proteinChain": {"sequence": sequences[second], "count": 1, "useStructureTemplate": True}},
            ],
            "dialect": "alphafoldserver",
            "version": 3,
        }]
        path = jobs_dir / f"{label}.json"
        path.write_text(json.dumps(payload, indent=1))
        written.append(str(path.relative_to(ROOT)))
    # domain-trimmed variant for the WRKY, which is largely disordered full length
    wrky = sequences["XP_004233015.1"]
    trimmed = wrky[170:234]
    payload = [{
        "name": "MAPK_plus_WRKY_domain_171_234",
        "modelSeeds": [],
        "sequences": [
            {"proteinChain": {"sequence": sequences["Solyc06g068990.4.1"], "count": 1, "useStructureTemplate": True}},
            {"proteinChain": {"sequence": trimmed, "count": 1, "useStructureTemplate": True}},
        ],
        "dialect": "alphafoldserver",
        "version": 3,
    }]
    path = jobs_dir / "MAPK_plus_WRKY_domain_171_234.json"
    path.write_text(json.dumps(payload, indent=1))
    written.append(str(path.relative_to(ROOT)))

    write_json(RESULTS / "v2_user_interaction_prep_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "purpose": "free-route inputs for the three-protein interaction question",
        "paid_service_used": False,
        "monomer_models_source": "AlphaFold Server results already present locally, produced 2026-09-12",
        "monomer_models_terms": "AlphaFold Server output is for non-commercial use; see terms_of_use.md in the source archive",
        "plddt_cutoff": PLDDT_CUTOFF,
        "minimum_segment_length": MIN_SEGMENT,
        "prepared_structures": prepared,
        "alphafold_server_job_files": written,
        "wrky_domain_trim": {"range": "171-234", "length": len(trimmed), "sequence": trimmed},
        "versions": {"python": sys.version, "numpy": np.__version__},
        "input_hashes": {
            str((INPUTS / "candidate_proteins.fasta").relative_to(ROOT)): sha256(
                INPUTS / "candidate_proteins.fasta"
            )
        },
    })
    print(json.dumps(rows, indent=1, default=str))


if __name__ == "__main__":
    main()
