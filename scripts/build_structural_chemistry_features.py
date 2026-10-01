"""Structural chemistry descriptors for every cysteine, from AlphaFold models.

The reactivity work left one information dimension untried.  The frozen feature
set carries only five structural columns (rel_sasa, plddt, pka, ncys8, nacid)
and about 20% of sites have none of them.  Everything else the project has
tested is sequence or digest information, and the conclusion there is now firm:
the annotation label is explained by peptide visibility, and a model that does
predict chemistry (ESM-2 on ABPP engagement, Spearman 0.4394) scores at chance
on the annotation benchmark.

This program computes real three-dimensional chemistry around each cysteine from
the AlphaFold v6 models of the cohort proteins: burial, packing at several
radii, the local electrostatic environment, hydrogen-bonding geometry to the
thiol sulfur, secondary-structure context and the distance to the nearest
partner cysteine in space rather than in sequence.

No label is read anywhere.  The output keeps the frozen site order so it can be
concatenated with `features/v2_features.npz`.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from common import INPUTS, ROOT, load_frozen, sha256, write_json

PDB_DIR = ROOT / "external/alphafold_tomato/pdb"
OUT = ROOT / "features"

# Kyte-Doolittle hydropathy and formal side-chain charge near pH 7
HYDROPATHY = {
    "ALA": 1.8, "CYS": 2.5, "ASP": -3.5, "GLU": -3.5, "PHE": 2.8, "GLY": -0.4,
    "HIS": -3.2, "ILE": 4.5, "LYS": -3.9, "LEU": 3.8, "MET": 1.9, "ASN": -3.5,
    "PRO": -1.6, "GLN": -3.5, "ARG": -4.5, "SER": -0.8, "THR": -0.7, "VAL": 4.2,
    "TRP": -0.9, "TYR": -1.3,
}
CHARGE = defaultdict(float, {"LYS": 1.0, "ARG": 1.0, "HIS": 0.1, "ASP": -1.0, "GLU": -1.0})
DONOR_ACCEPTOR = {
    "SER": ("OG",), "THR": ("OG1",), "TYR": ("OH",), "ASN": ("OD1", "ND2"),
    "GLN": ("OE1", "NE2"), "HIS": ("ND1", "NE2"), "LYS": ("NZ",),
    "ARG": ("NE", "NH1", "NH2"), "ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"),
    "TRP": ("NE1",), "CYS": ("SG",),
}
SHELLS = (4.0, 6.0, 8.0, 10.0, 12.0)
# Shrake-Rupley radii, coarse element mapping
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
PROBE = 1.40
SPHERE_POINTS = 128
# reference sidechain SASA of a free cysteine, used to make burial relative
CYS_REFERENCE_SASA = 110.0


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 8) // 2))
    parser.add_argument("--limit", type=int, default=0)
    return parser.parse_args()


def golden_sphere(n):
    index = np.arange(n, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * index / n)
    theta = math.pi * (1.0 + 5.0 ** 0.5) * index
    return np.stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1
    )


SPHERE = golden_sphere(SPHERE_POINTS)


def read_pdb(path):
    """Atoms, per-residue index and pLDDT from an AlphaFold model."""
    coords, elements, names, residues, seq_ids, plddt = [], [], [], [], [], []
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            element = line[76:78].strip() or line[12:16].strip()[0]
            if element == "H":
                continue
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            elements.append(element)
            names.append(line[12:16].strip())
            residues.append(line[17:20].strip())
            seq_ids.append(int(line[22:26]))
            plddt.append(float(line[60:66]))
    if not coords:
        return None
    return {
        "xyz": np.asarray(coords, dtype=np.float64),
        "element": np.asarray(elements),
        "name": np.asarray(names),
        "resname": np.asarray(residues),
        "resid": np.asarray(seq_ids, dtype=int),
        "plddt": np.asarray(plddt, dtype=float),
    }


def sidechain_sasa(structure, mask, neighbour_cut=15.0):
    """Shrake-Rupley accessible area of the atoms in ``mask``."""
    xyz = structure["xyz"]
    radii = np.asarray([VDW.get(e, 1.70) for e in structure["element"]]) + PROBE
    total = 0.0
    targets = np.flatnonzero(mask)
    for i in targets:
        centre = xyz[i]
        near = np.flatnonzero(
            (np.abs(xyz[:, 0] - centre[0]) < neighbour_cut)
            & (np.abs(xyz[:, 1] - centre[1]) < neighbour_cut)
            & (np.abs(xyz[:, 2] - centre[2]) < neighbour_cut)
        )
        near = near[near != i]
        points = centre + SPHERE * radii[i]
        if near.size:
            distances = np.linalg.norm(
                points[:, None, :] - xyz[near][None, :, :], axis=2
            )
            buried = (distances < radii[near][None, :]).any(axis=1)
        else:
            buried = np.zeros(points.shape[0], dtype=bool)
        fraction = 1.0 - buried.mean()
        total += 4.0 * math.pi * radii[i] ** 2 * fraction
    return total


def site_features(structure, resid):
    residue = structure["resid"] == resid
    if not residue.any():
        return None
    sulfur = residue & (structure["name"] == "SG")
    if not sulfur.any():
        return None
    xyz = structure["xyz"]
    centre = xyz[sulfur][0]
    backbone = np.isin(structure["name"], ("N", "CA", "C", "O"))
    sidechain = residue & ~backbone

    values = {}
    values["struct_plddt_site"] = float(structure["plddt"][residue].mean())
    values["struct_sidechain_sasa"] = sidechain_sasa(structure, sidechain)
    values["struct_relative_sasa"] = values["struct_sidechain_sasa"] / CYS_REFERENCE_SASA
    values["struct_sg_sasa"] = sidechain_sasa(structure, sulfur)

    distances = np.linalg.norm(xyz - centre, axis=1)
    other = ~residue
    for radius in SHELLS:
        shell = other & (distances <= radius)
        tag = f"{int(radius)}a"
        values[f"struct_atoms_{tag}"] = float(shell.sum())
        resnames = structure["resname"][shell]
        unique_res = len({(r, i) for r, i in zip(resnames, structure["resid"][shell])})
        values[f"struct_residues_{tag}"] = float(unique_res)
        if shell.any():
            values[f"struct_charge_{tag}"] = float(sum(CHARGE[r] for r in resnames)
                                                   / max(1, unique_res))
            values[f"struct_hydropathy_{tag}"] = float(
                np.mean([HYDROPATHY.get(r, 0.0) for r in resnames])
            )
            values[f"struct_frac_polar_{tag}"] = float(
                np.mean([r in ("SER", "THR", "ASN", "GLN", "TYR", "HIS", "CYS")
                         for r in resnames])
            )
            values[f"struct_frac_aromatic_{tag}"] = float(
                np.mean([r in ("PHE", "TRP", "TYR", "HIS") for r in resnames])
            )
        else:
            for key in ("charge", "hydropathy", "frac_polar", "frac_aromatic"):
                values[f"struct_{key}_{tag}"] = 0.0
        values[f"struct_plddt_{tag}"] = float(
            structure["plddt"][shell].mean()) if shell.any() else 0.0

    # positive and negative charge centres weighted by inverse distance
    for label, selector in (
        ("positive", ("LYS", "ARG")), ("negative", ("ASP", "GLU")), ("histidine", ("HIS",))
    ):
        hits = other & np.isin(structure["resname"], selector) & (structure["element"] == "N") \
            if label == "histidine" else other & np.isin(structure["resname"], selector)
        if hits.any():
            d = distances[hits]
            d = d[d > 0.1]
            values[f"struct_{label}_inverse_distance"] = float(np.sum(1.0 / d)) if d.size else 0.0
            values[f"struct_{label}_min_distance"] = float(d.min()) if d.size else 99.0
        else:
            values[f"struct_{label}_inverse_distance"] = 0.0
            values[f"struct_{label}_min_distance"] = 99.0

    # hydrogen-bond partners to the thiol sulfur
    partners = 0
    closest = 99.0
    for i in np.flatnonzero(other):
        name = structure["name"][i]
        resname = structure["resname"][i]
        if name in DONOR_ACCEPTOR.get(resname, ()) or name in ("N", "O"):
            d = distances[i]
            if d <= 4.0:
                partners += 1
            closest = min(closest, d)
    values["struct_hbond_partners_4a"] = float(partners)
    values["struct_nearest_polar_atom"] = float(closest)

    # nearest cysteine sulfur in space, and the disulfide-range flag
    sulfurs = (structure["name"] == "SG") & other
    if sulfurs.any():
        d = distances[sulfurs]
        values["struct_nearest_cys_sg"] = float(d.min())
        values["struct_cys_sg_within_8a"] = float((d <= 8.0).sum())
        values["struct_disulfide_range"] = 1.0 if d.min() <= 2.5 else 0.0
    else:
        values["struct_nearest_cys_sg"] = 99.0
        values["struct_cys_sg_within_8a"] = 0.0
        values["struct_disulfide_range"] = 0.0

    # depth: distance from the sulfur to the protein centroid, and to the surface
    centroid = xyz.mean(axis=0)
    values["struct_distance_to_centroid"] = float(np.linalg.norm(centre - centroid))
    values["struct_radius_of_gyration"] = float(
        np.sqrt(np.mean(np.sum((xyz - centroid) ** 2, axis=1)))
    )
    values["struct_relative_depth"] = (
        values["struct_distance_to_centroid"] / max(1e-6, values["struct_radius_of_gyration"])
    )
    return values


def main():
    args = parse_args()
    started = time.time()
    X, y, folds, proteins, components, positions, meta = load_frozen()
    del X, y, folds, components

    order = defaultdict(list)
    for index, (accession, position) in enumerate(zip(proteins, positions)):
        order[accession].append((index, int(position)))
    accessions = sorted(order)
    if args.limit:
        accessions = accessions[: args.limit]

    names = None
    matrix = None
    missing_structure = []
    missing_site = 0
    done = 0
    for accession in accessions:
        path = PDB_DIR / f"AF-{accession}-F1-model_v6.pdb"
        if not path.exists():
            missing_structure.append(accession)
            continue
        structure = read_pdb(path)
        if structure is None:
            missing_structure.append(accession)
            continue
        for index, position in order[accession]:
            values = site_features(structure, position)
            if values is None:
                missing_site += 1
                continue
            if names is None:
                names = sorted(values)
                matrix = np.full((len(proteins), len(names)), np.nan, dtype=np.float32)
            matrix[index] = [values[name] for name in names]
        done += 1
        if done % 200 == 0:
            print(json.dumps({"proteins_done": done, "of": len(accessions),
                              "minutes": round((time.time() - started) / 60, 1)}), flush=True)

    if names is None:
        raise SystemExit("No structures were parsed")
    covered = np.isfinite(matrix).any(axis=1)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT / "structural_chemistry_features.npz",
        X=matrix, names=np.asarray([f"S:{n}" for n in names]),
        proteins=proteins, positions=positions, covered=covered,
    )
    audit = {
        "completed": True,
        "elapsed_minutes": round((time.time() - started) / 60, 2),
        "host": platform.node(),
        "source": "AlphaFold DB v6 models of the cohort proteins",
        "n_features": len(names),
        "feature_names": [f"S:{n}" for n in names],
        "n_sites": int(len(proteins)),
        "n_sites_with_structure": int(covered.sum()),
        "coverage": round(float(covered.mean()), 4),
        "n_proteins_parsed": done,
        "n_proteins_without_model": len(missing_structure),
        "n_sites_without_a_cysteine_sulfur_in_the_model": missing_site,
        "labels_used": False,
        "sasa_method": f"Shrake-Rupley, {SPHERE_POINTS} sphere points, probe {PROBE} A",
        "shells_angstrom": list(SHELLS),
        "relative_sasa_reference": CYS_REFERENCE_SASA,
        "limits": [
            "AlphaFold models are predictions; pLDDT is carried per site and per shell so low-confidence regions stay visible",
            "single conformer, so no dynamics and no alternative rotamers",
            "no explicit pKa calculation; the electrostatic descriptors are geometric proxies",
            "a residue numbering mismatch between the model and the cohort sequence leaves the site uncovered rather than guessed",
        ],
        "versions": {"python": sys.version, "numpy": np.__version__},
        "input_hashes": {
            "inputs/primary_site_folds.csv": sha256(INPUTS / "primary_site_folds.csv"),
        },
    }
    write_json(OUT / "structural_chemistry_audit.json", audit)
    print(json.dumps({k: v for k, v in audit.items()
                      if k not in ("feature_names", "input_hashes", "limits")}, indent=2))


if __name__ == "__main__":
    main()
