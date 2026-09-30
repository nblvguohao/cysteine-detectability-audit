"""Minimum heavy-atom distance from each candidate Cys to the partner chain.

Reads the AF3 top-ranked WRKY x PP2A-B56 model and reports, for every Cys in
both chains, how far it sits from the other protein, next to the v2 chemistry
model's within-protein rank.  Pure standard library so it runs with the system
interpreter.  Output: results/v2_interface_proximity_vs_rank.csv
"""
from __future__ import annotations
import csv, json, math, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
PDB = ROOT / "interaction/af3_results/wrky_plus_pp2a_b56/best_model_0.pdb"
SITES = ROOT / "results/v2_user_answer_sites.csv"
OUT = ROOT / "results/v2_interface_proximity_vs_rank.csv"
AUDIT = ROOT / "results/v2_interface_proximity_audit.json"

CHAIN_TO_PROTEIN = {"A": "XP_004233015.1", "B": "Solyc12g006920.2.1"}


def load_chains(path):
    chains = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        index = int(line[22:26])
        name = line[17:20].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        residues = chains.setdefault(chain, {})
        residue = residues.setdefault(index, {"name": name, "atoms": []})
        residue["atoms"].append(xyz)
    return chains


def min_distance(residue, other_chain):
    best = math.inf
    for x in residue["atoms"]:
        for other in other_chain.values():
            for y in other["atoms"]:
                d = math.dist(x, y)
                if d < best:
                    best = d
    return best


def main():
    chains = load_chains(PDB)
    ranks = {}
    with SITES.open(encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            ranks[(row["protein"], int(row["position"]))] = row
    rows = []
    for chain, residues in chains.items():
        partner = chains["B" if chain == "A" else "A"]
        protein = CHAIN_TO_PROTEIN[chain]
        for index, residue in sorted(residues.items()):
            if residue["name"] != "CYS":
                continue
            site = ranks.get((protein, index), {})
            rows.append({
                "protein": protein,
                "chain": chain,
                "position": index,
                "min_heavy_atom_distance_to_partner_angstrom": round(min_distance(residue, partner), 2),
                "v2_chem_rank": site.get("v2_chem_rank", ""),
                "v2_chem_score": site.get("v2_chem_score", ""),
                "v2_full_rank": site.get("v2_full_rank", ""),
                "evidence_grade": site.get("evidence_grade", ""),
                "assayable_standard_trypsin": site.get("assayable_standard_trypsin", ""),
            })
    rows.sort(key=lambda r: (r["chain"], r["min_heavy_atom_distance_to_partner_angstrom"]))
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    AUDIT.write_text(json.dumps({
        "model": str(PDB.relative_to(ROOT)),
        "note": "distances come from the single top-ranked AF3 diffusion sample; "
                "they are not filtered by cross-chain PAE, so proximity outside "
                "the P13-E16 anchor region sits in the low-confidence part of the model",
        "n_cys_chain_A": sum(1 for r in rows if r["chain"] == "A"),
        "n_cys_chain_B": sum(1 for r in rows if r["chain"] == "B"),
    }, indent=2), encoding="utf-8")
    for row in rows:
        print(row["protein"], "C%d" % row["position"],
              row["min_heavy_atom_distance_to_partner_angstrom"], "A",
              "chem rank", row["v2_chem_rank"])


if __name__ == "__main__":
    main()
