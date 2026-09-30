"""Short-linear-motif scan for the interaction and phosphorylation questions.

This is a transparent regular-expression scan against published consensus
definitions.  It is not a trained predictor and it produces hypotheses to test,
not site assignments.  Every hit is reported with its sequence context so it can
be judged by eye.

Motifs scanned
--------------
* ``MAPK_phosphosite``   S/T followed by proline, the minimal MAPK consensus.
* ``MAPK_phosphosite_strong`` P-x-S/T-P, the stronger proline-directed form.
* ``MAPK_docking_D_motif`` a basic patch followed by a hydrophobic L/I/V-x-L/I/V,
  the docking groove ligand used by MAPK substrates and regulators.
* ``B56_LxxIxE``         the PP2A-B56 regulatory-subunit binding motif.
* ``SnRK2_like``         R-x-x-S/T, the basophilic consensus shared by SnRK2 and
  several calcium-dependent kinases.
* ``CDPK_like``          basic residue at -3 and hydrophobic at -5.
* ``CK2_like``           S/T-x-x-D/E, acidophilic casein-kinase-2 consensus.
* ``PKA_like``           R-R/K-x-S/T.
"""
from __future__ import annotations

import json
import re
import sys
import time

from common import INPUTS, RESULTS, ROOT, sha256, write_csv, write_json
import score_experimental_candidates as v1

# Which protein a motif is informative *in*.  A docking motif is carried by the
# ligand, not by the receptor that reads it, so a hit inside the receptor itself
# is not interpretable and is labelled as such.
MOTIF_ROLE = {
    "MAPK_docking_D_motif": "carried by MAPK substrates and regulators; a hit inside the kinase core is usually spurious",
    "B56_LxxIxE": "carried by PP2A-B56 substrates; a hit inside the B56 subunit itself is not interpretable",
    "MAPK_phosphosite": "site on a MAPK substrate",
    "MAPK_phosphosite_strong": "site on a MAPK substrate",
    "SnRK2_like": "site on a basophilic-kinase substrate",
    "CDPK_like": "site on a basophilic-kinase substrate",
    "CK2_like": "site on a casein-kinase-2 substrate",
    "PKA_like": "site on a basophilic-kinase substrate",
    "kinase_activation_loop_TxY": "the MAPK's own activation loop, phosphorylated by an upstream MAP2K",
}
RECEPTOR_PROTEIN = {
    "B56_LxxIxE": "Solyc12g006920.2.1",
}

MOTIFS = [
    ("MAPK_phosphosite", r"(?=([ST])P)", "S/T-P", "phosphorylation"),
    ("MAPK_phosphosite_strong", r"(?=P.([ST])P)", "P-x-S/T-P", "phosphorylation"),
    ("B56_LxxIxE", r"(?=(L..[IL].E))", "LxxIxE", "interaction"),
    ("SnRK2_like", r"(?=R..([ST]))", "R-x-x-S/T", "phosphorylation"),
    ("CDPK_like", r"(?=[ILVFM].[KR]..([ST]))", "hydrophobic-x-basic-x-x-S/T", "phosphorylation"),
    ("CK2_like", r"(?=([ST])..[DE])", "S/T-x-x-D/E", "phosphorylation"),
    ("PKA_like", r"(?=R[RK].([ST]))", "R-R/K-x-S/T", "phosphorylation"),
]
DOCKING = re.compile(r"(?=([KR]{2,3}.{1,6}[LIV].[LIV]))")
ACTIVATION_LOOP = re.compile(r"(?=(T[DE]Y))")
KINASE_CORE = re.compile(r"HRDLKP|HRDIKA|DFG[LM]AR")
TARGETS = ("XP_004233015.1", "Solyc06g068990.4.1", "Solyc12g006920.2.1")


def context(sequence, start, stop, flank=6):
    left = max(0, start - flank)
    right = min(len(sequence), stop + flank)
    return (
        sequence[left:start].lower()
        + sequence[start:stop]
        + sequence[stop:right].lower()
    )


def main():
    started = time.time()
    sequences = {
        name: value.strip().upper().rstrip("*")
        for name, value in v1.read_fasta(INPUTS / "candidate_proteins.fasta").items()
    }
    rows = []
    for name in TARGETS:
        sequence = sequences[name]
        for motif, pattern, consensus, purpose in MOTIFS:
            for match in re.finditer(pattern, sequence):
                group = match.group(1)
                start = match.start(1)
                residue = sequence[start]
                rows.append({
                    "protein": name,
                    "motif": motif,
                    "consensus": consensus,
                    "purpose": purpose,
                    "position": start + 1,
                    "residue": residue,
                    "match": group,
                    "context_13aa": context(sequence, start, start + len(group)),
                    "motif_role": MOTIF_ROLE.get(motif, ""),
                    "interpretable_in_this_protein": RECEPTOR_PROTEIN.get(motif) != name,
                })
        kinase_core = [m.start() for m in KINASE_CORE.finditer(sequence)]
        for match in DOCKING.finditer(sequence):
            group = match.group(1)
            start = match.start(1)
            inside_core = any(abs(start - c) < 120 for c in kinase_core)
            rows.append({
                "protein": name,
                "motif": "MAPK_docking_D_motif",
                "consensus": "(K/R)2-3-x(1-6)-L/I/V-x-L/I/V",
                "purpose": "interaction",
                "position": start + 1,
                "residue": sequence[start],
                "match": group,
                "context_13aa": context(sequence, start, start + len(group)),
                "motif_role": MOTIF_ROLE["MAPK_docking_D_motif"],
                "interpretable_in_this_protein": not inside_core,
            })
        for match in ACTIVATION_LOOP.finditer(sequence):
            group = match.group(1)
            start = match.start(1)
            if not any(abs(start - c) < 120 for c in kinase_core):
                continue
            rows.append({
                "protein": name,
                "motif": "kinase_activation_loop_TxY",
                "consensus": "T-D/E-Y inside the kinase domain",
                "purpose": "phosphorylation",
                "position": start + 1,
                "residue": sequence[start],
                "match": group,
                "context_13aa": context(sequence, start, start + len(group)),
                "motif_role": MOTIF_ROLE["kinase_activation_loop_TxY"],
                "interpretable_in_this_protein": True,
            })
    rows.sort(key=lambda r: (r["protein"], r["motif"], r["position"]))
    write_csv(RESULTS / "v2_user_motif_scan.csv", rows)

    counts = {}
    for row in rows:
        counts.setdefault(row["protein"], {}).setdefault(row["motif"], 0)
        counts[row["protein"]][row["motif"]] += 1
    write_json(RESULTS / "v2_user_motif_audit.json", {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "method": "regular-expression scan against published consensus definitions",
        "not_a_trained_predictor": True,
        "interpretation": "each hit is a hypothesis to test by mutagenesis or by a dedicated predictor; motif density alone does not rank sites",
        "counts": counts,
        "sequence_lengths": {name: len(sequences[name]) for name in TARGETS},
        "versions": {"python": sys.version},
        "input_hashes": {
            str((INPUTS / "candidate_proteins.fasta").relative_to(ROOT)): sha256(
                INPUTS / "candidate_proteins.fasta"
            )
        },
    })
    print(json.dumps(counts, indent=1))
    for name in TARGETS:
        print("="*70)
        print(name)
        for row in rows:
            if row["protein"] == name and row["motif"] in (
                "MAPK_docking_D_motif", "B56_LxxIxE", "MAPK_phosphosite_strong",
                "MAPK_phosphosite", "kinase_activation_loop_TxY",
            ) and row["interpretable_in_this_protein"]:
                print(f"  {row['motif']:24s} pos {row['position']:>4d} {row['match']:<12s} {row['context_13aa']}")


if __name__ == "__main__":
    main()
