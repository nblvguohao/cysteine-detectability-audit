"""Convert a MaxQuant ABE-style (acyl-biotin exchange / capture-and-alkylate) search into a Cys-Audit input.

Written for PXD063463 (Kim et al., Proteomics 2026; mouse macrophages, S-palmitoylation by ABE; NEM blocks free
thiols, captured thiols are alkylated with iodoacetamide, so Carbamidomethyl (C) is the site readout and
NEM (C) marks a cysteine seen unmodified). Works for any MaxQuant txt folder with the same variable mods.

Definitions (fixed before the conversion was run on data):
  arm        evidence rows whose Experiment name ends with f"_Hyd_{H}_{ENZYME}" (H = P or N, ENZYME as in the
             experiment names: Trypsin, AspN, CT, GluC)
  filters    Reverse != '+', Potential contaminant != '+', protein not starting with 'CON__' or 'REV__'
  protein    the peptide's Leading razor protein (from peptides_ABE.txt via Peptide ID); written as the UniProt
             accession (second |-field)
  observed   every cysteine of every peptide with >= 1 evidence row in the arm -> detected = 1
  positive   a cysteine carrying Carbamidomethyl (C) with localisation probability >= 0.75 in >= 1 evidence row
             of the arm (probability read from the 'Carbamidomethyl (C) Probabilities' string) -> label = 1
  n_cys_peptide  smallest cysteine count over the arm's peptides that observed the cysteine
  peptide_mass   'Mass' of the highest-Score evidence row that supports a positive
  abundance  proteinGroups_Global.txt 'Intensity' (> 0), filters as above; a protein in several groups gets the
             median of their intensities
Output columns: protein, position, label, detected, n_cys_peptide, abundance, peptide_mass
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys

LOC_MIN = 0.75
PROB = re.compile(r"C\(([0-9.eE+-]+)\)")


def acc(tok):
    parts = tok.split("|")
    return parts[1] if len(parts) >= 3 else tok


def bad_protein(tok):
    return (not tok) or tok.startswith("CON__") or tok.startswith("REV__")


def cys_probs(prob_string, sequence):
    """Return {0-based index in sequence: probability} for cysteines annotated in a MaxQuant probability string."""
    out = {}
    i = 0
    j = 0
    s = prob_string
    while j < len(s):
        ch = s[j]
        if ch == "(":
            k = s.index(")", j)
            if i > 0 and sequence[i - 1] == "C":
                out[i - 1] = float(s[j + 1:k])
            j = k + 1
            continue
        i += 1
        j += 1
    return out


def convert(data_dir, enzyme, hyd, out, global_groups=True):
    csv.field_size_limit(10 ** 9)
    peps = {}
    with open(f"{data_dir}/peptides_ABE.txt", encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Reverse"] == "+" or r["Potential contaminant"] == "+" or bad_protein(r["Leading razor protein"]):
                continue
            peps[r["id"]] = (r["Sequence"], int(r["Start position"]), acc(r["Leading razor protein"]))
    suffix = f"_Hyd_{hyd}_{enzyme}"
    observed, positive, ncys, best = set(), set(), {}, {}
    n_ev = 0
    with open(f"{data_dir}/evidence_ABE.txt", encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["Experiment"].endswith(suffix):
                continue
            if r["Reverse"] == "+" or r["Potential contaminant"] == "+":
                continue
            p = peps.get(r["Peptide ID"])
            if p is None:
                continue
            n_ev += 1
            seq, start, prot = p
            nc = seq.count("C")
            for i, ch in enumerate(seq):
                if ch == "C":
                    key = (prot, start + i)
                    observed.add(key)
                    ncys[key] = min(ncys.get(key, 99), nc)
            if r["Carbamidomethyl (C)"] not in ("", "0"):
                probs = cys_probs(r["Carbamidomethyl (C) Probabilities"], seq)
                score = float(r["Score"] or 0)
                for i, pr in probs.items():
                    if pr >= LOC_MIN:
                        key = (prot, start + i)
                        positive.add(key)
                        if key not in best or score > best[key][0]:
                            best[key] = (score, float(r["Mass"]))
    ab = {}
    if global_groups:
        tmp = {}
        with open(f"{data_dir}/proteinGroups_Global.txt", encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["Reverse"] == "+" or r["Potential contaminant"] == "+":
                    continue
                try:
                    inten = float(r["Intensity"])
                except ValueError:
                    continue
                if inten <= 0:
                    continue
                for t in r["Protein IDs"].split(";"):
                    if not bad_protein(t):
                        tmp.setdefault(acc(t), []).append(inten)
        ab = {k: statistics.median(v) for k, v in tmp.items()}
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["protein", "position", "label", "detected", "n_cys_peptide", "abundance", "peptide_mass"])
        for prot, pos in sorted(observed):
            a = ab.get(prot)
            m = best.get((prot, pos))
            w.writerow([prot, pos, int((prot, pos) in positive), 1, ncys[(prot, pos)],
                        "" if a is None else repr(a), "" if m is None else repr(m[1])])
    return {"evidence_rows_in_arm": n_ev, "observed": len(observed), "positive": len(positive & observed),
            "positive_not_observed": len(positive - observed), "proteins_with_abundance": len(ab)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--enzyme", required=True, choices=["Trypsin", "AspN", "CT", "GluC"])
    ap.add_argument("--hydroxylamine", default="P", choices=["P", "N"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    print(convert(a.data_dir, a.enzyme, a.hydroxylamine, a.out))


if __name__ == "__main__":
    sys.exit(main())
