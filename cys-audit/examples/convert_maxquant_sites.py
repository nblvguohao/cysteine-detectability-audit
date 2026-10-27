"""Convert a MaxQuant txt folder with a per-modification site table into a Cys-Audit input.

Reproduces the PRIMARY definition of the analysis tree's public Artefact 4 analysis
(scripts/analyse_public_artefact4_2026-09-21.py, PXD048216):
  observed   every cysteine of peptides.txt rows with Reverse != '+', Potential contaminant != '+', a non-empty
             Leading razor protein; position = Start position + index; protein = Leading razor protein (kept verbatim)
  positive   rows of '<MOD>Sites.txt' (default 'Carbamidomethyl (C)Sites.txt') with Localization prob >= 0.75,
             keyed by (Protein, Position)
Output columns: protein, position, label, detected
"""
from __future__ import annotations

import argparse
import csv
import sys

LOC_MIN = 0.75


def convert(txt_dir, site_file, out):
    csv.field_size_limit(10 ** 9)
    observed = set()
    with open(f"{txt_dir}/peptides.txt", encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("Reverse") == "+" or r.get("Potential contaminant") == "+" or not r["Leading razor protein"]:
                continue
            s0 = int(r["Start position"])
            for i, ch in enumerate(r["Sequence"]):
                if ch == "C":
                    observed.add((r["Leading razor protein"], s0 + i))
    sites = set()
    with open(f"{txt_dir}/{site_file}", encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                lp = float(r["Localization prob"])
            except (ValueError, KeyError):
                continue
            if lp >= LOC_MIN:
                sites.add((r["Protein"], int(float(r["Position"]))))
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["protein", "position", "label", "detected"])
        for prot, pos in sorted(observed):
            w.writerow([prot, pos, int((prot, pos) in sites), 1])
    return {"observed": len(observed), "positive": len(sites & observed), "sites_not_observed": len(sites - observed)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--txt-dir", required=True)
    ap.add_argument("--site-file", default="Carbamidomethyl (C)Sites.txt")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    print(convert(a.txt_dir, a.site_file, a.out))


if __name__ == "__main__":
    sys.exit(main())
