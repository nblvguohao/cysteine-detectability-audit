"""POST HOC revision analysis (2026-09-30), item C_pxd063463_specific - step 0: exploration only.

Loads the eight PXD063463 arm tables, checks every (protein, position) against the UniProt 2026_03 mouse FASTA and
prints overlap counts between HydP and HydN arms. Writes nothing except stdout.
"""
import gzip
import os
import sys

import pandas as pd

W = "/path/to/local/_cys_repo_work/public/revision_2026-09-30"
IN = f"{W}/inputs/phase4_inputs"
FA = f"{W}/external/UP000000589_10090.fasta.gz"
ARMS = ["Trypsin", "AspN", "CT", "GluC"]


def read_fasta_gz(path):
    seqs, cur, buf = {}, None, []

    def flush():
        if cur is not None:
            s = "".join(buf).upper().rstrip("*")
            for k in cur:
                seqs.setdefault(k, s)
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                flush()
                tok = line[1:].split()[0] if line[1:].split() else ""
                keys = [tok]
                parts = tok.split("|")
                if len(parts) >= 3:
                    keys.insert(0, parts[1])
                cur, buf = keys, []
            elif line:
                buf.append(line.strip())
        flush()
    return seqs


seqs = read_fasta_gz(FA)
print("fasta entries (keys)", len(seqs))
tabs = {}
for a in ARMS:
    for h in ("HydP", "HydN"):
        d = pd.read_csv(f"{IN}/pxd063463_{a}_{h}.tsv", sep="\t")
        d["in_fasta"] = d.protein.isin(seqs)
        d["is_c"] = [p in seqs and pos <= len(seqs[p]) and seqs[p][pos - 1] == "C" for p, pos in zip(d.protein, d.position)]
        tabs[(a, h)] = d
        print(a, h, "rows", len(d), "label1", int(d.label.sum()), "proteins", d.protein.nunique(),
              "absent", int((~d.in_fasta).sum()), "not_C", int((d.in_fasta & ~d.is_c).sum()),
              "abund_cov", round(d.abundance.notna().mean(), 3), "detected_all1", bool((d.detected == 1).all()))
for a in ARMS:
    P, N = tabs[(a, "HydP")], tabs[(a, "HydN")]
    kp = set(zip(P.protein, P.position))
    kn = set(zip(N.protein, N.position))
    camP = set(zip(P.protein[P.label == 1], P.position[P.label == 1]))
    camN = set(zip(N.protein[N.label == 1], N.position[N.label == 1]))
    nocamN = kn - camN
    print(f"{a}: HydP id {len(kp)} CAM {len(camP)} | HydN id {len(kn)} CAM {len(camN)} | HydN id not in HydP {len(kn - kp)}"
          f" | HydP CAM also CAM in HydN {len(camP & camN)} | strict (CAM P & id N without CAM) {len(camP & nocamN)}"
          f" | lenient (CAM P & not CAM N) {len(camP - camN)} | HydN CAM not CAM in HydP {len(camN - camP)}")
    print("   proteins HydP", P.protein.nunique(), "HydN", N.protein.nunique(), "shared", len(set(P.protein) & set(N.protein)))
# cross-arm overlaps (HydP identified)
ids = {a: set(zip(tabs[(a, 'HydP')].protein, tabs[(a, 'HydP')].position)) for a in ARMS}
prots = {a: set(tabs[(a, 'HydP')].protein) for a in ARMS}
for a in ARMS[1:]:
    print(a, "cys shared with trypsin HydP", len(ids[a] & ids['Trypsin']), "of", len(ids[a]),
          "| proteins shared", len(prots[a] & prots['Trypsin']), "of", len(prots[a]))
