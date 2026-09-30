"""POST-HOC (not pre-registered): fill species for the claim cohorts whose builders opened no reference-proteome
file (species_tested = 'undetermined' in results/claim_audit_2026-09-22.csv).

Rule (fixed before running): for each site-level cohort exported by the Phase 5 runner
(external/intake/phase5_claim_inputs_20260922/<claim>.tsv), take its distinct accessions and count membership in
each local reference proteome (external/proteomes: hsa, mmu, ath, osa as FASTA; pvu, tgo as phase2 TSV). Species =
the proteome holding >= 0.90 of the accessions; otherwise 'ambiguous' with the shares reported. Protein-level
cohorts were not exported and are left as they are. Output: results/claim_audit_species_posthoc_2026-09-22.csv
and _audit.json. The pre-registered matrix is not modified. Standard library only.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = os.path.join(ROOT, "external", "intake", "phase5_claim_inputs_20260922")
PROT = os.path.join(ROOT, "external", "proteomes")
MATRIX = os.path.join(ROOT, "results", "claim_audit_2026-09-22.csv")
OUT = os.path.join(ROOT, "results", "claim_audit_species_posthoc_2026-09-22.csv")
AUD = os.path.join(ROOT, "results", "claim_audit_species_posthoc_2026-09-22_audit.json")
NAMES = {"hsa": "Homo sapiens", "mmu": "Mus musculus", "ath": "Arabidopsis thaliana", "osa": "Oryza sativa",
         "pvu": "Phaseolus vulgaris", "tgo": "Toxoplasma gondii"}


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def accs_fasta(p):
    out = set()
    with gzip.open(p, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].split()[0].split("|")
                out.add(h[1] if len(h) > 2 else h[0])
    return out


def accs_tsv(p):
    with gzip.open(p, "rt", encoding="utf-8") as fh:
        return {r["Entry"].strip() for r in csv.DictReader(fh, delimiter="\t") if r.get("Entry")}


def main():
    for p in (OUT, AUD):
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    sets = {}
    for k in ("hsa", "ath", "osa"):
        sets[k] = accs_fasta(os.path.join(PROT, f"{k}.fasta.gz"))
    for k in ("mmu", "pvu", "tgo"):
        sets[k] = accs_tsv(os.path.join(PROT, f"phase2_{k}_uniprot.tsv.gz"))
    rows = []
    for r in csv.DictReader(open(MATRIX, encoding="utf-8")):
        cid = r["claim_id"]
        tsv = os.path.join(IN, f"{cid}.tsv")
        if not os.path.exists(tsv):
            rows.append({"claim_id": cid, "species_tested_prereg": r["species_tested"], "species_posthoc": "",
                         "n_accessions": "", "shares": "", "note": "not exported (protein-level or no baseline)"})
            continue
        accs = {x["protein"] for x in csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t")}
        shares = {k: round(len(accs & s) / len(accs), 4) for k, s in sets.items()}
        best = max(shares, key=shares.get)
        sp = NAMES[best] if shares[best] >= 0.90 else "ambiguous"
        rows.append({"claim_id": cid, "species_tested_prereg": r["species_tested"], "species_posthoc": sp,
                     "n_accessions": len(accs), "shares": json.dumps(shares), "note": "post-hoc accession membership"})
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    conflicts = [x["claim_id"] for x in rows if x["species_posthoc"] and x["species_tested_prereg"] not in ("undetermined", "", x["species_posthoc"])]
    audit = {"script": "scripts/posthoc_species_claim_cohorts_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
             "posthoc": True, "matrix_sha256": sha(MATRIX), "conflicts_with_prereg_species": conflicts,
             "outputs": {os.path.relpath(OUT, ROOT): sha(OUT)}}
    json.dump(audit, open(AUD, "w"), indent=1)
    for x in rows:
        if x["species_posthoc"]:
            print(x["claim_id"], x["species_tested_prereg"], "->", x["species_posthoc"], x["shares"])
    print("conflicts", conflicts)


if __name__ == "__main__":
    main()
