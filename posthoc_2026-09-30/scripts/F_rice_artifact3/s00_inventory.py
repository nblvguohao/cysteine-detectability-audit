"""Step 0 (POST HOC): inventory of accessions and of reference-proteome coverage.

Writes results/F_rice_artifact3/s00_inventory.json and s00_accessions_to_fetch.txt (accessions that
appear in the quantification report or the peptide table and are absent from UP000059680 2026_03).
"""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

pro = read_tsv(PRO)
pep = read_tsv(PEP)
ref = read_fasta_any(FASTA_REF)

lead_group = []           # leading accession of every group
all_group = set()
for p in pro:
    a = split_accs(p.get("PG.ProteinAccessions"))
    if a:
        lead_group.append(a[0])
        all_group.update(a)
lrp = {(r.get("Leading razor protein") or "").strip() for r in pep} - {""}
lrp_real = {a for a in lrp if not is_decoy_or_contaminant(a)}
prot_col = set()
for r in pep:
    for a in split_accs(r.get("Proteins")):
        if not is_decoy_or_contaminant(a):
            prot_col.add(a)

# cohort proteins of the stored self-audit site table (read-only repo product)
import pandas as pd
sc = pd.read_csv(REPO_SITE_SCORES, encoding="utf-8-sig", low_memory=False)
cohort = set(sc.loc[sc.dataset == "rice_PXD072089", "accession"])


def cov(s):
    s = set(s)
    return {"n": len(s), "in_UP000059680_2026_03": len(s & set(ref)),
            "fraction": round(len(s & set(ref)) / len(s), 4) if s else None}


inv = {
    "analysis_label": "POST HOC revision analysis 2026-09-30 (not registered)",
    "n_reference_entries": len(ref),
    "protein_groups_rows": len(pro),
    "coverage": {
        "group_leading_accessions": cov(lead_group),
        "all_group_accessions": cov(all_group),
        "peptide_table_leading_razor_real": cov(lrp_real),
        "peptide_table_proteins_column_real": cov(prot_col),
        "self_audit_cohort_640": cov(cohort),
    },
    "n_leading_razor_including_REV_CON": len(lrp),
}
need = sorted((set(lead_group) | all_group | lrp_real | prot_col | cohort) - set(ref))
inv["n_accessions_to_fetch"] = len(need)
(RES / "s00_accessions_to_fetch.txt").write_text("\n".join(need) + "\n", encoding="utf-8")
write_json(RES / "s00_inventory.json", inv)
print(json.dumps(inv, indent=1))
