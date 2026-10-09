"""Step 5 (POST HOC): reconcile the rice cohort sizes.

(a) EXPLORATORY rule grid over the deposit's peptide table: how many cysteine sites does each simple
    counting rule give, and which rules land near the source publication's 1,691? The grid is
    reported in full; a near-match identifies a plausible counting rule, it does not establish
    the authors' rule (their pipeline is not deposited).
(b) Reconstruction of the manuscript's rice cohort (640 proteins, 5,598 cysteines, 817 reported
    sites; 594 cysteines / 80 sites excluded -> 5,004 / 737). The cohort table was built on
    2026-09-12 by a second-tree script (24_build_persulfidation_table.py, not in the repository;
    see s05 report) as PXD072089_mapped_rice.json; its site set survives in the repository's
    external/public_cohort_site_scores.csv (rice rows). Candidate rules are compared with that
    stored set, site by site.

Outputs: s05_rule_grid.csv, s05_cohort_rules.csv, s05_reconcile.json
"""
import itertools

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

TARGET_PUB = 1691
seqs, src = load_all_sequences()
ref = read_fasta_any(FASTA_REF)
pep = pd.read_csv(PEP, sep="\t", low_memory=False, dtype=str)
pep["C Count"] = pep["C Count"].astype(int)
pep["lrp"] = pep["Leading razor protein"].fillna("")
pep["start"] = pd.to_numeric(pep["Start position"], errors="coerce")
id_cols = [c for c in pep.columns if c.startswith("Identification type ")]
pep["n_runs_any"] = sum((pep[c].fillna("") != "").astype(int) for c in id_cols)
pep["n_runs_msms"] = sum((pep[c].fillna("") == "By MS/MS").astype(int) for c in id_cols)
pep["is_rev"] = (pep["Reverse"].fillna("") == "+") | pep.lrp.str.startswith("REV__")
pep["is_con"] = (pep["Potential contaminant"].fillna("") == "+") | pep.lrp.str.startswith("CON__")
cys = pep[pep["C Count"] > 0].copy()


def sites_lrp(df):
    out = set()
    for s, a, st in zip(df.Sequence, df.lrp, df.start):
        if np.isnan(st):
            continue
        for i, ch in enumerate(s):
            if ch == "C":
                out.add((a, int(st) + i))
    return out


def sites_any(df):
    out = set()
    for s, prots in zip(df.Sequence, df.Proteins.fillna("")):
        for a in split_accs(prots):
            if is_decoy_or_contaminant(a) or a not in seqs:
                continue
            ps = seqs[a]; k = ps.find(s)
            while k >= 0:
                for i, ch in enumerate(s):
                    if ch == "C":
                        out.add((a, k + 1 + i))
                k = ps.find(s, k + 1)
    return out


# ---------------------------------------------------------------- (a) rule grid -----
grid = []
for assign, cyssel, decoy, uniq, evid in itertools.product(
        ("leading_razor", "any_protein"), ("all_Cys_peptides", "single_Cys_peptides"),
        ("REV_CON_removed", "REV_CON_kept"), ("all", "unique_groups", "unique_proteins"),
        ("any_run", "MSMS_in_>=1_run", ">=2_runs")):
    df = cys
    if cyssel == "single_Cys_peptides":
        df = df[df["C Count"] == 1]
    if decoy == "REV_CON_removed":
        df = df[~df.is_rev & ~df.is_con]
    if uniq == "unique_groups":
        df = df[df["Unique (Groups)"] == "yes"]
    elif uniq == "unique_proteins":
        df = df[df["Unique (Proteins)"] == "yes"]
    if evid == "MSMS_in_>=1_run":
        df = df[df.n_runs_msms >= 1]
    elif evid == ">=2_runs":
        df = df[df.n_runs_any >= 2]
    S = sites_lrp(df) if assign == "leading_razor" else sites_any(df)
    grid.append({"protein_assignment": assign, "peptides": cyssel, "decoys_contaminants": decoy,
                 "uniqueness": uniq, "evidence": evid, "n_peptides": int(len(df)), "n_sites": len(S),
                 "n_proteins": len({a for a, _ in S}),
                 "rel_diff_to_1691": (len(S) - TARGET_PUB) / TARGET_PUB})
GRID = pd.DataFrame(grid).sort_values(["protein_assignment", "peptides", "decoys_contaminants", "uniqueness", "evidence"])
GRID.to_csv(RES / "s05_rule_grid.csv", index=False)

# ---------------------------------------------------------------- (b) cohort reconstruction -----
sc = pd.read_csv(REPO_SITE_SCORES, encoding="utf-8-sig", low_memory=False)
co = sc[sc.dataset == "rice_PXD072089"][["accession", "position", "observed_in_retrospective_dataset"]]
stored_sites = set(map(tuple, co[co.observed_in_retrospective_dataset == 1][["accession", "position"]].itertuples(index=False)))
stored_prot = set(co.accession)
real = cys[~cys.is_rev & ~cys.is_con]
single = real[real["C Count"] == 1]


def ncys(a):
    return seqs[a].count("C") if a in seqs else None


def summarize(name, S, note=""):
    P = {a for a, _ in S}
    return {"rule": name, "n_sites": len(S), "n_proteins": len(P),
            "n_cys_in_proteins": int(sum(ncys(a) or 0 for a in P)),
            "proteins_without_sequence": int(sum(1 for a in P if a not in seqs)),
            "site_intersection_with_stored": len(S & stored_sites),
            "stored_not_in_rule": len(stored_sites - S), "rule_not_in_stored": len(S - stored_sites),
            "protein_intersection_with_stored": len(P & stored_prot),
            "jaccard_sites": len(S & stored_sites) / len(S | stored_sites) if S else 0.0, "note": note}


def ge2(S):
    return {(a, q) for a, q in S if (ncys(a) or 0) >= 2}


def mixed(S):
    by = {}
    for a, q in S:
        by.setdefault(a, set()).add(q)
    return {(a, q) for a, q in S if (ncys(a) or 0) >= 2 and len(by[a]) < (ncys(a) or 0)}


def first_protein_with_seq(df, pool):
    """map each peptide to the FIRST accession of its Proteins column present in `pool`."""
    out = set()
    for s, prots in zip(df.Sequence, df.Proteins.fillna("")):
        for a in split_accs(prots):
            if a in pool:
                k = pool[a].find(s)
                if k >= 0:
                    for i, ch in enumerate(s):
                        if ch == "C":
                            out.add((a, k + 1 + i))
                    break
    return out


L1 = sites_lrp(single)
A1 = sites_any(single)
rows = [
    summarize("LRP_single_Cys", L1),
    summarize("LRP_single_Cys_proteins_ge2_Cys", ge2(L1)),
    summarize("LRP_single_Cys_proteins_ge2_Cys_mixed", mixed(L1), "protein keeps >=1 unobserved Cys"),
    summarize("ANY_single_Cys", A1),
    summarize("ANY_single_Cys_proteins_ge2_Cys", ge2(A1)),
    summarize("FIRSTREF_single_Cys", first_protein_with_seq(single, ref), "first accession of 'Proteins' found in UP000059680 2026_03"),
    summarize("FIRSTREF_single_Cys_ge2", ge2(first_protein_with_seq(single, ref))),
    summarize("FIRSTANY_single_Cys", first_protein_with_seq(single, seqs), "first accession of 'Proteins' with any sequence"),
    summarize("FIRSTANY_single_Cys_ge2", ge2(first_protein_with_seq(single, seqs))),
    summarize("LRP_all_Cys", sites_lrp(real)),
    summarize("LRP_all_Cys_ge2", ge2(sites_lrp(real))),
]
# stored sites: which peptides cover them?
cover = {}
for s, prots, c in zip(real.Sequence, real.Proteins.fillna(""), real["C Count"]):
    for a in split_accs(prots):
        if a in seqs:
            k = seqs[a].find(s)
            while k >= 0:
                for i, ch in enumerate(s):
                    if ch == "C":
                        cover.setdefault((a, k + 1 + i), set()).add(int(c))
                k = seqs[a].find(s, k + 1)
st_cov = {"stored_sites_covered_by_single_Cys_peptide_any_protein": sum(1 for x in stored_sites if 1 in cover.get(x, set())),
          "stored_sites_covered_only_by_multi_Cys_peptides": sum(1 for x in stored_sites if x in cover and 1 not in cover[x]),
          "stored_sites_not_covered_by_any_identified_peptide": sum(1 for x in stored_sites if x not in cover)}
stored_lrp = {a for a, _ in stored_sites}
st_cov["stored_proteins_that_are_leading_razor_of_a_single_Cys_peptide"] = len(stored_lrp & set(single.lrp))
st_cov["stored_proteins_min_cys"] = int(co.groupby("accession").size().min())
st_cov["stored_proteins_all_cys_observed"] = int((co.groupby("accession").observed_in_retrospective_dataset.mean() == 1).sum())
tot = {a: ncys(a) for a in stored_prot}
st_cov["stored_cohort_total_cys_from_sequences"] = int(sum(v for v in tot.values() if v is not None))
st_cov["stored_cohort_proteins_with_sequence"] = int(sum(v is not None for v in tot.values()))
per = co.groupby("accession").size()
st_cov["stored_cohort_per_protein_cys_equal_sequence_count"] = int(sum(1 for a in stored_prot if tot[a] is not None and tot[a] == per[a]))
nr = {a for a in stored_prot if a not in ref}
st_cov["absent_from_UP000059680_2026_03"] = {"proteins": len(nr), "cysteines": int(per[list(nr)].sum()),
                                             "stored_sites": int(sum(1 for a, _ in stored_sites if a in nr))}
st_cov["manuscript_exclusion"] = {"cysteines": 594, "sites": 80, "remaining_cysteines": 5004, "remaining_sites": 737}
# ---------------------------------------------------------------- (c) the filter chain -----
# Found post hoc by inspecting (b): the stored proteins are the canonical accessions of the leading
# razor proteins of single-Cys peptides (isoform suffix '-n' stripped), with >=2 Cys, whose
# accession was still ACTIVE in UniProtKB when the table was built (sequence retrievable); the
# excluded ones were deleted in UniProt 2026_02 ('Not part of a reference proteome').
fs_rows = list(csv.DictReader(open(FETCHED / "fetched_sequences.tsv", encoding="utf-8"), delimiter="\t"))
active_2026_03 = set(ref) | {r["requested_accession"] for r in fs_rows if r["source"] == "uniprot_rest_active"}


def canonical_sites(df):
    out = set()
    for s, a in zip(df.Sequence, df.lrp):
        c = a.split("-")[0]
        if c not in seqs:
            continue
        ps = seqs[c]; k = ps.find(s)
        while k >= 0:
            for i, ch in enumerate(s):
                if ch == "C":
                    out.add((c, k + 1 + i))
            k = ps.find(s, k + 1)
    return out


chain = []
def step(name, S, note):
    P = {a for a, _ in S}
    chain.append({"step": name, "n_sites": len(S), "n_proteins": len(P),
                  "n_cys_in_proteins": int(sum(ncys(a) or 0 for a in P)),
                  "equals_stored_sites": S == stored_sites, "equals_stored_proteins": P == stored_prot,
                  "intersection_with_stored": len(S & stored_sites), "note": note})
    return S

s0 = step("0_all_Cys_peptides_every_C_LRP_incl_REV_CON", sites_lrp(cys), "every C of every identified Cys-containing peptide, on its leading razor protein, decoys/contaminants kept")
s1 = step("1_REV_CON_removed", sites_lrp(real), "decoys and contaminants removed")
s2 = step("2_single_Cys_peptides", sites_lrp(single), "peptides with exactly one C (site assignment unambiguous)")
s3 = step("3_isoform_LRP_to_canonical", canonical_sites(single), "isoform accessions (-n) replaced by the canonical accession, peptide re-mapped")
s4 = step("4_proteins_ge2_Cys", ge2(s3), "protein has >=2 cysteines (within-protein ranking needs >=2 candidates)")
s5 = step("5_accession_active_in_UniProtKB", {(a, q) for a, q in s4 if a in active_2026_03},
          "accession still active (sequence retrievable); the others were deleted in UniProt 2026_02")
s6 = step("6_protein_in_reference_proteome_2026_03", {(a, q) for a, q in s5 if a in ref},
          "protein present in UP000059680 (2026_03); the manuscript's exclusion used an earlier release (594 Cys / 80 sites excluded, 5,004 / 737 left)")
pd.DataFrame(chain).to_csv(RES / "s05_filter_chain.csv", index=False)
# also the all-Cys-peptide analogue of the stored cohort (what the same chain gives without step 2)
s3m = canonical_sites(real)
alt = {"all_Cys_peptides_canonical_ge2_active": len({(a, q) for a, q in ge2(s3m) if a in active_2026_03}),
       "proteins": len({a for a, q in ge2(s3m) if a in active_2026_03})}

# kept (active) vs lost (deleted in 2026_02) candidate proteins: are the lost ones different?
cand_prot = {a for a, _ in s4}
by_prot_sites = {}
for a, q in s4:
    by_prot_sites[a] = by_prot_sites.get(a, 0) + 1
G = pd.read_csv(RES / "s02_groups.csv")
acc2ibaq = {}
for accs_, ib in zip(G.accessions, G.ibaq_median):
    for a_ in accs_.split(";"):
        acc2ibaq.setdefault(a_, ib)
kv = []
for a in sorted(cand_prot):
    kv.append({"accession": a, "kept_in_stored_cohort": a in active_2026_03, "length": len(seqs[a]),
               "n_cys": ncys(a), "n_single_cys_sites": by_prot_sites[a],
               "quantified_in_PXD072035": a in acc2ibaq, "ibaq_median": acc2ibaq.get(a)})
KV = pd.DataFrame(kv)
KV.to_csv(RES / "s05_kept_vs_deleted_proteins.csv", index=False)
kv_sum = KV.groupby("kept_in_stored_cohort").agg(
    n=("accession", "size"), median_length=("length", "median"), median_n_cys=("n_cys", "median"),
    mean_sites_per_protein=("n_single_cys_sites", "mean"), frac_quantified=("quantified_in_PXD072035", "mean"),
    median_ibaq_if_quantified=("ibaq_median", "median")).reset_index()
kv_sum["site_fraction_of_cys"] = [KV[KV.kept_in_stored_cohort == k].n_single_cys_sites.sum() /
                                   KV[KV.kept_in_stored_cohort == k].n_cys.sum() for k in kv_sum.kept_in_stored_cohort]
kv_sum.to_csv(RES / "s05_kept_vs_deleted_summary.csv", index=False)
print(kv_sum.to_string())

R = pd.DataFrame(rows)
R.to_csv(RES / "s05_cohort_rules.csv", index=False)
write_json(RES / "s05_reconcile.json", {
    "analysis_label": "POST HOC revision analysis 2026-09-30 (not registered); rule grid is exploratory",
    "target_publication_sites": TARGET_PUB,
    "grid_rules_within_1pct_of_1691": GRID[GRID.rel_diff_to_1691.abs() <= 0.01].to_dict("records"),
    "grid_rules_within_5pct_of_1691": int((GRID.rel_diff_to_1691.abs() <= 0.05).sum()),
    "grid_size": int(len(GRID)),
    "stored_cohort": {"n_sites": len(stored_sites), "n_proteins": len(stored_prot), "n_cysteines": int(len(co))},
    "stored_site_coverage": st_cov,
    "filter_chain": chain,
    "same_chain_without_single_Cys_restriction": alt,
})
print(pd.DataFrame(chain).to_string())
print(alt)
pd.set_option("display.width", 250)
print(GRID.to_string())
print(R.to_string())
print(json.dumps(st_cov, indent=1))
