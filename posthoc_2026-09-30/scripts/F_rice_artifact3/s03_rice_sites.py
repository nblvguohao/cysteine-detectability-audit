"""Step 3 (POST HOC): derive rice persulfidation cysteine sites from the deposit's peptide table
(PXD072089 SS-all-peptides.tsv) under explicit rules, and reconcile them with the numbers used in
the manuscript (817 / 737 sites, 5,598 cysteines in 640 proteins) and with the source publication
(1,691 persulfidated Cys sites, Lin et al. PNAS 2026, abstract as quoted in the review and in
Supplemental Data 7 row pubmed:42479832).

Site = (protein accession, 1-based residue position of a C in an identified peptide).
Position on the leading razor protein = 'Start position' + index of C in 'Sequence' (MaxQuant
reports Start/End for the leading razor protein only); every such position is checked against the
sequence (reference proteome 2026_03, else the sequence fetched in s01).
Position on any protein of the 'Proteins' column = every exact occurrence of the peptide in that
protein's sequence (all occurrences counted).

Outputs (results/F_rice_artifact3/):
  s03_site_rules.csv              one row per rule: n sites, n proteins, notes
  s03_sites_leading_razor.csv     every (leading razor protein, position) with per-run evidence
  s03_cohort_comparison.json      exact comparison with the stored self-audit cohort (640 proteins)
  s03_per_condition.csv           sites per condition / run and per identification type
"""
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

pep = pd.read_csv(PEP, sep="\t", low_memory=False, dtype=str)
seqs, src = load_all_sequences()
ref = read_fasta_any(FASTA_REF)

id_cols = [c for c in pep.columns if c.startswith("Identification type ")]
runs = [c.replace("Identification type ", "") for c in id_cols]          # e.g. '0-1'
cond_of = {r: r.split("-")[0] for r in runs}
conds = sorted(set(cond_of.values()), key=lambda x: int(x))

pep["C Count"] = pep["C Count"].astype(int)
pep["is_rev"] = pep["Reverse"].fillna("").eq("+") | pep["Leading razor protein"].fillna("").str.startswith("REV__")
pep["is_con"] = pep["Potential contaminant"].fillna("").eq("+") | pep["Leading razor protein"].fillna("").str.startswith("CON__")
pep["lrp"] = pep["Leading razor protein"].fillna("")
pep["start"] = pd.to_numeric(pep["Start position"], errors="coerce")
for r in runs:
    pep[f"id_{r}"] = pep[f"Identification type {r}"].fillna("")
pep["n_runs_any"] = sum((pep[f"id_{r}"] != "").astype(int) for r in runs)
pep["n_runs_msms"] = sum((pep[f"id_{r}"] == "By MS/MS").astype(int) for r in runs)
pep["n_conditions_any"] = sum(
    np.logical_or.reduce([(pep[f"id_{r}"] != "") for r in runs if cond_of[r] == c]).astype(int) for c in conds)

cys = pep[pep["C Count"] > 0].copy()

# ------------------------------------------------ leading-razor-protein sites -----
recs, pos_check = [], {"ok": 0, "mismatch": 0, "no_sequence": 0, "no_start": 0}
for _, r in cys.iterrows():
    s = r["Sequence"]
    idx = [i for i, ch in enumerate(s) if ch == "C"]
    if np.isnan(r["start"]):
        pos_check["no_start"] += 1
        continue
    st = int(r["start"])
    prot_seq = seqs.get(r["lrp"])
    if prot_seq is None:
        pos_check["no_sequence"] += 1
        status = "no_sequence"
    elif prot_seq[st - 1: st - 1 + len(s)] == s:
        pos_check["ok"] += 1
        status = "ok"
    else:
        pos_check["mismatch"] += 1
        status = "mismatch"
    for i in idx:
        recs.append({"protein": r["lrp"], "position": st + i, "peptide": s, "c_count": r["C Count"],
                     "is_rev": bool(r["is_rev"]), "is_con": bool(r["is_con"]),
                     "unique_groups": r["Unique (Groups)"], "pep_status": status,
                     "n_runs_any": int(r["n_runs_any"]), "n_runs_msms": int(r["n_runs_msms"]),
                     **{f"id_{x}": r[f"id_{x}"] for x in runs}})
S = pd.DataFrame(recs)
real = S[~S.is_rev & ~S.is_con]


def site_count(df):
    keys = set(zip(df.protein, df.position))
    return len(keys), len({k[0] for k in keys})


rules = []


def add(rule, df, note=""):
    n, p = site_count(df)
    rules.append({"rule": rule, "n_sites": n, "n_proteins": p, "note": note})


add("LRP_all_Cys_peptides_incl_REV_CON", S, "every C of every Cys-containing peptide, leading razor protein, decoys and contaminants kept")
add("LRP_all_Cys_peptides", real, "every C of every Cys-containing peptide, leading razor protein; REV/CON removed")
add("LRP_single_Cys_peptides", real[real.c_count == 1], "only peptides with exactly one C")
add("LRP_multi_Cys_peptides_only", real[real.c_count >= 2], "only peptides with >=2 C (all their C counted)")
add("LRP_all_Cys_peptides_ByMSMS_in_>=1_run", real[real.n_runs_msms >= 1], "peptide identified 'By MS/MS' in at least one run")
add("LRP_single_Cys_ByMSMS_in_>=1_run", real[(real.c_count == 1) & (real.n_runs_msms >= 1)], "")
add("LRP_all_Cys_peptides_>=2_runs", real[real.n_runs_any >= 2], "identification type populated in >=2 of 15 runs")
add("LRP_all_Cys_peptides_unique_groups", real[real.unique_groups == "yes"], "peptide unique to one protein group")
add("LRP_single_Cys_unique_groups", real[(real.c_count == 1) & (real.unique_groups == "yes")], "")

# per condition / per run
pc = []
for c in conds:
    rr = [x for x in runs if cond_of[x] == c]
    m_any = np.logical_or.reduce([real[f"id_{x}"] != "" for x in rr])
    m_ms = np.logical_or.reduce([real[f"id_{x}"] == "By MS/MS" for x in rr])
    for lab, m in (("any_identification", m_any), ("By_MSMS_only", m_ms)):
        n, p = site_count(real[m]); n1, p1 = site_count(real[m & (real.c_count == 1)])
        pc.append({"condition": c, "runs": ",".join(rr), "evidence": lab, "n_sites_all_cys_peptides": n,
                   "n_proteins": p, "n_sites_single_cys": n1, "n_proteins_single_cys": p1})
for x in runs:
    for lab, m in (("any_identification", real[f"id_{x}"] != ""), ("By_MSMS_only", real[f"id_{x}"] == "By MS/MS")):
        n, p = site_count(real[m])
        pc.append({"condition": cond_of[x], "runs": x, "evidence": lab, "n_sites_all_cys_peptides": n,
                   "n_proteins": p, "n_sites_single_cys": site_count(real[m & (real.c_count == 1)])[0],
                   "n_proteins_single_cys": site_count(real[m & (real.c_count == 1)])[1]})
PC = pd.DataFrame(pc)
PC.to_csv(RES / "s03_per_condition.csv", index=False)
cond_tot = PC[(PC.runs.str.contains(","))]
rules.append({"rule": "LRP_all_Cys_sum_over_5_conditions_any", "n_sites": int(cond_tot[cond_tot.evidence == "any_identification"].n_sites_all_cys_peptides.sum()),
              "n_proteins": None, "note": "sum of per-condition site counts (a site seen in k conditions counted k times)"})
rules.append({"rule": "LRP_all_Cys_sum_over_5_conditions_MSMS", "n_sites": int(cond_tot[cond_tot.evidence == "By_MSMS_only"].n_sites_all_cys_peptides.sum()),
              "n_proteins": None, "note": "as above, By MS/MS only"})

# peptide-level counts (a 'site' counted per peptide or per C in peptide)
rules.append({"rule": "Cys_peptides_count_real", "n_sites": int(((cys.is_rev == False) & (cys.is_con == False)).sum()),
              "n_proteins": None, "note": "number of Cys-containing peptide sequences, REV/CON removed"})
rules.append({"rule": "Cys_peptides_count_incl_REV_CON", "n_sites": int(len(cys)), "n_proteins": None, "note": ""})
rules.append({"rule": "Cys_residues_summed_over_peptides_real", "n_sites": int(cys[~cys.is_rev & ~cys.is_con]["C Count"].sum()),
              "n_proteins": None, "note": "sum of C Count over Cys-containing peptides (overlapping peptides double-count)"})
rules.append({"rule": "Cys_residues_summed_over_peptides_incl_REV_CON", "n_sites": int(cys["C Count"].sum()), "n_proteins": None, "note": ""})

# ------------------------------------------------ any-protein sites -----
anyrec, missing_prot, not_found = [], set(), 0
for _, r in cys[~cys.is_rev & ~cys.is_con].iterrows():
    s = r["Sequence"]
    idx = [i for i, ch in enumerate(s) if ch == "C"]
    for a in split_accs(r["Proteins"]):
        if is_decoy_or_contaminant(a):
            continue
        ps = seqs.get(a)
        if ps is None:
            missing_prot.add(a); continue
        k = ps.find(s)
        if k < 0:
            not_found += 1; continue
        while k >= 0:
            for i in idx:
                anyrec.append({"protein": a, "position": k + 1 + i, "c_count": r["C Count"]})
            k = ps.find(s, k + 1)
A = pd.DataFrame(anyrec)
add("ANY_protein_all_Cys_peptides", A, "every C of every Cys-containing peptide mapped to every protein of the 'Proteins' column (all exact occurrences)")
add("ANY_protein_single_Cys_peptides", A[A.c_count == 1], "")
# collapse proteins with identical sequences (redundant entries from overlapping rice annotation sets)
A["seqhash"] = A.protein.map(lambda a: hashlib.md5(seqs[a].encode()).hexdigest())
n_seqcollapsed = len(set(zip(A.seqhash, A.position)))
rules.append({"rule": "ANY_protein_all_Cys_collapsed_identical_sequences", "n_sites": n_seqcollapsed,
              "n_proteins": int(A.seqhash.nunique()), "note": "(sequence, position) pairs; identical protein sequences merged"})

R = pd.DataFrame(rules)
R.to_csv(RES / "s03_site_rules.csv", index=False)

# per-site table (leading razor, real)
agg = real.groupby(["protein", "position"]).agg(
    n_peptides=("peptide", "nunique"), min_c_count=("c_count", "min"), max_c_count=("c_count", "max"),
    any_single_cys=("c_count", lambda v: bool((v == 1).any())),
    n_runs_any=("n_runs_any", "max"), n_runs_msms=("n_runs_msms", "max"),
    pep_status=("pep_status", lambda v: ";".join(sorted(set(v))))).reset_index()
agg["residue_is_C"] = [bool(seqs.get(p) and 0 < q <= len(seqs[p]) and seqs[p][q - 1] == "C") for p, q in zip(agg.protein, agg.position)]
agg["in_reference_proteome_2026_03"] = agg.protein.map(lambda a: a in ref)
agg["seq_source"] = agg.protein.map(lambda a: src.get(a, "none"))
agg.to_csv(RES / "s03_sites_leading_razor.csv", index=False)

# ------------------------------------------------ comparison with the stored cohort -----
sc = pd.read_csv(REPO_SITE_SCORES, encoding="utf-8-sig", low_memory=False)
co = sc[sc.dataset == "rice_PXD072089"][["accession", "position", "observed_in_retrospective_dataset"]].copy()
co_obs = set(map(tuple, co[co.observed_in_retrospective_dataset == 1][["accession", "position"]].itertuples(index=False)))
co_prot = set(co.accession)
lrp_all = set(zip(real.protein, real.position))
lrp_single = set(zip(real[real.c_count == 1].protein, real[real.c_count == 1].position))

# proteins' cysteine totals from our sequences
cys_tot = {a: seqs[a].count("C") for a in co_prot if a in seqs}
per_prot_n = co.groupby("accession").size()
agree_n = sum(1 for a in co_prot if a in seqs and cys_tot[a] == per_prot_n[a])
# 1-Cys proteins among single-Cys-peptide site proteins (candidate exclusion rule)
single_prot = {p for p, _ in lrp_single}
single_sites_by_prot = {}
for p, q in lrp_single:
    single_sites_by_prot.setdefault(p, set()).add(q)
cand = {p for p in single_prot if p in seqs and seqs[p].count("C") >= 2}
cand_sites = {(p, q) for p, q in lrp_single if p in cand}
not_ref = {a for a in co_prot if a not in ref}
comp = {
    "analysis_label": "POST HOC revision analysis 2026-09-30 (not registered)",
    "stored_cohort": {"n_cysteines": int(len(co)), "n_proteins": len(co_prot), "n_observed": len(co_obs)},
    "stored_observed_vs_LRP_single_Cys_sites": {
        "n_LRP_single": len(lrp_single), "intersection": len(co_obs & lrp_single),
        "stored_only": len(co_obs - lrp_single), "LRP_single_only": len(lrp_single - co_obs)},
    "stored_observed_vs_LRP_all_Cys_sites": {
        "n_LRP_all": len(lrp_all), "intersection": len(co_obs & lrp_all),
        "stored_only": len(co_obs - lrp_all), "LRP_all_only": len(lrp_all - co_obs)},
    "rule_single_Cys_sites_in_proteins_with_>=2_Cys": {"n_sites": len(cand_sites), "n_proteins": len(cand),
        "equals_stored_observed_set": cand_sites == co_obs,
        "equals_stored_protein_set": cand == co_prot},
    "stored_cohort_cys_total_matches_sequence_count": {"n_proteins_with_seq": len(cys_tot), "n_agree": agree_n},
    "stored_cohort_proteins_not_in_reference_2026_03": {"n_proteins": len(not_ref),
        "n_cysteines": int(co[co.accession.isin(not_ref)].shape[0]),
        "n_observed": int(co[co.accession.isin(not_ref) & (co.observed_in_retrospective_dataset == 1)].shape[0])},
    "leading_razor_position_check_per_Cys_peptide": pos_check,
    "any_protein_mapping": {"n_missing_protein_sequences": len(missing_prot), "n_peptide_protein_pairs_not_found": not_found},
}
# inspect the discrepancy, if any
diff_stored = sorted(co_obs - cand_sites)[:30]
diff_ours = sorted(cand_sites - co_obs)[:30]
comp["examples_stored_not_in_rule"] = [list(x) for x in diff_stored]
comp["examples_rule_not_in_stored"] = [list(x) for x in diff_ours]
write_json(RES / "s03_cohort_comparison.json", comp)
print(R.to_string())
print(json.dumps(comp, indent=1)[:4000])
