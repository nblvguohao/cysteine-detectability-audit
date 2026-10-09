"""Step 11 (POST HOC; revision round 2, 2026-09-30): the rice site-list filter chain with positions taken
from MaxQuant's reported start position, for the proposed Experimental Procedures text.

Why this step exists. s05's chain re-located every single-cysteine peptide on the canonical sequence by
string search and counted every exact occurrence. One peptide (IVCTLPR) occurs twice in B9F8C3, so s05
counts 1,378 sites at step 3, 1,341 at step 4 and 524 lost sites, whereas positions from the deposit's
'Start position' (the site definition used everywhere else in this item and in the verifier's rebuild)
give 1,377, 1,340 and 523. The final 817-site list is identical under both conventions (B9F8C3 is one of
the 400 deleted entries). This step recomputes the chain under the start-position convention, states
the one-site difference explicitly, and collects the cohort facts the proposed text quotes:
  - the deposit (PXD072089) contains a MaxQuant peptide table and no site list;
  - 817 sites / 640 proteins / 5,598 cysteines result from the chain below and equal the stored
    self-audit list site by site;
  - 400 proteins (523 sites) were lost at the last step, all UniProtKB entries deleted in release 2026_02;
  - the self-audit then excluded 594 cysteines (80 sites) in 59 proteins, leaving 5,004 cysteines,
    737 sites and 581 proteins (repository audit: rows_kept, positives, proteins_matched).
Isoform leading razor proteins ('-n' suffix) are mapped to the canonical entry by locating the peptide
in the canonical sequence; every such peptide occurs exactly once (asserted).
Round 3 (2026-09-30, after verifier round 3): the 22 decoy rows among the cysteine-containing peptides have no
'Start position'. Round 2 skipped them silently, so its step 0 ("decoys and contaminants kept", 1,769 sites /
1,170 proteins) in fact kept contaminants only. They are now keyed by (protein, peptide, offset of the C), as
in the verifier's rebuild: step 0 is 1,796 sites / 1,192 proteins, and a new step 0b (decoys removed,
contaminants kept) gives the former 1,769 / 1,170. Steps 1-5 are unchanged (asserted).
Outputs: s11_rice_chain.csv, s11_rice_chain.json
"""
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa

seqs, src = load_all_sequences()
ref = read_fasta_any(FASTA_REF)
pep = pd.read_csv(PEP, sep="\t", low_memory=False, dtype=str)
pep["ccount"] = pep["C Count"].astype(int)
pep["lrp"] = pep["Leading razor protein"].fillna("")
pep["start"] = pd.to_numeric(pep["Start position"], errors="coerce")
pep["is_rev"] = (pep["Reverse"].fillna("") == "+") | pep.lrp.str.startswith("REV__")
pep["is_con"] = (pep["Potential contaminant"].fillna("") == "+") | pep.lrp.str.startswith("CON__")
cys = pep[pep.ccount > 0].copy()
real = cys[~cys.is_rev & ~cys.is_con]
single = real[real.ccount == 1]
fs_rows = list(csv.DictReader(open(FETCHED / "fetched_sequences.tsv", encoding="utf-8"), delimiter="\t"))
active = set(ref) | {r["requested_accession"] for r in fs_rows if r["source"] == "uniprot_rest_active"}
inactive_reasons = json.loads((FETCHED / "inactive_reasons_400.json").read_text(encoding="utf-8"))


def ncys(a):
    return seqs[a].count("C") if a in seqs else 0


def start_sites(df):
    """(leading razor protein, Start position + index of each C); positions verified against the sequence.
    Rows without a start position (only decoys; asserted) are keyed by (protein, 'peptide:offset')."""
    out, bad = set(), 0
    for s, a, st, rev in zip(df.Sequence, df.lrp, df.start, df.is_rev):
        if np.isnan(st):
            assert rev, (a, s)
            for i, ch in enumerate(s):
                if ch == "C":
                    out.add((a, f"{s}:{i}"))
            continue
        st = int(st)
        if a in seqs and seqs[a][st - 1:st - 1 + len(s)] != s:
            bad += 1
        for i, ch in enumerate(s):
            if ch == "C":
                out.add((a, st + i))
    assert bad == 0
    return out


def canonical_start_sites(df):
    """As start_sites, but isoform accessions are replaced by the canonical accession; the peptide is
    located in the canonical sequence and must occur exactly once."""
    out, n_iso = set(), 0
    for s, a, st in zip(df.Sequence, df.lrp, df.start):
        if "-" in a:
            c = a.split("-")[0]
            k = seqs[c].find(s)
            assert k >= 0 and seqs[c].find(s, k + 1) < 0, (a, s)
            n_iso += 1
            for i, ch in enumerate(s):
                if ch == "C":
                    out.add((c, k + 1 + i))
        else:
            st = int(st)
            for i, ch in enumerate(s):
                if ch == "C":
                    out.add((a, st + i))
    return out, n_iso


sc = pd.read_csv(REPO_SITE_SCORES, encoding="utf-8-sig", low_memory=False)
co = sc[sc.dataset == "rice_PXD072089"][["accession", "position", "observed_in_retrospective_dataset"]]
stored = set(map(tuple, co[co.observed_in_retrospective_dataset == 1][["accession", "position"]].itertuples(index=False)))
stored_prot = set(co.accession)

chain = []


def step(name, S, note):
    P = {a for a, _ in S}
    chain.append({"step": name, "n_sites": len(S), "n_proteins": len(P), "n_cys_in_proteins": int(sum(ncys(a) for a in P)),
                  "equals_stored_sites": S == stored, "equals_stored_proteins": P == stored_prot,
                  "intersection_with_stored": len(S & stored), "note": note})
    return S


s0 = step("0_every_C_of_every_Cys_peptide_incl_REV_CON", start_sites(cys),
          "leading razor protein, Start position + index of C; the 22 decoy rows without a Start position are "
          "keyed by (protein, peptide, offset)")
s0b = step("0b_REV_removed_CON_kept", start_sites(cys[~cys.is_rev]), "decoys removed, contaminants kept")
s1 = step("1_REV_CON_removed", start_sites(real), "decoys and contaminants removed")
# round 3: step 0b must equal round 2's step 0, and step 1 must be unchanged
assert (len(s0b), len({a for a, _ in s0b})) == (1769, 1170)
assert (len(s1), len({a for a, _ in s1})) == (1753, 1164)
nan_rows = cys[cys.start.isna()]
decoy_info = {"cys_peptide_rows_without_start_position": int(len(nan_rows)),
              "all_of_them_decoys": bool(nan_rows.is_rev.all()),
              "decoy_sites": len(s0 - s0b), "decoy_proteins": len({a for a, _ in s0} - {a for a, _ in s0b}),
              "contaminant_sites": len(s0b - s1), "contaminant_proteins": len({a for a, _ in s0b} - {a for a, _ in s1}),
              "round2_step0_value": [1769, 1170],
              "round2_step0_note": "round 2 skipped the decoy rows without a start position, so its step 0 kept "
                                   "contaminants only; it equals step 0b here"}
s2 = step("2_single_Cys_peptides", start_sites(single), "peptides with exactly one C")
s3set, n_iso = canonical_start_sites(single)
s3 = step("3_isoforms_mapped_to_canonical", s3set, f"{n_iso} isoform peptides located (uniquely) in the canonical sequence")
s4 = step("4_proteins_with_ge2_Cys", {(a, q) for a, q in s3 if ncys(a) >= 2}, "protein has >= 2 cysteines")
s5 = step("5_UniProtKB_entry_active_2026_03", {(a, q) for a, q in s4 if a in active},
          "the 400 proteins lost here are all entries deleted in UniProtKB 2026_02")
lost_prot = {a for a, _ in s4} - {a for a, _ in s5}
lost_sites = s4 - s5
assert s5 == stored and {a for a, _ in s5} == stored_prot
CH = pd.DataFrame(chain)
CH.to_csv(RES / "s11_rice_chain.csv", index=False)

old = pd.read_csv(RES / "s05_filter_chain.csv")
audit = json.loads((REPO / "results" / "self_audit_public_cohorts_2026-09-20_audit.json").read_text(encoding="utf-8"))
cov = audit["gates"]["coverage_rice_PXD072089"]
out = {
    "analysis_label": "POST HOC revision analysis 2026-09-30, round 2 (not registered)",
    "deposit_has_site_list": False,
    "deposit_peptide_table": "SS-all-peptides.tsv (MaxQuant peptides); only site-level modification column: 'Oxidation (M) site IDs'",
    "chain_start_position_convention": chain,
    "chain_all_occurrence_convention_s05": old[["step", "n_sites", "n_proteins"]].to_dict("records"),
    "difference_between_conventions": "one site: peptide IVCTLPR occurs twice in B9F8C3 (a deleted entry); "
                                      "s05 counted both occurrences (1,378 / 1,341 / 524 lost), the start position gives one "
                                      "(1,377 / 1,340 / 523 lost); the final 817-site list is identical",
    "lost_at_step5": {"proteins": len(lost_prot), "sites": len(lost_sites),
                      "cysteines_in_lost_proteins": int(sum(ncys(a) for a in lost_prot)),
                      "all_lost_inactive_in_2026_03": all(a not in active for a in lost_prot),
                      "entry_records": inactive_reasons},
    "final_list": {"sites": len(s5), "proteins": len({a for a, _ in s5}), "cysteines": int(sum(ncys(a) for a in {a for a, _ in s5})),
                   "equals_stored_self_audit_list_site_by_site": s5 == stored},
    "self_audit_exclusion_repository_audit": {
        "rows_in_table": cov["rows_in_table"], "rows_kept": cov["rows_kept"], "proteins_in_table": cov["proteins_in_table"],
        "proteins_matched": cov["proteins_matched"], "positives_kept": cov["positives"],
        "excluded_cysteines": cov["rows_in_table"] - cov["rows_kept"], "excluded_sites": len(stored) - cov["positives"],
        "excluded_proteins": cov["proteins_in_table"] - cov["proteins_matched"]},
    "complete_list_without_step5": {"sites": len(s4), "proteins": len({a for a, _ in s4})},
    "round3_decoys_and_contaminants": decoy_info,
    # round 3: the saved UniProtKB entry records of the lost accessions carry no sequence
    "round3_lost_entry_records": {
        "n_lost": len(lost_prot),
        "n_records_saved": sum((FETCHED / "inactive_reasons" / f"{a}.json").exists() for a in lost_prot),
        "n_records_with_sequence_field": sum(
            "sequence" in json.loads((FETCHED / "inactive_reasons" / f"{a}.json").read_text(encoding="utf-8"))
            for a in lost_prot if (FETCHED / "inactive_reasons" / f"{a}.json").exists())},
}
write_json(RES / "s11_rice_chain.json", out)
print(CH.to_string())
print(json.dumps({k: v for k, v in out.items() if k not in ("chain_start_position_convention",)}, indent=1))
