"""VERIFIER round 4 (adversarial; POST HOC), item F_rice_artifact3: the rice site-list facts that the proposed
Experimental Procedures / Discussion / Supplemental Note texts state, rebuilt with the verifier's own code:
  - chain step 0 (decoys and contaminants kept, including the decoy rows without a start position), step 0b,
    step 1 and the later steps down to the stored 817-site list (repo public_cohort_site_scores.csv);
  - the 400 lost proteins, their sites and cysteines; kept vs lost descriptives;
  - the evidence cited for "UniProtKB returns no sequence for a deleted entry": the saved entry records
    (fetched by s01b) and the s01 batch query that requested sequences;
  - the self-audit exclusion (repo audit JSON) and the 1,692 'unique to one group' count.
No network. Output: verify_r4/rice_r4.json (+ stdout).
"""
import gzip
import json
import pathlib
import re

import numpy as np
import pandas as pd

W = pathlib.Path("/path/to/local/_cys_repo_work/public/revision_2026-09-30")
REPO = pathlib.Path("/path/to/local/_cys_repo_work/repo")
RES = W / "results" / "F_rice_artifact3"
OUT = RES / "verify_r4"
FET = RES / "fetched"


def read_fasta(p):
    seqs, acc, buf = {}, None, []
    with gzip.open(p, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(buf)
                acc = line[1:].split("|")[1]; buf = []
            else:
                buf.append(line.strip())
    if acc is not None:
        seqs[acc] = "".join(buf)
    return seqs


ref = read_fasta(W / "external" / "UP000059680_39947.fasta.gz")
fs = pd.read_csv(FET / "fetched_sequences.tsv", sep="\t", dtype=str, keep_default_na=False)
seqs = dict(ref)
for a, s in zip(fs.requested_accession, fs.sequence):
    if s and a not in seqs:
        seqs[a] = s
active = set(ref) | set(fs.loc[fs.source == "uniprot_rest_active", "requested_accession"])

pep = pd.read_csv(W / "external" / "SS-all-peptides.tsv", sep="\t", dtype=str, keep_default_na=False)
pep["cc"] = pep["C Count"].astype(int)
lrp = pep["Leading razor protein"]
pep["dec"] = (pep["Reverse"] == "+") | lrp.str.startswith("REV__")
pep["con"] = (pep["Potential contaminant"] == "+") | lrp.str.startswith("CON__")
cys = pep[pep.cc > 0]


def sites(df, key_missing_start=True):
    out = set()
    for s, a, st in zip(df.Sequence, df["Leading razor protein"], df["Start position"]):
        offs = [m.start() for m in re.finditer("C", s)]
        if st == "":
            for o in offs:
                out.add((a, f"{s}:{o}"))
        else:
            for o in offs:
                out.add((a, int(st) + o))
    return out


chain = {}
s0 = sites(cys); chain["0_all_incl_decoys_contaminants"] = s0
s0b = sites(cys[~cys.dec]); chain["0b_decoys_removed"] = s0b
real = cys[~cys.dec & ~cys.con]
s1 = sites(real); chain["1_decoys_contaminants_removed"] = s1
single = real[real.cc == 1]
s2 = sites(single); chain["2_single_cys"] = s2
# isoforms -> canonical by unique location in the canonical sequence
s3 = set(); iso = 0; iso_nonunique = 0
for s, a, st in zip(single.Sequence, single["Leading razor protein"], single["Start position"]):
    o = s.index("C")
    if "-" in a:
        c = a.split("-")[0]; q = seqs[c]
        k = q.find(s)
        iso += 1
        if k < 0 or q.find(s, k + 1) >= 0:
            iso_nonunique += 1
        s3.add((c, k + 1 + o))
    else:
        s3.add((a, int(st) + o))
chain["3_isoforms_to_canonical"] = s3
s4 = {(a, p) for a, p in s3 if a in seqs and seqs[a].count("C") >= 2}; chain["4_ge2_cys"] = s4
s5 = {(a, p) for a, p in s4 if a in active}; chain["5_active_2026_03"] = s5

sc = pd.read_csv(REPO / "external" / "public_cohort_site_scores.csv", encoding="utf-8-sig", low_memory=False)
rice = sc[sc.dataset == "rice_PXD072089"]
stored = set(zip(rice.loc[rice.observed_in_retrospective_dataset == 1, "accession"],
                 rice.loc[rice.observed_in_retrospective_dataset == 1, "position"].astype(int)))
res = {"chain": {k: {"sites": len(v), "proteins": len({a for a, _ in v}),
                     "cys_in_proteins": int(sum(seqs[a].count("C") for a in {a for a, _ in v} if a in seqs))}
                 for k, v in chain.items()},
       "rows_without_start_position": int((cys["Start position"] == "").sum()),
       "rows_without_start_all_decoys": bool(cys.loc[cys["Start position"] == "", "dec"].all()),
       "isoform_peptides": iso, "isoform_nonunique": iso_nonunique,
       "final_equals_stored": s5 == stored, "stored_sites": len(stored), "stored_rows": int(len(rice)),
       "stored_proteins": int(rice.accession.nunique())}
lost_p = {a for a, _ in s4} - {a for a, _ in s5}
lost_s = s4 - s5
res["lost"] = {"proteins": len(lost_p), "sites": len(lost_s), "cys": int(sum(seqs[a].count("C") for a in lost_p))}

# ---- evidence for "UniProtKB returns no sequence for a deleted entry"
recs = {a: json.loads((FET / "inactive_reasons" / f"{a}.json").read_text(encoding="utf-8")) for a in lost_p
        if (FET / "inactive_reasons" / f"{a}.json").exists()}
keys = sorted({k for r in recs.values() for k in r})
s01b_src = (W / "scripts" / "F_rice_artifact3" / "s01b_inactive_reasons.py").read_text(encoding="utf-8")
m = re.search(r'rest\.uniprot\.org/uniprotkb/\{a\}\?([^"]+)"', s01b_src)
requested = {a.strip() for a in (RES / "s00_accessions_to_fetch.txt").read_text(encoding="utf-8").split()}
active_batches = set()
for f in sorted(FET.glob("active_batch_*.tsv")):
    for ln in f.read_text(encoding="utf-8").splitlines()[1:]:
        if ln.strip():
            active_batches.add(ln.split("\t")[0])
s01_src = (W / "scripts" / "F_rice_artifact3" / "s01_fetch_sequences.py").read_text(encoding="utf-8")
res["deleted_entry_sequence_evidence"] = {
    "entry_records_saved": len(recs), "json_keys_in_saved_records": keys,
    "s01b_query_string": m.group(1) if m else None,
    "note": "the s01b request restricted the returned fields to 'accession', so no saved record could carry a "
            "sequence whatever UniProtKB holds; the 0/400 check is uninformative by construction",
    "lost_requested_in_s01_batch_query": len(lost_p & requested),
    "lost_returned_by_s01_batch_query_with_sequence_field": len(lost_p & active_batches),
    "s01_batch_query_fields": re.search(r'uniprotkb/accessions\?accessions=\{[^}]+\}&format=tsv"\s*f"&fields=([a-z_,]+)', s01_src).group(1),
    "lost_source_in_fetched_sequences": fs[fs.requested_accession.isin(lost_p)].source.value_counts().to_dict()}

# ---- kept vs lost descriptives
kept_p = {a for a, _ in s5}
def med(ps, f):
    return float(np.median([f(a) for a in ps]))
G = pd.read_csv(OUT / "groups_r4.csv")
quant = set()
for accs in G.accessions:
    quant.update(accs.split(";"))
rep = pd.read_csv(W / "external" / "20240326_061934_xyj_proteome_Report.tsv", sep="\t", dtype=str, keep_default_na=False)
allq = set()
for accs in rep["PG.ProteinAccessions"]:
    allq.update(a.strip() for a in accs.split(";") if a.strip())
spp = lambda ps, S: len(S) / len(ps)
res["kept_vs_lost"] = {
    "median_length": [med(kept_p, lambda a: len(seqs[a])), med(lost_p, lambda a: len(seqs[a]))],
    "median_cys": [med(kept_p, lambda a: seqs[a].count("C")), med(lost_p, lambda a: seqs[a].count("C"))],
    "sites_per_protein": [spp(kept_p, s5), spp(lost_p, lost_s)],
    "frac_in_quantified_groups_7692": [len(kept_p & quant) / len(kept_p), len(lost_p & quant) / len(lost_p)],
    "frac_in_any_report_row_7825": [len(kept_p & allq) / len(kept_p), len(lost_p & allq) / len(lost_p)]}

# ---- 1,692: peptides unique to one protein group
uq = real[real["Unique (Groups)"] == "yes"]
res["unique_groups_every_C"] = len(sites(uq))
# ---- self-audit exclusion
aud = json.loads((REPO / "results" / "self_audit_public_cohorts_2026-09-20_audit.json").read_text(encoding="utf-8"))
res["self_audit_gate"] = aud["gates"]["coverage_rice_PXD072089"]
# ---- the 1,691 in the census
census = pd.read_csv(pathlib.Path("/path/to/local/巯基化/MCP/supplemental/Supplemental_Data_7_dataset_census.csv"),
                     dtype=str, keep_default_na=False)
hit = census[census.apply(lambda r: "42479832" in " ".join(r.values), axis=1)]
res["census_rows_for_pubmed_42479832"] = hit.to_dict("records")
(OUT / "rice_r4.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
print(json.dumps({k: v for k, v in res.items() if k != "census_rows_for_pubmed_42479832"}, indent=1, default=str))
print("census:", json.dumps(res["census_rows_for_pubmed_42479832"], ensure_ascii=False)[:1500])
