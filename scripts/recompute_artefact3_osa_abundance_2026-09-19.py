"""Artefact 3 recomputed from public data: does protein abundance explain protein-level commonality in the rice persulfidome?

Replaces a set of numbers that entered the manuscript through a collaborator status table and that no
product in this tree could produce (CLAUDE.md 9.25.4, 9.37.1, 9.38). The underlying data is public.

INPUTS (both downloaded from PRIDE 2026-09-19, public depositions, no consent required)
  external/pride_osa_abundance_20260919/SS-all-peptides.tsv
      PXD072089, persulfidation peptides, rice leaf. PNAS 2026, doi 10.1073/pnas.2608150123.
  external/pride_osa_abundance_20260919/20240326_061934_xyj_proteome_Report.tsv
      PXD072035, the companion total-proteome quantification, 15 runs (5 timepoints x 3 replicates).
  external/proteomes/osa.fasta.gz   rice reference proteome, for length and cysteine count.

=========================== CRITERIA, WRITTEN BEFORE THE RUN ===========================

UNIT AND JOIN
  The unit is the protein group of the quantification table. Its identity is its FIRST accession in
  PG.ProteinAccessions. A group is labelled PERSULFIDATED when any of its accessions appears in the
  'Leading razor protein' column of the peptide table, and UNMODIFIED otherwise. Groups with no
  usable abundance are dropped from both arms and counted.

ABUNDANCE
  Primary: iBAQ. The PG.IBAQ cells are semicolon-separated per accession in the group, in the same
  order as PG.ProteinAccessions (verified on the file: 7,366 rows with a parseable value have equal
  segment and accession counts, 0 mismatched). The group's value in a run is the FIRST segment, i.e.
  the leading accession's. A group's abundance is the MEDIAN over the runs in which it has a value.
  Sensitivity: PG.Quantity, one value per group per run, same median rule.
  iBAQ is primary BECAUSE it is normalised by the number of theoretically observable peptides.
  Raw intensity rises with protein length, and protein length is one of the two outcomes here, so
  using it as the primary would build the answer into the covariate.

OUTCOME 1 - the abundance difference
  ratio = median(abundance | persulfidated) / median(abundance | unmodified)
  Interval: bootstrap over protein groups, 5,000 resamples, seed 20260919 (CLAUDE.md 2). The unit of
  resampling is the protein group because the unit of analysis is the protein group. Homology
  components are NOT available for rice in this tree, so the clustering CLAUDE.md 2 asks for is not
  applied here; this is declared, not silently skipped, and it makes the interval narrower than a
  homology-clustered one would be.

OUTCOME 2 - do length and cysteine count separate the two arms, before and after abundance matching
  Statistic: the Mann-Whitney U z (normal approximation, tie-corrected). "Before" is the whole
  comparison. "After" is the stratified (van Elteren) combination of the per-stratum z over abundance
  DECILES of the pooled distribution, equal weights, strata with fewer than 10 in either arm dropped
  and counted.

PRE-WRITTEN BRANCHES (the verdict is whichever fires; no other reading is permitted afterwards)
  B1 ratio interval excludes 1.0 AND both |z| fall below 2 after matching
     -> the original claim reproduces: the protein-level signal is abundance.
  B2 ratio interval excludes 1.0 AND either |z| stays at or above 2 after matching
     -> abundance does not account for it; the manuscript's claim would have to be weakened.
  B3 ratio interval includes 1.0
     -> there is no abundance difference to control for; the premise of Artefact 3 fails here.
  B4 fewer than 100 persulfidated groups carry an abundance
     -> underpowered; no ratio and no z are reported.

POSITIVE AND NEGATIVE CONTROLS
  PC1 both input sha256 recorded; row counts reported.
  PC2 the persulfidated accession set must be non-empty and its size reported.
  PC3 length and cysteine count must be recovered for at least 80% of analysed groups, else REFUSE.
  PC4 label shuffle: permute the persulfidated label over the analysed groups (same seed), recompute
      the ratio. Its interval MUST include 1.0, else the pipeline is broken and the run REFUSES.
  PC5 the two arms must partition the analysed set: n_pos + n_neg == n_analysed, exactly.

NOT ATTEMPTED THIS ROUND, AND THEREFORE NOT REPORTED
  The annotation-term and keyword parts of the original paragraph (twelve terms with |z|>3, 42
  keywords, four surviving). They need a rice keyword annotation source that this tree does not hold.
  No substitute is invented for them.

This script computes. It does not edit the manuscript.
"""
import csv, gzip, hashlib, json, pathlib, random, re, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
D = ROOT / "external/pride_osa_abundance_20260919"
PEP = D / "SS-all-peptides.tsv"
PRO = D / "20240326_061934_xyj_proteome_Report.tsv"
FASTA = ROOT / "external/proteomes/osa.fasta.gz"
OUT = ROOT / "results/artefact3_osa_abundance_recomputed_2026-09-19.csv"
AUD = ROOT / "results/artefact3_osa_abundance_recomputed_2026-09-19_audit.json"
SEED = 20260919
csv.field_size_limit(10**7)
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

# ---------- load ----------
with open(PEP) as f:
    peps = list(csv.DictReader(f, delimiter="\t"))
pos_accs = {(x.get("Leading razor protein") or "").strip() for x in peps}
pos_accs.discard("")
with open(PRO) as f:
    prot = list(csv.DictReader(f, delimiter="\t"))
IBAQ = [c for c in prot[0] if c.endswith("PG.IBAQ")]
QTY  = [c for c in prot[0] if c.endswith("PG.Quantity")]

def group_value(row, cols, first_segment):
    vals = []
    for c in cols:
        raw = (row.get(c) or "").strip()
        if not raw or raw == "NaN":
            continue
        tok = raw.split(";")[0] if first_segment else raw
        try:
            v = float(tok)
        except ValueError:
            continue
        if v > 0:
            vals.append(v)
    return float(np.median(vals)) if vals else None

# ---------- sequences ----------
seqs = {}
with gzip.open(FASTA, "rt") as fh:
    acc, buf = None, []
    for line in fh:
        if line.startswith(">"):
            if acc: seqs[acc] = "".join(buf)
            m = re.match(r">(?:\w+\|)?([A-Za-z0-9_.\-]+)", line)
            acc, buf = (m.group(1).split("|")[0] if m else None), []
        else:
            buf.append(line.strip())
    if acc: seqs[acc] = "".join(buf)

rows, n_no_abund = [], 0
for p in prot:
    accs = [a.strip() for a in (p.get("PG.ProteinAccessions") or "").split(";") if a.strip()]
    if not accs: continue
    ib = group_value(p, IBAQ, True)
    qt = group_value(p, QTY, False)
    if ib is None:
        n_no_abund += 1
        continue
    s = None
    for a in accs:
        if a in seqs: s = seqs[a]; break
    rows.append({"group_id": accs[0], "n_accessions": len(accs),
                 "persulfidated": int(any(a in pos_accs for a in accs)),
                 "ibaq_median": ib, "quantity_median": qt,
                 "length": len(s) if s else None,
                 "cys_count": s.count("C") if s else None})

n_analysed = len(rows)
pos = [r for r in rows if r["persulfidated"]]
neg = [r for r in rows if not r["persulfidated"]]
pc5 = (len(pos) + len(neg) == n_analysed)
with_seq = sum(1 for r in rows if r["length"] is not None)
pc3 = with_seq >= 0.80 * n_analysed
pc2 = len(pos_accs) > 0

def ratio_of_medians(a, b):
    return float(np.median(a) / np.median(b)) if len(a) and len(b) and np.median(b) > 0 else float("nan")

def boot_ratio(values, labels, n=5000, seed=SEED):
    rng = np.random.default_rng(seed)
    v, l = np.asarray(values, float), np.asarray(labels, int)
    idx = np.arange(len(v)); out = []
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        vs, ls = v[s], l[s]
        if ls.sum() == 0 or (1 - ls).sum() == 0: continue
        mb = np.median(vs[ls == 0])
        if mb > 0: out.append(np.median(vs[ls == 1]) / mb)
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))) if out else (float("nan"),) * 2

def mwu_z(x, y):
    """Mann-Whitney U z, normal approximation with tie correction."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2: return float("nan")
    allv = np.concatenate([x, y])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(len(allv), float)
    sorted_v = allv[order]
    i = 0
    while i < len(sorted_v):
        j = i
        while j + 1 < len(sorted_v) and sorted_v[j + 1] == sorted_v[i]: j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    R1 = ranks[:nx].sum()
    U = R1 - nx * (nx + 1) / 2.0
    mu = nx * ny / 2.0
    _, counts = np.unique(allv, return_counts=True)
    tie = (counts ** 3 - counts).sum()
    N = nx + ny
    sd = np.sqrt(nx * ny / 12.0 * ((N + 1) - tie / (N * (N - 1))))
    return float((U - mu) / sd) if sd > 0 else float("nan")

def stratified_z(rows_, key, nstrata=10):
    vals = np.array([r["ibaq_median"] for r in rows_], float)
    edges = np.percentile(vals, np.linspace(0, 100, nstrata + 1))
    edges[0] -= 1e-9; edges[-1] += 1e-9
    zs, used, dropped = [], 0, 0
    for k in range(nstrata):
        sel = [r for r in rows_ if edges[k] < r["ibaq_median"] <= edges[k + 1] and r[key] is not None]
        a = [r[key] for r in sel if r["persulfidated"]]
        b = [r[key] for r in sel if not r["persulfidated"]]
        if len(a) < 10 or len(b) < 10: dropped += 1; continue
        z = mwu_z(a, b)
        if not np.isnan(z): zs.append(z); used += 1
    if not zs: return float("nan"), used, dropped
    return float(np.sum(zs) / np.sqrt(len(zs))), used, dropped

verdict_notes = []
if len(pos) < 100:
    branch = "B4_underpowered"
    ratio = ci = None
else:
    vals = [r["ibaq_median"] for r in rows]
    labs = [r["persulfidated"] for r in rows]
    ratio = ratio_of_medians([r["ibaq_median"] for r in pos], [r["ibaq_median"] for r in neg])
    ci = boot_ratio(vals, labs)
    have = [r for r in rows if r["length"] is not None]
    z_len_before = mwu_z([r["length"] for r in have if r["persulfidated"]],
                         [r["length"] for r in have if not r["persulfidated"]])
    z_cys_before = mwu_z([r["cys_count"] for r in have if r["persulfidated"]],
                         [r["cys_count"] for r in have if not r["persulfidated"]])
    z_len_after, used_l, drop_l = stratified_z(have, "length")
    z_cys_after, used_c, drop_c = stratified_z(have, "cys_count")
    excludes_one = not (ci[0] <= 1.0 <= ci[1])
    if not excludes_one:
        branch = "B3_no_abundance_difference_to_control"
    elif abs(z_len_after) < 2 and abs(z_cys_after) < 2:
        branch = "B1_original_claim_reproduces"
    else:
        branch = "B2_abundance_does_not_account_for_it"

# PC4 label shuffle
rng = random.Random(SEED)
sh = [r["persulfidated"] for r in rows]; rng.shuffle(sh)
sh_ratio = ratio_of_medians([r["ibaq_median"] for r, s in zip(rows, sh) if s],
                            [r["ibaq_median"] for r, s in zip(rows, sh) if not s])
sh_ci = boot_ratio([r["ibaq_median"] for r in rows], sh)
pc4 = sh_ci[0] <= 1.0 <= sh_ci[1]

# sensitivity on Quantity
q = [r for r in rows if r["quantity_median"]]
q_ratio = ratio_of_medians([r["quantity_median"] for r in q if r["persulfidated"]],
                           [r["quantity_median"] for r in q if not r["persulfidated"]])

ok = pc2 and pc3 and pc4 and pc5
res = [
 ("n_protein_groups_in_quantification", len(prot), "", "rows in the PXD072035 report"),
 ("n_groups_dropped_no_usable_iBAQ", n_no_abund, "", ""),
 ("n_groups_analysed", n_analysed, "", ""),
 ("n_persulfidated", len(pos), "", "group has >=1 accession among the peptide table's leading razor proteins"),
 ("n_unmodified", len(neg), "", ""),
 ("n_unique_leading_razor_proteins_in_peptide_table", len(pos_accs), "", ""),
 ("n_groups_with_sequence", with_seq, "", "length and cysteine count recovered"),
]
if branch != "B4_underpowered":
    res += [
     ("median_iBAQ_ratio_persulfidated_over_unmodified", round(ratio, 4), f"[{ci[0]:.4f}, {ci[1]:.4f}]", "5,000 bootstrap over groups, seed 20260919"),
     ("z_protein_length_before_matching", round(z_len_before, 4), "", "Mann-Whitney U z, tie-corrected"),
     ("z_protein_length_after_abundance_matching", round(z_len_after, 4), "", f"stratified over iBAQ deciles; strata used {used_l}, dropped {drop_l}"),
     ("z_cysteine_count_before_matching", round(z_cys_before, 4), "", "Mann-Whitney U z, tie-corrected"),
     ("z_cysteine_count_after_abundance_matching", round(z_cys_after, 4), "", f"stratified over iBAQ deciles; strata used {used_c}, dropped {drop_c}"),
     ("sensitivity_median_Quantity_ratio", round(q_ratio, 4), "", "not length-normalised; declared sensitivity only"),
     ("negative_control_shuffled_label_ratio", round(sh_ratio, 4), f"[{sh_ci[0]:.4f}, {sh_ci[1]:.4f}]", "PC4; interval must include 1.0"),
    ]
with open(OUT, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["quantity", "value", "interval_95", "note"]); w.writerows(res)

AUD.write_text(json.dumps({
 "script": "scripts/recompute_artefact3_osa_abundance_2026-09-19.py", "script_sha256": sha(__file__),
 "interpreter": sys.version.split()[0], "numpy": np.__version__, "seed": SEED,
 "inputs": {str(PEP.relative_to(ROOT)): sha(PEP), str(PRO.relative_to(ROOT)): sha(PRO),
            str(FASTA.relative_to(ROOT)): sha(FASTA)},
 "input_rows": {"peptides": len(peps), "protein_groups": len(prot)},
 "branch_taken": branch,
 "controls": {"PC2_positive_set_nonempty": pc2, "PC3_sequences_ge_80pct": pc3,
              "PC4_shuffled_label_interval_includes_1": pc4, "PC5_arms_partition": pc5,
              "all_pass": ok},
 "not_attempted": ["annotation terms with |z|>3", "the 42-keyword screen"],
 "declared_deviation_from_house_rule_2": "the bootstrap resamples protein groups, not homology "
   "components; no rice homology clustering exists in this tree. The interval is therefore narrower "
   "than a component-clustered interval would be.",
 "what_this_does_not_establish": [
   "It does not reproduce the collaborator's pipeline; group definition, abundance summarisation and "
   "the z statistic are defined here and may differ from theirs. Agreement in magnitude would be "
   "evidence; disagreement does not by itself show either side is wrong.",
   "It says nothing about tomato, where the same comparison is reported to run the other way.",
   "iBAQ is itself a model of abundance, not a measurement of copy number."],
}, ensure_ascii=False, indent=1, sort_keys=True))

for r in res: print(f"  {r[0]:<52} {r[1]:>12} {r[2]}")
print("\nbranch:", branch)
print("controls:", {"PC2": pc2, "PC3": pc3, "PC4": pc4, "PC5": pc5})
sys.exit(0 if ok else 1)
