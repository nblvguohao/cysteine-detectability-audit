"""VERIFIER round 4 (adversarial; POST HOC), item F_rice_artifact3.

Independent rebuild of the Artifact 3 protein-group table from the raw public inputs, with the verifier's
own parsing, own tryptic digest and own label rules (written from the report's definitions, not imported
from the item's code). Compared, group by group, with the item's s08_group_features_windows.csv.

Inputs (read-only): W/external/20240326_061934_xyj_proteome_Report.tsv (PXD072035),
W/external/SS-all-peptides.tsv (PXD072089), W/external/UP000059680_39947.fasta.gz and the item's
fetched_sequences.tsv (the only record of the UniProt/UniParc sequences of accessions outside the
reference proteome; the verifier cannot re-fetch 13k accessions).
Outputs: W/results/F_rice_artifact3/verify_r4/groups_r4.csv, build_r4.json
Deterministic (no randomness).
"""
import csv
import gzip
import hashlib
import json
import pathlib

import numpy as np
import pandas as pd

W = pathlib.Path("C:/Users/admin/Desktop/小论文/_cys_repo_work/public/revision_2026-09-30")
RES = W / "results" / "F_rice_artifact3"
OUT = RES / "verify_r4"
OUT.mkdir(parents=True, exist_ok=True)
PRO = W / "external" / "20240326_061934_xyj_proteome_Report.tsv"
PEP = W / "external" / "SS-all-peptides.tsv"
FASTA = W / "external" / "UP000059680_39947.fasta.gz"
FETCHED = RES / "fetched" / "fetched_sequences.tsv"
csv.field_size_limit(10 ** 9)


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def read_fasta(p):
    seqs, acc, buf = {}, None, []
    with gzip.open(p, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(buf)
                acc = line[1:].split("|")[1]
                buf = []
            else:
                buf.append(line.strip())
    if acc is not None:
        seqs[acc] = "".join(buf)
    return seqs


def tryptic_detectable_mask(s, lo=7, hi=30, mc=2, block_p=True):
    """Boolean array: residue covered by >=1 fully tryptic peptide (cut after K/R, not before P if block_p)
    with <= mc missed cleavages and length in [lo, hi]."""
    n = len(s)
    cuts = [0]
    for i, ch in enumerate(s):
        if ch in "KR" and i + 1 < n and not (block_p and s[i + 1] == "P"):
            cuts.append(i + 1)
    cuts.append(n)
    cuts = sorted(set(cuts))
    cov = np.zeros(n + 1, int)
    for a in range(len(cuts) - 1):
        for k in range(mc + 1):
            b = a + 1 + k
            if b >= len(cuts):
                break
            L = cuts[b] - cuts[a]
            if lo <= L <= hi:
                cov[cuts[a]] += 1
                cov[cuts[b]] -= 1
    return np.cumsum(cov)[:n] > 0


def median_pos(vals):
    v = []
    for x in vals:
        if x is None:
            continue
        x = str(x).strip()
        if not x or x == "NaN" or x == "nan":
            continue
        try:
            f = float(x)
        except ValueError:
            continue
        if f > 0:
            v.append(f)
    return float(np.median(v)) if v else None


def main():
    # ---------------- sequences
    seqs = read_fasta(FASTA)
    n_ref = len(seqs)
    fetched = pd.read_csv(FETCHED, sep="\t", dtype=str, keep_default_na=False)
    for a, s in zip(fetched.requested_accession, fetched.sequence):
        if s and a not in seqs:
            seqs[a] = s
    # ---------------- peptides
    pep = pd.read_csv(PEP, sep="\t", dtype=str, keep_default_na=False)
    pep["cc"] = pep["C Count"].astype(int)
    lrp = pep["Leading razor protein"]
    decoy = (pep["Reverse"] == "+") | lrp.str.startswith("REV__")
    cont = (pep["Potential contaminant"] == "+") | lrp.str.startswith("CON__")
    lrp_all = set(lrp) - {""}
    real_cys = pep[(pep.cc > 0) & ~decoy & ~cont]
    lrp_cys = set(real_cys["Leading razor protein"])
    lrp_single = set(real_cys.loc[real_cys.cc == 1, "Leading razor protein"])
    peps_by_lrp = {}
    for s, a in zip(real_cys.Sequence, real_cys["Leading razor protein"]):
        peps_by_lrp.setdefault(a, []).append(s)
    # any accession appearing anywhere in the 'Proteins' column (razor or not), decoys/contaminants removed
    prot_any = set()
    for s in pep.loc[~decoy & ~cont, "Proteins"]:
        for a in s.split(";"):
            if a.strip():
                prot_any.add(a.strip())
    # ---------------- protein groups
    rep = pd.read_csv(PRO, sep="\t", dtype=str, keep_default_na=False)
    ib_cols = [c for c in rep.columns if c.endswith(".PG.IBAQ")]
    q_cols = [c for c in rep.columns if c.endswith(".PG.Quantity")]
    assert len(ib_cols) == 15 and len(q_cols) == 15
    rows, n_drop = [], 0
    alpha_sorted = 0
    for _, r in rep.iterrows():
        accs = [a.strip() for a in r["PG.ProteinAccessions"].split(";") if a.strip()]
        if not accs:
            continue
        alpha_sorted += int(accs == sorted(accs))
        ib = median_pos([str(r[c]).split(";")[0] for c in ib_cols])
        q = median_pos([r[c] for c in q_cols])
        if ib is None:
            n_drop += 1
            continue
        lead = accs[0]
        s = seqs.get(lead)
        rec = {"group_id": lead, "accessions": ";".join(accs), "n_acc": len(accs), "ibaq": ib, "quantity": q,
               "log10_ibaq": float(np.log10(ib)),
               "A": int(any(a in lrp_all for a in accs)),
               "S": int(any(a in lrp_cys for a in accs)),
               "S1": int(any(a in lrp_single for a in accs)),
               "in_proteins_col_any": int(any(a in prot_any for a in accs)),
               "has_seq": int(s is not None)}
        if s is not None:
            det = tryptic_detectable_mask(s)
            cpos = [i for i, ch in enumerate(s) if ch == "C"]
            sites = set()
            for a in accs:
                for ps in peps_by_lrp.get(a, []):
                    k = s.find(ps)
                    while k >= 0:
                        for i, ch in enumerate(ps):
                            if ch == "C":
                                sites.add(k + i)
                        k = s.find(ps, k + 1)
            rec.update({"length": len(s), "n_cys": len(cpos), "n_det": int(det[cpos].sum()) if cpos else 0,
                        "S_map": int(any(det[p] for p in sites))})
            # windows used in the report's robustness section
            for (lo, hi, mc, bp, tag) in [(6, 35, 2, True, "w6_35"), (7, 40, 3, True, "w7_40"),
                                          (5, 50, 3, True, "w5_50"), (7, 30, 2, False, "wTP")]:
                d2 = tryptic_detectable_mask(s, lo, hi, mc, bp)
                rec[f"ndet_{tag}"] = int(d2[cpos].sum()) if cpos else 0
                rec[f"smap_{tag}"] = int(any(d2[p] for p in sites))
        rows.append(rec)
    G = pd.DataFrame(rows)
    G["n_undet"] = G.n_cys - G.n_det
    G.to_csv(OUT / "groups_r4.csv", index=False, lineterminator="\n")

    # ---------------- compare with the item's table
    F = pd.read_csv(RES / "s08_group_features_windows.csv")
    M = G.merge(F[["group_id", "A", "S", "S1", "S_map", "log10_ibaq", "length", "n_cys", "n_det",
                   "ndet_w6_35_mc2", "ndet_w7_40_mc3", "ndet_w5_50_mc3", "ndet_w7_30_mc2_trypsinP", "quantity"]],
                on="group_id", how="outer", suffixes=("", "_item"), indicator=True)
    comp = {"merge": M["_merge"].value_counts().to_dict()}
    both = M[M["_merge"] == "both"]
    for c in ["A", "S", "S1", "S_map", "length", "n_cys", "n_det"]:
        comp[f"mismatch_{c}"] = int((both[c].astype(float) != both[c + "_item"].astype(float)).sum())
    comp["max_abs_diff_log10_ibaq"] = float((both.log10_ibaq - both.log10_ibaq_item).abs().max())
    comp["max_rel_diff_quantity"] = float(((both.quantity - both.quantity_item).abs() / both.quantity_item).max())
    for mine, theirs in [("ndet_w6_35", "ndet_w6_35_mc2"), ("ndet_w7_40", "ndet_w7_40_mc3"),
                         ("ndet_w5_50", "ndet_w5_50_mc3"), ("ndet_wTP", "ndet_w7_30_mc2_trypsinP")]:
        comp[f"mismatch_{mine}"] = int((both[mine] != both[theirs]).sum())
    F2 = F.set_index("group_id")
    for mine, theirs in [("smap_w6_35", "smap_w6_35_mc2"), ("smap_w7_40", "smap_w7_40_mc3"),
                         ("smap_w5_50", "smap_w5_50_mc3"), ("smap_wTP", "smap_w7_30_mc2_trypsinP")]:
        comp[f"mismatch_{mine}"] = int((G.set_index("group_id")[mine] != F2.loc[G.group_id, theirs].to_numpy()).sum())
    summ = {
        "label": "VERIFIER round 4 (POST HOC), independent rebuild",
        "inputs_sha256": {"report": sha(PRO), "peptides": sha(PEP), "fasta": sha(FASTA), "fetched": sha(FETCHED)},
        "n_ref_entries": n_ref,
        "report_rows": int(len(rep)), "groups_dropped_no_ibaq": n_drop, "groups": int(len(G)),
        "accession_lists_alphabetical": alpha_sorted,
        "labels": {k: int(G[k].sum()) for k in ("A", "S", "S1", "S_map")},
        "A_not_S": int(((G.A == 1) & (G.S == 0)).sum()),
        "with_sequence": int(G.has_seq.sum()),
        "n_det_ge1": int((G.n_det >= 1).sum()), "no_cys": int((G.n_cys == 0).sum()),
        "frac_cys_detectable": float(G.n_det.sum() / G.n_cys.sum()),
        "frac_groups_no_undet": float((G.n_undet == 0).mean()),
        "negatives": int((G.A == 0).sum()),
        "negatives_with_an_accession_in_Proteins_column": int(((G.A == 0) & (G.in_proteins_col_any == 1)).sum()),
        "comparison_with_item": comp,
    }
    (OUT / "build_r4.json").write_text(json.dumps(summ, indent=1), encoding="utf-8")
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
