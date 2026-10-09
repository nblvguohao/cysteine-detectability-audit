"""POST HOC adversarial verification, round 4, item J_provenance_species (2026-09-30).

Extra checks: (1) T:pep_cys_count (proline-blocked tryptic peptide, as in v2_features.py) as a
reader of rice non-coverability, pooled and within protein; (2) the rice step at the ceiling;
(3) every token of artifact12_number_registry.csv on its stated manuscript line;
(4) the rice protein-by-protein reference figures and the set they refer to.
Writes only to results/J_provenance_species/verify_r4/.
"""
import csv
import gzip
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
W = Path(r"/path/to/local/_cys_repo_work/public/revision_2026-09-30")
REPO = Path(r"/path/to/local/_cys_repo_work/repo")
MCP = Path(r"/path/to/local/巯基化/MCP")
OUT = W / "results" / "J_provenance_species" / "verify_r4"


def auc(y, s):
    y = np.asarray(y, int)
    s = np.asarray(s, float)
    pos, neg = s[y == 1], s[y == 0]
    allv = np.concatenate([pos, neg])
    order = np.argsort(allv, kind="mergesort")
    sv = allv[order]
    ranks = np.empty(sv.size)
    i = 0
    while i < sv.size:
        j = i
        while j + 1 < sv.size and sv[j + 1] == sv[i]:
            j += 1
        ranks[i:j + 1] = (i + j) / 2 + 1
        i = j + 1
    r = np.empty(sv.size)
    r[order] = ranks
    return float((r[:pos.size].sum() - pos.size * (pos.size + 1) / 2) / (pos.size * neg.size))


def load_rice():
    ref, acc, buf = {}, None, []
    with gzip.open(W / "external/UP000059680_39947.fasta.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    ref[acc] = "".join(buf)
                acc, buf = line[1:].split("|")[1], []
            else:
                buf.append(line.strip())
    ref[acc] = "".join(buf)
    with open(W / "results/F_rice_artifact3/fetched/fetched_sequences.tsv", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["source"] == "uniprot_rest_active" and r["requested_accession"] not in ref:
                ref[r["requested_accession"]] = r["sequence"].strip().upper()
    rows = []
    with open(REPO / "external/public_cohort_site_scores.csv", encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["dataset"] == "rice_PXD072089":
                rows.append((r["accession"], int(r["position"]), int(r["observed_in_retrospective_dataset"])))
    return ref, rows


def pb_boundaries(seq):
    n = len(seq)
    b = [0]
    for i, a in enumerate(seq):
        if a in "KR" and (i + 1 >= n or seq[i + 1] != "P"):
            b.append(i + 1)
    if b[-1] != n:
        b.append(n)
    return b


def minimal_nc(seq, i):
    start = 0
    for j in range(i - 1, -1, -1):
        if seq[j] in "KR":
            start = j + 1
            break
    end = len(seq) - 1
    for j in range(i + 1, len(seq)):
        if seq[j] in "KR":
            end = j
            break
    return seq[start:end + 1].count("C") >= 2


def main():
    import bisect
    res = {"label": "POST HOC adversarial verification round 4, extra checks"}
    ref, rows = load_rice()
    nc, pcc, lab, accs = [], [], [], []
    for a, p, y in rows:
        s = ref[a]
        i = p - 1
        b = pb_boundaries(s)
        k = bisect.bisect_right(b, i) - 1
        k = min(max(k, 0), len(b) - 2)
        pcc.append(s[b[k]:b[k + 1]].count("C"))
        nc.append(int(minimal_nc(s, i)))
        lab.append(y)
        accs.append(a)
    nc, pcc, lab, accs = np.array(nc), np.array(pcc), np.array(lab), np.array(accs)
    res["pep_cys_count"] = {
        "all_nc_have_ge2": bool((pcc[nc == 1] >= 2).all()),
        "coverable_with_ge2": int(((pcc >= 2) & (nc == 0)).sum()),
        "coverable_n": int((nc == 0).sum()),
        "labelled_with_ge2": int(((pcc >= 2) & (lab == 1)).sum()),
        "auc_for_nc": auc(nc, pcc),
        "auc_for_nc_among_unlabelled": auc(nc[lab == 0], pcc[lab == 0]),
    }
    # within-protein AUC of pep_cys_count (negated) for the label (what the detect-only set can reach from it alone)
    by = defaultdict(list)
    for k, a in enumerate(accs):
        by[a].append(k)
    per = []
    for a, idx in by.items():
        yy = lab[idx]
        if yy.min() == yy.max():
            continue
        per.append(auc(yy, -pcc[idx]))
    res["pep_cys_count"]["within_protein_auc_for_label_negated"] = float(np.mean(per))
    res["pep_cys_count"]["n_proteins"] = len(per)
    # (2) ceiling steps
    def lor(k, N, D):
        return math.log2(((k + .5) / (100 - k + .5)) / ((D - k + .5) / ((N - 100) - (D - k) + .5)))
    res["ceiling_step_99_100"] = {"human": round(lor(100, 18567, 14670) - lor(99, 18567, 14670), 4),
                                  "rice": round(lor(100, 5004, 3803) - lor(99, 5004, 3803), 4)}
    # (3) registry tokens on their lines
    ms = (MCP / "01_manuscript.tex").read_text(encoding="utf-8").split("\n")
    reg = list(csv.DictReader(open(W / "results/J_provenance_species/artifact12_number_registry.csv", encoding="utf-8-sig")))
    miss = []
    def norm(t):
        return t.replace("−", "$-$")
    for r in reg:
        line = ms[int(r["tex_line"]) - 1]
        tok = r["printed"]
        cands = {tok, norm(tok), tok.replace(",", "{,}"), tok.lstrip("+"), tok.replace("-", "$-$", 1) if tok.startswith("-") else tok}
        if not any(c in line for c in cands):
            miss.append((r["entry_id"], r["tex_line"], tok))
    res["registry"] = {"n": len(reg), "lines": sorted({int(r["tex_line"]) for r in reg}), "tokens_not_on_line": miss}
    # (4) rice protein-draw reference figures
    rc = json.loads((W / "results/J_provenance_species/ablation_top100_rice_cluster_sensitivity.json").read_text(encoding="utf-8"))
    res["rice_cluster_sensitivity_file"] = rc
    with open(OUT / "extra_verify_r4.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print(json.dumps(res, indent=1, ensure_ascii=False)[:8000])


if __name__ == "__main__":
    main()
