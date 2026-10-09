"""Step 9 (POST HOC, revision after verification, 2026-09-30): does Artifact 2's rice site count
(756 sites / 600 proteins asymmetric; 689 / 545 symmetric; stored by second-tree script 36 per the
repository's baseline registry) follow from the stored 640-protein / 817-site self-audit list?

Why: the verifier (round 1, minor 8) noted that conflict 12 asserted this dependence without checking.
The repository's rerun audit (baseline/rerun_correction_2026-09-21_audit.json and the docstring of
scripts/baseline_phase0_rerun_correction_2026-09-21.py) documents that script 36 reads script 24's
640-protein product and that its outputs change (687/624) when script 24 loses its per-accession
sequence file (581 proteins). This script tries to rebuild 756/600 and 689/545 from the stored
817-site list under a grid of detectability rules. It is exploratory: the original rule and
the 8-angstrom structure criterion (which needs a model per protein) are not available, so a miss does
not refute the dependence and a hit is reported with the rule that produced it.

Output: s09_artifact2_trace.csv, s09_artifact2_trace.json
"""
import itertools
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
sys.path.insert(0, str(CYS_AUDIT_SRC))
from cys_audit.proteases import digest

T = pd.read_csv(REPO_SITE_SCORES, usecols=["dataset", "accession", "position", "observed_in_retrospective_dataset"],
                encoding="utf-8-sig")
T = T[T.dataset == "rice_PXD072089"].copy()
seqs, _ = load_all_sequences()
prots = sorted(T.accession.unique())
sites = {a: set(g.loc[g.observed_in_retrospective_dataset == 1, "position"].astype(int))
         for a, g in T.groupby("accession")}
assert len(prots) == 640 and sum(len(v) for v in sites.values()) == 817

# the complete list (s05 step 4: single-Cys peptides, isoforms -> canonical, proteins with >=2 Cys,
# BEFORE the loss of UniProt-deleted accessions): 1,341 sites / 1,040 proteins
pep = pd.read_csv(PEP, sep="\t", low_memory=False, dtype=str)
pep["C Count"] = pep["C Count"].astype(int)
pep["lrp"] = pep["Leading razor protein"].fillna("")
bad = (pep["Reverse"].fillna("") == "+") | (pep["Potential contaminant"].fillna("") == "+") | \
      pep.lrp.str.startswith("REV__") | pep.lrp.str.startswith("CON__")
single = pep[~bad & (pep["C Count"] == 1)]
full = {}
for s_, a_ in zip(single.Sequence, single.lrp):
    c = a_.split("-")[0]
    if c not in seqs:
        continue
    ps = seqs[c]; k = ps.find(s_)
    while k >= 0:
        full.setdefault(c, set()).add(k + 1 + s_.index("C"))
        k = ps.find(s_, k + 1)
full = {a: v for a, v in full.items() if seqs[a].count("C") >= 2}
assert len(full) == 1040 and sum(len(v) for v in full.values()) == 1341
assert all(sites.get(a, set()) <= full.get(a, set()) for a in prots)


def count(site_map, enz, mc, lo, hi):
    a_s = a_p = s_s = s_p = a_s_all = 0
    for acc in sorted(site_map):
        s = seqs[acc]
        det = np.zeros(len(s), bool)
        reach1 = set()
        for x, y in digest(s, enz, missed=mc, min_len=lo, max_len=hi):
            det[x - 1:y] = True
            sub = s[x - 1:y]
            if sub.count("C") == 1:
                reach1.add(x + sub.index("C"))
        cpos = [i + 1 for i, ch in enumerate(s) if ch == "C"]
        st = site_map[acc]
        if sum(det[p - 1] for p in cpos) >= 2:
            k = sum(det[p - 1] for p in st)
            if k:
                a_s += k; a_p += 1
            a_s_all += len(st)
        if sum(p in reach1 for p in cpos) >= 2:
            k = sum(p in reach1 for p in st)
            if k:
                s_s += k; s_p += 1
    return a_s, a_p, a_s_all, s_s, s_p


TARGETS = {"asym_sites": 756, "asym_prots": 600, "sym_sites": 689, "sym_prots": 545}
complete_rows = []
for enz, mc, lo, hi in (("trypsin", 0, 6, 30), ("trypsin", 2, 7, 30)):
    for nm, mp in (("stored_640", {a: sites.get(a, set()) for a in prots}), ("complete_1040", full)):
        a_s, a_p, a_all, s_s, s_p = count(mp, enz, mc, lo, hi)
        complete_rows.append({"list": nm, "enzyme": enz, "missed": mc, "min_len": lo, "max_len": hi,
                              "asym_sites": a_s, "asym_prots": a_p, "sym_sites": s_s, "sym_prots": s_p})
pd.DataFrame(complete_rows).to_csv(RES / "s09_artifact2_stored_vs_complete.csv", index=False)
print(pd.DataFrame(complete_rows).to_string())
rows = []
for enz, mc, lo, hi in itertools.product(("trypsin", "trypsin/p"), (0, 1, 2), (6, 7), (25, 30, 35, 40, 50)):
    a_s = a_p = s_s = s_p = 0
    a_s_all = 0
    for acc in prots:
        s = seqs[acc]
        det = np.zeros(len(s), bool)
        reach1 = set()
        for x, y in digest(s, enz, missed=mc, min_len=lo, max_len=hi):
            det[x - 1:y] = True
            sub = s[x - 1:y]
            if sub.count("C") == 1:
                reach1.add(x + sub.index("C"))
        cpos = [i + 1 for i, ch in enumerate(s) if ch == "C"]
        ndet = sum(det[p - 1] for p in cpos)
        st = sites.get(acc, set())
        if ndet >= 2:
            k = sum(det[p - 1] for p in st)
            if k:
                a_s += k; a_p += 1
            a_s_all += len(st)
        n1 = sum(p in reach1 for p in cpos)
        if n1 >= 2:
            k = sum(p in reach1 for p in st)
            if k:
                s_s += k; s_p += 1
    rows.append({"enzyme": enz, "missed": mc, "min_len": lo, "max_len": hi,
                 "asym_sites_on_det_in_prots_ge2det": a_s, "asym_prots_with_site": a_p,
                 "all_sites_in_prots_ge2det": a_s_all,
                 "sym_sites_on_single_cys_reach_in_prots_ge2": s_s, "sym_prots_with_site": s_p})
R = pd.DataFrame(rows)
R["hits"] = [
    ";".join(k for k, v in (("756", r.asym_sites_on_det_in_prots_ge2det == 756 or r.all_sites_in_prots_ge2det == 756),
                            ("600", r.asym_prots_with_site == 600), ("689", r.sym_sites_on_single_cys_reach_in_prots_ge2 == 689),
                            ("545", r.sym_prots_with_site == 545)) if v)
    for r in R.itertuples()]
R.to_csv(RES / "s09_artifact2_trace.csv", index=False)
closest = R.assign(d=(R.asym_sites_on_det_in_prots_ge2det - 756).abs()).sort_values("d").head(5)
write_json(RES / "s09_artifact2_trace.json", {
    "analysis_label": "POST HOC revision analysis 2026-09-30 (after verification); exploratory rule grid",
    "targets": TARGETS, "n_rules": int(len(R)), "rules_with_any_hit": R[R.hits != ""].to_dict("records"),
    "closest_to_756": closest.drop(columns="d").to_dict("records"),
    "range_asym_sites": [int(R.asym_sites_on_det_in_prots_ge2det.min()), int(R.asym_sites_on_det_in_prots_ge2det.max())],
    "range_sym_sites": [int(R.sym_sites_on_single_cys_reach_in_prots_ge2.min()), int(R.sym_sites_on_single_cys_reach_in_prots_ge2.max())],
    "stored_vs_complete_list": complete_rows,
    "repository_evidence": {
        "baseline/baseline_metrics.csv rows A-047..A-049": "756 sites / 600 proteins (asymmetric) and 689 / 545 "
            "(symmetric) are stored by second-tree results/21_crossspecies_summary.json, 35_symmetric_ncys_test.json "
            "and 36_ncys_decompose.json",
        "scripts/baseline_phase0_rerun_correction_2026-09-21.py docstring": "without results/24_extra_sequences.json "
            "(sequences retrieved one by one by accession) script 24 wrote 24_mapped_osa_persulfidation.json with "
            "581 proteins instead of 640 and every downstream row of 36_ncys_decompose differed",
        "baseline/rerun_comparison_2026-09-21.csv rows tree2_36_ncys": "n sites 756 -> 687 and 689 -> 624 in that run",
        "baseline/rerun_correction_2026-09-21_audit.json": "with the file restored, 24 and 36 outputs IDENTICAL"},
    "inputs": {"site_table": str(REPO_SITE_SCORES), "sha256": sha256_file(REPO_SITE_SCORES),
               "peptides": str(PEP), "peptides_sha256": sha256_file(PEP)}})
pd.set_option("display.width", 250)
print(R.to_string())
