# -*- coding: utf-8 -*-
"""J7. Revision round 3 checks for item J_provenance_species.

POST HOC revision analysis (2026-09-30), revision round 3; not registered, not pre-specified.

A  Manuscript anchoring. The round-2 report quoted line numbers of an earlier version of 01_manuscript.tex
   (sha256 32242ded..., 475 lines). This part locates every target of the report in the manuscript on disk by
   a unique anchor phrase (exactly one hit is asserted), records the new line, and checks that every passage
   the report quotes as "current text" still occurs verbatim on that line. The same is done for the
   Supplemental Note and Supplemental Data pointers. A4 lists, read-only, the line pointers of the sibling
   items' reports that fall on the same lines, so the numbering can be compared.
B  Human window power. The human coverability check is window-based. Among cysteines whose 15-residue window
   holds the whole minimal Trypsin/P segment (j1 flag trypsinP_complete) the check is exact; this part
   counts non-coverable cysteines by label there and compares the completeness of the two label classes.
C  Rice: can the chemistry arms read which cysteines are never labelled? Three columns that every chemistry
   arm (including both deletion arms) carries are recomputed from full UniProt 2026_03 sequences with the
   definitions of repo/scripts/v2_features.py::_site_features (B:cys_dist_min, B:cys_in_4,
   B:win10_count_cys), and their AUC for the non-coverable flag of j5 is computed over the 5,598 released
   rice cysteines, with a protein-clustered bootstrap (5,000 replicates, seed 20260930 + 3).
D  Row sets of Supplemental Note 11 (rice). Reads the repository's refit analysis script, the post hoc paired-
   difference script and its output, and Note 11, and records which row set each rice number rests on: the arm
   tables use the rows with a sequence (the `keep` mask of analyse_public_refit_2026-09-21.py), the blend's
   paired differences read the out-of-fold files without that mask.

Outputs (results/J_provenance_species/): manuscript_anchor_map.csv, supplement_anchor_map.csv,
sibling_line_pointers.csv, round3_checks.json, _provenance_fragment_j7_round3_checks.json.
"""
from __future__ import annotations

import csv
import datetime
import glob
import os
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from jcommon import MCP, OUT, REPO, REPORTS, W, Inputs, sha256_file, sha256_lf, write_csv, write_fragment, write_json  # noqa: E402
from j3_ablation_statistic import cluster_weights, make_weighted_auc  # noqa: E402

SEED = 20260930 + 3
BOOT = 5000
TEX = MCP + "/01_manuscript.tex"
SUPP = MCP + "/supplemental"
OLD_TEX_SHA = "32242ded279cd1f5648dce59efa7b160dfc926e5dc6012d383c67a153d6991a1"
J1_SITES = OUT + "/species_peptide_check_sites.csv"
J5_SITES = OUT + "/singlecys_construction_sites.csv"
RICE_FASTA = W + "/external/UP000059680_39947.fasta.gz"
RICE_FETCHED = W + "/results/F_rice_artifact3/fetched/fetched_sequences.tsv"   # item F, read-only
F_CHAIN_R2 = W + "/results/F_rice_artifact3/s11_rice_chain.csv"                 # item F round 2, read-only (cited)

# ------------------------------------------------------------------------------------------------ part A
# (id, where the report uses it, line in the round-2 report, anchor, [verbatim quotes the report uses])
POINTERS = [
    ("abstract", "proposed text 6; conflicts", 27,
     "Deleting every suspect feature from a site ranker of our own left its cleavage-geometry ordering intact",
     ["Deleting every suspect feature from a site ranker of our own left its cleavage-geometry ordering intact"]),
    ("methods_permutation", "J2 registry; conflicts", 54,
     "the permutation P of 0.7435 agreed with the analytic value",
     ["the permutation P of 0.7435 agreed with the analytic value"]),
    ("ep_public_cohorts", "proposed text 1; conflicts", 62,
     "The self-audit uses two public persulfidation deposits",
     ["The self-audit uses two public persulfidation deposits: human (PXD044043; 18,567 cysteines in 1,519 proteins, "
      "3,206 of them reported in the deposit's enriched site list)",
      "The label is \\emph{reported in the deposit's site list}",
      "(a NEG\\_A negative in the tiers of the label-semantics table)",
      "Folds and bootstrap clusters are proteins, because no homology components exist for these cohorts; intervals are "
      "protein-clustered bootstrap (5,000 replicates, fixed seed) and are narrower than component clustering would give"]),
    ("ep_ablation", "proposed text 2; conflicts", 87,
     "For the three-arm ablation, the verdict statistic",
     ["the verdict statistic, the branches for reading it, and the requirement that the label-permuted control produce no "
      "enrichment were fixed in the script's documentation before any arm was fitted"]),
    ("edsr_revision_notes", "round 3: note placement", None,
     "The analyses added in revision (Supplemental Notes 14--20) are post hoc responses to review",
     ["The analyses added in revision (Supplemental Notes 14--20) are post hoc responses to review and are labelled as such."]),
    ("artifact1_raw", "J2 registry; results 2(d)", 130,
     "The raw enrichment reached z = +15.72",
     ["the largest contributors are PMID 37277371 and 25699590"]),
    ("artifact1_profile", "proposed text 8; conflicts", 132,
     "in a public Arabidopsis S-nitrosylation dataset the K/R fraction at position",
     ["in a public Arabidopsis S-nitrosylation dataset the K/R fraction at position $-$1 is 0.4576"]),
    ("artifact1_matching", "J2 registry", 134,
     "Restricting the background to cysteines that are themselves theoretically detectable", ["+5.59 $\\rightarrow$ +0.67"]),
    ("artifact1_protease", "J2 registry", 136,
     "Of the eight datasets analyzed in this study, six were documented", ["187, 728 and 884 identified cysteines against 5,328"]),
    ("artifact2_values", "conflicts; results 2(f)", 141,
     "The human persulfidation dataset (PXD044043 [48]) behaves the same way",
     ["The human persulfidation dataset (PXD044043 [48]) behaves the same way: $-$8.78",
      "from $-$0.98 to $-$11.09", "from +4.53 to $-$3.64"]),
    ("artifact2_remedy", "proposed text 9; conflicts; section 1b", 143,
     "The remedy must be applied symmetrically",
     ["\\textbf{The remedy must be applied symmetrically:} either both sets are filtered or the feature is not reported.",
      "On these grounds we do not report neighboring-cysteine values for two peptide-inferred datasets; datasets with "
      "instrument-assigned site localisation are unaffected and need no filter."]),
    ("cohorts_of_self_audit", "proposed text 6; conflicts", 218,
     "\\textbf{Cohorts of the self-audit.}",
     ["neither has homology clustering, so cross-validation is protein-disjoint but not verified homology-disjoint",
      "(PXD072089; the deposit's own site list reports 817 sites, before the exclusions of Experimental Procedures)",
      "Our ranker was trained on each of two public cohorts in turn", "but not verified homology-disjoint."]),
    ("self_audit_zero_shot", "proposed text 10; conflicts", 221,
     "We applied the same instruments to a site-ranking model",
     ["Within each stratum the label is constant, so this ordering cannot be inherited from it: a binary annotation can "
      "separate label classes, but it carries no information about which cysteine within a class lies nearer a cleavage site."]),
    ("ablation_intro", "conflicts", 229,
     "One explanation remained",
     ["One explanation remained", "the verdict statistic was fixed before any arm was fitted"]),
    ("ablation_result", "proposed text 3; conflicts", 231,
     "All six ablation intervals exclude zero, so the ordering survives deletion of the named columns",
     ["All six ablation intervals exclude zero, so the ordering survives deletion of the named columns",
      "The rice cohort reproduces the human direction in all three arms.",
      "A detectability bias therefore cannot be removed by deleting features whose names look detectability-related, "
      "because the signal remains reachable through the remaining columns; control has to be built into the evaluation design."]),
    ("detectability_only", "proposed text 11; conflicts", 234,
     "A \\textbf{25-feature set that contains no chemistry at all}",
     ["A high recovered share is what these labels predict, because they are close to ``was this residue reported’’; "
      "the result is a statement about the labels, not evidence that no chemical information exists.",
      "Detectability is not a contaminant in this task; it is most of the task"]),
    ("discussion_self_audit", "proposed text 6; conflicts", 243,
     "The self-audit carries the same lesson to predictive models",
     ["In our own ranker, a detectability-only feature set recovered almost all of the within-protein discrimination of the "
      "full model on the human cohort, while deleting features that looked detectability-related left the cleavage-geometry "
      "ordering in place;"]),
    ("discussion_residual", "J2 registry", 245,
     "The residual in Artifact 1 is the quantitative case for that ordering", ["a reduction of 48\\% and 53\\%"]),
    ("data_availability_pride", "proposed text 7", 259,
     "The PRIDE [45] deposits re-analyzed here are",
     ["PXD044043 [48,72] and PXD072089 [47,73] (Artifact 2 and self-audit)"]),
    ("data_availability_source", "proposed text 7; conflicts", 261,
     "Source data for Figures 1--7 are provided as Supplemental Source Data files",
     ["the analysis scripts are in the code repository (Data Availability)",
      "Source data for Figures 1--7 are provided"]),
    ("code_availability", "conflicts", 263,
     "The Cys-Audit tool (version 0.2.2, tag v3.0.0), with its 55-test suite",
     ["together with the analysis scripts, the pre-registration protocols and a SHA-256 manifest of every file"]),
    ("supplemental_notes_list", "round 3: note placement", None,
     "\\textbf{Supplemental Notes 1--13.}",
     ["the feature ablation (Note 11)",
      "and the non-null check of the normal approximation behind the binary-feature statistic (Note 13)"]),
    ("table2_header", "proposed text 4 (Table 2 block, first line)", 447,
     "cohort & arm & columns & top-100 distal log2 OR & within-protein AUC",
     ["cohort & arm & columns & top-100 distal log2 OR & within-protein AUC"]),
    ("table2_caption", "proposed text 4; conflicts (Table 2 block, last line)", 456,
     "\\textbf{Table 2 | }Three-arm feature ablation on the two public cohorts",
     ["Values are top-100 distal log2 odds ratios with cluster-bootstrap intervals and within-protein AUC, read through the "
      "LightGBM ranking member of each arm (Experimental Procedures)."]),
    ("fig7_legend", "proposed text 5; conflicts", 472,
     "\\noindent\\textbf{Figure 7.} Self-audit of the ranker on the two public cohorts. Color marks",
     ["read through the LightGBM ranking member refit within each cohort (Supplemental Note 11)",
      "clustered by protein because these cohorts have no homology grouping", "Source data are provided."]),
]

SUPP_POINTERS = [
    ("SN2", "Supplemental_Note_2_public_ranker_training.md", 3, "PRIDE PXD044043 (human) to PXD072089 (rice)",
     ["(a site reported in the deposit's site list; an unreported cysteine is not known to be unmodified)"]),
    ("SN2", "Supplemental_Note_2_public_ranker_training.md", 5, "3,206 reported sites",
     ["3,206 reported sites", "PXD072089: 5,598 residues in 640 proteins in the deposit"]),
    ("SN2", "Supplemental_Note_2_public_ranker_training.md", 7, "Neither cohort has homology clustering",
     ["Neither cohort has homology clustering", "the deposit's own site list reports 817 sites.",
      "Because each direction is audited on the other species, no homology leakage between training and audit cohorts "
      "is possible."]),
    ("SN9", "Supplemental_Note_9_public_cleavage_geometry_self_audit.md", 3, "Labels are 'reported in the deposit's site list'",
     ["Labels are 'reported in the deposit's site list' (NEG_A: unreported is not known to be unmodified)."]),
    ("SN9", "Supplemental_Note_9_public_cleavage_geometry_self_audit.md", 67, "because neither cohort has homology grouping",
     ["Bootstrap clusters are proteins, not homology components, because neither cohort has homology grouping. "
      "Intervals are narrower than component clustering would give and must be quoted with this note."]),
    ("SN10", "Supplemental_Note_10_public_five_protease_self_audit.md", 3, "Labels are 'reported in the deposit's site list'",
     ["Labels are 'reported in the deposit's site list' (NEG_A)."]),
    ("SN11", "Supplemental_Note_11_public_feature_ablation.md", 5, "The labels are 'reported in the deposit's site list'",
     ["The labels are 'reported in the deposit's site list' (a NEG_A negative: unreported is not known to be unmodified).",
      "Bootstrap clusters are proteins, not homology components, because none exist for these cohorts; intervals are "
      "therefore narrower than component clustering would give."]),
    ("SN11", "Supplemental_Note_11_public_feature_ablation.md", 26, "All six ablation intervals exclude zero in both cohorts",
     ["All six ablation intervals exclude zero in both cohorts", "and the label-permuted control crosses zero in both."]),
    ("SN12", "Supplemental_Note_12_public_detectability_only.md", 5, "The labels are 'reported in the deposit's site list'",
     ["The labels are 'reported in the deposit's site list' (a NEG_A negative: unreported is not known to be unmodified).",
      "Bootstrap clusters are proteins, not homology components, because none exist for these cohorts; intervals are "
      "therefore narrower than component clustering would give."]),
    ("SN12", "Supplemental_Note_12_public_detectability_only.md", 23, "a high recovery ratio is the expected result",
     ["a high recovery ratio is the expected result"]),
    ("SN13", "Supplemental_Note_13_nonnull_permutation_check.md", 8, "in Artifact 2; +15.72 to",
     ["(z = −2.61, +3.34, −8.78, +2.67 in Artifact 2; +15.72 to −0.98 in Artifact 1)",
      "+15.72 to −0.98 in Artifact 1"]),
    ("SD12", "Supplemental_Data_12_analysis_provenance.csv", 2, "A1_matching,", ["A1_matching,Artifact 1"]),
    ("SD12", "Supplemental_Data_12_analysis_provenance.csv", 3, "A2_filter,", ["A2_filter,Artifact 2"]),
    ("SD12", "Supplemental_Data_12_analysis_provenance.csv", 13, "Selfaudit,", ["Selfaudit,"]),
    ("SD12", "Supplemental_Data_12_analysis_provenance.csv", 14, "Ablation,", ["Ablation,Three-arm feature ablation"]),
    ("SD12", "Supplemental_Data_12_analysis_provenance.csv", 15, "Detect_only,", ["Detect_only,Detectability-only refit"]),
]


def find_unique(lines, anchor):
    hits = [i + 1 for i, t in enumerate(lines) if anchor in t]
    return hits


def part_a(inp, res):
    text = inp.read_text(TEX)
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines = lines[:-1]
    st = os.stat(TEX)
    res["A_manuscript"] = {
        "path": "MCP/01_manuscript.tex", "sha256": sha256_file(TEX), "sha256_lf": sha256_lf(TEX),
        "lines": len(lines), "bytes": st.st_size,
        "last_modified_local": datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M"),
        "round2_report_stated_sha256": OLD_TEX_SHA,
        "round2_report_stated_version_matches_disk": sha256_file(TEX) == OLD_TEX_SHA,
    }
    rows = []
    for pid, where, old, anchor, quotes in POINTERS:
        hits = find_unique(lines, anchor)
        assert len(hits) == 1, (pid, anchor, hits)
        new = hits[0]
        qok = [q in lines[new - 1] for q in quotes]
        rows.append({"pointer": pid, "used_in_report": where, "line_round2_report": old if old is not None else "",
                     "line_current": new, "offset": (new - old) if old is not None else "",
                     "anchor": anchor, "quotes_checked": len(quotes), "quotes_verbatim_on_line": sum(qok),
                     "all_quotes_verbatim": int(all(qok)),
                     "missing_quotes": " || ".join(q for q, ok in zip(quotes, qok) if not ok)})
    out1 = write_csv(OUT + "/manuscript_anchor_map.csv", rows)
    # Table 2 block range
    t2_begin = max(i + 1 for i, t in enumerate(lines[:rows[[r["pointer"] for r in rows].index("table2_header")]["line_current"]])
                   if t.startswith("\\begin{table}"))
    t2_end = min(i + 1 for i, t in enumerate(lines) if i + 1 > t2_begin and t.startswith("\\end{table}"))
    res["A_table2_block"] = {"begin_table": t2_begin, "end_table": t2_end,
                             "rows_human_rice": [i + 1 for i, t in enumerate(lines) if re.match(r"^(human|rice) & ", t)]}
    res["A_pointers"] = {"n": len(rows), "all_unique": True,
                         "all_quotes_verbatim": all(r["all_quotes_verbatim"] for r in rows),
                         "quotes_checked": sum(r["quotes_checked"] for r in rows),
                         "offsets": {r["pointer"]: r["offset"] for r in rows}}
    # supplements
    srows = []
    for note, fname, old, anchor, quotes in SUPP_POINTERS:
        p = SUPP + "/" + fname
        sl = inp.read_text(p).split("\n")
        hits = find_unique(sl, anchor)
        new = hits[0] if len(hits) == 1 else None
        qok = [new is not None and q in sl[new - 1] for q in quotes]
        srows.append({"file": fname, "note": note, "line_round2_report": old, "line_current": new if new else "",
                      "anchor_hits": len(hits), "all_quotes_verbatim": int(all(qok)), "sha256": sha256_file(p),
                      "missing_quotes": " || ".join(q for q, ok in zip(quotes, qok) if not ok)})
    out2 = write_csv(OUT + "/supplement_anchor_map.csv", srows)
    res["A_supplements"] = {"n": len(srows), "quotes_checked": sum(len(q) for *_, q in SUPP_POINTERS),
                            "all_unique": all(r["anchor_hits"] == 1 for r in srows),
                            "all_at_stated_line": all(r["line_current"] == r["line_round2_report"] for r in srows),
                            "all_quotes_verbatim": all(r["all_quotes_verbatim"] for r in srows)}
    # A4 sibling reports (read-only): pointers to the same lines
    targets = {r["line_current"] for r in rows}
    t2 = set(range(t2_begin, t2_end + 1))
    sib = []
    for f in sorted(glob.glob(REPORTS + "/*.md")):
        base = os.path.basename(f)
        if base.startswith("J_") or base == "DOWNLOADS.md":
            continue
        t = inp.read_text(f)
        stated = re.search(r"(79ce7159|32242ded|f79cc644)[0-9a-f]*", t)
        for m in re.finditer(r"\b[Ll]ines? (\d{2,3})\b", t):
            n = int(m.group(1))
            if n in targets or n in t2:
                ctx = t[max(0, m.start() - 70):m.end() + 90].replace("\n", " ").replace("|", "/")
                sib.append({"report": base, "manuscript_version_stated": stated.group(0)[:8] if stated else "",
                            "pointer": m.group(0), "line": n,
                            "this_report_target": ";".join(r["pointer"] for r in rows if r["line_current"] == n) or
                            ("table2_block" if n in t2 else ""),
                            "context": ctx})
    out3 = write_csv(OUT + "/sibling_line_pointers.csv", sib)
    res["A4_sibling_reports"] = {"pointers_on_shared_lines": len(sib),
                                 "reports": sorted({s["report"] for s in sib}),
                                 "versions_stated": {s["report"]: s["manuscript_version_stated"] for s in sib}}
    return [out1, out2, out3]


# ------------------------------------------------------------------------------------------------ part B
def part_b(inp, res):
    j1 = {(r["accession"], int(r["position"])): r for r in inp.read_csv(J1_SITES)}
    j5 = [r for r in inp.read_csv(J5_SITES) if r["dataset"] == "human_PXD044043"]
    assert len(j1) == len(j5) == 18567
    out = {}
    for cls, lab in (("labelled", "1"), ("unlabelled", "0")):
        sub = [r for r in j5 if r["label"] == lab]
        comp = [r for r in sub if j1[(r["accession"], int(r["position"]))]["trypsinP_complete"] == "1"]
        inc = [r for r in sub if j1[(r["accession"], int(r["position"]))]["trypsinP_complete"] == "0"]
        seglen = np.asarray([int(j1[(r["accession"], int(r["position"]))]["trypsinP_len"]) for r in comp])
        out[cls] = {"n": len(sub), "window_holds_whole_minimal_segment": len(comp),
                    "share_complete": len(comp) / len(sub),
                    "not_coverable_among_complete": sum(int(r["window_second_cys_trypsinP"]) for r in comp),
                    "not_coverable_among_incomplete_windows": sum(int(r["window_second_cys_trypsinP"]) for r in inc),
                    "median_complete_segment_length": float(np.median(seglen))}
        out[cls]["share_not_coverable_among_complete"] = out[cls]["not_coverable_among_complete"] / len(comp)
    out["reading"] = ("the window sees the whole minimal segment equally often in both label classes, so the contrast 0 of "
                      "3,206 against 26% is not an artefact of window truncation; where the window is complete the check is "
                      "exact")
    res["B_human_window_power"] = out


# ------------------------------------------------------------------------------------------------ part C
def read_fasta_gz(path):
    import gzip
    out, acc, buf = {}, None, []
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    out[acc] = "".join(buf)
                m = re.match(r">\w+\|([^|]+)\|", line)
                acc = m.group(1) if m else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        out[acc] = "".join(buf)
    return out


def site_cys_columns(seq, position):
    """B:cys_dist_min, B:cys_in_4 and B:win10_count_cys as in repo/scripts/v2_features.py::_site_features
    (0-based index = position - 1; other = every other C of the protein; win10 = residues index-10..index+10
    clipped to the protein, the site excluded)."""
    idx = position - 1
    n = len(seq)
    cpos = np.asarray([i for i, a in enumerate(seq) if a == "C"])
    other = cpos[cpos != idx]
    dist_min = float(np.abs(other - idx).min()) if other.size else float(n)
    in4 = float(np.sum(np.abs(other - idx) <= 4))
    low, high = max(0, idx - 10), min(n, idx + 11)
    win10 = float(sum(1 for j in range(low, high) if j != idx and seq[j] == "C"))
    return dist_min, in4, win10


def part_c(inp, res):
    ref = read_fasta_gz(inp.use(RICE_FASTA))
    fetched = {}
    for r in inp.read_csv(RICE_FETCHED, encoding="utf-8", delimiter="\t"):
        if r["source"] == "uniprot_rest_active":
            fetched.setdefault(r["requested_accession"], r["sequence"])
    rows = [r for r in inp.read_csv(J5_SITES) if r["dataset"] == "rice_PXD072089"]
    acc, y, lab, cols = [], [], [], {"B:cys_dist_min (negated)": [], "B:cys_in_4": [], "B:win10_count_cys": []}
    for r in rows:
        a, p = r["accession"], int(r["position"])
        s = ref.get(a) if r["full_sequence_source"] == "UP000059680_2026_03" else fetched.get(a)
        assert s is not None and s[p - 1] == "C", (a, p)
        d, i4, w10 = site_cys_columns(s, p)
        cols["B:cys_dist_min (negated)"].append(-d)
        cols["B:cys_in_4"].append(i4)
        cols["B:win10_count_cys"].append(w10)
        acc.append(a)
        y.append(int(float(r["full_second_cys_trypsinP"])))
        lab.append(int(r["label"]))
    acc, y, lab = np.asarray(acc), np.asarray(y), np.asarray(lab)
    rng = np.random.default_rng(SEED)
    W_, _, _ = cluster_weights(acc, BOOT, rng)
    out = {"n": int(y.size), "not_coverable": int(y.sum()), "proteins": int(np.unique(acc).size),
           "definition": "AUC of each column for the non-coverable flag (j5, Trypsin/P, full sequence) over all released rice "
                         "cysteines; protein-clustered bootstrap, 5,000 replicates, seed 20260930 + 3", "columns": {}}
    for name, v in cols.items():
        v = np.asarray(v)
        auc_w = make_weighted_auc(y, v)
        point = float(roc_auc_score(y, v))
        assert abs(point - auc_w(np.ones_like(v))) < 1e-12
        draws = np.asarray([auc_w(W_[b]) for b in range(BOOT)])
        draws = draws[np.isfinite(draws)]
        unl = lab == 0
        out["columns"][name] = {"auc": point, "ci95": [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))],
                                "auc_among_unlabelled": float(roc_auc_score(y[unl], v[unl])),
                                "median_not_coverable": float(np.median(v[y == 1])), "median_coverable": float(np.median(v[y == 0]))}
    out["reading"] = ("columns present in the chemistry set and in both deletion arms identify the cysteines that the "
                      "single-cysteine construction can never label; whether the refit models used them is not known")
    res["C_rice_coverability_readable"] = out


def part_d(inp, res):
    """Row sets behind the rice numbers of Supplemental Note 11 (read-only inspection of repository code and outputs)."""
    ana_p = REPO + "/scripts/analyse_public_refit_2026-09-21.py"
    cost_p = REPO + "/scripts/posthoc_public_refit_performance_cost_2026-09-21.py"
    cost_csv = REPO + "/results/public_refit_performance_cost_posthoc_2026-09-21.csv"
    ana = inp.read_text(ana_p).replace("\r\n", "\n").split("\n")
    cost = inp.read_text(cost_p).replace("\r\n", "\n")
    keep = [(i + 1, t.strip()) for i, t in enumerate(ana) if t.strip().startswith("keep = ") and "seqs" in t]
    wp = [(i + 1, t.strip()) for i, t in enumerate(ana) if ("protein_rows(" in t and "[m]" in t) or t.strip().endswith("m = d[\"y\"], d[\"score\"], keep")]
    per_protein = cost[cost.index("def per_protein"):cost.index("def main")]
    blend_rows = [r for r in inp.read_csv(cost_csv) if r["cohort"] == "rice_PXD072089"]
    note11 = inp.read_text(SUPP + "/Supplemental_Note_11_public_feature_ablation.md").split("\n")
    arm_rows = [t for t in note11 if t.startswith("| rice_PXD072089 |") and t.count("|") == 11]
    pair_rows = [t for t in note11 if t.startswith("| rice_PXD072089 |") and t.count("|") == 8]
    res["D_note11_rice_row_sets"] = {
        "analyse_public_refit_keep_mask": keep,
        "analyse_public_refit_within_protein_call": wp,
        "paired_difference_script_reads_oof_without_sequence_filter": ("_oof.csv" in per_protein and "seqs" not in per_protein
                                                                       and "keep" not in per_protein),
        "paired_difference_blend_csv_rice_n_proteins": sorted({int(r["n_proteins"]) for r in blend_rows}),
        "note11_arm_tables_rice_n_proteins_scored": sorted({int(t.split("|")[-2]) for t in arm_rows}),
        "note11_paired_tables_rice_n_proteins": sorted({int(t.split("|")[3]) for t in pair_rows}),
        "reading": ("the rice within-protein AUCs of the arm tables (and Table 2) use the rows with a sequence (5,004 rows, 578 "
                    "proteins with both classes); the paired differences report 637 proteins, the count for all 5,598 refit "
                    "rows; for the blend this follows from the script, which reads the out-of-fold files without the mask; the "
                    "ranking member's paired table reports the same 637, but its script is not on this machine"),
    }


def main():
    inp = Inputs()
    res = {"label": "POST HOC revision analysis (2026-09-30), item J_provenance_species, revision round 3; not registered, "
                    "not pre-specified"}
    outputs = part_a(inp, res)
    part_b(inp, res)
    part_c(inp, res)
    part_d(inp, res)
    # cited, not computed on: item F's corrected rice filter chain (1,769 -> 1,753 -> 1,377 -> 1,377 -> 1,340 -> 817)
    res["F_rice_chain_cited"] = [{"step": r["step"], "n_sites": int(r["n_sites"]), "n_proteins": int(r["n_proteins"])}
                                 for r in inp.read_csv(F_CHAIN_R2)]
    outputs.append(write_json(OUT + "/round3_checks.json", res))
    write_fragment(os.path.basename(__file__), inp, outputs, {"part C protein-clustered bootstrap": SEED},
                   {"A": "manuscript and supplement pointers located by unique anchor phrases; quotes checked verbatim",
                    "B": "human window power from j1 (trypsinP_complete) and j5 (window_second_cys_trypsinP) site tables",
                    "C": "three chemistry-arm columns recomputed with the v2_features.py definitions on UniProt 2026_03 sequences",
                    "D": "row sets behind the rice numbers of Supplemental Note 11, from the repository's scripts and outputs"})
    a = res["A_pointers"]
    print("manuscript:", res["A_manuscript"]["sha256"][:12], res["A_manuscript"]["lines"], "lines,",
          res["A_manuscript"]["last_modified_local"])
    print("pointers:", a["n"], "all quotes verbatim:", a["all_quotes_verbatim"], "offsets:", a["offsets"])
    print("table 2 block:", res["A_table2_block"])
    print("supplements:", res["A_supplements"])
    print("sibling pointers on shared lines:", res["A4_sibling_reports"])
    print("B:", {k: v for k, v in res["B_human_window_power"].items() if k != "reading"})
    print("C:", res["C_rice_coverability_readable"]["columns"])


if __name__ == "__main__":
    main()
