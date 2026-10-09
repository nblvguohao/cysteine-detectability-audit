# -*- coding: utf-8 -*-
"""J1. Species of PRIDE PXD044043: collect the internal evidence and run two local checks.

POST HOC revision analysis (2026-09-30), written in response to a pre-submission review; not
registered, not pre-specified. No network access; no downloads.

QUESTION
PRIDE lists PXD044043 as Mus musculus; the manuscript calls the cohort human. What evidence do our
own records hold for the human assignment, and what is missing?

PART A (records only; nothing recomputed): extract, from files already in the internal repository,
the evidence items (PRIDE organism field, the deposit's own MaxQuant parameters and experiment
names, species tokens in the input-arm proteinGroups, the cohort's accessions against the human
reference proteome file, the publication's census row) and write them with the sha256 of the file
each came from.

PART B (new, local): two checks on the released cohort table external/public_cohort_site_scores.csv
(18,567 human-cohort cysteines with accession, position, label and a 15-residue context window):
  B1 accession check: how many of the 1,519 cohort accessions occur in the mouse reference proteome
     UP000000589 (UniProt 2026_03, W/external)?
  B2 peptide-level check: for every cohort cysteine, the tryptic segment that contains it is cut out
     of its context window (primary rule Trypsin/P, i.e. cleavage after every K/R, which gives the
     shortest peptide and is therefore conservative for a "human-specific" call; sensitivity: the
     proline-blocked trypsin rule). A segment of at least 7 residues that occurs in no mouse
     reference-proteome sequence (I and L treated as identical, since MS cannot tell them apart) is
     "human-specific": any peptide containing it, with or without missed cleavages, cannot be
     produced by mouse material. A complete segment (both cleavage sites inside the window or at a
     protein terminus) of 7 or more residues that does occur in mouse is "conserved".
     Reading: if the labelled cysteines came from mouse material searched against a human database,
     labelled cysteines would be (almost) absent from human-specific segments, because only peptides
     shared by the two species could be identified. If the material is human, the labelled share in
     human-specific segments should be of the same order as in conserved segments. The ratio of the
     two labelled shares is reported with a protein-clustered bootstrap interval (5,000 replicates,
     seed 20260930). This check speaks to the origin of the peptides behind the labels; it cannot
     show which cell line was used, and it cannot replace a check on the deposit's own peptide table,
     which is not on this machine.
  B3 conservation-depth control (added after B2 had been run): the same comparison against the rice
     reference proteome UP000059680 (W/external), where no rice material can be present, to see how
     much of any labelled-share difference between conserved and non-conserved segments arises from
     sequence conservation itself.

REVISION AFTER VERIFICATION (round 1, 2026-09-30; still post hoc)
  * Window termini. Round 1 treated the window as starting at the protein N-terminus only when the
    cysteine was at position <= 7; for position 8 the 15-residue window also starts at residue 1, so
    the rule is now (position - index of C in window) == 1. The C-terminal side had the symmetric
    problem (a cysteine exactly 7 residues from the C-terminus has a full-length window that ends at the
    C-terminus); it is now decided with the protein length recorded in
    results/v3_human_homology_components.csv (sequence_length, from the same human proteome file),
    after checking that every window truncated on the right ends exactly at that length. Every
    cysteine whose class changes is listed in species_terminus_fix_changes.csv.
  * Clustering sensitivity. The labelled-share intervals are also computed with the k-mer homology
    components of that table (0.40 threshold) as bootstrap clusters instead of proteins.

OUTPUTS (results/J_provenance_species/)
  species_evidence_items.csv, species_accession_check.csv, species_peptide_check_sites.csv,
  species_peptide_check_summary.csv, species_cohort_structure.csv, species_terminus_fix_changes.csv
"""
from __future__ import annotations

import gzip
import os
import re
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jcommon import (TREE2, MCP, OUT, BACKUP_DIR, REPO, W, Inputs, write_csv,  # noqa: E402
                     write_fragment)

SEED = 20260930
BOOT = 5000
MOUSE_FASTA = W + "/external/UP000000589_10090.fasta.gz"
RICE_FASTA = W + "/external/UP000059680_39947.fasta.gz"   # used only for the conservation-depth control B3
SCORES = REPO + "/external/public_cohort_site_scores.csv"
COMP = REPO + "/results/v3_human_homology_components.csv"   # protein lengths and k-mer components (human)
SCREEN_PROTOCOL = REPO + "/protocols/phase1_dataset_screening_preregistration_2026-09-22.json"
MIN_LEN = 7
# Post hoc reading rule for B2, written before this rule was applied but after the first run had
# shown the counts: the check "supports human material" if more than 5% of the labelled cysteines
# lie in human-specific segments. 5% is a generous allowance for false identifications at a 1%
# peptide FDR (MaxQuant default; the deposit's own setting is not in our records) plus mouse
# isoforms that are absent from the reference proteome.
HS_ALLOWANCE = 0.05

_BIG = None  # mouse proteome, I->L, '|'-separated; set in worker initialiser


def read_mouse(path):
    accs, seqs = [], []
    acc, buf = None, []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    seqs.append("".join(buf))
                m = re.match(r">(?:sp|tr)\|([^|]+)\|", line)
                acc = m.group(1) if m else line[1:].split()[0]
                accs.append(acc)
                buf = []
            else:
                buf.append(line.strip())
    if acc is not None:
        seqs.append("".join(buf))
    return accs, seqs


def _init(big):
    global _BIG
    _BIG = big


def _found(seg):
    return seg in _BIG


def tryptic_segment(w, c_idx, left_is_nterm, right_is_cterm, block_proline):
    """Tryptic segment containing w[c_idx] inside the context window w.

    Returns (start, end, left_complete, right_complete); w[start:end] is always a substring of
    the true 0-missed-cleavage peptide, whether or not the peptide extends beyond the window.
    """
    n = len(w)
    s, left_complete = 0, bool(left_is_nterm)
    for j in range(c_idx - 1, -1, -1):
        if w[j] in "KR" and (not block_proline or w[j + 1] != "P"):
            s, left_complete = j + 1, True
            break
    e, right_complete = n, bool(right_is_cterm)
    for j in range(c_idx + 1, n):
        if w[j] in "KR":
            if j + 1 < n:
                if not block_proline or w[j + 1] != "P":
                    e, right_complete = j + 1, True
                    break
            else:  # K/R is the last residue of the window
                if right_is_cterm or not block_proline:
                    e, right_complete = j + 1, True
                else:  # next residue unknown: peptide ends here or continues; w[s:n] is a substring either way
                    e, right_complete = n, False
                break
    return s, e, left_complete, right_complete


def classify(seg_len, complete, found):
    if seg_len >= MIN_LEN and not found:
        return "human_specific"
    if seg_len >= MIN_LEN and complete and found:
        return "conserved"
    if complete and seg_len < MIN_LEN:
        return "short_complete"
    return "undetermined"


def cluster_boot_rates(groups, cls, y, classes, boot=BOOT, seed=SEED, chunk=250):
    """Protein-clustered bootstrap of the labelled share within each class and of their ratio."""
    uniq, inv = np.unique(groups, return_inverse=True)
    rng = np.random.default_rng(seed)
    masks = {c: (cls == c) for c in classes}
    # per-protein sums, so that each replicate is a weighted sum over proteins
    per = {}
    for c in classes:
        m = masks[c]
        per[c] = (np.bincount(inv, weights=(m & (y == 1)).astype(float), minlength=len(uniq)),
                  np.bincount(inv, weights=m.astype(float), minlength=len(uniq)))
    draws = {c: [] for c in classes}
    ratio = []
    remaining = boot
    while remaining > 0:
        take = min(chunk, remaining)
        wts = rng.multinomial(len(uniq), np.full(len(uniq), 1.0 / len(uniq)), size=take).astype(float)
        vals = {}
        for c in classes:
            pos, tot = per[c]
            num, den = wts @ pos, wts @ tot
            with np.errstate(invalid="ignore", divide="ignore"):
                vals[c] = num / den
            draws[c].extend(vals[c].tolist())
        if len(classes) == 2:
            with np.errstate(invalid="ignore", divide="ignore"):
                ratio.extend((vals[classes[0]] / vals[classes[1]]).tolist())
        remaining -= take

    def q(a):
        a = np.asarray(a, dtype=float)
        a = a[np.isfinite(a)]
        return (float(np.quantile(a, 0.025)), float(np.quantile(a, 0.975)), int(a.size)) if a.size else (np.nan, np.nan, 0)

    return {c: q(draws[c]) for c in classes}, (q(ratio) if ratio else None)


def mh_within_protein(groups, cls, y, a="human_specific", b="conserved", boot=BOOT, seed=SEED, chunk=250):
    """Mantel-Haenszel odds ratio of being labelled, class a vs class b, stratified by protein
    (only proteins carrying both classes contribute), with a protein-clustered bootstrap interval."""
    keep = (cls == a) | (cls == b)
    g, c, yy = groups[keep], cls[keep], y[keep]
    uniq, inv = np.unique(g, return_inverse=True)
    x = (c == a).astype(float)
    n11 = np.bincount(inv, weights=x * yy, minlength=len(uniq))          # a, labelled
    n10 = np.bincount(inv, weights=x * (1 - yy), minlength=len(uniq))    # a, unlabelled
    n01 = np.bincount(inv, weights=(1 - x) * yy, minlength=len(uniq))    # b, labelled
    n00 = np.bincount(inv, weights=(1 - x) * (1 - yy), minlength=len(uniq))
    tot = n11 + n10 + n01 + n00
    both = ((n11 + n10) > 0) & ((n01 + n00) > 0)
    num_p = np.where(both, n11 * n00 / np.where(tot > 0, tot, 1), 0.0)
    den_p = np.where(both, n10 * n01 / np.where(tot > 0, tot, 1), 0.0)
    point = float(num_p.sum() / den_p.sum()) if den_p.sum() > 0 else float("nan")
    rng = np.random.default_rng(seed)
    draws, remaining = [], boot
    while remaining > 0:
        take = min(chunk, remaining)
        wts = rng.multinomial(len(uniq), np.full(len(uniq), 1.0 / len(uniq)), size=take).astype(float)
        num, den = wts @ num_p, wts @ den_p
        with np.errstate(invalid="ignore", divide="ignore"):
            draws.extend((num / den).tolist())
        remaining -= take
    d = np.asarray(draws)
    d = d[np.isfinite(d)]
    return {"mh_or": point, "ci_low": float(np.quantile(d, 0.025)), "ci_high": float(np.quantile(d, 0.975)),
            "n_proteins_with_both_classes": int(both.sum()),
            "n_cysteines_in_those_proteins": int(tot[both].sum())}


def main():
    inp = Inputs()
    outputs = []

    # ------------------------------------------------------------------ PART A: records
    items = []

    def item(eid, evidence, value, source, supports, note):
        items.append({"id": eid, "evidence": evidence, "value": value,
                      "source_file": source.replace(REPO, "repo").replace(BACKUP_DIR, "BACKUP").replace(MCP, "MCP").replace(W, "W"),
                      "source_sha256": inp.files[os.path.normpath(source).replace("\\", "/")]["sha256"] if os.path.exists(source) else "",
                      "supports": supports, "note": note})

    cand = inp.read_csv(REPO + "/results/phase1_dataset_screening_candidates_2026-09-22.csv")
    r = [x for x in cand if x["accession"] == "PXD044043"][0]
    prot = inp.read_json(SCREEN_PROTOCOL)
    item("A1", "PRIDE project metadata: organism field and file listing",
         r["organisms"], REPO + "/results/phase1_dataset_screening_candidates_2026-09-22.csv", "mouse",
         f"PRIDE record retrieved 2026-09-22; submissionType {r['submissionType']}; files {r['n_files']} "
         f"(raw-named {r['n_raw']}, other {r['n_nonraw']}); search archives {r['search_archives']}; "
         f"site-table-named files: {r['site_table_files'] or 'none'}; peptide-table-named files: "
         f"{r['peptide_table_files'] or 'none'}; lab PI {r['lab_pis']}. The third non-raw file's name and type are "
         "not recorded: it matched none of the screening protocol's name patterns (site table "
         f"/{prot['regex']['SITE_TABLE_NAME']}/, peptide table /{prot['regex']['PEPTIDE_TABLE_NAME']}/, search "
         f"archive /{prot['regex']['SEARCH_ARCHIVE_NAME']}/), so a results table not named after sites or peptides "
         "cannot be excluded; check the PRIDE listing before submission")

    aud = inp.read_json(REPO + "/results/phase1_download_request_2026-09-22_audit.json")
    od = aud["on_disk_PXD044043"]
    for arm in ("MQ-txt BTA", "MQ-txt input"):
        a = od["arms"][arm]
        fpath = a["Fasta file"]
        year_dir = [seg for seg in fpath.split("\\") if re.fullmatch(r"20\d\d", seg)]
        item("A2" + ("a" if arm.endswith("BTA") else "b"),
             f"Deposit's own MaxQuant parameters.txt, arm '{arm}': FASTA searched",
             fpath.split("\\")[-1], REPO + "/results/phase1_download_request_2026-09-22_audit.json", "human",
             f"full path as written by the depositors: {fpath}; the release year is inferred only from the folder "
             f"names in that path ({', '.join(year_dir) or 'none'}), i.e. a {year_dir[0] if year_dir else 'n/a'} copy per the path, "
             f"not from a release number; MaxQuant {a['Version']}; fixed modification {a['Fixed modifications']}; "
             f"variable modifications {a['experiments'][0]['variable_modifications']} (read 2026-09-22 from the deposit's archive)")
        exps = "; ".join(e["experiment"] + " [" + e["raw_file"] + "]" for e in a["experiments"])
        item("A3" + ("a" if arm.endswith("BTA") else "b"),
             f"Deposit's own MaxQuant summary, arm '{arm}': experiment names", exps,
             REPO + "/results/phase1_download_request_2026-09-22_audit.json", "human",
             "all six run names of the arm contain 'THP1'; THP-1 is a human monocytic cell line. Whether the cells "
             "were differentiated to macrophages is not stated in the files we read (the publication title says "
             "'activated macrophages'), so 'THP-1 cells (run names)' is the most the deposit files support")
    pg = od["input_proteinGroups"]
    item("A4", "Input-arm proteinGroups.txt: species tokens in protein identifiers",
         f"rows {pg['rows']}; tokens {pg['species_tokens']}", REPO + "/results/phase1_download_request_2026-09-22_audit.json",
         "human", "the enriched (BTA) arm's proteinGroups.txt was not token-counted; archive sha256 "
         f"txt_BTA.rar {od['txt_BTA.rar']['sha256']}, txt_input.rar {od['txt_input.rar']['sha256']}")
    adj = inp.read_csv(REPO + "/results/phase1_adjudicated_candidates_2026-09-22.csv")
    ra = [x for x in adj if x["accession"] == "PXD044043"][0]
    item("A5", "Internal adjudication note (2026-09-22)", ra["why"],
         REPO + "/results/phase1_adjudicated_candidates_2026-09-22.csv", "human", "a statement, based on A2-A4")
    ass = inp.read_csv(REPO + "/results/pride_coincidence_assessed_2026-09-19.csv")
    rc = [x for x in ass if x["accession"] == "PXD044043"][0]
    members = {k: sorted(od[k]["members"]) for k in ("txt_BTA.rar", "txt_input.rar")}
    site_like = {k: [m for m in v if re.search(r"sites?\.txt$", m, re.I)] for k, v in members.items()}
    item("A6", "Contents of the two deposited MaxQuant archives: cysteine site table", rc["status"],
         REPO + "/results/pride_coincidence_assessed_2026-09-19.csv", "neutral",
         "neither MaxQuant archive has a cysteine-modification site table; their only site tables are "
         f"{site_like['txt_BTA.rar']} (BTA) and {site_like['txt_input.rar']} (input), i.e. methionine oxidation "
         "(members of txt_BTA.rar: " + ", ".join(m.split("/")[-1] for m in members["txt_BTA.rar"]) + "). "
         "Carbamidomethyl (C) is fixed and no cysteine modification is variable in either arm, so the cohort's "
         "labels were derived from the identified peptides by a mapping step whose script is not in the repository. "
         "This holds for the MaxQuant outputs; the third non-raw file of the PRIDE listing is unnamed in our records (A1)")

    comp = inp.read_csv(REPO + "/results/v3_human_homology_components.csv")
    comp_aud = inp.read_json(REPO + "/results/v3_human_homology_components_audit.json")
    px = [x for x in comp if "pxd044043" in x["cohorts"]]
    item("A7", "Cohort accessions found in the human reference proteome file used by the repository",
         f"{len(px)} of {comp_aud['scope']['per_cohort']['pxd044043_persulfidome']}",
         REPO + "/results/v3_human_homology_components.csv", "human",
         "proteome file external/proteomes/hsa.fasta.gz sha256 "
         f"{comp_aud['input_hashes']['external/proteomes/hsa.fasta.gz']}; cohort file inputs/PXD044043_mapped_human.json "
         f"sha256 {comp_aud['input_hashes']['inputs/PXD044043_mapped_human.json']} (both absent from this clone); "
         f"gene symbols are human (e.g. {', '.join(x['gene'] for x in px[:5])})")
    cen = inp.read_csv(REPO + "/supplemental/Supplemental_Data_7_dataset_census.csv")
    rp = [x for x in cen if x["accession_or_doi"] == "10.1016/j.redox.2024.103125"]
    item("A8", "Census row for the source publication (Salti et al. 2024, Redox Biol 72:103125)",
         rp[0]["species"] if rp else "", REPO + "/supplemental/Supplemental_Data_7_dataset_census.csv", "neutral",
         "the publication used mouse and human macrophages; PRIDE's organism field may describe the "
         "publication rather than the deposited proteomics")
    lit = inp.read_csv(REPO + "/results/parallel_label_literature_screen.csv")
    rl = [x for x in lit if x["doi"] == "10.1016/j.redox.2024.103125"]
    item("A9", "Abstract-level note on the source publication", rl[0]["screen_note"] if rl else "",
         REPO + "/results/parallel_label_literature_screen.csv", "neutral",
         "our screening note, written from the abstract, says '~800 proteins'; the publication's own wording, "
         "page and criteria for calling a protein persulfidated are not recorded in our records and must be "
         "quoted from the paper itself (the brief's '>800' is not in our records)")
    scr = inp.read_text(REPO + "/scripts/build_dataset_screening_tables_2026-09-22.py")
    ln = [ln for ln in scr.splitlines() if ln.strip().startswith('"PXD044043": "A3 PRIMARY')]
    item("A10", "Internal note on how the cohort label relates to the enriched arm",
         ln[0].strip() if ln else "", REPO + "/scripts/build_dataset_screening_tables_2026-09-22.py", "neutral",
         "translation: 'the BTA arm's CAM is a fixed modification; positive = detected holds by construction'")
    lab = inp.read_text(W + "/inputs/repo_reports/LABEL_PROVENANCE_AND_SEMANTICS.md")
    ll = [ln for ln in lab.splitlines() if "PXD044043" in ln]
    item("A11", "Internal note on how positives were inferred", ll[0].strip() if ll else "",
         W + "/inputs/repo_reports/LABEL_PROVENANCE_AND_SEMANTICS.md", "neutral",
         "translation: 'human PXD044043 and rice PXD072089 are site-level persulfidation data (positives inferred "
         "from single-cysteine peptides)'")
    ph1 = inp.read_text(BACKUP_DIR + "/PHASE_1_REPORT.md")
    lp = [ln for ln in ph1.splitlines() if "PRIDE 把它标成小鼠" in ln]
    item("A12", "Phase 1 report (2026-09-22) conclusion", lp[0].strip() if lp else "",
         BACKUP_DIR + "/PHASE_1_REPORT.md", "human",
         "translation: 'species: input arm proteinGroups 4,103 rows, 7,770 _HUMAN tokens, 0 _MOUSE; search database "
         "human UniProt UP000005640. PRIDE labels it mouse, which is wrong.'")

    # ------------------------------------------------------------------ PART B: cohort table
    rows = [x for x in inp.read_csv(SCORES) if x["dataset"] == "human_PXD044043"]
    acc = np.asarray([x["accession"] for x in rows])
    pos = np.asarray([int(x["position"]) for x in rows])
    y = np.asarray([int(float(x["observed_in_retrospective_dataset"])) for x in rows])
    win = [x["context_15aa_or_terminal_shorter"] for x in rows]
    prot_pos = Counter(a for a, yy in zip(acc, y) if yy == 1)
    prot_all = Counter(acc)
    structure = [
        {"quantity": "cysteines", "value": len(rows)},
        {"quantity": "labelled cysteines", "value": int(y.sum())},
        {"quantity": "proteins", "value": len(prot_all)},
        {"quantity": "proteins with >=1 labelled cysteine", "value": len(prot_pos)},
        {"quantity": "proteins with both labelled and unlabelled cysteines",
         "value": sum(1 for a in prot_all if 0 < prot_pos.get(a, 0) < prot_all[a])},
        {"quantity": "proteins in which every cysteine is labelled",
         "value": sum(1 for a in prot_all if prot_pos.get(a, 0) == prot_all[a])},
    ]
    outputs.append(write_csv(OUT + "/species_cohort_structure.csv", structure))

    # B1 accessions
    maccs, mseqs = read_mouse(inp.use(MOUSE_FASTA))
    mset = set(maccs)
    uacc = sorted(set(acc))
    in_mouse = [a for a in uacc if a in mset]
    base = [a.split("-")[0] for a in uacc]
    in_mouse_base = sorted({b for b in base if b in mset})
    acc_rows = [
        {"check": "cohort accessions (unique)", "n": len(uacc), "note": "human_PXD044043 rows of the released site table"},
        {"check": "found in mouse reference proteome UP000000589 (2026_03)", "n": len(in_mouse),
         "note": f"{len(maccs)} mouse entries; examples: {', '.join(in_mouse[:5])}"},
        {"check": "found after stripping isoform suffix", "n": len(in_mouse_base), "note": ""},
        {"check": "found in the human reference proteome file used by the repository (record A7)", "n": len(px),
         "note": "external/proteomes/hsa.fasta.gz, read by build_v3_human_homology_components.py"},
    ]
    outputs.append(write_csv(OUT + "/species_accession_check.csv", acc_rows))

    # B2 peptide-level check
    # protein lengths (and k-mer components, used for the clustering sensitivity) from the repository's
    # human component table, built from the same human proteome file as the cohort (record A7)
    seqlen = {x["accession"]: int(x["sequence_length"]) for x in comp if "pxd044043" in x["cohorts"]}
    component = {x["accession"]: x["component"] for x in comp if "pxd044043" in x["cohorts"]}
    big = "|".join(s.replace("I", "L") for s in mseqs)
    site_rows, segs = [], {}
    bad_windows = 0
    term = Counter()        # bookkeeping of the terminus rule (revision after verification)
    for i in range(len(rows)):
        w = win[i]
        p = int(pos[i])
        c_idx = p - 1 if p <= 7 else 7
        if c_idx >= len(w) or w[c_idx] != "C":
            bad_windows += 1
            continue
        right_len = len(w) - c_idx - 1
        L = seqlen.get(acc[i])
        # corrected rule: the window starts at residue 1 whenever position - c_idx == 1 (p <= 8), and ends
        # at the C-terminus whenever its last residue is the protein's last residue
        left_is_nterm = (p - c_idx) == 1
        if L is not None:
            right_is_cterm = (p + right_len) == L
            if right_len < 7:
                term["right_truncated_windows"] += 1
                term["right_truncated_consistent_with_length" if right_is_cterm else "right_truncated_INCONSISTENT"] += 1
            if p > L:
                term["position_beyond_length"] += 1
        else:
            right_is_cterm = right_len < 7
            term["no_length_fallback"] += 1
        # round-1 rule, kept only to list the cysteines whose class the correction changes
        left_r1, right_r1 = p <= 7, right_len < 7
        term["left_nterm_newly_true"] += int(left_is_nterm and not left_r1)
        term["right_cterm_newly_true"] += int(right_is_cterm and not right_r1)
        rec = {"accession": acc[i], "position": p, "label": int(y[i]), "window": w,
               "window_starts_at_residue_1": int(left_is_nterm), "window_ends_at_c_terminus": int(right_is_cterm)}
        for rule, block in (("trypsinP", False), ("trypsin", True)):
            s, e, lc, rc_ = tryptic_segment(w, c_idx, left_is_nterm, right_is_cterm, block)
            seg = w[s:e]
            rec[f"{rule}_segment"] = seg
            rec[f"{rule}_len"] = len(seg)
            rec[f"{rule}_complete"] = int(lc and rc_)
            s1, e1, lc1, rc1 = tryptic_segment(w, c_idx, left_r1, right_r1, block)
            rec[f"_{rule}_seg_round1"] = w[s1:e1]
            rec[f"_{rule}_complete_round1"] = int(lc1 and rc1)
            segs[seg.replace("I", "L")] = None
            segs[w[s1:e1].replace("I", "L")] = None
        site_rows.append(rec)
    assert term["right_truncated_INCONSISTENT"] == 0 and term["position_beyond_length"] == 0, term
    keys = sorted(segs)
    with Pool(4, initializer=_init, initargs=(big,)) as pool:
        found = pool.map(_found, keys, chunksize=200)
    fmap = dict(zip(keys, found))
    changes = []
    for rec in site_rows:
        for rule in ("trypsinP", "trypsin"):
            f = fmap[rec[f"{rule}_segment"].replace("I", "L")]
            rec[f"{rule}_in_mouse"] = int(f)
            rec[f"{rule}_class"] = classify(rec[f"{rule}_len"], bool(rec[f"{rule}_complete"]), f)
            seg1 = rec.pop(f"_{rule}_seg_round1")
            cpl1 = rec.pop(f"_{rule}_complete_round1")
            cls1 = classify(len(seg1), bool(cpl1), fmap[seg1.replace("I", "L")])
            if cls1 != rec[f"{rule}_class"]:
                changes.append({"accession": rec["accession"], "position": rec["position"], "label": rec["label"],
                                "rule": rule, "window": rec["window"],
                                "window_starts_at_residue_1": rec["window_starts_at_residue_1"],
                                "window_ends_at_c_terminus": rec["window_ends_at_c_terminus"],
                                "class_round1": cls1, "class_corrected": rec[f"{rule}_class"],
                                "segment_round1": seg1, "segment_corrected": rec[f"{rule}_segment"]})
    outputs.append(write_csv(OUT + "/species_terminus_fix_changes.csv", changes,
                             ["accession", "position", "label", "rule", "window", "window_starts_at_residue_1",
                              "window_ends_at_c_terminus", "class_round1", "class_corrected", "segment_round1",
                              "segment_corrected"]))

    summ = []
    g = np.asarray([r["accession"] for r in site_rows])
    gcomp = np.asarray([component[a] for a in g])     # k-mer homology components (clustering sensitivity)
    yy = np.asarray([r["label"] for r in site_rows])
    for rule in ("trypsinP", "trypsin"):
        cls = np.asarray([r[f"{rule}_class"] for r in site_rows])
        comp_ = np.asarray([r[f"{rule}_complete"] for r in site_rows])
        ln_ = np.asarray([r[f"{rule}_len"] for r in site_rows])
        subsets = {
            "all windows": np.ones(len(site_rows), bool),
            "complete segments of 7-30 residues only": (comp_ == 1) & (ln_ >= 7) & (ln_ <= 30),
        }
        for sname, smask in subsets.items():
            c2 = np.where(smask, cls, "excluded")
            ci, rci = cluster_boot_rates(g, c2, yy, ["human_specific", "conserved"])
            cic, rcic = cluster_boot_rates(gcomp, c2, yy, ["human_specific", "conserved"])
            n_prot_hs_pos = len({a for a, c, l in zip(g, c2, yy) if c == "human_specific" and l == 1})
            for c in ("human_specific", "conserved", "undetermined", "short_complete"):
                m = c2 == c
                rate = float(yy[m].mean()) if m.any() else float("nan")
                row = {"rule": rule, "subset": sname, "class": c, "n_cysteines": int(m.sum()),
                       "n_labelled": int(yy[m].sum()), "labelled_share": rate,
                       "ci_low": "", "ci_high": "", "n_proteins": len(set(g[m]))}
                if c in ci:
                    row["ci_low"], row["ci_high"] = ci[c][0], ci[c][1]
                    row["ci_low_component_clustered"], row["ci_high_component_clustered"] = cic[c][0], cic[c][1]
                if c == "human_specific":
                    row["n_proteins_with_labelled_human_specific_cysteine"] = n_prot_hs_pos
                summ.append(row)
            hs = c2 == "human_specific"
            cv = c2 == "conserved"
            ratio = float(yy[hs].mean() / yy[cv].mean()) if hs.any() and cv.any() else float("nan")
            summ.append({"rule": rule, "subset": sname, "class": "ratio human_specific / conserved",
                         "n_cysteines": int(hs.sum() + cv.sum()), "n_labelled": int(yy[hs].sum() + yy[cv].sum()),
                         "labelled_share": ratio, "ci_low": rci[0], "ci_high": rci[1], "n_proteins": "",
                         "ci_low_component_clustered": rcic[0], "ci_high_component_clustered": rcic[1]})
            mh = mh_within_protein(g, c2, yy)
            summ.append({"rule": rule, "subset": sname,
                         "class": "Mantel-Haenszel OR within protein, human_specific vs conserved",
                         "n_cysteines": mh["n_cysteines_in_those_proteins"], "n_labelled": "",
                         "labelled_share": mh["mh_or"], "ci_low": mh["ci_low"], "ci_high": mh["ci_high"],
                         "n_proteins": mh["n_proteins_with_both_classes"]})
    # B3 conservation-depth control: the same comparison against the rice reference proteome, where no
    # rice material can be present. If labelled shares also differ between rice-conserved and
    # non-conserved segments, a deficit in human-specific (vs mouse) segments reflects properties of
    # conserved sequence (e.g. abundance, detectability), not non-human material.
    raccs, rseqs = read_mouse(inp.use(RICE_FASTA))
    rbig = "|".join(sq.replace("I", "L") for sq in rseqs)
    rkeys = sorted({r["trypsinP_segment"].replace("I", "L") for r in site_rows})
    with Pool(4, initializer=_init, initargs=(rbig,)) as pool:
        rfound = pool.map(_found, rkeys, chunksize=200)
    rmap = dict(zip(rkeys, rfound))
    for rec in site_rows:
        f = rmap[rec["trypsinP_segment"].replace("I", "L")]
        rec["trypsinP_in_rice"] = int(f)
        c = classify(rec["trypsinP_len"], bool(rec["trypsinP_complete"]), f)
        rec["trypsinP_class_vs_rice"] = {"human_specific": "not_in_rice", "conserved": "conserved_in_rice"}.get(c, c)
    cls_r = np.asarray([r["trypsinP_class_vs_rice"] for r in site_rows])
    comp_ = np.asarray([r["trypsinP_complete"] for r in site_rows])
    ln_ = np.asarray([r["trypsinP_len"] for r in site_rows])
    for sname, smask in (("all windows", np.ones(len(site_rows), bool)),
                         ("complete segments of 7-30 residues only", (comp_ == 1) & (ln_ >= 7) & (ln_ <= 30))):
        c2 = np.where(smask, cls_r, "excluded")
        ci, rci = cluster_boot_rates(g, c2, yy, ["not_in_rice", "conserved_in_rice"])
        cic, rcic = cluster_boot_rates(gcomp, c2, yy, ["not_in_rice", "conserved_in_rice"])
        for c in ("not_in_rice", "conserved_in_rice"):
            m = c2 == c
            summ.append({"rule": "trypsinP vs rice proteome (control)", "subset": sname, "class": c,
                         "n_cysteines": int(m.sum()), "n_labelled": int(yy[m].sum()),
                         "labelled_share": float(yy[m].mean()) if m.any() else float("nan"),
                         "ci_low": ci[c][0], "ci_high": ci[c][1], "n_proteins": len(set(g[m])),
                         "ci_low_component_clustered": cic[c][0], "ci_high_component_clustered": cic[c][1]})
        a, b = c2 == "not_in_rice", c2 == "conserved_in_rice"
        summ.append({"rule": "trypsinP vs rice proteome (control)", "subset": sname, "class": "ratio not_in_rice / conserved_in_rice",
                     "n_cysteines": int(a.sum() + b.sum()), "n_labelled": int(yy[a].sum() + yy[b].sum()),
                     "labelled_share": float(yy[a].mean() / yy[b].mean()) if a.any() and b.any() else float("nan"),
                     "ci_low": rci[0], "ci_high": rci[1], "n_proteins": "",
                     "ci_low_component_clustered": rcic[0], "ci_high_component_clustered": rcic[1]})
        mh = mh_within_protein(g, c2, yy, a="not_in_rice", b="conserved_in_rice")
        summ.append({"rule": "trypsinP vs rice proteome (control)", "subset": sname,
                     "class": "Mantel-Haenszel OR within protein, not_in_rice vs conserved_in_rice",
                     "n_cysteines": mh["n_cysteines_in_those_proteins"], "n_labelled": "",
                     "labelled_share": mh["mh_or"], "ci_low": mh["ci_low"], "ci_high": mh["ci_high"],
                     "n_proteins": mh["n_proteins_with_both_classes"]})
    outputs.append(write_csv(OUT + "/species_peptide_check_sites.csv", site_rows))
    outputs.append(write_csv(OUT + "/species_peptide_check_summary.csv", summ))

    items.append({"id": "B1", "evidence": "Cohort accessions found in the mouse reference proteome (new local check)",
                  "value": f"{len(in_mouse)} of {len(uacc)}", "source_file": "W/external/UP000000589_10090.fasta.gz + repo/external/public_cohort_site_scores.csv",
                  "source_sha256": inp.files[os.path.normpath(MOUSE_FASTA).replace("\\", "/")]["sha256"],
                  "supports": "human" if not in_mouse else "mixed", "note": "UniProt accessions are species-specific entries"})
    hs_row = [s for s in summ if s["rule"] == "trypsinP" and s["subset"] == "all windows" and s["class"] == "human_specific"][0]
    cv_row = [s for s in summ if s["rule"] == "trypsinP" and s["subset"] == "all windows" and s["class"] == "conserved"][0]
    rt_row = [s for s in summ if s["rule"] == "trypsinP" and s["subset"] == "all windows" and s["class"].startswith("ratio")][0]
    rc_row = [s for s in summ if s["rule"] == "trypsinP" and s["subset"].startswith("complete") and s["class"].startswith("ratio")][0]
    mh_row = [s for s in summ if s["rule"] == "trypsinP" and s["subset"] == "all windows" and s["class"].startswith("Mantel")][0]
    items.append({"id": "B2", "evidence": "Labelled share in human-specific vs conserved tryptic segments (new local check)",
                  "value": (f"human-specific {hs_row['n_labelled']}/{hs_row['n_cysteines']} = {hs_row['labelled_share']:.4f} "
                            f"({hs_row['n_proteins_with_labelled_human_specific_cysteine']} proteins); "
                            f"conserved {cv_row['n_labelled']}/{cv_row['n_cysteines']} = {cv_row['labelled_share']:.4f}; "
                            f"ratio {rt_row['labelled_share']:.3f} [{rt_row['ci_low']:.3f}, {rt_row['ci_high']:.3f}] "
                            f"(homology-component clusters: [{rt_row['ci_low_component_clustered']:.3f}, "
                            f"{rt_row['ci_high_component_clustered']:.3f}]); "
                            f"complete 7-30-residue segments only: ratio {rc_row['labelled_share']:.3f} "
                            f"[{rc_row['ci_low']:.3f}, {rc_row['ci_high']:.3f}]; within-protein Mantel-Haenszel OR "
                            f"{mh_row['labelled_share']:.3f} [{mh_row['ci_low']:.3f}, {mh_row['ci_high']:.3f}]"),
                  "source_file": "repo/external/public_cohort_site_scores.csv + W/external/UP000000589_10090.fasta.gz",
                  "source_sha256": inp.files[os.path.normpath(SCORES).replace("\\", "/")]["sha256"],
                  "supports": "human" if hs_row["n_labelled"] > HS_ALLOWANCE * int(yy.sum()) else "not decisive",
                  "note": ("mouse material searched against a human database can yield a peptide absent from the mouse "
                           "one-per-gene reference proteome only through a false identification or through a mouse "
                           "sequence outside that proteome (isoforms, variants, unreviewed entries; not checked here: the "
                           "mouse additional-isoform file is not on this machine and no download was allowed). "
                           f"{hs_row['n_labelled']} of {int(yy.sum())} labelled cysteines "
                           f"({hs_row['n_labelled'] / int(yy.sum()):.4f}) lie in such segments, against a post hoc allowance of "
                           f"{HS_ALLOWANCE:.2f} (set after the first run) for false identifications (MaxQuant default 1% FDR; the "
                           "deposit's own FDR setting is not in our records) and for mouse sequences outside the reference "
                           "proteome; how large the second share is has not been measured. The ratio below one is not a "
                           "species-mixture estimate: conserved and human-specific segments differ in length, completeness "
                           f"and detectability. Windows whose centre is not C: {bad_windows}")})
    rr = {(s["subset"], s["class"]): s for s in summ if s["rule"].startswith("trypsinP vs rice")}
    ra = rr[("all windows", "ratio not_in_rice / conserved_in_rice")]
    rcpl = rr[("complete segments of 7-30 residues only", "ratio not_in_rice / conserved_in_rice")]
    rmh = rr[("all windows", "Mantel-Haenszel OR within protein, not_in_rice vs conserved_in_rice")]
    items.append({"id": "B3", "evidence": "Conservation-depth control: the same comparison against the rice reference proteome (no rice material can be present)",
                  "value": (f"ratio not-in-rice / conserved-in-rice {ra['labelled_share']:.3f} [{ra['ci_low']:.3f}, {ra['ci_high']:.3f}]; "
                            f"complete 7-30-residue segments {rcpl['labelled_share']:.3f} [{rcpl['ci_low']:.3f}, {rcpl['ci_high']:.3f}]; "
                            f"within-protein Mantel-Haenszel OR {rmh['labelled_share']:.3f} [{rmh['ci_low']:.3f}, {rmh['ci_high']:.3f}]"),
                  "source_file": "repo/external/public_cohort_site_scores.csv + W/external/UP000059680_39947.fasta.gz",
                  "source_sha256": inp.files[os.path.normpath(RICE_FASTA).replace("\\", "/")]["sha256"],
                  "supports": "interpretation of B2",
                  "note": ("a labelled-share deficit in non-conserved segments that also appears against a plant proteome, "
                           "where no plant material can contribute, shows that part of the B2 deficit reflects properties of "
                           "conserved sequence rather than non-human material; it cannot exclude a minor non-human share")})
    outputs.append(write_csv(OUT + "/species_evidence_items.csv", items))

    n_comp = len(set(component[a] for a in uacc))
    write_fragment(os.path.basename(__file__), inp, outputs, {"bootstrap_seed": SEED, "bootstrap_replicates": BOOT},
                   {"second_tree_path_recorded": TREE2, "bad_windows": bad_windows,
                    "mouse_entries": len(maccs), "rice_entries": len(raccs), "unique_segments_searched": len(keys),
                    "terminus_rule_bookkeeping (revision after verification)": dict(term),
                    "class_changes_vs_round1_rule": dict(Counter(f"{c['rule']}: {c['class_round1']} -> {c['class_corrected']}" for c in changes)),
                    "component_clusters_for_cohort": n_comp})
    print("done", len(site_rows), "sites;", len(keys), "segments;", "bad windows", bad_windows)
    print("terminus bookkeeping", dict(term))
    print("class changes", Counter(f"{c['rule']}: {c['class_round1']} -> {c['class_corrected']}" for c in changes))
    for s in summ:
        print(s)


if __name__ == "__main__":
    main()
