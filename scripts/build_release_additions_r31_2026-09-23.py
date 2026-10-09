#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R31 release additions: the material three referees asked for, staged from files this tree already holds.

WHY (round-5 review, results/manuscript_review_actions_round5_2026-09-23.csv)
  R5-02 (blocking for two referees) the per-paper survey coding table is in neither the package nor the
        released repository, although Code availability says it is. It is in this tree. It ships here as
        Supplementary Data 10, with its evidence quotations, so the survey's proportions can be checked.
  R5-01 (blocking for two referees) every Artefact 1 and Artefact 2 number is produced by scripts in a separate
        working tree. The scripts and their stored outputs are staged here with source hashes so that they can
        be deposited with the code. The PlantPTMViewer export they read is a third-party database download and
        is NOT copied; it is identified instead, by category code and by the PubMed identifiers its own `exps`
        column carries, so a reader can fetch the same rows.
  R5-12 / R5-22 two Source Data gaps: the per-spectrum precursor-error spans behind Figure 3b, and the stored
        protease-reach table behind the values quoted in Artefact 4.

WHAT THIS DOES (copies and extracts only; it computes no new statistic)
  1. supplementary_r31/Supplementary_Data_10_survey_coding_table.csv        byte copy of the 942-row table
  2. supplementary_r31/Supplementary_Data_10_survey_coding_table_README.md  column glossary and funnel counts,
     printed from the stored recount audit
  3. release_staging_artefact12_2026-09-23/{scripts,results}/              byte copies with a manifest
  4. source_data_r31/Source_Data_Fig3b_tie_precursor_error_spans.csv       one row per tied spectrum
  5. source_data_r31/Source_Data_text_protease_reach_public.csv            byte copy of the stored table
  6. source_data_r31/SOURCE_DATA_INDEX_R31.csv                             rewritten to list what is there

GATES (fixed before the run; any failure REFUSES and writes nothing new)
  G1 every byte copy has sha256(source) == sha256(copy)
  G2 the coding table has 942 rows and exactly 120 rows flagged as classified, of which the background classes
     are a 34, b 1, c 4, e 35 and 46 blank, which is what the manuscript's survey section rests on
  G3 the Figure 3b extract has exactly 2079 rows, one per tied spectrum, and its span column reproduces the
     maximum 18.829 that the manuscript prints
  G4 the staged Artefact 1 and 2 outputs contain the z values the manuscript prints (15.72, 10.03, 5.59, 8.22,
     4.71, 0.67, 6.04, 6.62) and the Artefact 2 values (-2.61, 3.34, -8.78, 2.67); the value 2.54, which has no
     stored table, must NOT be found, and that absence is recorded rather than repaired
  G5 the index lists exactly the files present in source_data_r31/
Interpreter: any Python 3.9+, standard library only. Deterministic; two runs are byte-identical.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
from collections import Counter, OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREE2 = "/path/to/project"
SUPP31 = os.path.join(ROOT, "supplementary_r31")
SD31 = os.path.join(ROOT, "source_data_r31")
STAGE = os.path.join(ROOT, "release_staging_artefact12_2026-09-23")
AUD = os.path.join(ROOT, "results", "release_additions_r31_2026-09-23_audit.json")
CODING = os.path.join(ROOT, "external", "lit_survey_2026-09-13", "lit_survey_background_controls.csv")
RECOUNT_AUD = os.path.join(ROOT, "results", "lit_survey_recount_audit.json")
TIES = os.path.join(ROOT, "results", "pxd015307_posthoc_score_ties.csv")
REACH = os.path.join(ROOT, "results", "protease_reach_public_2026-09-17.csv")
STAGE_SCRIPTS = ["09_map_large_scale.py", "12_detectability_matched.py", "21_crossspecies_test.py",
                 "36_ncys_decompose.py", "10_large_scale_motif.py"]
STAGE_RESULTS = ["12_matched_sno.json", "12_matched_so.json", "12_matched_ox.json",
                 "21_crossspecies_summary.json", "21_crossspecies_summary.csv", "35_symmetric_ncys_test.json"]
PTMVIEWER = {"sno": "S-nitrosylation", "so": "S-sulfenylation", "ox": "reversible cysteine oxidation"}


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def refuse(m):
    sys.stderr.write("REFUSE: " + m + "\n")
    raise SystemExit(2)


def copy(src, dst, record):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)
    if sha(src) != sha(dst):
        refuse("G1 copy mismatch " + src)
    record[os.path.relpath(dst, ROOT)] = {"source": src if not src.startswith(ROOT) else os.path.relpath(src, ROOT),
                                          "sha256": sha(dst)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    if os.path.exists(AUD):
        refuse(AUD + " exists")
    copied = OrderedDict()
    gates = OrderedDict()

    # ---- G2: the coding table
    with open(CODING, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    classified = [r for r in rows if r["in_classified_sample"] == "1"]
    bg = Counter(r["background_class"] for r in classified)
    ok2 = (len(rows) == 942 and len(classified) == 120 and bg.get("a") == 34 and bg.get("b") == 1
           and bg.get("c") == 4 and bg.get("e") == 35)
    gates["G2_coding_table"] = {"rows": len(rows), "classified": len(classified), "background_classes": dict(bg),
                               "pass": ok2}
    if not ok2:
        refuse("G2: coding table does not reconcile: %s" % dict(bg))

    # ---- G3: the Figure 3b extract
    with open(TIES, encoding="utf-8") as fh:
        ties = [r for r in csv.DictReader(fh)
                if "sulfide" in r["chemistries_at_best_xcorr"] and "dioxidation" in r["chemistries_at_best_xcorr"]]
    spans = [float(r["ppm_span"]) for r in ties]
    ok3 = len(ties) == 2079 and abs(max(spans) - 18.829) < 1e-9
    gates["G3_fig3b_extract"] = {"rows": len(ties), "max_span_ppm": max(spans),
                                 "zero_span_rows": sum(1 for s in spans if s == 0.0), "pass": ok3}
    if not ok3:
        refuse("G3: tie extract %d rows, max span %s" % (len(ties), max(spans)))

    # ---- G4: the staged Artefact 1 and 2 outputs
    staged_text = ""
    for name in STAGE_RESULTS:
        p = os.path.join(TREE2, "results", name)
        if not os.path.exists(p):
            refuse("G4: missing staged result " + p)
        with open(p, encoding="utf-8") as fh:
            staged_text += fh.read()
    want = ["15.72", "10.03", "5.59", "8.22", "4.71", "0.67", "6.04", "6.62", "-2.61", "3.34", "-8.78", "2.67"]
    missing = [w for w in want if w not in staged_text]
    has254 = re.search(r"(?<![\d.])2\.54(?![\d])", staged_text) is not None
    gates["G4_staged_outputs_carry_the_printed_values"] = {
        "values_not_found": missing, "value_2_54_found_in_any_stored_table": has254,
        "pass": not missing and not has254}
    if missing:
        refuse("G4: staged outputs lack %s" % missing)
    if has254:
        refuse("G4: 2.54 unexpectedly found; the review recorded it as having no stored table")

    if args.dry_run:
        print(json.dumps(gates, indent=1))
        return

    # ---- 1 and 2: Supplementary Data 10
    os.makedirs(SUPP31, exist_ok=True)
    copy(CODING, os.path.join(SUPP31, "Supplementary_Data_10_survey_coding_table.csv"), copied)
    ra = json.load(open(RECOUNT_AUD, encoding="utf-8"))
    funnel = ra["funnel_recomputed_from_the_copied_table"]
    mods = ra["modification_counts_over_the_classified_corpus"]
    glossary = [
        "# Supplementary Data 10. The survey coding table", "",
        "One row per record retained at title and abstract screening (%d rows). The columns are as coded; "
        "`in_classified_sample` marks the %d records drawn for full-text coding, and `background_class` is "
        "filled only for the %d that proved to be site-level feature analyses."
        % (len(rows), len(classified), sum(v for k, v in bg.items() if k)), "",
        "**Funnel, recomputed from this table**", "",
        "| step | records |", "|---|---:|"] + \
        ["| %s | %s |" % (k, v) for k, v in funnel.items()] + ["",
        "**Background classes over the classified corpus**", "",
        "| class | meaning | papers |", "|---|---|---:|",
        "| a | all residues of that type, unrestricted | %d |" % bg.get("a", 0),
        "| b | matched on theoretical detectability or abundance | %d |" % bg.get("b", 0),
        "| c | residues detected in the same experiment without a modification assignment | %d |" % bg.get("c", 0),
        "| e | other, or no statement found | %d |" % bg.get("e", 0), "",
        "There is no member of class d (an equal-size shuffled set) in the classified corpus. "
        "`equal_size_random_background` is a separate, non-exclusive flag and is set for %d papers."
        % sum(1 for r in classified if r["equal_size_random_background"] == "1"), "",
        "**Modification composition of the classified corpus**", "", "| modification | papers |", "|---|---:|"] + \
        ["| %s | %s |" % (k, v) for k, v in sorted(mods.items(), key=lambda kv: (-kv[1], kv[0]))] + ["",
        "**Coding limits.** One reader coded the sample; there is no second coder and no agreement estimate. "
        "Title and abstract screening was model-assisted and its false-negative rate was not measured. "
        "Evidence quotations are recorded per coded field in the `evidence_quote_*` columns, and `ambiguity` "
        "records where the coder judged the paper's statement unclear.", ""]
    gp = os.path.join(SUPP31, "Supplementary_Data_10_survey_coding_table_README.md")
    with open(gp, "w", encoding="utf-8") as fh:
        fh.write("\n".join(glossary))
    copied[os.path.relpath(gp, ROOT)] = {"source": "generated from " + os.path.relpath(RECOUNT_AUD, ROOT),
                                         "sha256": sha(gp)}

    # ---- 3: staged Artefact 1 and 2 material
    for name in STAGE_SCRIPTS:
        copy(os.path.join(TREE2, "scripts", name), os.path.join(STAGE, "scripts", name), copied)
    for name in STAGE_RESULTS:
        copy(os.path.join(TREE2, "results", name), os.path.join(STAGE, "results", name), copied)
    ptm_sources = {}
    for code, label in PTMVIEWER.items():
        p = os.path.join(TREE2, "data", "ptmviewer", "%s_allSPECIES.csv" % code)
        with open(p, encoding="utf-8") as fh:
            rr = [x for x in csv.DictReader(fh) if x["species"] == "ath" and x["modified_aa"] == "C"]
        pmids = Counter()
        for x in rr:
            for pm in (x["exps"] or "").replace(";", ",").split(","):
                if pm.strip():
                    pmids[pm.strip()] += 1
        ptm_sources[code] = {"modification": label, "arabidopsis_cysteine_rows": len(rr),
                             "top_source_pubmed_ids": [p for p, _ in pmids.most_common(4)],
                             "file_sha256_not_redistributed": sha(p)}
    manifest = OrderedDict([
        ("what_this_is", "The scripts and stored outputs behind Artefact 1 and Artefact 2, staged for deposit "
                         "with the released code. They were produced in a separate working tree; the copies here "
                         "carry the source hash so that the deposit can be checked against them."),
        ("source_tree", TREE2),
        ("files", {k: v for k, v in copied.items() if k.startswith(os.path.relpath(STAGE, ROOT))}),
        ("input_database", {"name": "PlantPTMViewer", "note": "third-party database export, identified rather "
                            "than redistributed; each row's `exps` column carries the PubMed identifiers of the "
                            "studies the site came from", "categories": ptm_sources}),
        ("value_without_a_stored_table", {"value": "2.54", "where_it_is_printed": "Artefact 2, the unfiltered "
                                          "neighbouring-cysteine z with P 0.011 on 980 sites",
                                          "status": "no stored table produces it; R31 removes it from the text"}),
    ])
    mp = os.path.join(STAGE, "MANIFEST.json")
    with open(mp, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=1)
    copied[os.path.relpath(mp, ROOT)] = {"source": "generated", "sha256": sha(mp)}

    # ---- 4 and 5: the two Source Data additions
    os.makedirs(SD31, exist_ok=True)
    fp = os.path.join(SD31, "Source_Data_Fig3b_tie_precursor_error_spans.csv")
    with open(fp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["figure", "panel", "file", "scan", "best_xcorr", "expectation_value",
                                           "ppm_span", "source_table"])
        w.writeheader()
        for r in ties:
            w.writerow({"figure": "Fig3", "panel": "b", "file": r["arm"], "scan": r["scan"],
                        "best_xcorr": r["best_xcorr"], "expectation_value": r["e_value_best"],
                        "ppm_span": r["ppm_span"], "source_table": "pxd015307_posthoc_score_ties.csv"})
    copied[os.path.relpath(fp, ROOT)] = {"source": os.path.relpath(TIES, ROOT), "sha256": sha(fp)}
    copy(REACH, os.path.join(SD31, "Source_Data_text_protease_reach_public.csv"), copied)

    # ---- 6: rewrite the index
    index = os.path.join(SD31, "SOURCE_DATA_INDEX_R31.csv")
    present = sorted(f for f in os.listdir(SD31) if f != os.path.basename(index))
    with open(index, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["source_data_file", "role", "figure", "panel"])
        w.writeheader()
        for name in present:
            m = re.search(r"Fig(\d+)([a-d])?", name)
            w.writerow({"source_data_file": name,
                        "role": "figure_underlying_values" if m else "text-cited table",
                        "figure": "Figure %s" % m.group(1) if m else "",
                        "panel": m.group(2) if (m and m.group(2)) else ""})
    gates["G5_index_matches_directory"] = {"files": len(present), "pass": True}

    audit = OrderedDict()
    audit["script"] = os.path.relpath(os.path.abspath(__file__), ROOT)
    audit["script_sha256"] = sha(os.path.abspath(__file__))
    audit["date"] = "2026-09-23"
    audit["gates"] = gates
    audit["copied"] = copied
    audit["index"] = {os.path.relpath(index, ROOT): sha(index)}
    audit["what_this_does_not_do"] = [
        "It computes no statistic and changes no value.",
        "It does not push anything to the public repository; depositing the staged material is an author action.",
        "It does not redistribute the PlantPTMViewer export; that database is identified, not copied.",
    ]
    with open(AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1)
    print(json.dumps({"gates": {k: v.get("pass") for k, v in gates.items()}, "files_written": len(copied)}, indent=1))


if __name__ == "__main__":
    main()
