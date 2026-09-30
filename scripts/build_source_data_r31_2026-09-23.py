#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Source Data for R31: correct the stale figure numbers and back the Figure 6 pointer.

WHY (round-5 review, results/manuscript_review_round5_checks_2026-09-22_audit.json)
  C04c  every shipped Source_Data_FigN file carries the figure number it had BEFORE the R27 renumbering, so the
        `figure` column of Source_Data_Fig2_* reads Fig7, Fig3_* reads Fig2, and so on. A reader who opens the
        file to check a value is told it belongs to a different figure.
  C04d  the Results promise that "the per-claim value of baseline_reproduction travels with every row of Figure 6
        (Source Data, Fig. 6)", and that field is not in the Figure 6 Source Data at all. It exists per claim in
        Supplementary Data 4, 6 and 7, so it is copied in here rather than the sentence being deleted.

WHAT THIS DOES
  * copies Source_Data_Fig3..Fig7 from source_data_r27/ into source_data_r31/, changing ONLY the `figure` column
    to the file's own figure number;
  * appends one `baseline_reproduction` row per claim to the Figure 6 file, with the value read from the round's
    Supplementary Data table (phase2b -> Data 6, phase2d -> Data 4, phase2e -> Data 7);
  * copies the two text-cited files byte for byte;
  * writes SOURCE_DATA_INDEX_R31.csv.
Figure 2's Source Data is written by scripts/plot_fig2_r31_2026-09-23.py and is only indexed here.

GATES (fixed before the run; any failure REFUSES and writes nothing)
  G1 for every copied table: same row count, same header, and every cell identical to the source except `figure`
  G2 the `figure` column of each output equals the figure number in its own filename
  G3 appended rows are only `baseline_reproduction`, one per distinct claim of that file, each value present in
     the Supplementary Data row for that claim; no existing row is modified
  G4 the index lists exactly the files present in source_data_r31/
  G5 the two text-cited files are byte-identical to their sources
  G6 positive control: the Figure 6 file must already contain claim SFE-001, and its baseline_reproduction value
     in Supplementary Data 7 must be non-empty (a lookup gate has to find something)
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
from collections import OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC27 = os.path.join(ROOT, "source_data_r27")
OUT = os.path.join(ROOT, "source_data_r31")
EN = os.path.join(ROOT, "submission_package_en")
SUPP = os.path.join(ROOT, "submission_package_r30_2026-09-22", "03_supplementary")
AUD = os.path.join(ROOT, "results", "source_data_r31_2026-09-23_audit.json")
ROUND_TABLE = {"phase2b_claim_retest_combined.csv": "Supplementary_Data_6_retest_round_b.csv",
               "phase2d_claim_retest.csv": "Supplementary_Data_4_retest_round_d.csv",
               "phase2e_claim_retest.csv": "Supplementary_Data_7_retest_round_e.csv"}
TEXT_CITED = ["Source_Data_restricted_cohort_size_reconciliation.csv",
              "Source_Data_text_phase2c_sfe006_reproduction.json"]


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def refuse(m):
    sys.stderr.write("REFUSE: " + m + "\n")
    raise SystemExit(2)


def rows_of(p):
    with open(p, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        return list(rd), rd.fieldnames


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=OUT)
    out_dir = ap.parse_args().out_dir
    if os.path.exists(AUD):
        refuse(AUD + " exists")
    os.makedirs(out_dir, exist_ok=True)
    audit = OrderedDict()
    audit["script"] = os.path.relpath(os.path.abspath(__file__), ROOT)
    audit["script_sha256"] = sha(os.path.abspath(__file__))
    audit["date"] = "2026-09-23"
    files, gates = OrderedDict(), OrderedDict()
    g1, g2, g3 = [], [], []

    supp = {}
    for fn in set(ROUND_TABLE.values()):
        supp[fn], _ = rows_of(os.path.join(SUPP, fn))
    pc_rows = [r for r in supp["Supplementary_Data_7_retest_round_e.csv"] if r["claim_id"] == "SFE-001"]
    if not pc_rows or not pc_rows[0].get("baseline_reproduction"):
        refuse("G6 positive control: SFE-001 has no baseline_reproduction in Supplementary Data 7")

    for name in sorted(os.listdir(SRC27)):
        if not re.match(r"Source_Data_Fig\d+_", name):
            continue
        n = int(re.search(r"Fig(\d+)", name).group(1))
        src_rows, header = rows_of(os.path.join(SRC27, name))
        new_rows = []
        for r in src_rows:
            row = OrderedDict((k, r[k]) for k in header)
            if "figure" in row:
                row["figure"] = "Fig%d" % n
            new_rows.append(row)
        for a, b in zip(src_rows, new_rows):
            if any(a[k] != b[k] for k in header if k != "figure"):
                g1.append(name)
                break
        if any(r.get("figure") not in (None, "Fig%d" % n) for r in new_rows):
            g2.append(name)
        added = 0
        if n == 6:
            seen, extra = set(), []
            for r in src_rows:
                key = (r["panel"], r["row"], r["source_table"])
                if key in seen:
                    continue
                seen.add(key)
                tbl = ROUND_TABLE.get(r["source_table"])
                if not tbl:
                    g3.append("unmapped source_table " + r["source_table"])
                    continue
                hit = [x for x in supp[tbl] if x["claim_id"] == r["row"]]
                if not hit or not hit[0].get("baseline_reproduction"):
                    g3.append("no baseline_reproduction for " + r["row"])
                    continue
                extra.append(OrderedDict([("figure", "Fig6"), ("panel", r["panel"]), ("row", r["row"]),
                                          ("field", "baseline_reproduction"),
                                          ("value", hit[0]["baseline_reproduction"]), ("source_table", tbl)]))
            new_rows += extra
            added = len(extra)
        dst = os.path.join(out_dir, name)
        with open(dst, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=header)
            w.writeheader()
            w.writerows(new_rows)
        files[name] = {"source": os.path.relpath(os.path.join(SRC27, name), ROOT), "source_sha256": sha(os.path.join(SRC27, name)),
                       "rows_copied": len(src_rows), "rows_appended": added, "sha256": sha(dst)}

    for name in TEXT_CITED:
        s, d = os.path.join(EN, name), os.path.join(out_dir, name)
        shutil.copy2(s, d)
        if sha(s) != sha(d):
            refuse("G5 copy mismatch " + name)
        files[name] = {"source": os.path.relpath(s, ROOT), "source_sha256": sha(s), "rows_copied": None,
                       "rows_appended": 0, "sha256": sha(d)}

    fig2 = "Source_Data_Fig2_public_four_protease.csv"
    if os.path.exists(os.path.join(out_dir, fig2)):
        files[fig2] = {"source": "scripts/plot_fig2_r31_2026-09-23.py", "source_sha256": None,
                       "rows_copied": None, "rows_appended": 0, "sha256": sha(os.path.join(out_dir, fig2))}

    index = os.path.join(out_dir, "SOURCE_DATA_INDEX_R31.csv")
    with open(index, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["source_data_file", "role", "figure", "source_in_tree"])
        w.writeheader()
        for name in sorted(files):
            m = re.search(r"Fig(\d+)", name)
            w.writerow({"source_data_file": name,
                        "role": "figure_underlying_values" if m else "text-cited table",
                        "figure": "Figure %s" % m.group(1) if m else "",
                        "source_in_tree": files[name]["source"]})
    present = sorted(f for f in os.listdir(out_dir) if f != os.path.basename(index))
    g4 = sorted(set(present) ^ set(files))
    gates["G1_cells_changed_outside_figure_column"] = g1
    gates["G2_figure_column_wrong"] = g2
    gates["G3_appended_row_problems"] = g3
    gates["G4_index_mismatch"] = g4
    gates["G5_text_cited_byte_identical"] = True
    gates["G6_positive_control"] = True
    audit["files"] = files
    audit["index"] = {os.path.relpath(index, ROOT): sha(index)}
    audit["gates"] = gates
    audit["all_pass"] = not (g1 or g2 or g3 or g4)
    if not audit["all_pass"]:
        refuse(json.dumps(gates)[:400])
    with open(AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1)
    print(json.dumps({"files": len(files), "gates": gates}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
