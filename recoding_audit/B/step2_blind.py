import csv, os

import os  # v3.1.1: paths made repository-relative (were absolute Windows paths of the authoring machine)
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
SRC = os.path.join(_ROOT, 'supplemental', 'Supplemental_Data_9_survey_coding_table.csv')
OUTDIR = os.path.join(_ROOT, 'recoding_audit', 'B')
OUT = os.path.join(OUTDIR, "blinded_74.csv")

KEEP_COLS = ["uid", "doi", "pmid", "year", "journal", "title", "modification",
             "study_type", "fulltext_stage_verdict", "evidence_quote_background",
             "evidence_quote_detectability", "ambiguity", "note", "verification"]

# columns that must NOT appear in the blinded extract
FORBIDDEN = ["background_class", "background_class_label", "equal_size_random_background",
             "abundance_normalisation", "discusses_detectability_bias",
             "reports_identified_fraction", "kr_basic_as_biological_motif",
             "evidence_quote_identified_fraction", "evidence_quote_kr",
             "evidence_quote_equal_size_bg", "sources", "queries"]

assert not (set(KEEP_COLS) & set(FORBIDDEN)), "leak in KEEP_COLS"

with open(SRC, newline="", encoding="utf-8-sig") as f:
    r = csv.DictReader(f)
    rows = [row for row in r]

classified = [row for row in rows if (row.get("background_class") or "").strip() in ("a", "b", "c", "e")]

os.makedirs(OUTDIR, exist_ok=True)
with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=KEEP_COLS)
    w.writeheader()
    for row in classified:
        w.writerow({k: row.get(k, "") for k in KEEP_COLS})

print("TOTAL_ROWS_IN_SOURCE:", len(rows))
print("CLASSIFIED_ROWS_WRITTEN:", len(classified))
print("OUT:", OUT)
print("COLUMNS_WRITTEN:", KEEP_COLS)

# integrity check: re-read the blinded file and confirm no forbidden column present
with open(OUT, newline="", encoding="utf-8") as f:
    hdr = next(csv.reader(f))
leak = [c for c in hdr if c in FORBIDDEN]
print("FORBIDDEN_COLUMNS_PRESENT:", leak)
print("BLINDED_COLCOUNT:", len(hdr))
