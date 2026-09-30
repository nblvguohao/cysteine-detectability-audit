"""Step 1 of the independent recoding protocol (coder A).

Reads the survey coding table, keeps the 74 rows classified after full-text
reading (background_class in {a,b,c,e}), and writes a BLINDED extract that
omits every column that could reveal the original code.

Nothing in this script prints or persists the original codes in a form the
coder can see while coding.
"""

import csv
import os
import sys

SRC = r"C:\Users\admin\Desktop\MCP\supplemental\Supplemental_Data_9_survey_coding_table.csv"
OUT = r"C:\Users\admin\Desktop\MCP\_recoding\A\blinded_74.csv"

KEEP_CLASSES = {"a", "b", "c", "e"}

# Columns permitted in the blinded file (exact task list, order preserved).
BLIND_COLS = [
    "uid",
    "doi",
    "pmid",
    "year",
    "journal",
    "title",
    "modification",
    "study_type",
    "fulltext_stage_verdict",
    "evidence_quote_background",
    "evidence_quote_detectability",
    "ambiguity",
    "note",
    "verification",
]

# Columns that must never appear (asserted at the end as a safety check).
FORBIDDEN = {
    "background_class",
    "background_class_label",
    "equal_size_random_background",
    "abundance_normalisation",
    "discusses_detectability_bias",
    "reports_identified_fraction",
    "kr_basic_as_biological_motif",
    "evidence_quote_identified_fraction",
    "evidence_quote_kr",
    "evidence_quote_equal_size_bg",
}


def main():
    with open(SRC, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        header = reader.fieldnames
        missing = [c for c in BLIND_COLS if c not in header]
        if missing:
            sys.exit("missing expected columns: %s" % missing)
        rows = [r for r in reader]

    print("total rows in source table: %d" % len(rows))

    classified = [r for r in rows if (r.get("background_class") or "").strip() in KEEP_CLASSES]
    print("rows with background_class in {a,b,c,e}: %d" % len(classified))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=BLIND_COLS, extrasaction="ignore")
        w.writeheader()
        for r in classified:
            w.writerow({c: r.get(c, "") for c in BLIND_COLS})

    # Read back and verify.
    with open(OUT, encoding="utf-8", newline="") as fh:
        back = list(csv.DictReader(fh))
    print("blinded_74.csv row count: %d" % len(back))
    leaked = FORBIDDEN.intersection(back[0].keys()) if back else set()
    if leaked:
        sys.exit("LEAK: forbidden columns present: %s" % sorted(leaked))
    print("blinded_74.csv contains no forbidden columns. OK")

    # uid uniqueness sanity check
    uids = [r["uid"] for r in back]
    print("unique uids: %d" % len(set(uids)))

    # Blank-evidence census (does NOT reveal codes) so we know how much is
    # genuinely unclassifiable from the record.
    blank_bg = sum(1 for r in back if not (r["evidence_quote_background"] or "").strip())
    blank_det = sum(1 for r in back if not (r["evidence_quote_detectability"] or "").strip())
    print("rows with EMPTY evidence_quote_background: %d" % blank_bg)
    print("rows with EMPTY evidence_quote_detectability: %d" % blank_det)


if __name__ == "__main__":
    main()
