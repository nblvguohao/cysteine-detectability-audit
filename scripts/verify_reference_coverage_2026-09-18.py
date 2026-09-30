#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check that the R2 manuscript's reference list and its in-text citations agree. 2026-09-18.

Gate N02 of the 2026-09-18 submission gates reports "method_citations_still_needed: 15" against the
R2 manuscript. That number is read from results/nc_reference_candidates_2026-09-18.csv, a diagnostic
table written BEFORE any reference existed, which lists what needed doing; it is not a measurement of
the manuscript. This script measures the manuscript instead, so the gate's stale count cannot be
mistaken for a live defect -- and so that, if the reference list really is incomplete, that is caught
rather than excused.

Output  results/reference_coverage_r2_2026-09-18.csv
        results/verify_reference_coverage_2026-09-18_audit.json

=========================== Criteria fixed before running (house rule 1) ===========================

Four checks, all fail-closed, all computed from the manuscript and the verified reference list:

  C1  Every reference in the list is cited at least once in the body. An uncited reference is a
      defect in Nature-family style and, more practically, a sign the numbering drifted.
  C2  Every in-text citation number resolves to a reference that exists. A citation to [38] when the
      list ends at 37 is the failure this check exists for.
  C3  The citation numbers, read in order of first appearance through the document, are
      non-decreasing and start at 1 with no gaps. Nature numbers references in order of first
      citation; a gap means a reference is listed but never cited, and an out-of-order first
      appearance means the list was not renumbered after an edit.
  C4  Each of the fifteen method citations that results/nc_reference_candidates_2026-09-18.csv
      listed as "identifier deliberately blank, still needed" is either present in the reference
      list, or explicitly recorded here as still absent with the reason. This is the check that
      answers the gate's stale 15.

C3 is deliberately the strict one. The block citation "refs 4-26" makes references 4 through 26 all
first-appear at the same point, which is correct Nature practice for a corpus citation; the check
therefore treats a range as first-appearing at its own position, not as 23 separate out-of-order
appearances.
"""

import csv
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timedelta, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(ROOT, "reports", "MANUSCRIPT_SUBMISSION_R2_2026-09-18.md")
REFS = os.path.join(ROOT, "results", "nc_references_verified_2026-09-18.json")
CANDIDATES = os.path.join(ROOT, "results", "nc_reference_candidates_2026-09-18.csv")
OUT_CSV = os.path.join(ROOT, "results", "reference_coverage_r2_2026-09-18.csv")
OUT_AUDIT = os.path.join(ROOT, "results", "verify_reference_coverage_2026-09-18_audit.json")

# The fifteen method-citation needs from the candidates table, mapped to a distinctive string that
# must appear in the reference list entry if that need has been met. Fixed before the run.
METHOD_NEEDS = {
    "relative solvent accessibility reference values": "Maximum allowed solvent accessibilities",
    "PROPKA3 thiol pKa prediction": "empirical calculation and rationalization of pKa",
    "DSSP hydrogen-bond secondary structure assignment": "Dictionary of protein secondary structure",
    "Shrake-Rupley solvent accessibility algorithm": "Environment and exposure to solvent",
    "AlphaFold DB monomer models": "AlphaFold Protein Structure Database",
    "AlphaFold prediction method": "Highly accurate protein structure prediction",
    "scikit-learn HistGradientBoostingClassifier": "Scikit-learn",
    "LightGBM": "LightGBM",
    "Europe PMC REST API": "Europe PMC",
    "OpenAlex API": "OpenAlex",
    "Benjamini-Hochberg procedure": "Controlling the false discovery rate",
    "Wilson binomial confidence interval": "Probable inference, the law of succession",
    "Comet search engine": "Comet: an open-source",
    "target-decoy FDR estimation": "Comet: an open-source",
    "monoisotopic masses of sulfide and cysteine dioxidation": "Unimod",
}


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    text = open(MS, encoding="utf-8").read()
    refs = json.load(open(REFS, encoding="utf-8"))
    n_refs = len(refs)

    body, _, _ = text.partition("## References")
    if not _:
        print("REFUSED: no References section found.", file=sys.stderr)
        return 1

    # citation markers in the body, expanded: [4-26] or [4–26] is a range, [28,29] a list
    cited_order, cited = [], set()
    for m in re.finditer(r"\[(\d+(?:\s*[,–-]\s*\d+)*)\]|refs (\d+)\s*[–-]\s*(\d+)", body):
        if m.group(1):
            inner = m.group(1)
            nums = set()
            for part in re.split(r"\s*,\s*", inner):
                rng = re.match(r"(\d+)\s*[–-]\s*(\d+)$", part)
                if rng:
                    nums |= set(range(int(rng.group(1)), int(rng.group(2)) + 1))
                else:
                    nums.add(int(part))
        else:
            nums = set(range(int(m.group(2)), int(m.group(3)) + 1))
        for n in sorted(nums):
            if n not in cited:
                cited_order.append(n)
                cited.add(n)

    failures, rows = [], []

    uncited = [n for n in range(1, n_refs + 1) if n not in cited]
    if uncited:
        failures.append("C1: {} references are never cited in the body: {}".format(len(uncited), uncited))
    dangling = sorted(n for n in cited if n < 1 or n > n_refs)
    if dangling:
        failures.append("C2: {} citations resolve to no reference: {}".format(len(dangling), dangling))
    if cited_order != sorted(cited_order):
        firstbad = next(i for i in range(1, len(cited_order)) if cited_order[i] < cited_order[i - 1])
        failures.append("C3: citation numbers are not in order of first appearance; {} appears "
                        "after {}".format(cited_order[firstbad], cited_order[firstbad - 1]))
    if cited_order and cited_order[0] != 1:
        failures.append("C3: the first citation in the body is [{}], not [1]".format(cited_order[0]))

    reflist = "\n".join(refs)
    for need, probe in METHOD_NEEDS.items():
        present = probe in reflist
        num = ""
        for i, r in enumerate(refs, 1):
            if probe in r:
                num = str(i); break
        rows.append({"method_citation_need": need, "probe_string": probe,
                     "present_in_reference_list": "yes" if present else "no",
                     "reference_number": num,
                     "cited_in_body": "yes" if (num and int(num) in cited) else "no"})
        if not present:
            failures.append("C4: method citation still absent from the reference list: {}".format(need))
        elif num and int(num) not in cited:
            failures.append("C4: method citation {} is in the list as [{}] but never cited".format(need, num))

    if failures:
        print("REFUSED: {} checks failed. Nothing written.\n".format(len(failures)), file=sys.stderr)
        print("\n".join("  " + f for f in failures), file=sys.stderr)
        return 1

    tmp = OUT_CSV + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["method_citation_need", "probe_string",
                                          "present_in_reference_list", "reference_number",
                                          "cited_in_body"])
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, OUT_CSV)

    tz = timezone(timedelta(hours=8))
    audit = {
        "generated_at": datetime.now(tz).isoformat(),
        "script": "scripts/verify_reference_coverage_2026-09-18.py",
        "script_sha256": sha256_of(os.path.abspath(__file__)),
        "manuscript": {"path": "reports/MANUSCRIPT_SUBMISSION_R2_2026-09-18.md",
                       "sha256": sha256_of(MS)},
        "reference_list": {"path": "results/nc_references_verified_2026-09-18.json",
                           "sha256": sha256_of(REFS), "entries": n_refs},
        "checks": {
            "C1_every_reference_cited": "pass",
            "C2_every_citation_resolves": "pass",
            "C3_numbering_in_order_of_first_appearance": "pass",
            "C4_method_citations_present_and_cited": "pass ({} of {} method needs met)".format(
                sum(1 for r in rows if r["present_in_reference_list"] == "yes"), len(rows)),
        },
        "citations_in_body": {"distinct_references_cited": len(cited),
                              "first_appearance_order": cited_order},
        "supersedes": {
            "file": "results/nc_reference_candidates_2026-09-18.csv",
            "why": "that table lists what needed doing before any reference existed; gate N02 reads "
                   "its row count and therefore still reports 'method_citations_still_needed: 15' "
                   "against a manuscript in which all fifteen are now present and cited. This audit "
                   "is the live measurement; the candidates table is a historical diagnostic.",
        },
        "limitations": [
            "This checks that a citation exists and resolves, not that it is the right paper for the "
            "sentence it is attached to; that judgement is the corresponding author's.",
            "Reference content was verified when the list was built "
            "(results/nc_reference_crossref_pubmed_crosscheck_2026-09-18.json); this script does not "
            "re-fetch anything.",
        ],
    }
    tmp = OUT_AUDIT + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(audit, f, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
    os.replace(tmp, OUT_AUDIT)

    print("references: {}   distinct cited: {}   method needs met: {}/{}".format(
        n_refs, len(cited), sum(1 for r in rows if r["present_in_reference_list"] == "yes"), len(rows)))
    print("C1-C4 all pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
