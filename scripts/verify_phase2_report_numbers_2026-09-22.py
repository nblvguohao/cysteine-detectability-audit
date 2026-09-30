"""House rule 10 back-check for PHASE_2_REPORT.md: find every number token in the stored products.

Same token/normalisation rules as scripts/verify_phase1_report_numbers_2026-09-22.py (whole-token match,
thousands commas stripped, accession/DOI-style digits excluded by the left-neighbour guard).
DERIVED entries are computed here from the same stored CSVs/JSON the report cites, with the formula given
alongside each one, so a mismatch is a real discrepancy, not a rounding choice made silently in prose.
CORRECTION C1 (run 2; run-1 kept as scripts/verify_phase2_report_numbers_2026-09-22_run1_superseded.py, its
outputs renamed *_run1_superseded.*). Run 1 compared token strings literally, so '1.00' in the report against
a stored '1.0', or '0.9960' against a stored '0.996', were flagged NO_PRODUCT even though they are the same
number. Run 1 was executed BEFORE a wording pass on the report and so also caught two real issues that the
wording pass had already fixed by the time run 1's tail was read (a stale-check risk noted elsewhere in this
tree). Run 2 normalises by stripping trailing fractional zeros (keeping at least one digit) before comparing,
which is a numeric-equality fix, not a loosening of the check -- a token whose STRIPPED form is not in the
corpus is still NO_PRODUCT.

Outputs (new names): results/phase2_report_number_check_2026-09-22.csv / _audit.json.
Interpreter: project venv (numpy, for the two derived ratio/z-approximation numbers).
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = os.path.join(ROOT, "PHASE_2_REPORT.md")
OUT = os.path.join(ROOT, "results", "phase2_report_number_check_2026-09-22.csv")
AUD = os.path.join(ROOT, "results", "phase2_report_number_check_2026-09-22_audit.json")
CORPUS = ["results/phase2_synthetic_benchmark_summary_2026-09-22.csv",
          "results/phase2_synthetic_benchmark_2026-09-22_audit.json",
          "results/phase2_calibration_precision_2026-09-22.csv",
          "results/phase2_calibration_precision_2026-09-22_audit.json",
          "protocols/phase2_synthetic_benchmark_preregistration_2026-09-22.json"]
TOK = re.compile(r"(?<![\w.\-/_])\d{1,3}(?:,\d{3})+(?:\.\d+)?(?![\d])|(?<![\w.\-/_])\d+(?:\.\d+)?(?![\d])")
STRUCTURAL = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "2026", "09", "22", "10", "4",
              "0.1", "0.05"}  # section numbers, dates, and the two thresholds already spelled out in prose


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def norm(n):
    """numeric-equality normalisation: strip thousands commas (done by the caller) and trailing
    fractional zeros, keeping at least one digit after the point if a point remains."""
    if "." in n:
        n = n.rstrip("0")
        if n.endswith("."):
            n = n[:-1]
    return n


def derived():
    with open(os.path.join(ROOT, "results/phase2_synthetic_benchmark_replicates_2026-09-22.csv"), newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for bench, nullv in [("B1_cleavage", 0.0), ("B2_abundance", 0.5), ("B3_multicys", 0.0), ("B4_depletion", 0.5)]:
        sub = [r for r in rows if r["benchmark"] == bench and float(r["strength"]) == 0.0]
        stats = np.array([float(r["statistic"]) for r in sub])
        hw = (np.array([float(r["ci_hi"]) for r in sub]) - np.array([float(r["ci_lo"]) for r in sub])) / 2
        sd = stats.std(ddof=1)
        ratio = hw.mean() / (1.96 * sd)
        out[f"ratio_{bench}"] = round(ratio, 3)
        z_eff = 1.96 * ratio
        fpr_pred = 2 * (1 - 0.5 * (1 + math.erf(z_eff / math.sqrt(2))))
        out[f"fpr_pred_{bench}"] = round(fpr_pred * 100, 1)  # as a percentage, matching the prose (%)
    return out


def main():
    for p in (OUT, AUD):
        if os.path.exists(p):
            sys.exit(f"REFUSE: {p} exists")
    corpus_text = []
    for rel in CORPUS:
        with open(os.path.join(ROOT, rel), encoding="utf-8", errors="replace") as fh:
            corpus_text.append(fh.read())
    big = "\n".join(corpus_text)
    corp_tokens = {norm(m.group(0).replace(",", "")) for m in TOK.finditer(big)}
    D = {k: norm(str(v)) for k, v in derived().items()}
    with open(REPORT, encoding="utf-8") as fh:
        rep = fh.read()
    rows = []
    for m in TOK.finditer(rep):
        t = m.group(0)
        n = norm(t.replace(",", ""))
        line = rep.count("\n", 0, m.start()) + 1
        if n in STRUCTURAL:
            state, note = "STRUCTURAL", ""
        elif n in corp_tokens:
            sig = len(n.replace(".", "").lstrip("0"))
            state = "ON_DISK_WEAK" if sig <= 2 else "ON_DISK"
            note = ""
        elif n in D.values():
            state, note = "DERIVED", next(k for k, v in D.items() if v == n)
        else:
            state, note = "NO_PRODUCT", ""
        rows.append({"line": line, "token": t, "normalised": n, "state": state, "note": note})
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    cnt = Counter(r["state"] for r in rows)
    bad = [r for r in rows if r["state"] == "NO_PRODUCT"]
    audit = {"script": "scripts/verify_phase2_report_numbers_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
             "report_sha256": sha(REPORT), "corpus": {rel: sha(os.path.join(ROOT, rel)) for rel in CORPUS},
             "derived_values": D, "counts": dict(cnt), "problems": bad,
             "outputs": {os.path.relpath(OUT, ROOT): sha(OUT)}}
    with open(AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1)
    print(json.dumps({"counts": dict(cnt), "problems": bad, "derived_values": D}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
