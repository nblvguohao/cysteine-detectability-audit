"""House rule 10 back-check for PHASE_3_REPORT.md: find every number token in the stored products.

Same token rules as scripts/verify_phase2_report_numbers_2026-09-22.py: whole-token match (accession/file-
name digits excluded by the left-neighbour guard), thousands commas stripped, and numeric-equality
normalisation (trailing fractional zeros stripped) so '0.9349' and a stored '0.9349' match regardless of
how each was formatted, and a genuinely different number does not.

CORRECTION C1 (found and fixed before the first real run produced a wrong verdict -- caught by manually
checking a handful of coefficient values that were visibly present in the corpus by eye but flagged
NO_PRODUCT). The regex's left-neighbour guard excluded a literal ASCII '-', intending only to stop a
sub-token match inside a compound like a date or accession. But every negative number in this tree's JSON
products is written as '-1.1734' style (ASCII minus, no space) -- so the guard ALSO permanently blocked the
corpus side from ever tokenising ANY negative number, for every report this whole session has back-checked,
not just this one. The report side was unaffected here only because its markdown tables use a typographic
minus (U+2212), a different character the guard doesn't match -- an accident of formatting, not a real
difference. Fixed by letting the token optionally consume a leading sign (ASCII '-' or U+2212) captured
right after the lookbehind (so the position immediately after it is no longer 'preceded by -' for the next
match attempt), and stripping the sign during normalisation (this convention has always compared magnitudes,
never signs, in every verifier this session wrote). Generalises to: A WHOLE-TOKEN NUMBER BACK-CHECK MUST
TOKENISE SIGNED NUMBERS ON BOTH SIDES, OR IT SILENTLY CANNOT VERIFY ANY NEGATIVE VALUE AT ALL.

Outputs (new names): results/phase3_report_number_check_2026-09-22.csv / _audit.json.
Interpreter: standard library only.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = os.path.join(ROOT, "PHASE_3_REPORT.md")
OUT = os.path.join(ROOT, "results", "phase3_report_number_check_2026-09-22.csv")
AUD = os.path.join(ROOT, "results", "phase3_report_number_check_2026-09-22_audit.json")
CORPUS = ["results/phase3_empirical_detectability_summary_2026-09-22.json",
          "results/phase3_empirical_detectability_2026-09-22_audit.json",
          "results/phase3_solver_crosscheck_2026-09-22.csv",
          "results/phase3_solver_crosscheck_2026-09-22_audit.json",
          "protocols/phase3_empirical_detectability_preregistration_2026-09-22.json"]
TOK = re.compile(r"(?<![\w.\-/_])[-−]?\d{1,3}(?:,\d{3})+(?:\.\d+)?(?![\d])"
                 r"|(?<![\w.\-/_])[-−]?\d+(?:\.\d+)?(?![\d])")
STRUCTURAL = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "2026", "09", "22", "0.05", "0.03", "0.5",
              "4600", "30", "9.17", "28"}  # section/list numbers, dates, thresholds, and a task-brief phase
              # name (Phase 5 = '28-claim audit') already spelled out in prose, not this tree's own data


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def norm(n):
    n = n.lstrip("-−")  # compare magnitudes, not signs (this session's convention throughout)
    if "." in n:
        n = n.rstrip("0")
        if n.endswith("."):
            n = n[:-1]
    return n


def derived(root):
    with open(os.path.join(root, "results/phase3_empirical_detectability_summary_2026-09-22.json"), encoding="utf-8") as fh:
        s = json.load(fh)
    coef = s["M3_full_model_coefficients"]
    ratio = coef["log_abundance"] / coef["length"]  # largest vs second-largest by |coefficient|
    return {"log_abundance_over_length_ratio": round(ratio, 1)}


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
    D = {k: norm(str(v)) for k, v in derived(ROOT).items()}
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
    cnt = Counter(r["state"] for r in rows)
    bad = [r for r in rows if r["state"] == "NO_PRODUCT"]
    audit = {"script": "scripts/verify_phase3_report_numbers_2026-09-22.py", "script_sha256": sha(os.path.abspath(__file__)),
             "report_sha256": sha(REPORT), "corpus": {rel: sha(os.path.join(ROOT, rel)) for rel in CORPUS},
             "counts": dict(cnt), "problems": bad, "outputs": {os.path.relpath(OUT, ROOT): sha(OUT)}}
    with open(AUD, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=1)
    print(json.dumps({"counts": dict(cnt), "problems": bad}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
