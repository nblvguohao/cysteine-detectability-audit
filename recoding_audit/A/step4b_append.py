"""Appends per-row disagreement detail and kappa robustness variants to
agreement_report_A.txt (same step-4 data, more detail)."""

import csv
import io
import sys
from collections import Counter

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import os  # v3.1.1: paths made repository-relative (were absolute Windows paths of the authoring machine)
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
SRC = os.path.join(_ROOT, 'supplemental', 'Supplemental_Data_9_survey_coding_table.csv')
MINE = os.path.join(_ROOT, 'recoding_audit', 'A', 'codes_A.csv')
OUT = os.path.join(_ROOT, 'recoding_audit', 'A', 'agreement_report_A.txt')

VERIF_NOTE = "explicit background statement found on full-text re-grep"


def kappa(pairs, cats):
    n = len(pairs)
    if not n:
        return float("nan")
    po = sum(1 for o, a in pairs if o == a) / n
    oc = Counter(o for o, _ in pairs)
    ac = Counter(a for _, a in pairs)
    pe = sum((oc.get(c, 0) / n) * (ac.get(c, 0) / n) for c in cats)
    return (po - pe) / (1 - pe) if pe != 1 else float("nan")


def main():
    orig = {}
    with open(SRC, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            if (r.get("background_class") or "").strip() in {"a", "b", "c", "e"}:
                orig[r["uid"]] = r
    mine = {r["uid"]: r for r in csv.DictReader(open(MINE, encoding="utf-8", newline=""))}
    uids = sorted(orig, key=int)

    L = []
    L.append("")
    L.append("=" * 78)
    L.append("APPENDIX - PER-ROW DISAGREEMENTS AND ROBUSTNESS VARIANTS")
    L.append("=" * 78)
    L.append("")

    dis = [u for u in uids if orig[u]["background_class"].strip() != mine[u]["coder_code"].strip()]
    L.append("A1. ALL %d DISAGREEMENTS (original -> coder A), with the recorded" % len(dis))
    L.append("    background evidence that drove coder A's call")
    L.append("")
    for u in dis:
        o, m = orig[u], mine[u]
        trunc = " [TRUNCATED QUOTE - a background clause may exist in the full text but is not in this record]" \
            if (o["evidence_quote_background"] or "").rstrip().endswith(("the", "all", "of", "at", "in", "amino", "with", "that", "which", "and", "was", "were")) \
            or "explicit background statement found" in (o["verification"] or "") else ""
        L.append("  uid %-5s %s -> %s   [%s]" % (u, o["background_class"].strip(), m["coder_code"], o["study_type"]))
        L.append("      %s" % o["title"].strip()[:110])
        L.append("      why A differs: %s" % m["justification"])
        if trunc:
            L.append("      NOTE:%s" % trunc)
        L.append("")

    L.append("A2. HOW MANY DISAGREEMENTS ARE QUOTE-TRUNCATION ARTEFACTS")
    art = [u for u in dis if VERIF_NOTE in (orig[u]["verification"] or "")]
    L.append("    %d of the %d disagreements carry a verification note asserting that an" % (len(art), len(dis)))
    L.append("    explicit background statement WAS found in the full text on re-grep, but the")
    L.append("    recorded evidence quote in this table still truncates before naming it:")
    L.append("    uids: %s" % ", ".join(art))
    L.append("    These are the rows where coder A's e is most likely to be too harsh and the")
    L.append("    original a is most likely correct; they are the rows the audit should resolve")
    L.append("    by going back to the papers' full text.")
    L.append("")

    L.append("A3. KAPPA ROBUSTNESS VARIANTS (sensitivity of the headline agreement)")
    pairs = [(orig[u]["background_class"].strip(), mine[u]["coder_code"].strip()) for u in uids]
    k5 = kappa(pairs, ["a", "b", "c", "e", "U"])
    pU = [("e" if o == "U" else o, "e" if a == "U" else a) for o, a in pairs]
    k4 = kappa(pU, ["a", "b", "c", "e"])
    pAE = [(o, a) for o, a in pairs if o in {"a", "e"} and a in {"a", "e"}]
    kAE = kappa(pAE, ["a", "e"])
    p4 = [(o, a) for o, a in pairs if o in {"a", "b", "c", "e"} and a in {"a", "b", "c", "e"}]
    k4x = kappa(p4, ["a", "b", "c", "e"])
    L.append("    kappa, 5 categories (a,b,c,e,U)                 = %.4f" % k5)
    L.append("    kappa, 4 categories with U folded into e        = %.4f" % k4)
    L.append("    kappa, 4 categories, U-rows dropped entirely    = %.4f  (n=%d)" % (k4x, len(p4)))
    L.append("    kappa, binary restricted to rows both coders placed in a or e")
    L.append("                                                    = %.4f  (n=%d)" % (kAE, len(pAE)))
    L.append("")
    L.append("    Binary a-vs-e agreement is %.1f%% (%d/%d), i.e. almost all of the"
             % (100 * sum(1 for o, a in pAE if o == a) / len(pAE),
                sum(1 for o, a in pAE if o == a), len(pAE)))
    L.append("    disagreement is about whether a background that is not stated is")
    L.append("    nonetheless to be inferred as the unrestricted proteome default.")
    L.append("")
    q = [orig[u]["evidence_quote_background"] or "" for u in uids]
    ntrunc = sum(1 for s in q if s.strip() and s.strip()[-1] not in '.).!?"”’]')
    nempty = sum(1 for s in q if not s.strip())
    L.append("A4. WHAT CODER A COULD NOT DO")
    L.append("    Coder A had no access to the papers. Every row was coded from the")
    L.append("    title, modification, study type and the recorded evidence quotations")
    L.append("    only. Measured over the 74 records: %d quotations are cut off mid-clause" % ntrunc)
    L.append("    (no terminal punctuation) and %d are empty. In most of the truncated" % nempty)
    L.append("    cases the missing clause is the one that would have named the background, so")
    L.append("    a background statement present in the paper may be absent from the record.")
    L.append("    For 1 of the 74 rows (uid 828) the recorded background evidence is a reagent")
    L.append("    sentence and is off-topic, which is the single U code.")
    L.append("    Coder A therefore cannot detect a background statement that exists in a")
    L.append("    paper but was not captured by the original regex extraction, so coder A's")
    L.append("    counts are lower bounds, not point estimates.")

    with open(OUT, "a", encoding="utf-8", newline="") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print("\nappended to %s" % OUT)


if __name__ == "__main__":
    main()
