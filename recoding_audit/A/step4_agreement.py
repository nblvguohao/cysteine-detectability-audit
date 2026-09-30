"""Step 4-5 of the independent recoding protocol (coder A).

Joins coder A's blinded codes back to the original coding table on uid and
reports raw agreement, Cohen's kappa, a confusion matrix, the disposition of
the rows the original coder called b or c, and the strict/lenient counts.
"""

import csv
import io
import sys
from collections import Counter, defaultdict

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

SRC = r"C:\Users\admin\Desktop\MCP\supplemental\Supplemental_Data_9_survey_coding_table.csv"
MINE = r"C:\Users\admin\Desktop\MCP\_recoding\A\codes_A.csv"
OUT = r"C:\Users\admin\Desktop\MCP\_recoding\A\agreement_report_A.txt"

CATS = ["a", "b", "c", "e", "U"]

buf = []
def p(s=""):
    buf.append(s)
    print(s)


def main():
    orig = {}
    with open(SRC, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            if (r.get("background_class") or "").strip() in {"a", "b", "c", "e"}:
                orig[r["uid"]] = r

    mine = {}
    with open(MINE, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            mine[r["uid"]] = r

    uids = sorted(set(orig) & set(mine), key=int)
    n = len(uids)

    p("=" * 78)
    p("CODER A INDEPENDENT RECODING - AGREEMENT REPORT")
    p("=" * 78)
    p("")
    p("Records joined: %d (original classified %d, coder A %d)" % (n, len(orig), len(mine)))
    p("")

    pairs = [(orig[u]["background_class"].strip(), mine[u]["coder_code"].strip()) for u in uids]
    agree = sum(1 for o, a in pairs if o == a)
    po = agree / n

    oc = Counter(o for o, _ in pairs)
    ac = Counter(a for _, a in pairs)
    cats_used = [c for c in CATS if oc.get(c, 0) or ac.get(c, 0)]
    pe = sum((oc.get(c, 0) / n) * (ac.get(c, 0) / n) for c in cats_used)
    kappa = (po - pe) / (1 - pe) if pe != 1 else float("nan")

    p("(i) RAW AGREEMENT")
    p("    %d/%d = %.3f (%.1f%%)" % (agree, n, po, 100 * po))
    p("")
    p("(ii) COHEN'S KAPPA (categories %s; U treated as its own category)" % ", ".join(cats_used))
    p("    po = %.4f" % po)
    p("    pe = %.4f" % pe)
    p("    kappa = %.4f" % kappa)
    p("    Landis & Koch band: " + ("slight" if kappa < 0.21 else "fair" if kappa < 0.41 else "moderate" if kappa < 0.61 else "substantial" if kappa < 0.81 else "almost perfect"))
    p("")
    p("    Original code distribution: %s" % dict(sorted(oc.items())))
    p("    Coder A  code distribution: %s" % dict(sorted(ac.items())))
    p("")
    p("(iii) CONFUSION MATRIX (rows = ORIGINAL, cols = CODER A)")
    hdr = "    %-10s" % "orig\\A" + "".join("%8s" % c for c in cats_used) + "%8s" % "total"
    p(hdr)
    for o in cats_used:
        row = [sum(1 for oo, aa in pairs if oo == o and aa == a) for a in cats_used]
        p("    %-10s" % o + "".join("%8d" % v for v in row) + "%8d" % sum(row))
    p("    %-10s" % "total" + "".join("%8d" % ac.get(a, 0) for a in cats_used) + "%8d" % n)
    p("")

    p("(iv) ROWS THE ORIGINAL CODER CALLED b OR c")
    flagged = [u for u in uids if orig[u]["background_class"].strip() in {"b", "c"}]
    p("     total: %d" % len(flagged))
    p("")
    agree_bc = 0
    for u in flagged:
        o = orig[u]; m = mine[u]
        same = o["background_class"].strip() == m["coder_code"].strip()
        if same:
            agree_bc += 1
        p("     --- uid %s ---" % u)
        p("       original code = %s (%s)" % (o["background_class"].strip(), (o["background_class_label"] or "").strip()))
        p("       coder A code  = %s   [%s confidence]" % (m["coder_code"], m["confidence"]))
        p("       agree         = %s" % ("YES" if same else "NO"))
        p("       title         = %s" % o["title"].strip())
        p("       doi           = %s" % (o["doi"] or "(none)"))
        p("       modification  = %s | study_type = %s" % (o["modification"], o["study_type"]))
        p("       original ambiguity note: %s" % ((o["ambiguity"] or "(empty)").strip()))
        p("       coder A justification: %s" % m["justification"])
        p("")
    p("     agreement on these %d rows: %d/%d" % (len(flagged), agree_bc, len(flagged)))
    p("")

    p("(v) ALL DISAGREEMENTS")
    dis = [u for u in uids if orig[u]["background_class"].strip() != mine[u]["coder_code"].strip()]
    p("     %d of %d rows disagree (%.0f%%)" % (len(dis), n, 100 * len(dis) / n))
    p("     direction summary (original -> coder A):")
    dirs = Counter((orig[u]["background_class"].strip(), mine[u]["coder_code"].strip()) for u in dis)
    for k, v in sorted(dirs.items(), key=lambda kv: -kv[1]):
        p("       %s -> %s : %d" % (k[0], k[1], v))
    p("")

    p("(5) REVISED COUNTS")
    strict = sum(1 for _, a in pairs if a in {"b", "c"})
    p("     coder A strict count of b-or-c : %d of %d" % (strict, n))
    p("     original count of b-or-c       : %d of %d" % (len(flagged), n))
    p("")
    p("     Coder A lenient count: rows coder A would defend as a genuine")
    p("     matched background (same-experiment background of detected but")
    p("     unassigned residues, or a background matched on detectability/")
    p("     abundance): uid 441, uid 480 (strict c) plus uid 1111, whose")
    p("     background is the full set of phosphosites identified in the same")
    p("     experiment - a same-experiment, detectability-matched reference,")
    p("     but one composed of modified sites, which is why it is not c under")
    p("     the strictness rule.")
    p("     lenient count = 3 of %d" % n)
    p("")
    p("     Coder A b-codes: %d (none). Coder A found no record whose recorded" % ac.get("b", 0))
    p("     evidence states the background was matched on detectability or")
    p("     abundance.")
    p("")

    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        fh.write("\n".join(buf) + "\n")
    print("written: %s" % OUT)


if __name__ == "__main__":
    main()
