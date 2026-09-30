import csv, os
from collections import Counter, defaultdict

SRC = r"C:\Users\admin\Desktop\MCP\supplemental\Supplemental_Data_9_survey_coding_table.csv"
CODES_B = r"C:\Users\admin\Desktop\MCP\_recoding\B\codes_B.csv"
OUT = r"C:\Users\admin\Desktop\MCP\_recoding\B\agreement_report_B.txt"

with open(SRC, newline="", encoding="utf-8-sig") as f:
    orig = {int(r["uid"]): r for r in csv.DictReader(f)}
with open(CODES_B, newline="", encoding="utf-8") as f:
    B = {int(r["uid"]): r for r in csv.DictReader(f)}

CATS = ["a", "b", "c", "e", "U"]
n = len(B)
pairs = []
for uid in sorted(B):
    o = (orig[uid].get("background_class") or "").strip()
    bcode = B[uid]["coder_code"].strip()
    pairs.append((uid, o, bcode))

agree = sum(1 for _, o, b in pairs if o == b)
po = agree / n

# Cohen's kappa, U treated as its own category
o_counts = Counter(o for _, o, _ in pairs)
b_counts = Counter(b for _, _, b in pairs)
pe = sum((o_counts[c] / n) * (b_counts[c] / n) for c in CATS)
kappa = (po - pe) / (1 - pe) if pe != 1 else float("nan")

conf = defaultdict(int)
for _, o, b in pairs:
    conf[(o, b)] += 1

lines = []
def P(s=""):
    lines.append(s)
    print(s)

P("=" * 78)
P("CODER B - INDEPENDENT SECOND CODING PASS: AGREEMENT REPORT")
P("=" * 78)
P("")
P("Source table : Supplemental_Data_9_survey_coding_table.csv (942 rows)")
P("Subset coded : 74 rows with background_class in {a,b,c,e}")
P("Coder B codes: codes_B.csv (blinded pass; U = unclassifiable from retrieved evidence)")
P("")
P("-" * 78)
P("(i) RAW AGREEMENT")
P("-" * 78)
P("Rows compared            : %d" % n)
P("Exact agreements         : %d" % agree)
P("Raw agreement (Po)       : %.4f  (%.1f%%)" % (po, 100 * po))
P("")
P("-" * 78)
P("(ii) COHEN'S KAPPA (categories a/b/c/e/U, U as its own category)")
P("-" * 78)
P("Observed agreement  Po   : %.4f" % po)
P("Expected agreement  Pe   : %.4f" % pe)
P("Cohen's kappa            : %.4f" % kappa)
P("")
P("Marginal distributions:")
P("  original code : " + ", ".join("%s=%d" % (c, o_counts.get(c, 0)) for c in CATS))
P("  coder B code  : " + ", ".join("%s=%d" % (c, b_counts.get(c, 0)) for c in CATS))
P("")
P("-" * 78)
P("(iii) CONFUSION MATRIX  (rows = ORIGINAL code, cols = CODER B code)")
P("-" * 78)
hdr = "orig\\B " + "".join("%6s" % c for c in CATS) + "%8s" % "total"
P(hdr)
for o in CATS:
    row = "".join("%6d" % conf.get((o, c), 0) for c in CATS)
    P("%-7s" % o + row + "%8d" % o_counts.get(o, 0))
P("%-7s" % "total" + "".join("%6d" % b_counts.get(c, 0) for c in CATS) + "%8d" % n)
P("")
P("-" * 78)
P("(iv) ROWS WHERE THE ORIGINAL CODE IS b OR c")
P("-" * 78)
bc_rows = [(uid, o, b) for uid, o, b in pairs if o in ("b", "c")]
P("Count of original b/c rows: %d" % len(bc_rows))
P("")
for i, (uid, o, b) in enumerate(bc_rows, 1):
    r = orig[uid]
    P("--- original %s/b/c row %d ---" % (o, i))
    P("  uid              : %s" % uid)
    P("  original code    : %s   (label: %s)" % (o, r.get("background_class_label", "")))
    P("  coder B code     : %s" % b)
    P("  agree?           : %s" % ("YES" if o == b else "NO"))
    P("  doi              : %s" % r.get("doi", ""))
    P("  title            : %s" % r.get("title", ""))
    P("  original ambiguity field: %s" % (r.get("ambiguity", "") or "(empty)"))
    P("  coder B justification   : %s" % B[uid]["justification"])
    P("")
P("-" * 78)
P("(v) REVISED COUNT OF MATCHED-BACKGROUND PAPERS OUT OF 74")
P("-" * 78)
b_cnt = b_counts.get("b", 0)
c_cnt = b_counts.get("c", 0)
u_cnt = b_counts.get("U", 0)
P("Coder B, STRICT (a/b/c/e as coded; U excluded):")
P("  b : %d" % b_cnt)
P("  c : %d" % c_cnt)
P("  b or c TOTAL : %d of 74" % (b_cnt + c_cnt))
P("  unclassifiable from retrieved evidence (U) : %d of 74" % u_cnt)
P("  (remaining: a=%d, e=%d)" % (b_counts.get("a", 0), b_counts.get("e", 0)))
P("")
P("Coder B, LENIENT (rows I would defend as a genuine matched background in front")
P("of a reviewer, including low-confidence calls):")
P("  b or c TOTAL : %d of 74" % (b_cnt + c_cnt))
P("")
P("Sensitivity: if uid 480 (native-MS detected set, also readable as 'all cysteines")
P("of the analysed proteins' = a) is read as a rather than c, the count falls to 1.")
P("")
P("-" * 78)
P("(vi) DISAGREEMENTS WITH THE ORIGINAL CODE")
P("-" * 78)
dis = [(uid, o, b) for uid, o, b in pairs if o != b]
P("Disagreeing rows: %d of %d" % (len(dis), n))
P("")
for uid, o, b in dis:
    r = orig[uid]
    P("  uid %-5s orig=%s  B=%s | %s" % (uid, o, b, (r.get("title") or "")[:78]))
    P("        reason: %s" % B[uid]["justification"])
P("")
P("-" * 78)
P("(vii) LIMITATIONS OF THIS CODING PASS")
P("-" * 78)
P("* Coding used ONLY the recorded/truncated evidence quotations in the blinded")
P("  extract. The papers' full texts were not consulted.")
P("* Many evidence_quote_background fields are cut off mid-clause at a fixed length,")
P("  in several cases immediately before the background description. Those rows are")
P("  coded U and cannot be independently verified from the record.")
P("* 12 rows carry the note 'explicit background statement found on full-text")
P("  re-grep (two-pass verification)' but the statement itself was not recorded;")
P("  11 of these are U. The audit trail is incomplete for those rows.")
P("* %d of 74 rows (%.0f%%) are U: the class cannot be reproduced from the evidence" % (u_cnt, 100 * u_cnt / n))
P("  that was retained. Any agreement statistic should be read alongside this.")

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print("\nwrote", OUT)
