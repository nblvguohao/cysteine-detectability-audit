"""Check the final claim tally, the claim estimates and the class-level retention tables against the public
supplemental tables of this repository (Supplemental Notes 14 and 15).

Run from anywhere:  python check_final_tables.py
Inputs (paths relative to the repository root):
  supplemental/Supplemental_Data_2_claim_verdict_tally.csv    final_verdict column
  supplemental/Supplemental_Data_4_retest_round_b.csv, _5_retest_round_d.csv, _6_retest_round_e.csv
  supplemental/Supplemental_Data_10_claim_independence.csv   unit (site or protein level)
  source_data_submitted/Source_Data_text_phase2c_sfe006_reproduction.json   the SFE-006 re-run on its authors' data
Checks:
  1. final_tally.csv equals the final_verdict column of Supplemental Data 2, and final_tally_counts.json its recount;
  2. every estimate in claim_estimates.csv equals the row of the re-test table it names;
  3. every retained fraction is the matched (or random-control) estimate over the baseline estimate;
  4. every class summary is the minimum, median and maximum of the per-claim retained fractions.
"""
import csv, io, json, pathlib, statistics, sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SUP, SRC = ROOT / "supplemental", ROOT / "source_data_submitted"


def read(p):
    return list(csv.DictReader(io.open(p, encoding="utf-8-sig")))


def num(x):
    return None if x in ("", None) else float(x)


fails = []


def check(ok, msg):
    if not ok:
        fails.append(msg)


# 1. verdicts
tally = {r["claim_id"]: r for r in read(HERE / "final_tally.csv")}
sd2 = {r["claim_id"]: r for r in read(SUP / "Supplemental_Data_2_claim_verdict_tally.csv")}
unit = {r["claim_id"]: r.get("unit", "") for r in read(SUP / "Supplemental_Data_10_claim_independence.csv")}
for cid, r in tally.items():
    check(sd2.get(cid, {}).get("final_verdict") == r["final_verdict"], f"{cid}: final verdict differs from Supplemental Data 2")
    check(unit.get(cid) == r["level"], f"{cid}: level differs from Supplemental Data 10")
counts = json.load(io.open(HERE / "final_tally_counts.json", encoding="utf-8"))
recount = {}
for r in tally.values():
    recount[r["final_verdict"]] = recount.get(r["final_verdict"], 0) + 1
check(recount == counts["overall"], f"recount {recount} differs from final_tally_counts.json {counts['overall']}")

# 2. estimates
sd4 = {r["claim_id"]: r for r in read(SUP / "Supplemental_Data_4_retest_round_b.csv")}
by_spec = {}
for t, fn in (("SD5", "Supplemental_Data_5_retest_round_d.csv"), ("SD6", "Supplemental_Data_6_retest_round_e.csv")):
    for r in read(SUP / fn):
        by_spec[(t, r["claim_id"], r["specification"])] = r
own = json.load(io.open(SRC / "Source_Data_text_phase2c_sfe006_reproduction.json", encoding="utf-8"))
own_row = {"baseline_log2_or": own["baseline_log2_or"], "baseline_ci_low": own["baseline_interval"][0],
           "baseline_ci_high": own["baseline_interval"][1], "matched_log2_or": own["matched_log2_or"],
           "matched_ci_low": own["matched_interval"][0], "matched_ci_high": own["matched_interval"][1],
           "random_control_log2_or": own["random_control_log2_or"], "random_control_ci_low": own["random_control_interval"][0],
           "random_control_ci_high": own["random_control_interval"][1]}
FIELDS = ["baseline_log2_or", "baseline_ci_low", "baseline_ci_high", "matched_log2_or", "matched_ci_low", "matched_ci_high",
          "random_control_log2_or", "random_control_ci_low", "random_control_ci_high"]
est = read(HERE / "claim_estimates.csv")
check(len(est) == 28, "claim_estimates.csv does not have 28 claims")
for r in est:
    cid = r["claim_id"]
    check(r["final_verdict"] == tally[cid]["final_verdict"], f"{cid}: verdict differs from final_tally.csv")
    if r["source_table"] == "SD4":
        src = sd4[cid]
        check(src["round"] == r["round"], f"{cid}: round differs from Supplemental Data 4")
    elif r["source_table"] in ("SD5", "SD6"):
        src = by_spec[(r["source_table"], cid, r["specification_used"])]
    else:
        src = own_row
    for f in FIELDS:
        a, b = num(r[f]), num(src[f])
        check(a is not None and b is not None and abs(a - b) < 1e-9, f"{cid}: {f} {a} differs from its source {b}")
    b, m, rc = num(r["baseline_log2_or"]), num(r["matched_log2_or"]), num(r["random_control_log2_or"])
    if r["retention_matched_over_baseline"]:
        check(abs(num(r["retention_matched_over_baseline"]) - m / b) < 1e-9, f"{cid}: matched retention is not matched/baseline")
        check(abs(num(r["retention_random_over_baseline"]) - rc / b) < 1e-9, f"{cid}: random retention is not random/baseline")

# 3-4. class-level retention
pc = read(HERE / "class_retention_per_claim.csv")
for r in pc:
    b = num(r["baseline_log2_or"])
    for col, est_col in (("retention_matched", "matched_log2_or"), ("retention_random", "random_control_log2_or")):
        if r[col]:
            check(abs(num(r[col]) - num(r[est_col]) / b) < 1e-9, f"{r['set']} {r['claim_id']}: {col} is not {est_col}/baseline")
for s in read(HERE / "class_retention_summary.csv"):
    rows = [r for r in pc if r["set"] == s["set"] and r["attribute_class"] == s["attribute_class"]]
    check(len(rows) == int(s["n_claims"]), f"{s['set']} {s['attribute_class']}: n_claims")
    for col in ("retention_matched", "retention_random"):
        vals = [num(r[col]) for r in rows if r[col]]
        check(len(vals) == int(s["n_estimable"]), f"{s['set']} {s['attribute_class']}: n_estimable")
        if vals:
            for stat, fn in (("min", min), ("max", max), ("median", statistics.median)):
                check(abs(fn(vals) - num(s[f"{stat}_{col}"])) < 1e-9, f"{s['set']} {s['attribute_class']}: {stat}_{col}")

if fails:
    print("FAILED:\n  " + "\n  ".join(fails))
    sys.exit(1)
print(f"all checks pass: {len(tally)} verdicts, {len(est)} claims x {len(FIELDS)} estimates, {len(pc)} retention rows")
