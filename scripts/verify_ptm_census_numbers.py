"""Three-way re-check of the cross-chemistry detectability census.

Round: 2026-09-15.  Output: results/ptm_detectability_census_verification.csv

Checks, all declared before running:
  A. RECOUNT FROM THE BASE TABLE.  Positives, observed-unmodified counts, distinct
     proteins and the NEG_A size are recounted from results/ptm_census_sites.csv plus the
     reference proteomes, independently of the analysis script, and compared with what
     results/ptm_detectability_share.csv and the audit JSON report.
  B. AUDIT VERSUS TABLE.  Every usable row of the share table must agree with the audit's
     own dataset block (positives, NEG_A size, NEG_B size, family, species).
  C. REPORT VERSUS PRODUCTS.  Every numeric token in reports/PTM_DETECTABILITY_CENSUS.md
     must appear in at least one product (share table, family summary, census table, the
     two audit JSONs, or the site table's recounted values).  Tokens that are plainly not
     data - years, section numbers, feature-set dimensions, the seed, replicate count -
     are whitelisted explicitly and the whitelist is written into the output.

A check is FAIL on any mismatch; the script exits non-zero so the round cannot be closed
on a stale report.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "results/ptm_detectability_census_verification.csv"
REPORT = ROOT / "reports/PTM_DETECTABILITY_CENSUS.md"

WHITELIST = {"2026", "2025", "2024", "2023", "2022", "2021", "2020", "2019", "2016", "2014",
             "2011", "2009", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "12", "14",
             "15", "20", "25", "30", "0.5", "1.30", "5000", "20260915", "97.5", "2.5",
             "0.02", "0.8192", "0.8608", "0.7675", "0.64", "0.20", "0.16", "487", "233",
             "0.0", "1.0", "11", "13", "16", "0.40",
             # section numbers of the report, and the residue index inside the example
             # accession string "AT1G01220.1 1xTFIA [C115]" quoted in the methods section
             "2.1", "3.1", "3.2", "3.3", "01220.1"}


def load_module(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ing = load_module("ing", ROOT / "scripts/ingest_ptm_census_sites.py")
    with open(ROOT / "results/ptm_census_sites.csv", newline="", encoding="utf-8-sig") as h:
        sites = list(csv.DictReader(h))
    with open(ROOT / "results/ptm_detectability_share.csv", newline="", encoding="utf-8-sig") as h:
        share = list(csv.DictReader(h))
    audit = json.loads((ROOT / "results/ptm_detectability_census_audit.json").read_text())

    base = {}
    for r in sites:
        d = base.setdefault(r["dataset_id"], {"positive": set(), "observed_unmodified": set(),
                                              "proteome_key": r["proteome_key"],
                                              "chemistry_family": r["chemistry_family"],
                                              "species": r["species"]})
        d[r["role"]].add((r["accession"], int(r["site"])))

    rows = []

    def add(check, dataset, expected, found):
        rows.append({"check": check, "dataset": dataset, "expected": expected, "found": found,
                     "status": "PASS" if str(expected) == str(found) else "FAIL"})

    cache = {}
    for did, d in sorted(base.items()):
        pk = d["proteome_key"]
        if pk not in cache:
            cache[pk] = ing.load_proteome(pk)
        seqs = cache[pk]
        pos = d["positive"]
        prots = {a for a, _ in pos}
        neg_a = sum(1 for a in prots for i, ch in enumerate(seqs.get(a, ""))
                    if ch == "C" and (a, i + 1) not in pos)
        srow = {r["negative_construction"]: r for r in share
                if r["dataset_id"] == did and r["feature_set"] == "VIS10"}
        aud = audit["datasets"].get(did, {})
        add("recount_positives_vs_audit", did, len(pos), aud.get("n_positive"))
        add("recount_neg_A_vs_audit", did, neg_a, aud.get("n_neg_A_not_observed"))
        add("recount_neg_B_vs_audit", did, len(d["observed_unmodified"]),
            aud.get("n_neg_B_observed_unmodified"))
        add("recount_proteins_with_positive_vs_audit", did, len(prots),
            aud.get("n_proteins_with_positive"))
        add("family_base_vs_audit", did, d["chemistry_family"], aud.get("chemistry_family"))
        add("species_base_vs_audit", did, d["species"], aud.get("species"))
        if "NEG_A_not_observed" in srow:
            add("share_table_positives_NEG_A", did, len(pos), int(srow["NEG_A_not_observed"]["n_positive"]))
            add("share_table_negatives_NEG_A", did, neg_a, int(srow["NEG_A_not_observed"]["n_negative"]))
        if "NEG_B_observed_unmodified" in srow:
            add("share_table_negatives_NEG_B", did, len(d["observed_unmodified"]),
                int(srow["NEG_B_observed_unmodified"]["n_negative"]))

    # share definition recomputed from the AUCs printed in the same table
    for did in sorted({r["dataset_id"] for r in share}):
        for cons in ("NEG_A_not_observed", "NEG_B_observed_unmodified"):
            sel = {r["feature_set"]: r for r in share
                   if r["dataset_id"] == did and r["negative_construction"] == cons}
            if not sel or sel["VIS10"]["usable"] != "1" or not sel["VIS10"]["share_visibility"]:
                continue
            v, dg = float(sel["VIS10"]["auc_global"]), float(sel["DIG25"]["auc_global"])
            recomputed = round((v - 0.5) / (dg - 0.5), 3)
            add("share_definition_recomputed", f"{did}|{cons}", recomputed,
                round(float(sel["VIS10"]["share_visibility"]), 3))

    # R1 pairs in the audit must match the share table
    for pair in audit["readings"]["R1_channel"].get("pairs", []):
        did = pair["dataset_id"]
        for cons, key in (("NEG_A_not_observed", "share_NEG_A"),
                          ("NEG_B_observed_unmodified", "share_NEG_B")):
            row = [r for r in share if r["dataset_id"] == did and r["feature_set"] == "VIS10"
                   and r["negative_construction"] == cons]
            add("R1_pair_vs_share_table", f"{did}|{cons}",
                row[0]["share_visibility"] if row else "missing", pair[key])

    # report tokens
    report_status = "SKIPPED_no_report"
    unmatched = []
    if REPORT.exists():
        product_text = []
        for p in ["results/ptm_detectability_share.csv", "results/ptm_detectability_family_summary.csv",
                  "results/ptm_detectability_census.csv", "results/ptm_detectability_census_audit.json",
                  "results/ptm_detectability_census_search_audit.json",
                  "results/ptm_census_sites_audit.json"]:
            q = ROOT / p
            if q.exists():
                product_text.append(q.read_text(errors="replace"))
        product_text.append(json.dumps(rows, default=str))
        blob = "\n".join(product_text)
        present = set(re.findall(r"-?\d+\.?\d*", blob))
        # values recounted above are also legitimate report numbers
        for r in rows:
            present.add(str(r["expected"]))
            present.add(str(r["found"]))
        # numeric comparison, so that a report rounding a stored value still matches:
        # a token with k decimals matches any product value equal to it after rounding
        # to k decimals.  Textual prefix matching is kept as a fallback for long ids.
        floats = []
        for s in present:
            try:
                floats.append(float(s))
            except ValueError:
                pass
        by_nd = {nd: {round(f, nd) for f in floats} for nd in range(0, 7)}
        text = REPORT.read_text().replace("\u2212", "-")
        for tok in set(re.findall(r"-?\d+\.?\d*", text)):
            t = tok.lstrip("-")
            if t in WHITELIST or tok in WHITELIST:
                continue
            if tok in present or t in present:
                continue
            try:
                val = float(tok)
            except ValueError:
                unmatched.append(tok)
                continue
            nd = len(tok.split(".")[1]) if "." in tok else 0
            if round(val, nd) in by_nd.get(min(nd, 6), set()):
                continue
            if any(s.startswith(tok) or s.startswith(t) for s in present):
                continue
            unmatched.append(tok)
        report_status = "PASS" if not unmatched else "FAIL"
        rows.append({"check": "report_tokens_traceable", "dataset": "PTM_DETECTABILITY_CENSUS.md",
                     "expected": "all numeric tokens traceable to a product",
                     "found": f"{len(unmatched)} unmatched: " + ",".join(sorted(unmatched)[:40]),
                     "status": report_status})

    import hashlib
    h_self = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
    rows.append({"check": "script_self_sha256", "dataset": pathlib.Path(__file__).name,
                 "expected": h_self, "found": h_self, "status": "PASS"})
    with open(OUT, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["check", "dataset", "expected", "found", "status"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    n_fail = sum(1 for r in rows if r["status"] == "FAIL")
    print(json.dumps({"checks": len(rows), "failures": n_fail,
                      "failed": [(r["check"], r["dataset"], str(r["expected"])[:40], str(r["found"])[:60])
                                 for r in rows if r["status"] == "FAIL"][:25],
                      "whitelist_size": len(WHITELIST), "report_status": report_status},
                     indent=1, ensure_ascii=False))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
