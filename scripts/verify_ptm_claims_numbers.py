#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三方回核：从底层表重数真值，再分别比对报告与审计 JSON 两侧（纪律第 4 条）。

真值来源（**只从底层表重数，不读任何中间打印**）：
  results/ptm_site_preference_claims.csv            逐论断一行
  results/ptm_claims_literature_search.csv          逐论文一行（持硫化/SNO/亚磺酰化轮）
  results/sulfenylation_claims_literature_search.csv 逐论文一行（范围纠正后的 sulfenylation 轮）
  results/ptm_claims_candidate_sentences.csv        候选句（定位用）

比对对象：
  A 侧 reports/PTM_SITE_PREFERENCE_CLAIMS.md 中的"数字对照表"块，
     每行形如  `- KEY = VALUE`，KEY 取自下面 KEYS 列表。
  B 侧 results/ptm_claims_extraction_audit.json 的 counts_recounted_from_written_table
     与 search_stage.counts。

任一侧与真值不一致即 exit 1 并逐条列出差异；一致则写
results/ptm_claims_number_verification.json（含本脚本 sha256 与三方数值）。
"""

import csv
import hashlib
import json
import os
import re
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLAIMS = os.path.join(ROOT, "results", "ptm_site_preference_claims.csv")
SEARCH = os.path.join(ROOT, "results", "ptm_claims_literature_search.csv")
SEARCH_SFE = os.path.join(ROOT, "results", "sulfenylation_claims_literature_search.csv")
SFE_AUDIT = os.path.join(ROOT, "results", "sulfenylation_claims_search_audit.json")
CAND = os.path.join(ROOT, "results", "ptm_claims_candidate_sentences.csv")
AUDIT = os.path.join(ROOT, "results", "ptm_claims_extraction_audit.json")
REPORT = os.path.join(ROOT, "reports", "PTM_SITE_PREFERENCE_CLAIMS.md")
OUT = os.path.join(ROOT, "results", "ptm_claims_number_verification.json")

KEYS = [
    "SEARCH_UNIQUE_RECORDS", "SEARCH_SCREEN_PASS", "SEARCH_FULLTEXT_OBTAINED",
    "SEARCH_SFE_UNIQUE_RECORDS", "SEARCH_SFE_SCREEN_PASS", "SEARCH_SFE_FULLTEXT_OBTAINED",
    "SEARCH_RECORDS_ALL_ROUNDS",
    "PAPERS_ENTERED_EXTRACTION", "CLAIMS_TOTAL", "CLAIMS_RETESTABLE", "CLAIMS_RETESTABLE_NOW",
    "CLAIMS_TIER_A", "CLAIMS_TIER_B", "CLAIMS_TIER_C",
    "CLAIMS_PERSULFIDATION", "CLAIMS_SNITROSYLATION", "CLAIMS_SULFENYLATION", "CLAIMS_SULFINYLATION",
    "CLAIMS_ABSTRACT_ONLY", "CLAIMS_FULLTEXT_METHODS", "CANDIDATE_SENTENCES",
]


def sha256_file(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def main():
    claims, search, cand = read_csv(CLAIMS), read_csv(SEARCH), read_csv(CAND)
    sfe = read_csv(SEARCH_SFE) if os.path.exists(SEARCH_SFE) else []
    truth = {
        "SEARCH_UNIQUE_RECORDS": len(search),
        "SEARCH_SCREEN_PASS": sum(1 for r in search if r["screen_pass"] == "True"),
        "SEARCH_FULLTEXT_OBTAINED": sum(1 for r in search if r["fulltext_cached"] == "True"),
        "SEARCH_SFE_UNIQUE_RECORDS": len(sfe),
        "SEARCH_SFE_SCREEN_PASS": sum(1 for r in sfe if r["screen_pass"] == "True"),
        "SEARCH_SFE_FULLTEXT_OBTAINED": sum(1 for r in sfe if r["fulltext_cached"] == "True"),
        "SEARCH_RECORDS_ALL_ROUNDS": len(search) + len(sfe),
        "CLAIMS_RETESTABLE_NOW": sum(1 for r in claims if r["retestable_now"] == "True"),
        "CLAIMS_SULFENYLATION": sum(1 for r in claims if r["family"] == "sulfenylation"),
        "PAPERS_ENTERED_EXTRACTION": len({r["doi"] or r["pmid"] or r["paper_short"] for r in claims}),
        "CLAIMS_TOTAL": len(claims),
        "CLAIMS_RETESTABLE": sum(1 for r in claims if r["retestable"] == "True"),
        "CLAIMS_TIER_A": sum(1 for r in claims if r["priority_tier"] == "A"),
        "CLAIMS_TIER_B": sum(1 for r in claims if r["priority_tier"] == "B"),
        "CLAIMS_TIER_C": sum(1 for r in claims if r["priority_tier"] == "C"),
        "CLAIMS_PERSULFIDATION": sum(1 for r in claims if r["family"] == "persulfidation"),
        "CLAIMS_SNITROSYLATION": sum(1 for r in claims if r["family"] == "s_nitrosylation"),
        "CLAIMS_SULFINYLATION": sum(1 for r in claims if r["family"] == "sulfinylation"),
        "CLAIMS_ABSTRACT_ONLY": sum(1 for r in claims if r["evidence_level"] == "abstract_only"),
        "CLAIMS_FULLTEXT_METHODS": sum(1 for r in claims if r["evidence_level"] == "fulltext_methods"),
        "CANDIDATE_SENTENCES": len(cand),
    }

    audit = json.load(open(AUDIT, encoding="utf-8"))
    c = audit["counts_recounted_from_written_table"]
    sc = (audit.get("search_stage") or {}).get("counts") or {}
    sfe_audit = json.load(open(SFE_AUDIT, encoding="utf-8")) if os.path.exists(SFE_AUDIT) else {}
    sfc = sfe_audit.get("counts") or {}
    audit_side = {
        "SEARCH_UNIQUE_RECORDS": sc.get("unique_records"),
        "SEARCH_SCREEN_PASS": sc.get("screen_pass"),
        "SEARCH_FULLTEXT_OBTAINED": sc.get("fulltext_obtained"),
        "SEARCH_SFE_UNIQUE_RECORDS": sfc.get("unique_records"),
        "SEARCH_SFE_SCREEN_PASS": sfc.get("screen_pass"),
        "SEARCH_SFE_FULLTEXT_OBTAINED": sfc.get("fulltext_obtained"),
        "SEARCH_RECORDS_ALL_ROUNDS": (sc.get("unique_records") or 0) + (sfc.get("unique_records") or 0),
        "CLAIMS_RETESTABLE_NOW": c.get("retestable_now_true"),
        "CLAIMS_SULFENYLATION": (c.get("claims_by_family") or {}).get("sulfenylation"),
        "PAPERS_ENTERED_EXTRACTION": c.get("papers_entered_extraction"),
        "CLAIMS_TOTAL": c.get("claims_total"),
        "CLAIMS_RETESTABLE": c.get("retestable_true"),
        "CLAIMS_TIER_A": (c.get("by_tier") or {}).get("A"),
        "CLAIMS_TIER_B": (c.get("by_tier") or {}).get("B"),
        "CLAIMS_TIER_C": (c.get("by_tier") or {}).get("C"),
        "CLAIMS_PERSULFIDATION": (c.get("claims_by_family") or {}).get("persulfidation"),
        "CLAIMS_SNITROSYLATION": (c.get("claims_by_family") or {}).get("s_nitrosylation"),
        "CLAIMS_SULFINYLATION": (c.get("claims_by_family") or {}).get("sulfinylation"),
        "CLAIMS_ABSTRACT_ONLY": (c.get("by_evidence_level") or {}).get("abstract_only"),
        "CLAIMS_FULLTEXT_METHODS": (c.get("by_evidence_level") or {}).get("fulltext_methods"),
        "CANDIDATE_SENTENCES": (audit.get("candidate_sentence_stage") or {}).get("counts", {}).get("candidate_sentences"),
    }

    report_side, report_text = {}, ""
    if os.path.exists(REPORT):
        report_text = open(REPORT, encoding="utf-8").read()
        for k in KEYS:
            m = re.search(r"^-\s*" + re.escape(k) + r"\s*=\s*([0-9]+)\s*$", report_text, flags=re.M)
            report_side[k] = int(m.group(1)) if m else None

    diffs = []
    for k in KEYS:
        if audit_side.get(k) != truth[k]:
            diffs.append("AUDIT  %s: audit=%r truth=%r" % (k, audit_side.get(k), truth[k]))
        if report_side.get(k) != truth[k]:
            diffs.append("REPORT %s: report=%r truth=%r" % (k, report_side.get(k), truth[k]))

    payload = {
        "script": os.path.relpath(os.path.abspath(__file__), ROOT),
        "script_sha256": sha256_file(os.path.abspath(__file__)),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "truth_from_tables": truth,
        "audit_side": audit_side,
        "report_side": report_side,
        "differences": diffs,
        "all_three_sides_agree": not diffs,
        "input_sha256": {os.path.relpath(p, ROOT): sha256_file(p)
                         for p in (CLAIMS, SEARCH, SEARCH_SFE, CAND, AUDIT, SFE_AUDIT) if os.path.exists(p)},
    }
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    print(json.dumps({"all_three_sides_agree": not diffs, "n_differences": len(diffs),
                      "truth": truth}, ensure_ascii=False, indent=1))
    if diffs:
        print("\n".join(diffs), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
