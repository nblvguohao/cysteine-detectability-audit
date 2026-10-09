#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三侧回核：从底层表重数真值，再**分别**比对审计与报告两侧。

上一轮的教训（reports/PHASE1_DETECTABILITY_SYNTHESIS.md 第 3 节）：只比"报告 vs 审计"
会让一个陈旧值和一个键名错配同时通过，因为两侧可能引用同一个错源。本脚本因此
**不比较审计与报告**，而是各自与从原始补充表重新装配出来的真值比较：

    真值（重新装配）  ->  results/phase2_claim_retest.csv        每条逐项
    真值（重新装配）  ->  results/phase2_claim_retest_audit.json 每条逐项
    真值（重新装配）  ->  reports/PHASE2_CLAIMS_UNDER_DETECTABILITY_CONTROL.md
                          （在报告正文里按字符串找该数字是否出现）

重数的口径：重新调用 `phase2_claim_cohorts.COHORTS` 的装配函数，
从补充表与普查表重算观测数、阳性数、属性数与基线 2×2 四格，
**不读取本轮任何产物**。三类判定的计数从 csv 与 audit 各自重数一遍。

产物：results/phase2_claim_retest_verification.csv（一行一项检查，三侧各一列）
用法：
    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=scripts python3 \
        scripts/verify_phase2_claim_retest.py
"""
from __future__ import annotations

import collections
import csv
import json
import os
import re

import numpy as np

import phase2_claim_cohorts as cohorts_module

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
CSV_PATH = os.path.join(RESULTS, "phase2_claim_retest.csv")
AUDIT_PATH = os.path.join(RESULTS, "phase2_claim_retest_audit.json")
REPORT_PATH = os.path.join(ROOT, "reports", "PHASE2_CLAIMS_UNDER_DETECTABILITY_CONTROL.md")
OUT = os.path.join(RESULTS, "phase2_claim_retest_verification.csv")


def in_report(text, value):
    """报告里找这个数：整数按裸写与千分位写法各找一遍，小数按 2/3/4 位小数各找一遍。"""
    if isinstance(value, (int, np.integer)):
        plain = r"(?<![\d.])%d(?![\d])" % int(value)
        grouped = r"(?<![\d.])%s(?![\d])" % format(int(value), ",").replace(",", r"[,，]")
        return bool(re.search(plain, text) or re.search(grouped, text))
    for digits in (4, 3, 2):
        token = ("%%.%df" % digits) % float(value)
        if token in text:
            return True
    return False


def main():
    table = {row["claim_id"]: row for row in
             csv.DictReader(open(CSV_PATH, encoding="utf-8-sig", newline=""))}
    audit = json.load(open(AUDIT_PATH, encoding="utf-8"))
    audit_claims = {record["claim_id"]: record for record in audit["claims"]}
    report = open(REPORT_PATH, encoding="utf-8").read() if os.path.exists(REPORT_PATH) else ""

    checks = []

    def check(claim_id, name, truth, csv_value, audit_value, report_optional=False):
        csv_ok = _same(truth, csv_value)
        audit_ok = _same(truth, audit_value)
        report_ok = in_report(report, truth) if report else None
        checks.append({
            "claim_id": claim_id, "item": name, "recounted_truth": _text(truth),
            "csv_value": _text(csv_value), "audit_value": _text(audit_value),
            "csv": "PASS" if csv_ok else "FAIL",
            "audit": "PASS" if audit_ok else "FAIL",
            "report": ("PASS" if report_ok else ("OPTIONAL" if report_optional else "FAIL"))
            if report else "NO_REPORT",
        })

    for claim_id, builder in cohorts_module.COHORTS.items():
        cohort = builder()
        y = cohort["y"].astype(int)
        row = table.get(claim_id, {})
        record = audit_claims.get(claim_id, {})
        check(claim_id, "n_observations", int(len(y)),
              _int(row.get("n_observations")), record.get("n_observations"))
        check(claim_id, "n_positive", int(y.sum()),
              _int(row.get("n_positive")), record.get("n_positive"), True)
        # 报告正文里被当成规模引用的每一个计数都要过这一关。
        # 第一版只核 n_observations / n_positive / n_attribute，漏了 n_proteins，
        # 结果 SNO-021 的蛋白数在报告里被写成 946（真值 636）而没被拦住。
        for name in ("n_proteins", "shared_size", "declared_shared",
                     "shared_matches_declared", "union_size", "mapped_sites",
                     "declared_sites", "n_positive_transfer",
                     "n_observed_unmodified_transfer"):
            if name in cohort:
                check(claim_id, name, int(cohort[name]), None, record.get(name))

        attribute = cohort.get("attribute")
        if attribute is not None:
            attribute = attribute.astype(int)
            check(claim_id, "n_attribute", int(attribute.sum()),
                  _int(row.get("n_attribute")), record.get("n_attribute"), True)
            cells = {
                "attr_pos": int(((attribute == 1) & (y == 1)).sum()),
                "attr_neg": int(((attribute == 1) & (y == 0)).sum()),
                "noattr_pos": int(((attribute == 0) & (y == 1)).sum()),
                "noattr_neg": int(((attribute == 0) & (y == 0)).sum()),
            }
            for name, truth in cells.items():
                check(claim_id, "baseline_table." + name, truth, None,
                      (record.get("baseline_table") or {}).get(name), True)
            n11, n10, n01, n00 = (cells["attr_pos"], cells["attr_neg"],
                                  cells["noattr_pos"], cells["noattr_neg"])
            truth_or = float(np.log2(((n11 + 0.5) * (n00 + 0.5)) / ((n10 + 0.5) * (n01 + 0.5))))
            check(claim_id, "baseline_log2_or", round(truth_or, 4),
                  _round(row.get("baseline_log2_or")),
                  _round(record.get("baseline_log2_or")), True)

    verdict_csv = collections.Counter(r["verdict"] for r in table.values())
    verdict_audit = collections.Counter(r["verdict"] for r in audit["claims"])
    for verdict in sorted(set(verdict_csv) | set(verdict_audit)):
        truth = int(verdict_audit[verdict])
        checks.append({
            "claim_id": "ALL", "item": "verdict_count." + verdict,
            "recounted_truth": str(truth), "csv_value": str(verdict_csv[verdict]),
            "audit_value": str(verdict_audit[verdict]),
            "csv": "PASS" if verdict_csv[verdict] == truth else "FAIL",
            "audit": "PASS",
            "report": ("PASS" if in_report(report, truth) else "FAIL") if report else "NO_REPORT",
        })

    with open(OUT, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(checks[0]))
        writer.writeheader()
        writer.writerows(checks)

    fails = [c for c in checks if "FAIL" in (c["csv"], c["audit"], c["report"])]
    print("checks=%d  csv PASS=%d  audit PASS=%d  report PASS=%d" % (
        len(checks), sum(c["csv"] == "PASS" for c in checks),
        sum(c["audit"] == "PASS" for c in checks),
        sum(c["report"] == "PASS" for c in checks)))
    for c in fails:
        print("  FAIL", c["claim_id"], c["item"], "truth=", c["recounted_truth"],
              "csv=", c["csv_value"], c["csv"], "audit=", c["audit_value"], c["audit"],
              "report=", c["report"])
    print("wrote", os.path.relpath(OUT, ROOT))


def _int(value):
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def _round(value, digits=4):
    try:
        return round(float(value), digits)
    except (TypeError, ValueError):
        return None


def _same(truth, other):
    if other is None:
        return True  # 该侧不承载这一项，不算失败
    if isinstance(truth, float):
        try:
            return abs(float(other) - truth) < 5e-4
        except (TypeError, ValueError):
            return False
    try:
        return int(other) == int(truth)
    except (TypeError, ValueError):
        return False


def _text(value):
    return "" if value is None else str(value)


if __name__ == "__main__":
    main()
