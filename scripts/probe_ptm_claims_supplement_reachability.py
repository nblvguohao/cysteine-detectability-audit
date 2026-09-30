#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""实测每篇被提取论文的补充位点表**在本环境里能不能真的拿到**，结果写成表供论断表机械引用。

为什么要单独测：可再检验性里"位点表是否可得"若靠出版社政策的印象来填，就是猜。
本脚本只记录实测结果，三种取值，判读规则运行前写定：

  in_hand                     本树 external/intake/ 或 external/ 下已持有该文的位点表/正负集
                              （逐篇写死在 IN_HAND 里，指向具体目录）
  verified_downloadable_here  Europe PMC supplementaryFiles 端点返回 zip（PK 魔数）且 > 50 kB
  not_reachable_from_sandbox  端点 404、或返回非 zip（通常是 164 字节的空响应 XML），
                              即需要出版社域名（沙箱未放行）才能取

**不做**的事：不判断"出版社是否允许公开下载"（那是政策，不是实测），不下载任何补充材料。
付费墙论文一律不下载、不缓存、不分发。

产物（文件名独立）：
  results/ptm_claims_supplement_reachability.csv
  results/ptm_claims_supplement_reachability_audit.json（含本脚本 sha256、逐篇 HTTP 状态与字节数）
"""

import csv
import hashlib
import json
import os
import sys
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_JSON = os.path.join(ROOT, "inputs", "ptm_claims_extraction_input_2026-09-15.json")
OUT_CSV = os.path.join(ROOT, "results", "ptm_claims_supplement_reachability.csv")
OUT_AUDIT = os.path.join(ROOT, "results", "ptm_claims_supplement_reachability_audit.json")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
MIN_ZIP_BYTES = 50000

# 本树已持有位点表/正负集的论文 -> 具体目录（相对 ROOT）
IN_HAND = {
    "10.1089/ars.2019.7777": "external/intake/ars2020_qtrp",
    "10.1038/srep29808": "external/intake/scirep2016_qpers_sid",
    "10.1093/bioinformatics/btaf078": "external/sul_bertgru",
}


def sha256_of_self():
    with open(os.path.abspath(__file__), "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def probe(pmcid):
    if not pmcid:
        return {"http_status": None, "bytes": 0, "magic": "", "verdict": "no_pmcid"}
    try:
        with urllib.request.urlopen(EPMC + "/" + pmcid + "/supplementaryFiles", timeout=90) as fh:
            body = fh.read()
            status = fh.status
    except Exception as exc:  # noqa: BLE001
        return {"http_status": getattr(exc, "code", None), "bytes": 0, "magic": "",
                "verdict": "http_error", "error": str(exc)[:120]}
    magic = body[:2].decode("latin-1", "replace")
    ok = magic == "PK" and len(body) > MIN_ZIP_BYTES
    return {"http_status": status, "bytes": len(body), "magic": magic,
            "verdict": "zip_ok" if ok else "not_zip_or_too_small"}


def main():
    payload = json.load(open(IN_JSON, encoding="utf-8"))
    papers = {}
    for c in payload["claims"]:
        papers.setdefault(c["doi"], {"paper_short": c["paper_short"], "pmid": c["pmid"],
                                     "pmcid": c["pmcid"], "family": c["family"],
                                     "declared_availability": c["sites_table_availability"],
                                     "claims": []})
        papers[c["doi"]]["claims"].append(c["claim_id"])

    rows, raw = [], {}
    for doi, meta in sorted(papers.items()):
        res = probe(meta["pmcid"])
        raw[doi] = res
        if doi in IN_HAND and os.path.isdir(os.path.join(ROOT, IN_HAND[doi])):
            reach = "in_hand"
            locator = IN_HAND[doi]
        elif res["verdict"] == "zip_ok":
            reach = "verified_downloadable_here"
            locator = "Europe PMC %s/supplementaryFiles (%d bytes)" % (meta["pmcid"], res["bytes"])
        else:
            reach = "not_reachable_from_sandbox"
            locator = "需出版社域名（沙箱未放行）；Europe PMC 端点结果：%s" % res["verdict"]
        rows.append({"doi": doi, "paper_short": meta["paper_short"], "family": meta["family"],
                     "pmid": meta["pmid"], "pmcid": meta["pmcid"],
                     "declared_availability": meta["declared_availability"],
                     "sites_table_reachability": reach, "reachability_locator": locator,
                     "http_status": res["http_status"], "bytes": res["bytes"],
                     "n_claims": len(meta["claims"]), "claim_ids": ";".join(sorted(meta["claims"]))})
        time.sleep(0.4)

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    counts = {}
    for r in rows:
        counts[r["sites_table_reachability"]] = counts.get(r["sites_table_reachability"], 0) + 1
    audit = {
        "script": os.path.relpath(os.path.abspath(__file__), ROOT),
        "script_sha256": sha256_of_self(),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "rules": {"in_hand": sorted(IN_HAND), "min_zip_bytes": MIN_ZIP_BYTES,
                  "endpoint": EPMC + "/{pmcid}/supplementaryFiles",
                  "note": "只记实测，不判断出版社政策；不下载任何补充材料"},
        "papers_probed": len(rows), "counts_by_reachability": counts, "raw_probe": raw,
        "outputs": {"csv": os.path.relpath(OUT_CSV, ROOT)},
    }
    with open(OUT_AUDIT, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    print(json.dumps({"papers": len(rows), "by_reachability": counts}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
