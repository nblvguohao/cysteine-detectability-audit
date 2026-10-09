#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 OA 全文缓存里机械抽取"可能承载位点偏好论断且自带统计量"的句子，供人工逐条阅读定位。

本脚本**不做**论断判定，只做定位。判读规则在运行前写定：

一、输入（两个缓存目录都扫，规则完全相同）
  external/ptm_claims_fulltext_cache_2026-09-15/*.xml       持硫化 / S-亚硝基化 / 亚磺酰化 轮次
  external/ptm_claims_sulfenylation_cache_2026-09-15/*.xml  范围纠正后补的 sulfenylation 轮次
  results/ptm_claims_literature_search.csv 与
  results/sulfenylation_claims_literature_search.csv （提供 pmcid -> family/pmid/doi/title 映射）

二、切句与分节
  用正则剥去 XML 标签，按 <sec> 标题归属到 section_title；句子切分用 (?<=[.!?])\\s+(?=[A-Z(])。
  句子长度 < 40 或 > 700 字符的丢弃（过短无信息，过长多为参考文献或表格串）。

三、入选规则（两个条件必须同时成立）
  (a) 句子命中至少一个"位点偏好关键词族"（CLAIM_KEYWORD_FAMILIES，与检索脚本同一份定义）；
  (b) 句子命中至少一个"统计量模式"（STAT_PATTERNS）：p 值、n=、百分比、倍数、AUC/AUROC、
      相关系数、富集倍数、位点/蛋白计数、平均±、置信区间。
  两条都满足 → 写入候选表，记 section_title、claim_families、stat_hits。

四、排序
  candidate_score = 命中关键词族数 * 2 + 命中统计量模式数；降序。
  节标题含 result/discussion/conclusion 的加 1 分（作者在这些节里下论断）。
  **这个分数只用于阅读顺序，不参与任何判定。**

五、产物（文件名独立）
  results/ptm_claims_candidate_sentences.csv
  results/ptm_claims_candidate_sentences_audit.json （含本脚本 sha256、句子计数、每篇命中数）

引用纪律：候选表里保存的是原文句子，**仅作内部定位用，不分发**；最终报告与论断表里
只写自己的话（20 词内）加短引语定位，不整段复制。
"""

import csv
import hashlib
import json
import os
import re
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_DIRS = [os.path.join(ROOT, "external", "ptm_claims_fulltext_cache_2026-09-15"),
              os.path.join(ROOT, "external", "ptm_claims_sulfenylation_cache_2026-09-15")]
CACHE_DIR = CACHE_DIRS[0]  # 向后兼容：其他脚本按名导入本模块的函数时仍可用
SEARCH_CSVS = [os.path.join(ROOT, "results", "ptm_claims_literature_search.csv"),
               os.path.join(ROOT, "results", "sulfenylation_claims_literature_search.csv")]
SEARCH_CSV = SEARCH_CSVS[0]
OUT_CSV = os.path.join(ROOT, "results", "ptm_claims_candidate_sentences.csv")
OUT_AUDIT = os.path.join(ROOT, "results", "ptm_claims_candidate_sentences_audit.json")

CLAIM_KEYWORD_FAMILIES = {
    "seq_motif": ["motif", "consensus sequence", "consensus motif", "flanking", "sequence context",
                  "sequence preference", "amino acid composition", "neighboring residue",
                  "neighbouring residue", "adjacent residue", "sequence logo", "linear motif"],
    "structure": ["solvent accessib", "secondary structure", "surface exposed", "surface-exposed",
                  "buried", "burial", "alpha-helix", "beta-sheet", "structural context", "plddt",
                  "intrinsically disordered", "disorder", "b-factor", "loop region", "relative accessibility"],
    "electrostatic": ["pka", "electrostatic", "acid-base", "acid-based", "acidic residue", "basic residue",
                      "charge", "nucleophilic", "nucleophilicity", "hydrophobic", "hydrophobicity"],
    "function": ["gene ontology", "go enrichment", "go term", "functional enrichment",
                 "pathway enrichment", "overrepresent", "over-represent", "enriched in",
                 "kegg", "functional categor"],
    "cooccurrence": ["co-occur", "cooccur", "crosstalk", "cross-talk", "overlap with",
                     "phosphorylation site", "other modification", "competing modification",
                     "shared site", "co-modif"],
}

STAT_PATTERNS = {
    "p_value": r"\b[Pp]\s*[=<>]\s*0?\.\d+|\b[Pp]\s*[=<>]\s*\d+(\.\d+)?[eE]-\d+|\bp-value",
    "n_equals": r"\bn\s*=\s*\d",
    "percent": r"\d+(\.\d+)?\s?%",
    "fold": r"\d+(\.\d+)?[-\s]?fold",
    "auc": r"\bAU[CR]O?C?\b|\bAUC\b",
    "corr": r"\br\s*=\s*-?0?\.\d+|\bR2\b|\bR\^?2\s*=",
    "counts": r"\b\d{2,6}\s+(sites?|cysteines?|proteins?|peptides?|residues?)\b",
    "mean_sd": r"\d+(\.\d+)?\s*[±]\s*\d+(\.\d+)?",
    "ci": r"95%\s*(confidence|CI)",
    "enrichment": r"enrichment\s+(score|factor|of)\s*[:=]?\s*\d|odds ratio",
}

SEC_BONUS = ("result", "discussion", "conclusion")
MIN_LEN, MAX_LEN = 40, 700


def sha256_of_self():
    with open(os.path.abspath(__file__), "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def load_meta():
    meta = {}
    for path in SEARCH_CSVS:
        if not os.path.exists(path):
            continue
        with open(path, newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                if row.get("pmcid"):
                    meta.setdefault(row["pmcid"], row)
    return meta


def sectionize(xml_text):
    """Yield (section_title, plain_text) chunks."""
    body = xml_text
    m = re.search(r"<body\b.*?>(.*)</body>", body, flags=re.S)
    if m:
        body = m.group(1)
    body = re.sub(r"<(table-wrap|fig|ref-list|xref|inline-formula|disp-formula)\b.*?</\1>", " ", body, flags=re.S)
    parts = re.split(r"<sec\b[^>]*>", body)
    out = []
    for part in parts:
        tm = re.search(r"<title>(.*?)</title>", part, flags=re.S)
        title = re.sub(r"<[^>]+>", " ", tm.group(1)) if tm else ""
        title = re.sub(r"\s+", " ", title).strip()[:120]
        txt = re.sub(r"<[^>]+>", " ", part)
        txt = re.sub(r"&[a-zA-Z#0-9]+;", " ", txt)
        txt = re.sub(r"\s+", " ", txt).strip()
        if txt:
            out.append((title, txt))
    return out


def sentences(text):
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(])", text)]


def main():
    meta = load_meta()
    rows, per_paper, n_sent_total = [], {}, 0
    files = []
    for cdir in CACHE_DIRS:
        if os.path.isdir(cdir):
            files += [(cdir, f) for f in sorted(os.listdir(cdir)) if f.endswith(".xml")]
    for cdir, fn in files:
        pmcid = fn[:-4]
        with open(os.path.join(cdir, fn), encoding="utf-8", errors="replace") as fh:
            xml = fh.read()
        m = meta.get(pmcid, {})
        hits = 0
        for sec_title, txt in sectionize(xml):
            for sent in sentences(txt):
                n_sent_total += 1
                if not (MIN_LEN <= len(sent) <= MAX_LEN):
                    continue
                low = sent.lower()
                fams = sorted(f for f, kws in CLAIM_KEYWORD_FAMILIES.items() if any(k in low for k in kws))
                if not fams:
                    continue
                stats = sorted(s for s, pat in STAT_PATTERNS.items() if re.search(pat, sent))
                if not stats:
                    continue
                score = 2 * len(fams) + len(stats) + (1 if any(b in sec_title.lower() for b in SEC_BONUS) else 0)
                rows.append({
                    "pmcid": pmcid, "pmid": m.get("pmid", ""), "doi": m.get("doi", ""),
                    "family": m.get("family", ""), "year": m.get("year", ""),
                    "title": m.get("title", ""), "section_title": sec_title,
                    "claim_families": ";".join(fams), "stat_hits": ";".join(stats),
                    "candidate_score": score, "sentence": sent,
                })
                hits += 1
        per_paper[pmcid] = hits
    rows.sort(key=lambda r: (-r["candidate_score"], r["pmcid"]))
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    audit = {
        "script": os.path.relpath(os.path.abspath(__file__), ROOT),
        "script_sha256": sha256_of_self(),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "inputs": {"fulltext_caches": [os.path.relpath(c, ROOT) for c in CACHE_DIRS],
                   "search_csvs": [os.path.relpath(c, ROOT) for c in SEARCH_CSVS],
                   "xml_files": len(files)},
        "rules": {"min_len": MIN_LEN, "max_len": MAX_LEN,
                  "require_claim_keyword_family": True, "require_stat_pattern": True,
                  "score": "2*n_claim_families + n_stat_patterns + 1 if section in result/discussion/conclusion",
                  "score_use": "reading order only; never a verdict"},
        "counts": {"sentences_scanned": n_sent_total, "candidate_sentences": len(rows),
                   "papers_with_candidates": sum(1 for v in per_paper.values() if v > 0),
                   "papers_scanned": len(files)},
        "candidates_per_paper_top20": sorted(per_paper.items(), key=lambda kv: -kv[1])[:20],
        "outputs": {"candidate_csv": os.path.relpath(OUT_CSV, ROOT)},
        "distribution": "候选句表含原文句子，仅内部定位用，不分发。",
    }
    with open(OUT_AUDIT, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)
    print(json.dumps(audit["counts"], ensure_ascii=False))


if __name__ == "__main__":
    sys.exit(main())
