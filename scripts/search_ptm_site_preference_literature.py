#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""检索三类 cysteine switch 化学中"位点偏好论断"的文献底库（第一阶段：只检索与筛选，不提取、不再检验）。

本脚本的全部判据在运行前写定，运行时不得改动。产物文件名独立，不覆盖任何既有产物。

===============================================================================
0. 范围
===============================================================================
三类修饰（"family"）：
  persulfidation   —— 持硫化 / S-sulfhydration / polysulfidation（Cys-SSH）
  s_nitrosylation  —— S-亚硝基化 / S-nitrosation（Cys-SNO）
  sulfinylation    —— 亚磺酰化 / cysteine sulfinic acid（Cys-SO2H）
只做"作者声称的位点偏好规律"这一类论断；不收集单点机理研究（单个 Cys 的功能报道）。

数据源：Europe PMC REST（https://www.ebi.ac.uk/europepmc/webservices/rest/search），
resultType=core，默认相关度排序。持硫化一支另外并入本树已有的
results/parallel_label_literature_screen.csv（233 条，2026-09-15 上一轮产物，只读）。

===============================================================================
1. 检索式与召回验证（本项目纪律第 6 条）
===============================================================================
每一类修饰都预先指定两个"已知阳性对照"（PC），在检索前写死在 POSITIVE_CONTROLS 里。
检索跑完后逐个 PC 检查是否被任一检索式召回，并写入审计 JSON。

  persulfidation:
    PC-A  Longen 2016, Sci Rep, qPerS-SID, PMID 27411966, doi 10.1038/srep29808
          （上一轮的教训：它标题写 "persulfide site identification"，不含 "persulfidation"）
    PC-B  Fu 2020, Antioxid Redox Signal, low-pH QTRP, PMID 31411056, doi 10.1089/ars.2019.7777
          （本项目在用的训练标签来源）
  s_nitrosylation:
    PC-A  Doulias 2010, PNAS, PMID 20837516, PMC2947911, doi 10.1073/pnas.1008036107
          （内源 SNO 位点结构画像，位点偏好论断的经典来源）
          注：本轮最初凭记忆把它的 doi 写成 10.1073/pnas.1005633107，那是同年 PNAS 的一篇钠通道
          beta4 论文（pmid 20566860），与 SNO 无关；已按标题检索改正，并加了 PC 标题自检
          （verify_pc_identities，核对不符即中止），同类错误以后不会静默通过。
    PC-B  Marino & Gladyshev 2010, J Mol Biol, PMID 19854201, doi 10.1016/j.jmb.2009.10.042
          （modified acid-base motif）
  sulfinylation:
    PC-A  Akter 2018, Nat Chem Biol, DiaAlk, PMID 30177848, doi 10.1038/s41589-018-0116-2
    PC-B  Wood, Poole & Karplus 2003, Science, doi 10.1126/science.1080405
          （GGLG / YF 结构基序决定 Prx 过氧化敏感性 —— 亚磺酰化位点偏好的经典结构论断）

判读规则：
  - PC 被至少一个检索式返回 → recalled=True。
  - 任一 PC 未被召回 → 该类修饰的检索结论只能写成"候选"，**不得写成"不存在"**；
    审计 JSON 的 recall_verified 记为 False，并在报告里明写。
  - 允许为未召回的 PC 追加一个"补救检索式"（rescue query），但必须单独标注 is_rescue=True，
    且原始式子的召回失败必须照实保留在审计里，不得抹掉。

===============================================================================
2. 抓取深度与预先声明的升级规则（超出即为"未穷尽"）
===============================================================================
本轮采用两级深度，两级都在运行前写定：
  TIER1_DEPTH = 200   一级（原始声明的上限）：每个检索式取相关度前 200 条。
  TIER2_DEPTH = 1000  二级（升级）：只要某类修饰有任一 PC 在一级深度内未被召回，
                      该类的全部检索式扩到相关度前 1000 条。
实现上只做一次抓取（抓到 TIER2_DEPTH），并给每条记录记下它在该检索式里的相关度**排名 rank**；
于是"深度 200 的召回"与"深度 1000 的召回"都能从同一次抓取里算出来，两者都写进审计。
命中数超过实际抓取数的检索式记 truncated=True。**本轮不做穷尽检索。**

第一次运行（script_sha256 见 results/ptm_claims_search_audit_pass1.json 的同名字段，
原始运行日志见本文件版本历史）在一级深度下 6 个 PC 只召回 3 个：
persulfidation PC-B、s_nitrosylation PC-A/PC-B、sulfinylation PC-B 全部落空，
13/15 个检索式在 200 条处被截断。这一事实按纪律第 6 条原样保留在审计的
recall_verification.at_depth_200 里，不得被二级深度的结果覆盖。

===============================================================================
3. 摘要级机械筛选（screen）—— 只用于排序，不用于"存在/不存在"判定
===============================================================================
对 title+abstract 做五类关键词族匹配（大小写不敏感，见 CLAIM_KEYWORD_FAMILIES）：
  seq_motif    序列基序 / 邻近残基组成
  structure    二级结构 / 溶剂可及性 / 埋藏度 / 无序度 / pLDDT
  electrostatic pKa / 静电 / 酸碱基序 / 亲核性
  function     GO / 通路 / 功能富集
  cooccurrence 与其它修饰共现 / crosstalk
规则：
  screen_score      = 命中的关键词族数（0–5）
  family_term_hit   = title+abstract 是否含该类修饰自身的词（见 FAMILY_TERMS）
  screen_pass       = (screen_score >= 1) AND family_term_hit
screen_pass 只是"进入全文阅读的候选池"，**不是**论断存在的判定；
最终"进入提取"的论文由第二阶段（全文阅读）决定，以 results/ptm_site_preference_claims.csv 为准。

补救检索式（rescue，预先声明、单独标注 is_rescue=True）：针对一级深度漏掉的三个 PC，
按"概念层面更窄"的方式各加一式，用来检验漏召回到底是深度造成的还是用词造成的。
Marino & Gladyshev 那一篇的标题写的是 "acid-**based** motif"（不是 acid-base），
因此补救式里两种拼法都放进去——这与上一轮 qPerS-SID 的教训是同一类错误。
补救检索式仍然可能失败；失败也照实记录，不得因此把任何论断写成"不存在"。

另外，三个漏召回的 PC 另按 DOI 直接入库（PC_DIRECT_INGEST），在表里标
matched_queries=pc_direct_ingest_<family>，以便它们能进入提取池。
**这一步不计入召回率**：召回率只看概念检索式，DOI 直取是补救入库，不是检索成功。

===============================================================================
4. 全文获取
===============================================================================
全文预算 FULLTEXT_FETCH_MAX = 400（预先声明）。按 (screen_score 降, 年份降) 排序后取前 400 条尝试，
外加三个 DOI 直取的 PC 无条件尝试。预算之外的 screen_pass 记录标 fulltext_skipped_budget=True，
它们只能停在摘要级，**不得**当作"读过方法"来提取。
对 screen_pass 且有 pmcid 且 isOpenAccess=Y 的记录，抓 Europe PMC fullTextXML 存入
external/ptm_claims_fulltext_cache_2026-09-15/。该缓存仅供本树内复现，**不分发**；
付费墙论文一律不抓、不存、不分发，只保留摘要级判定并在产物里标 evidence_level=abstract_only。
全文获取成功率写入审计 JSON（fulltext_attempted / fulltext_obtained）。

===============================================================================
5. 产物（文件名独立，不覆盖他人产物）
===============================================================================
  results/ptm_claims_literature_search.csv    逐论文一行（含 family、检索式命中、筛选标记）
  results/ptm_claims_search_audit.json        检索式、命中数、召回验证、上限、全文成功率、本脚本 sha256
  external/ptm_claims_fulltext_cache_2026-09-15/*.xml   OA 全文缓存（不分发）
  external/ptm_claims_fulltext_cache_2026-09-15/_search_raw.json  原始检索返回（可复现）

种子：本脚本无随机过程，不需要种子（自助法在第二阶段的再检验轮次才使用，种子 20260915）。
"""

import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_DIR = os.path.join(ROOT, "external", "ptm_claims_fulltext_cache_2026-09-15")
OUT_CSV = os.path.join(ROOT, "results", "ptm_claims_literature_search.csv")
OUT_AUDIT = os.path.join(ROOT, "results", "ptm_claims_search_audit.json")
PRIOR_SCREEN = os.path.join(ROOT, "results", "parallel_label_literature_screen.csv")

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
TIER1_DEPTH = 200
TIER2_DEPTH = 1000
FULLTEXT_FETCH_MAX = 400
PAGE_SIZE = 100
SLEEP = 0.35

POSITIVE_CONTROLS = {
    "persulfidation": [
        {"id": "PC-A", "label": "Longen 2016 qPerS-SID", "pmid": "27411966", "doi": "10.1038/srep29808",
         "title_must_contain": "persulfide site identification",
         "id_resolved_via": "NCBI ID converter; title re-verified via Europe PMC EXT_ID, 2026-09-15"},
        {"id": "PC-B", "label": "Fu 2020 low-pH QTRP", "pmid": "31411056", "doi": "10.1089/ars.2019.7777",
         "title_must_contain": "direct proteomic mapping of cysteine persulfidation",
         "id_resolved_via": "Europe PMC author+year probe; title re-verified via EXT_ID, 2026-09-15"},
    ],
    "s_nitrosylation": [
        {"id": "PC-A", "label": "Doulias 2010 PNAS SNO structural profiling", "pmid": "20837516", "doi": "10.1073/pnas.1008036107",
         "title_must_contain": "structural profiling of endogenous s-nitrosocysteine residues",
         "id_resolved_via": "**更正**：本轮最初把 PC-A 的 doi 写成 10.1073/pnas.1005633107（凭记忆），"
                            "该 doi 实际指向同年 PNAS 的一篇钠通道 beta4 论文（pmid 20566860），与 SNO 无关；"
                            "按标题检索改正为 doi 10.1073/pnas.1008036107 / pmid 20837516 / PMC2947911。"
                            "第二次运行（深度 1000）用的是错的 PC-A，其产物已作废，见 audit 的 discarded_runs。"},
        {"id": "PC-B", "label": "Marino & Gladyshev 2010 JMB acid-based motif", "pmid": "19854201", "doi": "10.1016/j.jmb.2009.10.042",
         "title_must_contain": "structural analysis of cysteine s-nitrosylation",
         "id_resolved_via": "Europe PMC title probe（凭记忆写的 doi 10.1016/j.jmb.2009.12.005 也是错的，"
                            "已改正为 .10.042）；title 再经 EXT_ID 复核，2026-09-15"},
    ],
    "sulfinylation": [
        {"id": "PC-A", "label": "Akter 2018 Nat Chem Biol DiaAlk", "pmid": "30177848", "doi": "10.1038/s41589-018-0116-2",
         "title_must_contain": "cysteine sulfinic acid reductase",
         "id_resolved_via": "Europe PMC title+author probe；title 再经 EXT_ID 复核，2026-09-15"},
        {"id": "PC-B", "label": "Wood 2003 Science Prx GGLG/YF hyperoxidation", "pmid": "12714747", "doi": "10.1126/science.1080405",
         "title_must_contain": "peroxiredoxin evolution and the regulation of hydrogen peroxide signaling",
         "id_resolved_via": "第一次探针因瞬时 TLS 错误**什么都没解析出来**（当时误称已解析，事后改正）；"
                            "改用 Europe PMC DOI+TITLE 复核得 pmid 12714747。NCBI ID converter 对该 Science DOI 不返回 PMID。"
                            "付费墙（isOpenAccess=N、hasSuppl=N）→ 只能摘要级。",
         "paywalled": True},
    ],
}

QUERIES = [
    # ---- persulfidation ----
    ("persulfidation", "P1_motif", '(persulfidation OR "S-sulfhydration" OR persulfidome OR polysulfidation) AND (motif OR "sequence context" OR "flanking residue" OR "site preference" OR "site selectivity")', False),
    ("persulfidation", "P2_site_id_synonyms", '("persulfide site identification" OR "persulfide proteome" OR "persulfidated proteome" OR sulfhydrome OR "persulfide site")', False),
    ("persulfidation", "P3_predictor", '(persulfidation OR "S-sulfhydration") AND (predictor OR "machine learning" OR "deep learning" OR "prediction of") AND cysteine', False),
    ("persulfidation", "P4_structure", '(persulfidation OR "S-sulfhydration") AND (pKa OR "solvent accessibility" OR "secondary structure" OR "intrinsically disordered" OR pLDDT OR "surface exposed")', False),
    ("persulfidation", "P5_function", '(persulfidation OR "S-sulfhydration" OR persulfidome) AND ("gene ontology" OR "functional enrichment" OR "pathway enrichment" OR "enriched in")', False),
    # ---- S-nitrosylation ----
    ("s_nitrosylation", "N1_motif", '("S-nitrosylation" OR "S-nitrosation" OR "S-nitrosoproteome") AND (motif OR "sequence context" OR "flanking residue" OR "consensus motif")', False),
    ("s_nitrosylation", "N2_site_proteomics", '("S-nitrosoproteome" OR "S-nitrosylation site" OR "S-nitrosylated site" OR SNOSID OR "SNO-RAC") AND (proteome OR proteomic OR "site-specific")', False),
    ("s_nitrosylation", "N3_structure", '("S-nitrosylation" OR "S-nitrosation") AND ("acid-base motif" OR "solvent accessibility" OR pKa OR hydrophobicity OR "structural determinants" OR "structural analysis")', False),
    ("s_nitrosylation", "N4_predictor", '("S-nitrosylation" OR "S-nitrosation") AND (predictor OR "machine learning" OR "deep learning" OR "prediction of") AND site', False),
    ("s_nitrosylation", "N5_function", '("S-nitrosylation" OR "S-nitrosoproteome") AND ("gene ontology" OR "functional enrichment" OR "pathway enrichment") AND (site OR cysteine)', False),
    # ---- sulfinylation ----
    ("sulfinylation", "S1_chemoproteomics", '(sulfinylation OR "S-sulfinylation" OR "sulfinic acid") AND cysteine AND (chemoproteomic OR proteomic OR "site-specific" OR "site-centric")', False),
    ("sulfinylation", "S2_motif", '("cysteine sulfinic acid" OR "sulfinic acid" OR sulfinylation) AND (motif OR "sequence context" OR "structural determinants" OR preference OR selectivity)', False),
    ("sulfinylation", "S3_hyperoxidation", '(hyperoxidation OR overoxidation OR "over-oxidation") AND (peroxiredoxin OR cysteine) AND (motif OR GGLG OR sensitivity OR structural)', False),
    ("sulfinylation", "S4_probe_sulfiredoxin", '(DiaAlk OR "sulfinic acid reductase" OR sulfiredoxin OR "BTD probe") AND (proteomic OR proteome OR chemoproteomic OR "site-specific")', False),
    ("sulfinylation", "S5_predictor", '(sulfinylation OR "cysteine sulfinic acid") AND (predictor OR "machine learning" OR "deep learning" OR prediction)', False),
    # ---- rescue queries (declared after pass-1 recall failure at depth 200; is_rescue=True) ----
    ("persulfidation", "R1_pers_site_mapping", 'TITLE:(persulfidation OR persulfide OR "S-sulfhydration") AND ("proteomic mapping" OR "direct proteomic" OR "site-level" OR "quantitative site" OR "mapping of cysteine")', True),
    ("s_nitrosylation", "R2_sno_structural", 'TITLE:("S-nitrosocysteine" OR "S-nitrosylation" OR "S-nitrosation") AND ("structural profiling" OR "structural analysis" OR "acid-based motif" OR "acid-base motif" OR "structural feature" OR "endogenous")', True),
    ("sulfinylation", "R3_prx_hyperoxidation_motif", '(peroxiredoxin OR Prx) AND (GGLG OR "YF motif" OR "C-terminal helix") AND (hyperoxidation OR "sensitivity to" OR inactivation OR "peroxide signaling")', True),
]

# 三个在一级深度漏召回的 PC 按 DOI 直接入库（不计入召回率，仅为让它们进入提取池）
PC_DIRECT_INGEST = [
    ("persulfidation", "10.1089/ars.2019.7777"),
    ("s_nitrosylation", "10.1073/pnas.1008036107"),
    ("s_nitrosylation", "10.1016/j.jmb.2009.10.042"),
    ("sulfinylation", "10.1126/science.1080405"),
]

FAMILY_TERMS = {
    "persulfidation": ["persulfid", "sulfhydrat", "polysulfid", "sulfhydrome", "ssh", "hydrogen sulfide", "h2s"],
    "s_nitrosylation": ["nitrosyl", "nitrosat", "nitroso", "sno"],
    "sulfinylation": ["sulfinyl", "sulfinic", "so2h", "hyperoxid", "overoxid", "over-oxid", "sulfiredoxin"],
}

CLAIM_KEYWORD_FAMILIES = {
    "seq_motif": ["motif", "consensus sequence", "consensus motif", "flanking", "sequence context",
                  "sequence preference", "amino acid composition", "neighboring residue",
                  "neighbouring residue", "adjacent residue", "sequence logo", "linear motif"],
    "structure": ["solvent accessib", "secondary structure", "surface exposed", "surface-exposed",
                  "buried", "burial", "alpha-helix", "beta-sheet", "structural context", "plddt",
                  "intrinsically disordered", "disorder", "b-factor", "loop region", "relative accessibility"],
    "electrostatic": ["pka", "electrostatic", "acid-base", "acidic residue", "basic residue",
                      "charge", "nucleophilic", "nucleophilicity", "hydrophobic", "hydrophobicity"],
    "function": ["gene ontology", "go enrichment", "go term", "functional enrichment",
                 "pathway enrichment", "overrepresent", "over-represent", "enriched in",
                 "kegg", "functional categor"],
    "cooccurrence": ["co-occur", "cooccur", "crosstalk", "cross-talk", "overlap with",
                     "phosphorylation site", "other modification", "competing modification",
                     "shared site", "co-modif"],
}


def sha256_of_self():
    with open(os.path.abspath(__file__), "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def get_json(url, tries=4):
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return json.load(fh)
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 + 3 * i)
    raise RuntimeError("GET failed: %s -> %s" % (url[:120], last))


def epmc_search(query, max_records=TIER2_DEPTH):
    """Return (hit_count, [records]) with ``_rank`` (1-based relevance rank) stamped on each record."""
    got, cursor, hit = [], "*", None
    while len(got) < max_records:
        url = (EPMC + "/search?query=" + urllib.parse.quote(query)
               + "&format=json&resultType=core&pageSize=%d&cursorMark=%s" % (PAGE_SIZE, urllib.parse.quote(cursor)))
        data = get_json(url)
        if hit is None:
            hit = data.get("hitCount", 0)
        batch = data.get("resultList", {}).get("result", [])
        if not batch:
            break
        for rec in batch:
            got.append(rec)
            rec["_rank"] = len(got)
        nxt = data.get("nextCursorMark")
        if not nxt or nxt == cursor:
            break
        cursor = nxt
        time.sleep(SLEEP)
    return hit or 0, got[:max_records]


def epmc_by_doi(doi):
    hit, recs = epmc_search('DOI:"%s"' % doi, max_records=PAGE_SIZE)
    return recs[0] if recs else None


def norm_text(rec):
    return ((rec.get("title") or "") + " " + (rec.get("abstractText") or "")).lower()


def screen(rec, family):
    txt = norm_text(rec)
    flags = sorted(fam for fam, kws in CLAIM_KEYWORD_FAMILIES.items() if any(k in txt for k in kws))
    fam_hit = any(t in txt for t in FAMILY_TERMS[family])
    return flags, len(flags), fam_hit, (len(flags) >= 1 and fam_hit)


def key_of(rec):
    for f in ("pmid", "doi", "pmcid", "id"):
        v = rec.get(f)
        if v:
            return f + ":" + str(v).lower()
    return "unknown:" + str(id(rec))


def pc_matches(rec, pc):
    return ((pc.get("pmid") and str(rec.get("pmid")) == pc["pmid"])
            or (pc.get("doi") and str(rec.get("doi", "")).lower() == pc["doi"].lower()))


def verify_pc_identities():
    """按 DOI 取回每个 PC 的记录，核对标题里必须出现 title_must_contain；不符即中止。

    这一步是本轮踩过的坑的硬化：PC 的 doi/pmid 若凭记忆写错，召回率会变成对着一篇
    无关论文测，而且不会报错。核对失败必须中止，不得降级继续。
    """
    checked = []
    for family, pcs in POSITIVE_CONTROLS.items():
        for pc in pcs:
            rec = epmc_by_doi(pc["doi"])
            title = (rec or {}).get("title", "") or ""
            need = pc["title_must_contain"].lower()
            ok = need in title.lower()
            checked.append({"family": family, "id": pc["id"], "doi": pc["doi"],
                            "pmid_declared": pc["pmid"], "pmid_returned": (rec or {}).get("pmid"),
                            "title_returned": re.sub(r"\s+", " ", title).strip(),
                            "title_must_contain": pc["title_must_contain"], "identity_ok": ok})
            if not ok:
                raise SystemExit("PC identity check FAILED: %s %s doi=%s -> title=%r (needs %r)"
                                 % (family, pc["id"], pc["doi"], title, pc["title_must_contain"]))
            if (rec or {}).get("pmid") and str(rec["pmid"]) != str(pc["pmid"]):
                raise SystemExit("PC pmid mismatch: %s %s declared=%s returned=%s"
                                 % (family, pc["id"], pc["pmid"], rec["pmid"]))
            time.sleep(SLEEP)
    return checked


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    pc_identity = verify_pc_identities()
    records, query_audit, raw_dump = {}, [], {}

    for family, qid, query, is_rescue in QUERIES:
        hit, batch = epmc_search(query)
        query_audit.append({
            "family": family, "query_id": qid, "query": query, "is_rescue": bool(is_rescue),
            "hit_count": hit, "retrieved": len(batch),
            "truncated_at_tier1_depth200": bool(hit > TIER1_DEPTH),
            "truncated_at_tier2": bool(hit > len(batch)),
            "tier1_depth": TIER1_DEPTH, "tier2_depth": TIER2_DEPTH,
        })
        raw_dump[qid] = [{k: r.get(k) for k in ("pmid", "pmcid", "doi", "title", "pubYear",
                                                "isOpenAccess", "_rank")} for r in batch]
        for rec in batch:
            k = key_of(rec)
            if k not in records:
                records[k] = {"rec": rec, "families": set(), "queries": set(), "best_rank": {}}
            records[k]["families"].add(family)
            records[k]["queries"].add(qid)
            prev = records[k]["best_rank"].get(qid)
            records[k]["best_rank"][qid] = min(prev, rec["_rank"]) if prev else rec["_rank"]
        time.sleep(SLEEP)

    # ---- recall verification at BOTH depths (discipline rule 6) ----
    def recall_block(depth, concept_only=True):
        block = {}
        for family, pcs in POSITIVE_CONTROLS.items():
            entries = []
            for pc in pcs:
                hits = []
                for v in records.values():
                    if not pc_matches(v["rec"], pc):
                        continue
                    for qid, rank in v["best_rank"].items():
                        is_rescue = next(q[3] for q in QUERIES if q[1] == qid)
                        if concept_only and is_rescue:
                            continue
                        if rank <= depth:
                            hits.append({"query_id": qid, "rank": rank, "is_rescue": is_rescue})
                entries.append({**{k: pc[k] for k in ("id", "label", "pmid", "doi")},
                                "recalled": bool(hits), "recalled_by": sorted(hits, key=lambda h: h["rank"])})
            block[family] = {"controls": entries,
                             "recall_verified": all(e["recalled"] for e in entries),
                             "n_recalled": sum(1 for e in entries if e["recalled"]),
                             "n_controls": len(entries)}
        return block

    recall_t1 = recall_block(TIER1_DEPTH, concept_only=True)
    recall_t2 = recall_block(TIER2_DEPTH, concept_only=True)
    recall_t2_rescue = recall_block(TIER2_DEPTH, concept_only=False)

    # ---- direct ingest of PCs still missing (NOT counted as recall) ----
    direct_ingested = []
    for family, doi in PC_DIRECT_INGEST:
        already = any(str(v["rec"].get("doi", "")).lower() == doi.lower() for v in records.values())
        if already:
            direct_ingested.append({"family": family, "doi": doi, "ingested": False,
                                    "reason": "already retrieved by a query"})
            continue
        rec = epmc_by_doi(doi)
        if rec is None:
            direct_ingested.append({"family": family, "doi": doi, "ingested": False,
                                    "reason": "Europe PMC DOI lookup returned nothing"})
            continue
        k = key_of(rec)
        records.setdefault(k, {"rec": rec, "families": set(), "queries": set(), "best_rank": {}})
        records[k]["families"].add(family)
        records[k]["queries"].add("pc_direct_ingest_" + family)
        direct_ingested.append({"family": family, "doi": doi, "ingested": True,
                                "pmid": rec.get("pmid"), "pmcid": rec.get("pmcid"),
                                "is_open_access": rec.get("isOpenAccess")})
        time.sleep(SLEEP)

    seen_pmids = {str(v["rec"].get("pmid")) for v in records.values() if v["rec"].get("pmid")}
    seen_dois = {str(v["rec"].get("doi", "")).lower() for v in records.values() if v["rec"].get("doi")}

    # ---- merge prior persulfidation screen (read-only) ----
    prior_n, prior_new = 0, 0
    if os.path.exists(PRIOR_SCREEN):
        with open(PRIOR_SCREEN, newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                prior_n += 1
                pmid = (row.get("pmid") or "").strip()
                doi = (row.get("doi") or "").strip().lower()
                if (pmid and pmid in seen_pmids) or (doi and doi in seen_dois):
                    continue
                prior_new += 1
                k = "prior:" + (pmid or doi or row.get("title", "")[:40])
                records[k] = {
                    "rec": {"pmid": pmid, "doi": doi, "pmcid": (row.get("pmcid") or "").strip(),
                            "title": row.get("title", ""), "pubYear": row.get("year", ""),
                            "abstractText": "", "isOpenAccess": row.get("open_access", ""),
                            "journalInfo": {"journal": {"title": row.get("journal", "")}}},
                    "families": {"persulfidation"}, "queries": {"prior_round_233"}, "best_rank": {},
                }

    # ---- build rows, then spend the full-text budget ----
    rows = []
    for k, val in sorted(records.items()):
        rec, fams = val["rec"], sorted(val["families"])
        fam_primary = fams[0]
        flags, score, fam_hit, spass = screen(rec, fam_primary)
        qids = sorted(val["queries"])
        ranks = val["best_rank"]
        in_tier1 = any(r <= TIER1_DEPTH for r in ranks.values()) if ranks else False
        rows.append({
            "record_key": k, "family": fam_primary, "families_all": ";".join(fams),
            "pmid": rec.get("pmid") or "", "pmcid": (rec.get("pmcid") or "").strip(),
            "doi": rec.get("doi") or "", "year": rec.get("pubYear") or "",
            "journal": (rec.get("journalInfo") or {}).get("journal", {}).get("title", "") or "",
            "title": re.sub(r"\s+", " ", rec.get("title") or "").strip(),
            "is_open_access": "Y" if str(rec.get("isOpenAccess") or "").upper().startswith("Y") else "N",
            "matched_queries": ";".join(qids),
            "best_rank_json": json.dumps(ranks, sort_keys=True),
            "in_tier1_depth200": in_tier1,
            "is_rescue_only": bool(qids) and all(q.startswith("R") or q.startswith("pc_direct") for q in qids),
            "screen_flags": ";".join(flags), "screen_score": score,
            "family_term_hit": fam_hit, "screen_pass": spass,
            "fulltext_cached": False, "fulltext_path": "", "fulltext_skipped_budget": False,
        })

    def fetch_priority(r):
        try:
            yr = int(str(r["year"])[:4])
        except Exception:  # noqa: BLE001
            yr = 0
        return (-r["screen_score"], -yr)

    force = [r for r in rows if any(q.startswith("pc_direct") for q in r["matched_queries"].split(";"))
             or any(pc_matches({"pmid": r["pmid"], "doi": r["doi"]}, pc)
                    for pcs in POSITIVE_CONTROLS.values() for pc in pcs)]
    cands = sorted([r for r in rows if r["screen_pass"]], key=fetch_priority)
    budget_set = {id(r) for r in cands[:FULLTEXT_FETCH_MAX]} | {id(r) for r in force}
    attempted = obtained = 0
    for r in rows:
        eligible = r["pmcid"] and r["is_open_access"] == "Y" and (r["screen_pass"] or id(r) in {id(x) for x in force})
        if not eligible:
            continue
        if id(r) not in budget_set:
            r["fulltext_skipped_budget"] = True
            continue
        attempted += 1
        dest = os.path.join(CACHE_DIR, r["pmcid"] + ".xml")
        if os.path.exists(dest) and os.path.getsize(dest) > 2000:
            r["fulltext_cached"], r["fulltext_path"] = True, os.path.relpath(dest, ROOT)
        else:
            try:
                with urllib.request.urlopen(EPMC + "/" + r["pmcid"] + "/fullTextXML", timeout=90) as fh:
                    body = fh.read()
                if len(body) > 2000:
                    with open(dest, "wb") as out:
                        out.write(body)
                    r["fulltext_cached"], r["fulltext_path"] = True, os.path.relpath(dest, ROOT)
            except Exception:  # noqa: BLE001
                pass
            time.sleep(SLEEP)
        if r["fulltext_cached"]:
            obtained += 1

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(CACHE_DIR, "_search_raw.json"), "w", encoding="utf-8") as fh:
        json.dump(raw_dump, fh, ensure_ascii=False)

    n_pass = sum(1 for r in rows if r["screen_pass"])
    audit = {
        "script": os.path.relpath(os.path.abspath(__file__), ROOT),
        "script_sha256": sha256_of_self(),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "stage": "phase1_search_and_screen_only__no_retest",
        "source": "Europe PMC REST resultType=core, relevance order",
        "depths": {"tier1": TIER1_DEPTH, "tier2": TIER2_DEPTH,
                   "escalation_rule": "any PC missed at tier1 -> that family's queries go to tier2; "
                                      "implemented as one retrieval to tier2 with per-record relevance rank"},
        "fulltext_budget": FULLTEXT_FETCH_MAX,
        "exhaustive": False,
        "queries": query_audit,
        "positive_control_identity_check": pc_identity,
        "discarded_runs": [
            {"run": "pass2_depth1000_first_attempt",
             "why_discarded": "该次运行用的 s_nitrosylation PC-A 标识符是错的"
                              "（doi 10.1073/pnas.1005633107 / pmid 20566860 实为同年 PNAS 一篇钠通道 beta4 论文），"
                              "召回率相当于对着一篇无关论文测，故作废、不出表；"
                              "其填好的 OA 全文缓存被本次运行复用。",
             "outputs_kept": "none (只保留 external/ptm_claims_fulltext_cache_2026-09-15 里的全文缓存)"},
        ],
        "recall_verification": {
            "at_depth_200_concept_queries_only": recall_t1,
            "at_depth_1000_concept_queries_only": recall_t2,
            "at_depth_1000_including_rescue_queries": recall_t2_rescue,
            "pc_direct_ingest_not_counted_as_recall": direct_ingested,
        },
        "recall_summary": {
            "depth200_concept": {f: [v["n_recalled"], v["n_controls"]] for f, v in recall_t1.items()},
            "depth1000_concept": {f: [v["n_recalled"], v["n_controls"]] for f, v in recall_t2.items()},
            "depth1000_with_rescue": {f: [v["n_recalled"], v["n_controls"]] for f, v in recall_t2_rescue.items()},
            "all_verified_depth200": all(v["recall_verified"] for v in recall_t1.values()),
            "all_verified_depth1000_with_rescue": all(v["recall_verified"] for v in recall_t2_rescue.values()),
        },
        "pass1_record": {
            "audit": "results/ptm_claims_search_audit_pass1.json",
            "table": "results/ptm_claims_literature_search_pass1.csv",
            "note": "第一次运行（深度 200、无补救式）6 个 PC 只召回 3 个，原样保留，不得覆盖。",
        },
        "prior_round_table": {"path": os.path.relpath(PRIOR_SCREEN, ROOT), "rows": prior_n,
                              "added_as_new_records": prior_new},
        "counts": {
            "unique_records": len(rows),
            "unique_records_in_tier1_depth200": sum(1 for r in rows if r["in_tier1_depth200"]),
            "by_family": {f: sum(1 for r in rows if r["family"] == f) for f in FAMILY_TERMS},
            "screen_pass": n_pass,
            "screen_pass_by_family": {f: sum(1 for r in rows if r["family"] == f and r["screen_pass"]) for f in FAMILY_TERMS},
            "fulltext_attempted": attempted,
            "fulltext_obtained": obtained,
            "fulltext_success_rate": round(obtained / attempted, 4) if attempted else None,
            "fulltext_skipped_budget": sum(1 for r in rows if r["fulltext_skipped_budget"]),
        },
        "outputs": {
            "search_csv": os.path.relpath(OUT_CSV, ROOT),
            "fulltext_cache_dir": os.path.relpath(CACHE_DIR, ROOT),
            "cache_distribution": "internal reproducibility only; NOT for redistribution",
        },
        "notes": [
            "screen_pass 是候选池标记，不是论断存在性判定；最终进入提取的论文以 results/ptm_site_preference_claims.csv 为准。",
            "付费墙论文不抓全文、不缓存、不分发，只保留摘要级判定。",
            "pc_direct_ingest 不计入召回率。",
        ],
    }
    with open(OUT_AUDIT, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)

    print(json.dumps({"unique_records": len(rows), "screen_pass": n_pass,
                      "fulltext": [attempted, obtained],
                      "recall": audit["recall_summary"]}, ensure_ascii=False))


if __name__ == "__main__":
    sys.exit(main())
