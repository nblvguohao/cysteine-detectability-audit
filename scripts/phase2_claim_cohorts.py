#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 A 层八条已发表论断各自装配成一张可再检验的表，一条论断一个函数。

本模块**只装配数据，不做检验**，检验在
`scripts/run_phase2_claims_under_detectability_control.py`。分成两个文件是因为
"按论文自己的口径"这件事全部落在装配上：每条论断的阳性集、阴性集（或背景集）、
被主张的属性，都必须照论文自述取，任何偏离都写在该函数的 docstring 与
返回结构的 `reproduction_notes` 里，由审计 JSON 逐条带出。

装配后的统一结构（一条论断一个 dict）
-------------------------------------
unit          protein | site          论断自身的统计单位
keys          list                    蛋白登录号，或 (登录号, 位点)
y             0/1 向量                论文的阳性指示（被检出/属某簇/四器官共有）
attribute     0/1 向量                被主张的属性（属某 GO 类/侧翼带 K）
groups        字符串向量              自助抽样聚类单位（见下）
covariates    n×d 实数矩阵            可检出性协变量
covariate_set VIS10 | DIG25 | PROT8   协变量口径名
observed_depth 实数向量或 None        论文自己的检出深度（肽段数/谱图计数），仅阳性侧可得时为 None

聚类单位（运行前写定）
--------------------
* **蛋白层论断（PERS-005/006/010/011/012）：单位就是蛋白**，因此聚类自助按蛋白
  重抽。项目纪律允许"按蛋白或同源组件"，这里取前者，因为一个蛋白就是一个观测，
  同源性带来的依赖只在"同一族多个旁系同源都进正集"时抬高精度；该限制写进报告，
  不假装已经控制。
* **位点层论断（SFE-006 转移检验、SNO-021）：单位是位点，聚类按蛋白整体重抽**，
  并另跑一版按同源组件（k=7 含有率≥0.40，与
  `scripts/build_v3_human_homology_components.py` 同一构造）作为敏感性；
  该构造只在这两条论断自己的蛋白集上跑（637 与数百个蛋白），不在整个参考蛋白组上跑。

可检出性协变量（运行前写定，定义**直接从探针脚本 import**，不重写）
------------------------------------------------------------------
* **VIS10** = `run_cross_protease_detectability_probe.VISIBILITY_ONLY` 的十个纯可见性
  特征，胰酶切割规则，位点层。
* **DIG25** = 同脚本 `FEATURE_NAMES` 全部二十五个消化特征，位点层。
* **PROT8** = 蛋白层聚合，本模块新声明（探针脚本没有蛋白层口径），八维：
  1 序列长度；2 log1p(胰酶肽段数)；3 半胱氨酸数；
  4 落在可检出窗口内的半胱氨酸数（VIS10 的 pep_detectable_both）；
  5 该比例；6 蛋白内可检出肽段数；7 含半胱氨酸肽段的平均 GRAVY；
  8 允许漏切后仍可检出的半胱氨酸比例。
  **聚合口径先写死再跑**：4/5/8 只对半胱氨酸位点聚合，1/2/6 是全蛋白量。

术语：sulfenylation = 次磺酰化（-SOH），sulfinylation = 亚磺酰化（-SO2H）。

用法：被 run_phase2_claims_under_detectability_control.py import；
单独运行时打印每条论断装配后的规模，供人工核对。
"""
from __future__ import annotations

import collections
import csv
import gzip
import os
import re

import json
import urllib.parse
import urllib.request

import numpy as np
import openpyxl

from run_cross_protease_detectability_probe import (
    FEATURE_NAMES, PROTEASES, VISIBILITY_ONLY, boundaries_for, features,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUPP = os.path.join(ROOT, "external", "intake", "phase2_claim_supp", "unpacked")
PROTEOMES = os.path.join(ROOT, "external", "proteomes")
RESULTS = os.path.join(ROOT, "results")

TRYPSIN = PROTEASES["Trypsin"]
VIS10_INDEX = [i for i, name in enumerate(FEATURE_NAMES) if name in VISIBILITY_ONLY]
DIG25_INDEX = list(range(len(FEATURE_NAMES)))
PROT8_NAMES = ["length", "log1p_n_peptides", "n_cys", "n_cys_detectable",
               "frac_cys_detectable", "n_detectable_peptides", "mean_cys_pep_gravy",
               "frac_cys_detectable_missed"]

# 被主张属性的判读，全部运行前写定 -------------------------------------------
# PERS-005 作者原文：参考拟南芥库里 GO:0006096 注释 57 个蛋白，30 个（52.6%）被检出
GO_GLYCOLYSIS = "GO:0006096"
# PERS-006 作者用 MapMan 本体的"初级代谢"（TCA、糖酵解、卡尔文循环）。
# 本环境没有 MapMan/Mercator4，**改用 GO 等价集并声明为替代口径**，不当作同一本体。
GO_PRIMARY_METABOLISM = ["GO:0006096", "GO:0006099", "GO:0019253", "GO:0015979", "GO:0006094"]
# PERS-011 作者原文：cluster 1 在 thylakoid 12.57 倍、photosynthetic 7.55 倍富集
GO_THYLAKOID = "GO:0009579"
GO_PHOTOSYNTHESIS = "GO:0015979"
# PERS-012 作者 Table 1 首行：GO:0019752 羧酸代谢，8.6% 对 2.3%，FDR 2.92e-16
GO_CARBOXYLIC = "GO:0019752"
MATAMOROS_TABLE1 = [("GO:0019752", 8.6, 2.3, 2.92e-16), ("GO:0006520", 7.0, 1.5, 2.92e-16),
                    ("GO:0006412", 8.9, 3.1, 9.10e-12), ("GO:0055086", 4.7, 1.2, 1.74e-9),
                    ("GO:0006096", 1.9, 0.3, 8.43e-6)]
# SNO-021 作者七个 Motif-X 基序里五个以 K 为核心：CXXXXXXXXXK(+10)、KXXXXXXXXXC(-10)、
# CXXXXXK(+6)、CXXXXK(+5)、CK(+1)；另两个以 I 为核心：IXXXXXXXXXC(-10)、CXI(+2)
SNO021_K_OFFSETS = [10, -10, 6, 5, 1]
SNO021_I_OFFSETS = [-10, 2]
# SFE-006 作者称最显著特征是 K/R 在 −10、−8~−6、−4、−2 与 +4~+8 位富集
SFE006_KR_OFFSETS = [-10, -8, -7, -6, -4, -2, 4, 5, 6, 7, 8]

K_MER = 7
MAX_POSTING = 50
COMPONENT_THRESHOLD = 0.40

# UniProt 登录号的官方正则（两种前缀形态，后者可带 4 位延长段）。
# 第一版只写了 [A-NR-Z] 那一支，把 O/P/Q 开头的登录号全丢了——
# 小鼠队列因此只剩 27 个蛋白，已修正，过程写在报告的更正节。
ACCESSION_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
                          r"[A-NR-Z][0-9](?:[A-Z0-9]{3}[0-9]){1,2})$")

# GO 富集工具（agriGO / DAVID / g:Profiler）一律把注释沿 GO 有向无环图上推到祖先项，
# UniProt 的 go_id 字段只给**直接注释**。因此属性判读必须"带任一后代项即算"。
# 后代集由 EBI QuickGO 取回（沙箱已放行 ebi.ac.uk）并缓存到
# external/intake/phase2_claim_supp/phase2_go_descendants.json。
# 这一条是第一次运行后**发现的缺陷修正**，不是按结果调口径：
# 直接注释下 GO:0019752 在菜豆参考蛋白组里一个都没有，属性列全零，
# 基线优势比 5.0087 完全是 Haldane 修正的产物。修正前后的值都记在报告里。
QUICKGO = ("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s/descendants"
           "?relations=is_a,part_of")
GO_CACHE = os.path.join(ROOT, "external", "intake", "phase2_claim_supp",
                        "phase2_go_descendants.json")


# ---------------------------------------------------------------- 基础读取

def go_with_descendants(term_ids):
    """{term: {term} | 其全部 is_a/part_of 后代}，取自 QuickGO，结果缓存到盘上。"""
    term_ids = sorted(set(term_ids))
    cache = {}
    if os.path.exists(GO_CACHE):
        cache = json.load(open(GO_CACHE, encoding="utf-8"))
    missing = [t for t in term_ids if t not in cache]
    for term in missing:
        url = QUICKGO % urllib.parse.quote(term)
        request = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(request, timeout=120) as handle:
            payload = json.loads(handle.read().decode("utf-8"))
        found = {term}
        for result in payload.get("results", []):
            found.update(result.get("descendants") or [])
        cache[term] = sorted(found)
    if missing:
        os.makedirs(os.path.dirname(GO_CACHE), exist_ok=True)
        with open(GO_CACHE, "w", encoding="utf-8") as handle:
            json.dump(cache, handle, indent=2)
    return {term: set(cache[term]) for term in term_ids}


def read_proteome_tsv(key):
    """返回 {accession: {"seq","go","gene","genes","length"}}，来自 phase2_*_uniprot.tsv.gz。"""
    path = os.path.join(PROTEOMES, "phase2_%s_uniprot.tsv.gz" % key)
    out = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = (row.get("Entry") or "").strip()
            if not accession:
                continue
            out[accession] = {
                "seq": (row.get("Sequence") or "").strip().upper(),
                "go": {g.strip() for g in (row.get("Gene Ontology IDs") or "").split(";") if g.strip()},
                "gene": (row.get("Gene Names (primary)") or "").strip(),
                "genes": (row.get("Gene Names") or "").strip(),
            }
    return out


def read_ath_tair_map():
    """AGI -> UniProt accession（多对多，取全部）。"""
    path = os.path.join(PROTEOMES, "ath_uniprot_tair_map.tsv.gz")
    agi_to_acc = collections.defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline()
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            accession = parts[0].strip()
            for agi in parts[1].split(";"):
                agi = agi.strip().upper()
                if agi:
                    agi_to_acc[agi].add(accession)
    return agi_to_acc


def ath_tair_sequences():
    """普查轮 `ath_tair` 视图的**逐字复刻**：TAIR 基因座 -> 最长的 UniProt 序列，
    取自 ath_uniprot_tair_map.tsv.gz 自己的 Sequence 列。

    必须用同一个视图，因为 `results/ptm_census_sites.csv` 里 fps2020_ath_sulfenyl 的
    位点编号就是按它编的；换成 UP000006548 的序列会让位点对不上残基。
    构造与 `scripts/ingest_ptm_census_sites.py` 的 load_proteome("ath_tair") 相同。
    """
    path = os.path.join(PROTEOMES, "ath_uniprot_tair_map.tsv.gz")
    sequences = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequence = (row.get("Sequence") or "").strip().upper()
            for locus in (row.get("TAIR") or "").split(";"):
                locus = locus.strip().upper()
                if not locus or not sequence:
                    continue
                if len(sequence) > len(sequences.get(locus, "")):
                    sequences[locus] = sequence
    return sequences


def sheet_rows(path, sheet=None):
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    worksheet = book[sheet] if sheet else book.worksheets[0]
    rows = [list(r) for r in worksheet.iter_rows(values_only=True)]
    book.close()
    return rows


# ---------------------------------------------------------------- 协变量

def site_feature_matrix(sequences, keys, subset):
    """对 (accession, position) 列表算位点层特征；position 为 1-based。"""
    index = VIS10_INDEX if subset == "VIS10" else DIG25_INDEX
    bounds_cache = {}
    rows = []
    for accession, position in keys:
        sequence = sequences[accession]
        if accession not in bounds_cache:
            bounds_cache[accession] = boundaries_for(sequence, TRYPSIN)
        row = features(sequence, position - 1, bounds_cache[accession], TRYPSIN)
        rows.append([row[i] for i in index])
    return np.asarray(rows, dtype=float)


def protein_feature_matrix(sequences, accessions):
    """PROT8：对每个蛋白的半胱氨酸位点聚合 VIS10/DIG25，口径见模块 docstring。"""
    name_index = {name: i for i, name in enumerate(FEATURE_NAMES)}
    rows = []
    for accession in accessions:
        sequence = sequences[accession]
        bounds = boundaries_for(sequence, TRYPSIN)
        n_peptides = len(bounds) - 1
        cys = [i for i, residue in enumerate(sequence) if residue == "C"]
        detectable, gravy, missed = [], [], []
        for i in cys:
            row = features(sequence, i, bounds, TRYPSIN)
            detectable.append(row[name_index["pep_detectable_both"]])
            gravy.append(row[name_index["pep_gravy"]])
            missed.append(row[name_index["pep_detectable_any_missed_cleavage"]])
        detectable_peptides = 0
        for k in range(n_peptides):
            length = bounds[k + 1] - bounds[k]
            if 7 <= length <= 30:
                detectable_peptides += 1
        rows.append([
            float(len(sequence)),
            float(np.log1p(n_peptides)),
            float(len(cys)),
            float(sum(detectable)),
            float(np.mean(detectable)) if cys else 0.0,
            float(detectable_peptides),
            float(np.mean(gravy)) if cys else 0.0,
            float(np.mean(missed)) if cys else 0.0,
        ])
    return np.asarray(rows, dtype=float)


def homology_components(sequences, accessions, threshold=COMPONENT_THRESHOLD):
    """k=7 含有率成分，与 build_v3_human_homology_components.py 同一构造。"""
    present = [a for a in accessions if sequences.get(a)]
    kmers = {a: {sequences[a][i:i + K_MER] for i in range(len(sequences[a]) - K_MER + 1)}
             for a in present}
    posting = collections.defaultdict(list)
    for accession, block in kmers.items():
        for kmer in block:
            posting[kmer].append(accession)
    shared = collections.Counter()
    for holders in posting.values():
        if not 2 <= len(holders) <= MAX_POSTING:
            continue
        for i in range(len(holders)):
            for j in range(i + 1, len(holders)):
                shared[(holders[i], holders[j])] += 1
    parent = {a: a for a in present}

    def find(item):
        root = item
        while parent[root] != root:
            root = parent[root]
        while parent[item] != root:
            parent[item], item = root, parent[item]
        return root

    for (left, right), count in shared.items():
        smaller = min(len(kmers[left]), len(kmers[right]))
        if smaller and count / smaller >= threshold:
            a, b = find(left), find(right)
            if a != b:
                parent[max(a, b)] = min(a, b)
    return {a: find(a) for a in present}


# ---------------------------------------------------------------- 各条论断

def ath_locus_view():
    """**以 TAIR 基因座（AGI）为单位**的拟南芥视图，而不是以 UniProt 登录号为单位。

    三条拟南芥论断（PERS-005/006/011）的作者单位都是基因座——PERS-005 原文写
    "57 proteins annotated in the reference A. thaliana database"，指的是 TAIR 基因模型，
    PERS-012 原文写 "% of genes in each GO category"。按登录号装配会因为一个基因座
    对应多条 UniProt 记录（同工型、未审阅条目）把分子和分母同时放大：
    实测 UP000006548 有 39272 条记录，而基因座只有约两万七千个。
    因此本视图做三件事，口径运行前写定：

    * 基因座集 = 出现在 UP000006548 里的全部 AGI 基因座
    * 基因座的 GO = 其名下所有登录号 GO 的并集
    * 基因座的序列 = 其名下最长的那条登录号序列（可检出性协变量的输入）
    """
    sequences = read_proteome_tsv("ath")
    agi_map = read_ath_tair_map()
    per_locus = collections.defaultdict(lambda: {"go": set(), "seq": "", "accessions": [],
                                                 "chosen": ""})
    # **必须按登录号排序后再遍历**：agi_map 的值是 set，集合的遍历顺序在不同进程里
    # 会因字符串哈希随机化而变，长度相同的同工型就会被随机挑中一条，
    # 于是协变量、倾向得分、匹配结果整条链都不可复现。
    # 并列规则写定：先取更长的序列，长度相同时取字典序最小的登录号。
    for agi in sorted(agi_map):
        for accession in sorted(agi_map[agi]):
            entry = sequences.get(accession)
            if not entry or not entry["seq"]:
                continue
            slot = per_locus[agi]
            slot["go"] |= entry["go"]
            slot["accessions"].append(accession)
            if (len(entry["seq"]) > len(slot["seq"])
                    or (len(entry["seq"]) == len(slot["seq"])
                        and (not slot["chosen"] or accession < slot["chosen"]))):
                slot["seq"] = entry["seq"]
                slot["chosen"] = accession
    return {agi: slot for agi, slot in per_locus.items() if slot["seq"]}


def cohort_pers005():
    """PERS-005 Aroca 2017：GO:0006096 注释 57 个蛋白中 30 个（52.6%）被检出持硫化。

    作者口径：单位 = TAIR 基因座（原文 "57 proteins annotated in the reference
    A. thaliana database"）；阳性 = 三个 WT 重复里 FDR<1% 至少一条肽段检出的持硫化蛋白
    （Dataset S1 工作表 "WT Proteins and locus"）；统计量 = 带 GO:0006096 注释的那 57 个
    基因座里被检出的比例 30/57。

    **装配上的队列范围与作者报的分母不是一回事，写清楚以免误读：**
    作者的 57 是他那张表"带该注释"那一行的规模；本函数返回的队列（keys）是
    UP000006548 的**全部**基因座，GO:0006096（含后代项）的 54 个基因座进 attribute==1，
    其余进 attribute==0——优势比需要这另一行才有分母。
    作者报的比例 30/57 在本装配里对应 attribute==1 行内的 y 均值。
    """
    view = ath_locus_view()
    rows = sheet_rows(os.path.join(SUPP, "aroca2017_jxb_PMC5853657",
                                   "erx294_suppl_supplementary_data_set_s1.xlsx"),
                      "WT Proteins and locus")
    positive_agis, declared_rows = set(), 0
    for row in rows[2:]:
        agi_cell = row[0]
        if agi_cell is None:
            continue
        declared_rows += 1
        for agi in re.split(r"[;,]", str(agi_cell)):
            agi = agi.strip().upper()
            if re.fullmatch(r"AT[1-5CM]G\d{5}", agi):
                positive_agis.add(agi)
    keys = sorted(view)
    y = np.asarray([1 if a in positive_agis else 0 for a in keys], dtype=int)
    glycolysis = go_with_descendants([GO_GLYCOLYSIS])[GO_GLYCOLYSIS]
    attribute = np.asarray([1 if view[a]["go"] & glycolysis else 0 for a in keys], dtype=int)
    sequence_of = {a: view[a]["seq"] for a in keys}
    return {
        "claim_id": "PERS-005", "unit": "protein", "keys": keys, "y": y,
        "attribute": attribute, "groups": np.asarray(keys),
        "covariates": protein_feature_matrix(sequence_of, keys),
        "covariate_set": "PROT8", "observed_depth": _aroca_depth(keys),
        "attribute_label": "GO:0006096 糖酵解",
        "author_statistic": {"kind": "proportion_in_annotated_category",
                             "value": 30 / 57, "numerator": 30, "denominator": 57,
                             "text": "30/57 = 52.6%"},
        "n_positive_set": int(len(positive_agis)), "n_positive_agi": int(len(positive_agis)),
        "declared_positive_rows": declared_rows,
        "sequences": sequence_of,
        "reproduction_notes": [
            "单位取 TAIR 基因座，与作者一致；按 UniProt 登录号装配会把分子分母同时放大"
            "（UP000006548 有 39272 条记录，基因座约两万七千个）",
            "阳性集按作者工作表原样取，多基因座行（如 AT1G07770; AT5G59850）拆开计入",
            "**队列（keys）是 UP000006548 的全部基因座，不是 GO:0006096 子集**——"
            "作者报的 57 是他那 2×2 表的属性行规模，不是队列规模；"
            "本装配把那 54 个带 GO:0006096（含后代项）的基因座放在 attribute==1 一行，"
            "其余基因座在 attribute==0 一行，优势比才有分母",
            "作者用的是其 2017 年的参考拟南芥库，GO 版本与注释来源不同，"
            "属性行规模 54 不等于 57",
        ],
    }


def _aroca_depth(keys):
    """Aroca 自己的检出深度：三个 WT 重复的 Num. Pept 之和，按 AGI 基因座取最大值；
    未检出的基因座记 0。"""
    rows = sheet_rows(os.path.join(SUPP, "aroca2017_jxb_PMC5853657",
                                   "erx294_suppl_supplementary_data_set_s1.xlsx"),
                      "Identification WT")
    per_locus = collections.defaultdict(float)
    for row in rows[3:]:
        agi_cell = row[1]
        if agi_cell is None:
            continue
        total = 0.0
        for column in (6, 9, 12):
            try:
                total += float(row[column] or 0)
            except (TypeError, ValueError, IndexError):
                pass
        for agi in re.split(r"[;,]", str(agi_cell)):
            agi = agi.strip().upper()
            if re.fullmatch(r"AT[1-5CM]G\d{5}", agi):
                per_locus[agi] = max(per_locus[agi], total)
    return np.asarray([per_locus.get(a, 0.0) for a in keys], dtype=float)


def cohort_pers006():
    """PERS-006 Aroca 2017：过表达的持硫化蛋白中 58.8% 属初级代谢（MapMan）。

    作者口径：统计单位 = Table S7 里 des1/wt 定量比较中**过表达**的那批基因座
    （log2 倍数>0 且 q≤0.05，作者数 80 个）；属性 = MapMan 本体的初级代谢
    （TCA、糖酵解、卡尔文循环）；分母是这 80 个而不是全蛋白组。
    作者未报"检出蛋白背景下的对照比例"，所以本条的背景是欠定义的——
    再检验时把背景取作同一次实验检出的全部持硫化蛋白，
    这是**本模块声明的补全**，不是作者口径。
    """
    view = ath_locus_view()
    rows = sheet_rows(os.path.join(SUPP, "aroca2017_jxb_PMC5853657",
                                   "erx294_suppl_supplementary_table_s7.xlsx"))
    over, under = set(), set()
    for row in rows[2:]:
        agi = str(row[0] or "").strip().upper()
        raw = row[3]
        if raw is None or not re.fullmatch(r"AT[1-5CM]G\d{5}", agi):
            continue
        try:
            value = float(raw)
        except (TypeError, ValueError):
            continue
        # 该表的 log2 值写成了整数（-1058 表示 -1.058），只用符号，不用数值
        (over if value > 0 else under).add(agi)
    base = cohort_pers005()
    detected = {a for a, flag in zip(base["keys"], base["y"]) if flag}
    keys = sorted((detected | over | under) & set(view))
    y = np.asarray([1 if a in over else 0 for a in keys], dtype=int)
    primary = set()
    for members in go_with_descendants(GO_PRIMARY_METABOLISM).values():
        primary |= members
    attribute = np.asarray([1 if view[a]["go"] & primary else 0 for a in keys], dtype=int)
    sequence_of = {a: view[a]["seq"] for a in keys}
    return {
        "claim_id": "PERS-006", "unit": "protein", "keys": keys, "y": y,
        "attribute": attribute, "groups": np.asarray(keys),
        "covariates": protein_feature_matrix(sequence_of, keys),
        "covariate_set": "PROT8", "observed_depth": _aroca_depth(keys),
        "attribute_label": "初级代谢 GO 等价集（MapMan 的替代口径）",
        "author_statistic": {"kind": "proportion_within_positive_set", "value": 0.588,
                             "numerator": None, "denominator": 80, "text": "58.8%"},
        "n_over": int(len(over)), "n_under": int(len(under)),
        "sequences": sequence_of,
        "reproduction_notes": [
            "Table S7 的 log2 倍数列写成整数（-1058 即 -1.058），本装配只用其符号",
            "MapMan 本体在本环境不可得，属性改用 GO 等价集 "
            + "/".join(GO_PRIMARY_METABOLISM) + "，**声明为替代口径**，"
            "因此比例不应与 58.8% 逐位对齐",
            "作者未定义背景集，再检验的背景由本模块补全为同实验检出的持硫化蛋白",
        ],
    }


def cohort_pers011():
    """PERS-011 Jurado-Flores 2023：cluster 1 在类囊体 12.57 倍、光合蛋白 7.55 倍富集。

    作者口径：阳性 = Dataset S9 列出的 cluster 1（1405 个基因座，严重碳饥饿下持硫化
    水平下降的一组）；背景 = 拟南芥全注释蛋白组（DAVID 默认全基因组背景）；
    统计量 = 倍数富集 (k/n)/(K/N) 与 EASE p 值。
    """
    view = ath_locus_view()
    rows = sheet_rows(os.path.join(SUPP, "juradoflores2023_antiox_PMC10135009",
                                   "inner", "Dataset S9.xlsx"))
    cluster_agis, declared_rows = set(), 0
    for row in rows[2:]:
        agi = str(row[0] or "").strip().upper()
        if not agi or agi == "NONE":
            continue
        declared_rows += 1
        if re.fullmatch(r"AT[1-5CM]G\d{5}", agi):
            cluster_agis.add(agi)
    keys = sorted(view)
    y = np.asarray([1 if a in cluster_agis else 0 for a in keys], dtype=int)
    expanded = go_with_descendants([GO_THYLAKOID, GO_PHOTOSYNTHESIS])
    attribute = np.asarray([1 if view[a]["go"] & expanded[GO_THYLAKOID] else 0
                            for a in keys], dtype=int)
    attribute_second = np.asarray([1 if view[a]["go"] & expanded[GO_PHOTOSYNTHESIS] else 0
                                   for a in keys], dtype=int)
    sequence_of = {a: view[a]["seq"] for a in keys}
    return {
        "claim_id": "PERS-011", "unit": "protein", "keys": keys, "y": y,
        "attribute": attribute, "attribute_secondary": attribute_second,
        "groups": np.asarray(keys),
        "covariates": protein_feature_matrix(sequence_of, keys),
        "covariate_set": "PROT8", "observed_depth": None,
        "attribute_label": "GO:0009579 类囊体",
        "attribute_secondary_label": "GO:0015979 光合作用",
        "author_statistic": {"kind": "fold_enrichment", "value": 12.57, "secondary": 7.55,
                             "p_value": 7.6e-17, "secondary_p": 1.9e-7,
                             "text": "thylakoid 12.57-fold p=7.6e-17；photosynthetic 7.55-fold p=1.9e-7"},
        "n_cluster": int(len(cluster_agis)), "n_cluster_agi": int(len(cluster_agis)),
        "declared_cluster_rows": declared_rows,
        "sequences": sequence_of,
        "reproduction_notes": [
            "cluster 1 按 Dataset S9 的 AGI 原样取，单位为 TAIR 基因座",
            "作者用 DAVID 的 thylakoid 与 photosynthetic 词条；本装配取 GO:0009579 与 "
            "GO:0015979，倍数因此依赖 GO 版本",
        ],
    }


def cohort_pers012():
    """PERS-012 Matamoros 2024：根瘤持硫化蛋白在生物过程域 55 个 GO 条目显著过表达。

    作者口径：阳性 = Tables S1–S7 列出的持硫化**植物**蛋白（菜豆）；
    背景 = AgriGO 单一富集分析的菜豆全注释背景；判据 FDR<0.001。
    Table 1 给出五个具名条目，首行 GO:0019752 羧酸代谢 8.6% 对 2.3%，FDR 2.92e-16，
    本装配以该条目为头条效应量（运行前写定），另把 55 这个计数作为复现统计量。
    """
    sequences = read_proteome_tsv("pvu")
    positives = set()
    rows_all = []
    path = os.path.join(SUPP, "matamoros2024_jxb_PMC11103110",
                        "erad436_suppl_supplementary_tables_s1-s9.xlsx")
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for index, worksheet in enumerate(book.worksheets):
        if index >= 7:  # Hoja8/Hoja9 是细菌（bacteroid）蛋白，论断针对植物侧
            continue
        rows = [list(r) for r in worksheet.iter_rows(values_only=True)]
        rows_all.append((worksheet.title, len(rows)))
        for row in rows[2:]:
            for cell in row[:3]:
                token = str(cell or "").strip()
                if ACCESSION_RE.match(token) and token in sequences:
                    positives.add(token)
    book.close()
    keys = sorted(a for a, entry in sequences.items() if entry["seq"])
    y = np.asarray([1 if a in positives else 0 for a in keys], dtype=int)
    carboxylic = go_with_descendants([GO_CARBOXYLIC])[GO_CARBOXYLIC]
    attribute = np.asarray([1 if sequences[a]["go"] & carboxylic else 0 for a in keys], dtype=int)
    return {
        "claim_id": "PERS-012", "unit": "protein", "keys": keys, "y": y,
        "attribute": attribute, "groups": np.asarray(keys),
        "covariates": protein_feature_matrix({a: sequences[a]["seq"] for a in keys}, keys),
        "covariate_set": "PROT8", "observed_depth": None,
        "attribute_label": "GO:0019752 羧酸代谢",
        "go_sets": {a: sequences[a]["go"] for a in keys},
        "author_statistic": {"kind": "enriched_term_count", "value": 55, "fdr": 1e-3,
                             "headline_term": GO_CARBOXYLIC, "headline_pct": 8.6,
                             "headline_ref_pct": 2.3, "headline_fdr": 2.92e-16,
                             "text": "55 GO terms FDR<0.001；头条 GO:0019752 8.6% vs 2.3%"},
        "table1": MATAMOROS_TABLE1, "sheets": rows_all,
        "sequences": {a: sequences[a]["seq"] for a in keys},
        "reproduction_notes": [
            "阳性集从 Tables S1–S7（植物侧）的 UniProt 列取，S8/S9 的细菌蛋白不计入",
            "背景为 UniProt 参考蛋白组 UP000000226；作者用 AgriGO 自带的菜豆注释背景，"
            "GO 版本与背景来源都不同，55 这个计数不应要求逐个对齐",
        ],
    }


def cohort_pers010():
    """PERS-010 Bithi 2021：限食后四器官共 1854 个持硫化蛋白，209 个（11.3%）为四器官共有。

    作者口径：四个器官各一张谱图计数表（Supplementary Data 1–4，肝/肾/肌/脑），
    阳性 = 四器官**都**检出的蛋白（Supplementary Data 5 列了 209 个）；
    背景 = 四器官并集 1854 个；统计量 = 209/1854 = 11.3%。
    """
    sequences = read_proteome_tsv("mmu")
    base = os.path.join(SUPP, "bithi2021_natcomm_PMC7979915")
    organs = {"liver": ("41467_2021_22001_MOESM4_ESM.xlsx", "WT Liver Prot ID Spectral Count"),
              "kidney": ("41467_2021_22001_MOESM5_ESM.xlsx", "WT Kidney Prot ID SpectralCount"),
              "muscle": ("41467_2021_22001_MOESM6_ESM.xlsx", "WT Muscle Protein ID Spectral"),
              "brain": ("41467_2021_22001_MOESM7_ESM.xlsx", "WT Brain Protein ID Spectral")}
    per_organ, spectral = {}, collections.defaultdict(float)
    for organ, (filename, sheet) in organs.items():
        rows = sheet_rows(os.path.join(base, filename), sheet)
        found = set()
        for row in rows[2:]:
            accession = str(row[1] or "").strip()
            if not ACCESSION_RE.match(accession):
                continue
            found.add(accession)
            try:
                spectral[accession] += abs(float(row[5] or 0))
            except (TypeError, ValueError):
                pass
        per_organ[organ] = found
    union = set().union(*per_organ.values())
    shared = set.intersection(*per_organ.values())
    declared = {str(r[0] or "").strip() for r in
                sheet_rows(os.path.join(base, "41467_2021_22001_MOESM8_ESM.xlsx"),
                           "209 shared prot Spectral Count")[1:]
                if ACCESSION_RE.match(str(r[0] or "").strip())}
    keys = sorted(a for a in union if a in sequences and sequences[a]["seq"])
    y = np.asarray([1 if a in shared else 0 for a in keys], dtype=int)
    n_organs = np.asarray([sum(a in per_organ[o] for o in organs) for a in keys], dtype=float)
    return {
        "claim_id": "PERS-010", "unit": "protein", "keys": keys, "y": y,
        "attribute": None, "groups": np.asarray(keys),
        "covariates": protein_feature_matrix({a: sequences[a]["seq"] for a in keys}, keys),
        "covariate_set": "PROT8",
        "observed_depth": np.asarray([spectral.get(a, 0.0) for a in keys], dtype=float),
        "n_organs": n_organs,
        "attribute_label": "四器官共有（属性即阳性本身，检验改为可检出性能否解释共有性）",
        "author_statistic": {"kind": "shared_fraction", "value": 209 / 1854,
                             "numerator": 209, "denominator": 1854,
                             "text": "209/1854 = 11.3%"},
        "per_organ": {o: len(v) for o, v in per_organ.items()},
        "union_size": int(len(union)), "shared_size": int(len(shared)),
        "declared_shared": int(len(declared)),
        "shared_matches_declared": int(len(shared & declared)),
        "sequences": {a: sequences[a]["seq"] for a in keys},
        "reproduction_notes": [
            "四器官表按 Supplementary Data 1–4 的 Accession Number 列取，谱图计数工作表",
            "共有集按四表交集自算，并与 Supplementary Data 5 声明的 209 个逐个比对",
        ],
    }


def cohort_sno021():
    """SNO-021 Wang 2023：七个 Motif-X 基序里五个以赖氨酸为核心。

    作者口径：阳性 = Table S1 的 983 个非冗余 S-亚硝基化位点（637 个蛋白）；
    背景 = Motif-X 以弓形虫全蛋白组序列为背景（即所有半胱氨酸都是候选）；
    统计量 = 七个基序及其匹配位点数（76/72/61/58/43/41/38）。
    本装配把"以 K 为核心"读成：五个 K 基序声明的偏移 +10/−10/+6/+5/+1 上出现 K。
    """
    sequences = read_proteome_tsv("tgo")
    gene_to_accession = {}
    # 同样按登录号排序，保证一个基因号映到哪条登录号不随进程变化
    for accession in sorted(sequences):
        for token in re.split(r"[;\s]+", sequences[accession]["genes"]):
            token = token.strip()
            if token:
                gene_to_accession.setdefault(token, accession)
    rows = sheet_rows(os.path.join(SUPP, "wang2023_molecules_PMC10649196", "inner", "Table S1.xlsx"))
    declared, mapped, unmapped, mismatched = 0, {}, [], []
    for row in rows[1:]:
        gene = str(row[0] or "").strip()
        try:
            position = int(row[1])
        except (TypeError, ValueError):
            continue
        declared += 1
        accession = gene_to_accession.get(gene)
        if accession is None:
            unmapped.append(gene)
            continue
        sequence = sequences[accession]["seq"]
        if not (1 <= position <= len(sequence)) or sequence[position - 1] != "C":
            mismatched.append((gene, position))
            continue
        mapped[(accession, position)] = True
    protein_set = sorted({a for a, _ in mapped})
    keys = []
    for accession in protein_set:
        sequence = sequences[accession]["seq"]
        for i, residue in enumerate(sequence):
            if residue == "C":
                keys.append((accession, i + 1))
    y = np.asarray([1 if k in mapped else 0 for k in keys], dtype=int)
    sequence_of = {a: sequences[a]["seq"] for a in protein_set}
    attribute = _flank_residue_flag(sequence_of, keys, SNO021_K_OFFSETS, "K")
    return {
        "claim_id": "SNO-021", "unit": "site", "keys": keys, "y": y,
        "attribute": attribute, "groups": np.asarray([a for a, _ in keys]),
        "covariates": site_feature_matrix(sequence_of, keys, "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "attribute_label": "五个 K 基序偏移 (+10,−10,+6,+5,+1) 上出现 K",
        "author_statistic": {"kind": "motif_count", "value": 5, "of": 7,
                             "matched_sites": [76, 72, 61, 58, 43, 41, 38],
                             "text": "7 motifs (5 K, 2 I)"},
        "k_offsets": SNO021_K_OFFSETS, "i_offsets": SNO021_I_OFFSETS,
        "declared_sites": declared, "mapped_sites": int(len(mapped)),
        "unmapped_genes": len(set(unmapped)), "position_mismatch": len(mismatched),
        "n_proteins": int(len(protein_set)),
        "sequences": sequence_of,
        "reproduction_notes": [
            "Table S1 的登录号是 ToxoDB 基因号（TGME49_*），经 UniProt UP000001529 的"
            "基因名列映射到登录号；未映射与位点不是 C 的条目逐个计数并报出",
            "阴性集按作者的 Motif-X 背景取作同一批蛋白里其余全部半胱氨酸",
        ],
    }


def _flank_residue_flag(sequence_of, keys, offsets, residue):
    flags = []
    for accession, position in keys:
        sequence = sequence_of[accession]
        hit = 0
        for offset in offsets:
            index = position - 1 + offset
            if 0 <= index < len(sequence) and sequence[index] == residue:
                hit = 1
                break
        flags.append(hit)
    return np.asarray(flags, dtype=int)


def flank_residue_flag(sequence_of, keys, offsets, residues):
    flags = []
    for accession, position in keys:
        sequence = sequence_of[accession]
        hit = 0
        for offset in offsets:
            index = position - 1 + offset
            if 0 <= index < len(sequence) and sequence[index] in residues:
                hit = 1
                break
        flags.append(hit)
    return np.asarray(flags, dtype=int)


def cohort_sfe006_transfer():
    """SFE-006 Bui 2016 的**转移检验**：次磺酰化位点侧翼富集 K/R。

    作者原数据取不到：1443 阳性 / 10521 阴性来自 SOHSite 网站与 Yang 2014 补充表，
    两者在本环境都不可达（记录见 external/intake/phase2_claim_supp/phase2_supp_fetch_log.json）。
    因此**本条不宣称复现基线**，改在本树已入库的另一个同化学队列上检验同一条论断：
    `results/ptm_census_sites.csv` 的 fps2020_ath_sulfenyl（拟南芥次磺酰化，1745 阳性 /
    411 同运行"被看见未修饰"）。属性按作者自述的偏移
    −10/−8~−6/−4/−2/+4~+8 上出现 K 或 R。

    两种阴性都装配：
      author_class  作者那一类——同一批已鉴定蛋白里其余全部半胱氨酸（未检出即阴性）
      same_run      同一次实验里被看见但未修饰的半胱氨酸
    """
    sequences = ath_tair_sequences()
    with open(os.path.join(RESULTS, "ptm_census_sites.csv"), encoding="utf-8-sig", newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r["dataset_id"] == "fps2020_ath_sulfenyl"]
    positives, observed_unmodified, off_residue = set(), set(), 0
    for row in rows:
        locus, site = row["accession"].strip().upper(), int(row["site"])
        sequence = sequences.get(locus, "")
        if not sequence or not (1 <= site <= len(sequence)) or sequence[site - 1] != "C":
            off_residue += 1
            continue
        (positives if row["role"] == "positive" else observed_unmodified).add((locus, site))
    protein_set = sorted({a for a, _ in positives | observed_unmodified})
    sequence_of = {a: sequences[a] for a in protein_set}
    keys = []
    for locus in protein_set:
        for i, residue in enumerate(sequence_of[locus]):
            if residue == "C":
                keys.append((locus, i + 1))
    y = np.asarray([1 if k in positives else 0 for k in keys], dtype=int)
    same_run_mask = np.asarray([1 if (k in positives or k in observed_unmodified) else 0
                                for k in keys], dtype=int)
    attribute = flank_residue_flag(sequence_of, keys, SFE006_KR_OFFSETS, {"K", "R"})
    return {
        "claim_id": "SFE-006", "unit": "site", "keys": keys, "y": y,
        "attribute": attribute, "groups": np.asarray([a for a, _ in keys]),
        "covariates": site_feature_matrix(sequence_of, keys, "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "same_run_mask": same_run_mask,
        "attribute_label": "侧翼 −10/−8~−6/−4/−2/+4~+8 上出现 K 或 R",
        "author_statistic": {"kind": "count_only", "value": None,
                             "n_positive": 1443, "n_negative": 10521, "n_proteins": 987,
                             "text": "1443 positives vs 10521 negatives (987 proteins)"},
        "kr_offsets": SFE006_KR_OFFSETS,
        "transfer_cohort": "fps2020_ath_sulfenyl",
        "n_positive_transfer": int(len(positives)),
        "n_observed_unmodified_transfer": int(len(observed_unmodified)),
        "n_proteins": int(len(protein_set)),
        "off_residue_rows": off_residue,
        "sequences": sequence_of,
        "reproduction_notes": [
            "作者原阳性集不可达（SOHSite 网站与 Yang 2014 补充表都取不到），"
            "**本条不报复现基线**",
            "转移到 fps2020_ath_sulfenyl 队列，物种与实验都不同，"
            "因此检验的是论断的内容而不是作者的数字",
        ],
    }


COHORTS = {
    "PERS-005": cohort_pers005, "PERS-006": cohort_pers006, "PERS-010": cohort_pers010,
    "PERS-011": cohort_pers011, "PERS-012": cohort_pers012, "SNO-021": cohort_sno021,
    "SFE-006": cohort_sfe006_transfer,
}
# SFI-004 不在此处装配：见 run_phase2_claims_under_detectability_control.py 的
# SFI004_EXCLUSION，原因是观测通道不匹配且位点级数据不在补充材料里。


def main():
    for claim_id, builder in COHORTS.items():
        cohort = builder()
        attribute = cohort.get("attribute")
        print(claim_id, cohort["unit"], "n=%d" % len(cohort["keys"]),
              "positives=%d" % int(cohort["y"].sum()),
              ("attribute=%d" % int(attribute.sum())) if attribute is not None else "attribute=n/a")


if __name__ == "__main__":
    main()
