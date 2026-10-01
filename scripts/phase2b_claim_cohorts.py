#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 B 层十四条已发表论断里**可装配的那些**各自做成一张可再检验的表。

本模块是 `scripts/phase2_claim_cohorts.py` 的 B 层对应件：**沿用第二阶段的
全部口径，不另立**。协变量定义（VIS10、PROT8）、聚类单位、装配后的统一结构
都直接 import 那个模块，只新写每条论断自己的装配函数——上一轮自述"装配函数是
唯一需新写的部分"，本模块就是那一部分。

装配后的统一结构与第二阶段完全一致（unit / keys / y / attribute / groups /
covariates / covariate_set / observed_depth / author_statistic / reproduction_notes）。

================================================================================
一 侧翼组成类属性的判读约定（**运行前写定，一条规则覆盖全部**）
================================================================================
第二阶段的两条侧翼论断（SFE-006、SNO-021）的作者都**点名了偏移**，所以属性是
"在任一点名偏移上出现该残基"。B 层有三条论断只说"邻域富集某几类残基"而没点名
偏移（PERS-003、PERS-009、SNO-016），因此需要一条把"邻域富集"变成二值属性的
约定。该约定在任何数字被看到之前写定，且**对这三条一视同仁**：

    FLANK_WINDOW  = 5      侧翼窗口 ±5，共 10 个位置（不含位点自身）
    FLANK_FRACTION = 0.30  属性 = 1 当且仅当这 10 个位置里**至少 3 个**
                           （= ceil(0.30 x 10)）属于该论断点名的残基集合

为什么取固定份额而不取"超过全蛋白组期望"：后者要先估背景频率，等于让数据参与
定阈值。固定份额 0.30 是约定，不是估计量；三条论断的残基集合大小不同
（3 类、4 类、7 类），同一份额下各自的属性流行度自然不同，这一点在报告里如实写。
敏感性口径同样写定：窗口改 ±10（20 个位置，阈值 ≥6），份额不变。

作者点名了偏移的（PERS-002 的 CXXC/CXC、SFE-007 的 Glu、SFE-008 的 K/R）
一律按点名偏移判读，不套上面的份额规则。

疏水残基集合写定为 {A, V, L, I, M, F, W}：非极性侧链，**不含 C**（含 C 会与
PERS-002 的 CXXC/CXC 属性以及多半胱氨酸肽段的可检出性耦合，那是另一条通道），
不含 P、G（构象残基，疏水性判读不一致）。

================================================================================
二 转移检验（**不宣称复现原文数字**）
================================================================================
四条论断的作者原始阳性集在本环境取不到（见
`external/intake/phase2b_claim_supp/phase2b_supp_fetch_log.json`：本轮所有外网
出口不可用），但它们主张的**内容**在本树已入库的同化学队列上可以检验。与第二阶段
SFE-006 的处理完全同型：**物种/实验不同，检的是论断的内容，不是作者的那张表。**

  SFE-007  Bui 2016 人源次磺酰化        -> fps2020_ath_sulfenyl（拟南芥次磺酰化）
  SFE-008  Bui 2016 跨修饰基序区分      -> fps2020_ath_sulfenyl 对 拟南芥 SNO 并集
  SNO-016  Lee 2011 人源 SNO（SNOSite） -> cysboost2019_human_sno_hela（人源 SNO）
  SFE-011  Yu 2022 甘蓝型油菜次磺酰化   -> fps2020_ath_sulfenyl 蛋白层（同为十字花科）

两条论断的材料在手，按作者口径直接装配，不是转移检验：
  PERS-002 / PERS-003  Longen 2016 qPerS-SID（results/qpers_sid_sites_normalised.csv）
  PERS-009             Wei 2025 Sul-BertGRU（external/sul_bertgru/）

================================================================================
三 检出深度（论文自报的那一种）
================================================================================
PERS-002 / PERS-003 的 qPerS-SID 有四个供体臂（GYY4137 / Na2S / Na2S4 / NaSH），
**深度 = 该位点在几个供体臂的 tier-B 表里出现过**，取值只有 1..4，因此深度分层
按臂计数原值分 4 层，**不取四分位**（四分位在只有四个取值时会退化）。这条也写在
运行前。其余论断的队列没有论文自报深度，记 None。

术语：sulfenylation = 次磺酰化（-SOH），sulfinylation = 亚磺酰化（-SO2H）。

用法：被 run_phase2b_claims_backfill.py import；单独运行时打印每条装配后的规模。
"""
from __future__ import annotations

import collections
import csv
import gzip
import os
import re

import numpy as np

from phase2_claim_cohorts import (
    RESULTS, ath_tair_sequences, ath_locus_view, flank_residue_flag,
    go_with_descendants, protein_feature_matrix, site_feature_matrix,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTEOMES = os.path.join(ROOT, "external", "proteomes")
SUL_BERTGRU = os.path.join(ROOT, "external", "sul_bertgru")
CENSUS = os.path.join(RESULTS, "ptm_census_sites.csv")
QPERS = os.path.join(RESULTS, "qpers_sid_sites_normalised.csv")

# ---- 运行前写定的属性判读常量 -------------------------------------------------
FLANK_WINDOW = 5
FLANK_FRACTION = 0.30
FLANK_MIN_HITS = 3          # ceil(0.30 * 10)
FLANK_WINDOW_SENS = 10
FLANK_MIN_HITS_SENS = 6     # ceil(0.30 * 20)
HYDROPHOBIC = frozenset("AVLIMFW")
# PERS-002 Longen 2016：CXXC 与 CXC。CXC = 另一个半胱氨酸在 ±2，CXXC = 在 ±3。
PERS002_CYS_OFFSETS = [-3, -2, 2, 3]
# PERS-009 Wei 2025：位点附近富集 A/K/R/V（阴性侧富集 C/R/S）
PERS009_RESIDUES = frozenset("AKRV")
# SNO-016 Lee 2011：MDD 子群在某一特定位置保守带正电残基 K/R/H
SNO016_RESIDUES = frozenset("KRH")
# SFE-007 Bui 2016：Glu 在 −3、+1、+3、+4 位明显增多（作者点名的偏移）
SFE007_GLU_OFFSETS = [-3, 1, 3, 4]
# SFE-007 的另一半：−1、+1、+2 位缺正电残基（作为声明的次属性一并报出）
SFE007_NEAR_KR_OFFSETS = [-1, 1, 2]
# SFE-008 Bui 2016：作者称次磺酰化最显著的特征是 K/R 在这些偏移上富集，
# 跨修饰区分检验就以这一特征为属性（与 SFE-006 同一偏移集，直接可比）
SFE008_KR_OFFSETS = [-10, -8, -7, -6, -4, -2, 4, 5, 6, 7, 8]
# SFE-011 Yu 2022：次磺酰化升高的蛋白富集于碳利用、光合作用、糖酵解
GO_PHOTOSYNTHESIS = "GO:0015979"
GO_GLYCOLYSIS = "GO:0006096"


# ---------------------------------------------------------------- 基础读取

def read_fasta_gz(path):
    """UniProt FASTA（gzip）-> {accession: sequence}，登录号取 sp|ACC| 的中段。"""
    sequences, name, buffer = {}, None, []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(buffer).upper()
                head = line[1:].split()[0]
                parts = head.split("|")
                name = parts[1] if len(parts) > 2 else head
                buffer = []
            else:
                buffer.append(line.strip())
    if name:
        sequences[name] = "".join(buffer).upper()
    return sequences


def read_fasta(path):
    """同上，非压缩。"""
    sequences, name, buffer = {}, None, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(buffer).upper()
                head = line[1:].split()[0]
                parts = head.split("|")
                name = parts[1] if len(parts) > 2 else head
                buffer = []
            else:
                buffer.append(line.strip())
    if name:
        sequences[name] = "".join(buffer).upper()
    return sequences


def census_rows(dataset_id):
    with open(CENSUS, encoding="utf-8-sig", newline="") as handle:
        return [r for r in csv.DictReader(handle) if r["dataset_id"] == dataset_id]


def flank_count_flag(sequence_of, keys, residues, window=FLANK_WINDOW,
                     min_hits=FLANK_MIN_HITS):
    """§一 的份额规则：±window 的 2*window 个侧翼位置里至少 min_hits 个属于 residues。

    窗口越界的位置不计入命中，也不缩小阈值——序列末端的位点因此更难满足该属性，
    这一点是规则的一部分（末端效应本身就是可检出性通道的一部分，由协变量承接）。
    """
    flags = []
    for accession, position in keys:
        sequence = sequence_of[accession]
        centre = position - 1
        hits = 0
        for offset in range(-window, window + 1):
            if offset == 0:
                continue
            index = centre + offset
            if 0 <= index < len(sequence) and sequence[index] in residues:
                hits += 1
        flags.append(1 if hits >= min_hits else 0)
    return np.asarray(flags, dtype=int)


def _all_cys_keys(sequence_of, accessions):
    keys = []
    for accession in sorted(accessions):
        sequence = sequence_of[accession]
        for i, residue in enumerate(sequence):
            if residue == "C":
                keys.append((accession, i + 1))
    return keys


# ---------------------------------------------------------------- qPerS-SID 队列

def _qpers_cohort_base():
    """Longen 2016 的同一次运行队列，PERS-002 与 PERS-003 共用。

    作者口径（`reports/QPERS_SID_LABEL_CRITERIA.md` 与
    `scripts/ingest_ptm_census_sites.py` 的 qpers_sid_tierB 规则）：
    阳性 = tier-B 富集判据（ratio_average >= 1.30，有 Significant 列时要求 '+'）；
    阴性 = **同一次运行**里同表 Elution 行低于该阈值者，加珠上留存的探针标记巯基，
    减去阳性。这正是论断表记的 same_run_mutually_exclusive 口径，
    因此本装配的队列 = 阳性 ∪ 同运行被看见未修饰，**不把从未被检出的半胱氨酸拉进来**。

    深度 = 该位点在几个供体臂（GYY4137/Na2S/Na2S4/NaSH）的 tier-B 表里出现过。
    """
    sequences = read_fasta_gz(os.path.join(PROTEOMES, "hsa.fasta.gz"))
    rows = census_rows("qpers_sid_tierB")
    positives, unmodified, off_residue = set(), set(), 0
    for row in rows:
        accession, site = row["accession"].strip(), int(row["site"])
        sequence = sequences.get(accession, "")
        if not sequence or not (1 <= site <= len(sequence)) or sequence[site - 1] != "C":
            off_residue += 1
            continue
        (positives if row["role"] == "positive" else unmodified).add((accession, site))
    unmodified -= positives

    arms = collections.defaultdict(set)
    with open(QPERS, encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["evidence_tier"] != "B_enrichment_based":
                continue
            donor = row["donor"].strip()
            if donor in ("GYY4137", "Na2S", "Na2S4", "NaSH"):
                try:
                    arms[(row["accession"].strip(), int(row["site"]))].add(donor)
                except (TypeError, ValueError):
                    continue

    keys = sorted(positives | unmodified)
    y = np.asarray([1 if k in positives else 0 for k in keys], dtype=int)
    depth = np.asarray([len(arms.get(k, ())) for k in keys], dtype=float)
    accessions = sorted({a for a, _ in keys})
    sequence_of = {a: sequences[a] for a in accessions}
    return {"keys": keys, "y": y, "depth": depth, "sequence_of": sequence_of,
            "n_positive": int(len(positives)), "n_unmodified": int(len(unmodified)),
            "off_residue_rows": off_residue, "n_proteins": len(accessions)}


def cohort_pers002():
    """PERS-002 Longen 2016：CXXC/CXC 基序**不是**持硫化的优先靶点（阴性论断）。

    作者原文 "such motifs are obviously not a preferential target for persulfide
    formation"，未给检验统计量。属性按点名基序判读：位点的 ±2（CXC）或 ±3（CXXC）
    位上另有一个半胱氨酸。阳性/阴性见 `_qpers_cohort_base` 的 docstring。

    **这是一条阴性论断**，判定按驱动脚本 §五 的阴性论断读法，不套正向论断的
    "保留一半以上"规则。
    """
    base = _qpers_cohort_base()
    attribute = flank_residue_flag(base["sequence_of"], base["keys"],
                                   PERS002_CYS_OFFSETS, {"C"})
    return {
        "claim_id": "PERS-002", "unit": "site", "keys": base["keys"], "y": base["y"],
        "attribute": attribute,
        "groups": np.asarray([a for a, _ in base["keys"]]),
        "covariates": site_feature_matrix(base["sequence_of"], base["keys"], "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "observed_depth_declared": base["depth"],
        "observed_depth_label": "供体臂计数（GYY4137/Na2S/Na2S4/NaSH 里出现过几个）",
        "claim_direction": "null_no_preference",
        "attribute_label": "±2（CXC）或 ±3（CXXC）位另有一个半胱氨酸",
        "author_statistic": {"kind": "none_reported", "value": None,
                             "n_positive_paper": 742, "n_protein_paper": 540,
                             "text": "no significant enrichment（作者未给检验）"},
        "n_positive_assembled": base["n_positive"],
        "n_observed_unmodified": base["n_unmodified"],
        "n_proteins": base["n_proteins"], "off_residue_rows": base["off_residue_rows"],
        "sequences": base["sequence_of"], "is_transfer": False,
        "reproduction_notes": [
            "材料在手（Sci Rep 开放获取补充表，已入库为 results/qpers_sid_sites_normalised.csv）",
            "阳性与阴性按作者自己的 1.30 富集阈值与同运行口径取，不换阴性定义",
            "作者未给该论断的检验统计量（none_reported），因此'复现基线'只能核对"
            "队列规模与方向，不能核对 p 值——本行记为基线可复现当且仅当队列规模与"
            "作者自报的 742 肽 / 540 蛋白同量级",
            "属性判读用点名基序（CXC/CXXC），不套 §一 的份额规则",
        ],
    }


def cohort_pers003():
    """PERS-003 Longen 2016：持硫化 Cys 周围有疏水残基（但频率低于 SNO）。

    作者原文 "We also detected hydrophobic amino acids surrounding the persulfide
    forming cysteine"，纯定性，未给频率表或检验，且未点名偏移——因此属性按 §一 的
    份额规则判读：±5 的十个侧翼位置里至少 3 个属于 {A,V,L,I,M,F,W}。

    **论断的后半句"频率低于 SNO"本轮不检验**：它的参照是 Doulias 的 SNO 组成表，
    那张表在本环境不可达（第二阶段 §2 已记录 Doulias 2010 的九条全部卡在此处）。
    只检验前半句，这一限制逐条写进报告。
    """
    base = _qpers_cohort_base()
    attribute = flank_count_flag(base["sequence_of"], base["keys"], HYDROPHOBIC)
    attribute_sens = flank_count_flag(base["sequence_of"], base["keys"], HYDROPHOBIC,
                                      FLANK_WINDOW_SENS, FLANK_MIN_HITS_SENS)
    return {
        "claim_id": "PERS-003", "unit": "site", "keys": base["keys"], "y": base["y"],
        "attribute": attribute, "attribute_secondary": attribute_sens,
        "groups": np.asarray([a for a, _ in base["keys"]]),
        "covariates": site_feature_matrix(base["sequence_of"], base["keys"], "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "observed_depth_declared": base["depth"],
        "observed_depth_label": "供体臂计数（GYY4137/Na2S/Na2S4/NaSH 里出现过几个）",
        "claim_direction": "positive_preference",
        "attribute_label": "±5 的十个侧翼位置里 ≥3 个属于 {A,V,L,I,M,F,W}",
        "attribute_secondary_label": "±10 的二十个侧翼位置里 ≥6 个属于同一集合（敏感性）",
        "author_statistic": {"kind": "none_reported", "value": None,
                             "n_positive_paper": 742,
                             "text": "定性描述（with less frequency），未给频率表或检验"},
        "n_positive_assembled": base["n_positive"],
        "n_observed_unmodified": base["n_unmodified"],
        "n_proteins": base["n_proteins"],
        "sequences": base["sequence_of"], "is_transfer": False,
        "reproduction_notes": [
            "材料在手；阳性与阴性同 PERS-002",
            "作者未点名偏移，属性按 §一 的 ±5/≥3 份额规则判读，±10/≥6 作敏感性",
            "**只检验前半句**：后半句'频率低于 SNO'的参照（Doulias SNO 组成）不可达",
        ],
    }


# ---------------------------------------------------------------- Sul-BertGRU

def cohort_pers009():
    """PERS-009 Wei 2025 Sul-BertGRU：持硫化位点邻域富集 A/K/R/V。

    作者口径：阳性/阴性划分沿用其发布的 positive1.txt / negative1.txt（已入库为
    external/sul_bertgru/mapped_sites.json，label 1/0），**阴性是他人数据集里
    未被标注为持硫化的半胱氨酸**（curated_negatives_other_db），这正是本项目所
    审计的偏倚结构，所以按作者原样用，不换。

    属性按 §一 的份额规则：±5 的十个侧翼位置里至少 3 个属于 {A,K,R,V}。
    另报两个声明的次属性：只看 {K,R}（胰酶切点子集）与 ±10 窗口。
    """
    import json
    sequences = read_fasta(os.path.join(SUL_BERTGRU, "uniprot.fasta"))
    with open(os.path.join(SUL_BERTGRU, "mapped_sites.json"), encoding="utf-8") as handle:
        entries = json.load(handle)
    keys, labels, off_residue, missing = [], [], 0, 0
    seen = set()
    for entry in entries:
        accession = str(entry["accession"]).strip()
        position = int(entry["position"])
        sequence = sequences.get(accession)
        if sequence is None:
            missing += 1
            continue
        if not (1 <= position <= len(sequence)) or sequence[position - 1] != "C":
            off_residue += 1
            continue
        key = (accession, position)
        if key in seen:
            continue
        seen.add(key)
        keys.append(key)
        labels.append(int(entry["label"]))
    order = sorted(range(len(keys)), key=lambda i: keys[i])
    keys = [keys[i] for i in order]
    y = np.asarray([labels[i] for i in order], dtype=int)
    accessions = sorted({a for a, _ in keys})
    sequence_of = {a: sequences[a] for a in accessions}
    attribute = flank_count_flag(sequence_of, keys, PERS009_RESIDUES)
    attribute_kr = flank_count_flag(sequence_of, keys, frozenset("KR"), FLANK_WINDOW, 2)
    return {
        "claim_id": "PERS-009", "unit": "site", "keys": keys, "y": y,
        "attribute": attribute, "attribute_secondary": attribute_kr,
        "groups": np.asarray([a for a, _ in keys]),
        "covariates": site_feature_matrix(sequence_of, keys, "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "claim_direction": "positive_preference",
        "attribute_label": "±5 的十个侧翼位置里 ≥3 个属于 {A,K,R,V}",
        "attribute_secondary_label": "±5 的十个侧翼位置里 ≥2 个属于 {K,R}（胰酶切点子集）",
        "author_statistic": {"kind": "none_reported", "value": None,
                             "text": "序列保守性图谱的定性描述；未给位置频率表或检验"},
        "n_proteins": len(accessions), "off_residue_rows": off_residue,
        "missing_accessions": missing,
        "sequences": sequence_of, "is_transfer": False,
        "reproduction_notes": [
            "材料在手（作者发布的正负集，本树已在用）",
            "阴性按作者原样取（他人数据集的未检出半胱氨酸），不换定义",
            "作者未点名偏移，属性按 §一 的 ±5/≥3 份额规则判读",
        ],
    }


# ---------------------------------------------------------------- 转移检验

def _ath_sulfenyl_site_cohort():
    """fps2020_ath_sulfenyl 的位点层队列，与第二阶段 cohort_sfe006_transfer 同构造。

    队列 = 该数据集涉及蛋白（TAIR 基因座）里的**全部**半胱氨酸，阳性 = 作者报的
    次磺酰化位点。这是 Bui 2016 自己的阴性类（同一批已鉴定蛋白里其余全部 Cys），
    因此与第二阶段 SFE-006 的基线逐位可比。
    """
    sequences = ath_tair_sequences()
    rows = census_rows("fps2020_ath_sulfenyl")
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
    keys = _all_cys_keys(sequence_of, protein_set)
    y = np.asarray([1 if k in positives else 0 for k in keys], dtype=int)
    same_run = np.asarray([1 if (k in positives or k in observed_unmodified) else 0
                           for k in keys], dtype=int)
    return {"keys": keys, "y": y, "sequence_of": sequence_of, "same_run_mask": same_run,
            "n_positive": len(positives), "n_observed_unmodified": len(observed_unmodified),
            "n_proteins": len(protein_set), "off_residue_rows": off_residue,
            "positives": positives}


def cohort_sfe007_transfer():
    """SFE-007 Bui 2016 的**转移检验**：紧邻位缺正电残基而 −3/+1/+3/+4 富集谷氨酸。

    作者原数据（1443 阳性 / 10521 阴性，来自 SOHSite 网站与 Yang 2014 补充表）
    与 SFE-006 同源，本环境同样不可达，**因此本条不宣称复现基线**，
    转移到 fps2020_ath_sulfenyl（拟南芥次磺酰化）检验同一条论断的内容。

    主属性 = 作者点名的 Glu 偏移 −3/+1/+3/+4 上出现 E。
    次属性 = 紧邻 −1/+1/+2 位**不出现** K 或 R（论断的另一半，"缺正电残基"），
    一并报出但不作主判定——两半是同一张图的两面，主判定取有方向且有点名残基的那半。
    """
    base = _ath_sulfenyl_site_cohort()
    attribute = flank_residue_flag(base["sequence_of"], base["keys"],
                                   SFE007_GLU_OFFSETS, {"E"})
    near_kr = flank_residue_flag(base["sequence_of"], base["keys"],
                                 SFE007_NEAR_KR_OFFSETS, {"K", "R"})
    return {
        "claim_id": "SFE-007", "unit": "site", "keys": base["keys"], "y": base["y"],
        "attribute": attribute, "attribute_secondary": 1 - near_kr,
        "groups": np.asarray([a for a, _ in base["keys"]]),
        "covariates": site_feature_matrix(base["sequence_of"], base["keys"], "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "same_run_mask": base["same_run_mask"],
        "claim_direction": "positive_preference",
        "attribute_label": "作者点名偏移 −3/+1/+3/+4 上出现 Glu",
        "attribute_secondary_label": "紧邻 −1/+1/+2 位不出现 K 或 R（论断的另一半）",
        "author_statistic": {"kind": "none_reported", "value": None,
                             "n_positive": 1443, "n_negative": 10521,
                             "text": "正文为图示描述，未给 p 值"},
        "transfer_cohort": "fps2020_ath_sulfenyl",
        "n_positive_transfer": base["n_positive"],
        "n_observed_unmodified_transfer": base["n_observed_unmodified"],
        "n_proteins": base["n_proteins"],
        "sequences": base["sequence_of"], "is_transfer": True,
        "reproduction_notes": [
            "作者原阳性集不可达（与 SFE-006 同源），**本条不报复现基线**",
            "转移到 fps2020_ath_sulfenyl；物种与实验都不同，检的是论断的内容",
            "阴性类与作者一致：同一批已鉴定蛋白里其余全部 Cys",
        ],
    }


def cohort_sfe008_transfer():
    """SFE-008 Bui 2016 的**转移检验**：次磺酰化基序可与 SNO 基序区分。

    作者用 TwoSampleLogo 对比三类位点（次磺酰化 / SNO / 谷胱甘肽化），未给区分度
    统计量。本轮把它化成第二阶段同一个 2x2：

        y = 1 该位点属次磺酰化类，y = 0 属 S-亚硝基化类
        属性 = 作者称次磺酰化最显著的那个特征（K/R 在 SFE-006 偏移集上出现）

    **谷胱甘肽化那一路不检验**：本树没有入库的谷胱甘肽化位点队列，逐条写进报告。

    两类都取拟南芥，以去掉物种差异：次磺酰化 = fps2020_ath_sulfenyl 阳性，
    SNO = natcomm2023_ath_sno 与 abiotech2025_ath_sno 阳性的**并集**（运行前写定
    取并集而不是挑一个）；两类都出现的位点按定义无法归类，整体剔除并计数。

    这条转移检验的限制比其余三条更重，必须随数字一起引用：**y 的含义是
    "这个位点被写进了哪一类实验的表"**，所以两臂的检出通道本来就不同；
    VIS10 只承接两臂共有的胰酶可见性，不承接两种富集化学的差别。
    """
    sequences = ath_tair_sequences()
    def positives_of(dataset_id):
        found = set()
        for row in census_rows(dataset_id):
            locus, site = row["accession"].strip().upper(), int(row["site"])
            if row["role"] != "positive":
                continue
            sequence = sequences.get(locus, "")
            if not sequence or not (1 <= site <= len(sequence)) or sequence[site - 1] != "C":
                continue
            found.add((locus, site))
        return found

    sulfenyl = positives_of("fps2020_ath_sulfenyl")
    sno = positives_of("natcomm2023_ath_sno") | positives_of("abiotech2025_ath_sno")
    overlap = sulfenyl & sno
    sulfenyl_only, sno_only = sulfenyl - overlap, sno - overlap
    keys = sorted(sulfenyl_only | sno_only)
    y = np.asarray([1 if k in sulfenyl_only else 0 for k in keys], dtype=int)
    protein_set = sorted({a for a, _ in keys})
    sequence_of = {a: sequences[a] for a in protein_set}
    attribute = flank_residue_flag(sequence_of, keys, SFE008_KR_OFFSETS, {"K", "R"})
    return {
        "claim_id": "SFE-008", "unit": "site", "keys": keys, "y": y,
        "attribute": attribute,
        "groups": np.asarray([a for a, _ in keys]),
        "covariates": site_feature_matrix(sequence_of, keys, "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "claim_direction": "positive_preference",
        "attribute_label": "K 或 R 出现在 SFE-006 偏移集（−10/−8~−6/−4/−2/+4~+8）",
        "author_statistic": {"kind": "none_reported", "value": None,
                             "n_positive": 1443, "n_negative": 10521,
                             "text": "TwoSampleLogo 对比；正文未给区分度统计量"},
        "transfer_cohort": "fps2020_ath_sulfenyl vs natcomm2023_ath_sno+abiotech2025_ath_sno",
        "n_sulfenyl_only": len(sulfenyl_only), "n_sno_only": len(sno_only),
        "n_overlap_dropped": len(overlap), "n_proteins": len(protein_set),
        "sequences": sequence_of, "is_transfer": True,
        "reproduction_notes": [
            "作者原三类数据不可达，**本条不报复现基线**",
            "转移为拟南芥内的两类对比（次磺酰化 对 SNO 并集），谷胱甘肽化那一路无队列，不检验",
            "两类都出现的 %d 个位点整体剔除（无法归类）" % len(overlap),
            "y 的含义是'被写进了哪一类实验的表'，两臂检出通道本来不同，"
            "VIS10 只承接共有的胰酶可见性——这条限制比其余转移检验更重",
        ],
    }


def cohort_sno016_transfer():
    """SNO-016 Lee 2011 的**转移检验**：多数 MDD 子群在特定位置保守带正电残基。

    作者原数据（SNOSite 训练集与 PLoS One Supporting Information）本轮取不到——
    所有外网出口不可用，记录在 phase2b_supp_fetch_log.json，**因此不宣称复现基线**。
    转移到本树已入库的同物种同化学队列 cysboost2019_human_sno_hela
    （人源 SNO，阳性 8184 / 同运行被看见未修饰 15457）。
    运行前写定取 HeLa 臂作主队列（两臂里规模更大的那个），SH-SY5Y 臂留作敏感性。

    阴性类按作者的类别取（all_cys_in_identified_proteins：数据集蛋白里未被报道为
    S-亚硝基化的 Cys），所以队列 = 这些蛋白的全部半胱氨酸。
    属性按 §一 的份额规则：±5 的十个侧翼位置里至少 3 个属于 {K,R,H}。
    """
    sequences = read_fasta_gz(os.path.join(PROTEOMES, "hsa.fasta.gz"))
    rows = census_rows("cysboost2019_human_sno_hela")
    positives, observed_unmodified, off_residue = set(), set(), 0
    for row in rows:
        accession, site = row["accession"].strip(), int(row["site"])
        sequence = sequences.get(accession, "")
        if not sequence or not (1 <= site <= len(sequence)) or sequence[site - 1] != "C":
            off_residue += 1
            continue
        (positives if row["role"] == "positive" else observed_unmodified).add((accession, site))
    protein_set = sorted({a for a, _ in positives | observed_unmodified})
    sequence_of = {a: sequences[a] for a in protein_set}
    keys = _all_cys_keys(sequence_of, protein_set)
    y = np.asarray([1 if k in positives else 0 for k in keys], dtype=int)
    same_run = np.asarray([1 if (k in positives or k in observed_unmodified) else 0
                           for k in keys], dtype=int)
    attribute = flank_count_flag(sequence_of, keys, SNO016_RESIDUES)
    attribute_sens = flank_count_flag(sequence_of, keys, SNO016_RESIDUES,
                                      FLANK_WINDOW_SENS, FLANK_MIN_HITS_SENS)
    return {
        "claim_id": "SNO-016", "unit": "site", "keys": keys, "y": y,
        "attribute": attribute, "attribute_secondary": attribute_sens,
        "groups": np.asarray([a for a, _ in keys]),
        "covariates": site_feature_matrix(sequence_of, keys, "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "same_run_mask": same_run,
        "claim_direction": "positive_preference",
        "attribute_label": "±5 的十个侧翼位置里 ≥3 个属于 {K,R,H}",
        "attribute_secondary_label": "±10 的二十个侧翼位置里 ≥6 个属于同一集合（敏感性）",
        "author_statistic": {"kind": "p_value_only", "value": 10,
                             "text": "10 个 MDD 子群有 K/R/H 保守基序（卡方，值未在正文给出）"},
        "transfer_cohort": "cysboost2019_human_sno_hela",
        "n_positive_transfer": len(positives),
        "n_observed_unmodified_transfer": len(observed_unmodified),
        "n_proteins": len(protein_set), "off_residue_rows": off_residue,
        "sequences": sequence_of, "is_transfer": True,
        "reproduction_notes": [
            "作者原数据本轮取不到（外网出口全不可用），**本条不报复现基线**",
            "转移到 cysboost2019_human_sno_hela：同物种（人）、同化学（SNO），实验不同",
            "阴性类与作者一致：数据集蛋白里其余全部 Cys",
            "作者未点名偏移，属性按 §一 的 ±5/≥3 份额规则判读",
        ],
    }


def cohort_sfe011_transfer():
    """SFE-011 Yu 2022 的**转移检验**：次磺酰化升高的蛋白富集于光合与糖酵解。

    作者原数据（Front Plant Sci 补充材料的 TMT 位点表）与甘蓝型油菜参考蛋白组
    本轮都取不到（外网出口全不可用），**因此不宣称复现基线**。转移到
    fps2020_ath_sulfenyl 的**蛋白层**：同为十字花科，同为次磺酰化化学，
    背景类与作者一致（全注释参考蛋白组，不是同次实验检出的蛋白）。

    单位 = TAIR 基因座，与第二阶段三条拟南芥论断一致（`ath_locus_view`）。
    阳性 = fps2020_ath_sulfenyl 里有次磺酰化位点被检出的基因座。
    主属性 = GO:0015979 光合作用（含 is_a/part_of 后代项，后代集取自第二阶段
    已缓存的 phase2_go_descendants.json，**不需要联网**）；
    次属性 = GO:0006096 糖酵解。协变量 PROT8。

    与作者论断的差别要写清：作者的阳性是"冻害下次磺酰化**升高**的蛋白"（定量比较），
    本转移队列的阳性是"被检出次磺酰化的蛋白"（定性检出）。差的是定量对比那一层，
    这一点必须随数字引用。
    """
    view = ath_locus_view()
    rows = census_rows("fps2020_ath_sulfenyl")
    positive_loci = set()
    for row in rows:
        if row["role"] != "positive":
            continue
        locus = row["accession"].strip().upper()
        if locus in view:
            positive_loci.add(locus)
    keys = sorted(view)
    y = np.asarray([1 if a in positive_loci else 0 for a in keys], dtype=int)
    expanded = go_with_descendants([GO_PHOTOSYNTHESIS, GO_GLYCOLYSIS])
    attribute = np.asarray([1 if view[a]["go"] & expanded[GO_PHOTOSYNTHESIS] else 0
                            for a in keys], dtype=int)
    attribute_second = np.asarray([1 if view[a]["go"] & expanded[GO_GLYCOLYSIS] else 0
                                   for a in keys], dtype=int)
    sequence_of = {a: view[a]["seq"] for a in keys}
    return {
        "claim_id": "SFE-011", "unit": "protein", "keys": keys, "y": y,
        "attribute": attribute, "attribute_secondary": attribute_second,
        "groups": np.asarray(keys),
        "covariates": protein_feature_matrix(sequence_of, keys),
        "covariate_set": "PROT8", "observed_depth": None,
        "claim_direction": "positive_preference",
        "attribute_label": "GO:0015979 光合作用（含后代项）",
        "attribute_secondary_label": "GO:0006096 糖酵解（含后代项）",
        "author_statistic": {"kind": "none_reported", "value": None,
                             "text": "GO 富集：碳利用、光合作用、糖酵解；正文未给 p 值"},
        "transfer_cohort": "fps2020_ath_sulfenyl（蛋白层）",
        "n_positive_loci": int(len(positive_loci)),
        "sequences": sequence_of, "is_transfer": True,
        "reproduction_notes": [
            "作者原位点表与甘蓝型油菜参考蛋白组本轮都取不到，**本条不报复现基线**",
            "转移到拟南芥（同科、同化学）蛋白层；背景类与作者一致（全注释参考蛋白组）",
            "阳性口径有差别：作者是'冻害下次磺酰化升高'的定量比较，"
            "本队列是'被检出次磺酰化'的定性检出，差的是定量对比那一层",
            "GO 后代项取自第二阶段已缓存的 phase2_go_descendants.json（本轮无外网）",
        ],
    }


COHORTS = {
    "PERS-002": cohort_pers002,
    "PERS-003": cohort_pers003,
    "PERS-009": cohort_pers009,
    "SNO-016": cohort_sno016_transfer,
    "SFE-007": cohort_sfe007_transfer,
    "SFE-008": cohort_sfe008_transfer,
    "SFE-011": cohort_sfe011_transfer,
}
# 其余七条不在此处装配，逐条理由写在 run_phase2b_claims_backfill.py 的
# OUT_OF_SCOPE 与 BLOCKED_MATERIAL 里。


def main():
    for claim_id, builder in COHORTS.items():
        cohort = builder()
        attribute = cohort.get("attribute")
        print(claim_id, cohort["unit"], "n=%d" % len(cohort["keys"]),
              "positives=%d" % int(cohort["y"].sum()),
              "attribute=%d" % int(attribute.sum()),
              "transfer=%s" % cohort["is_transfer"], flush=True)


if __name__ == "__main__":
    main()
