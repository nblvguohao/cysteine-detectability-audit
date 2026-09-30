#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第二阶段：把 A 层已发表位点/蛋白偏好论断放到可检出性控制下重做。

第一阶段已经证明两件事（reports/PHASE1_DETECTABILITY_SYNTHESIS.md）：可检出性主导
是**阴性定义**的性质；51 条已发表论断里 30 条的阴性集正是"未被检出的半胱氨酸"。
本轮回答剩下的那一问：**那些论断在可检出性控制下还成立吗？**

本文件里没有任何模型被称为 persulfidation site predictor；这里根本不拟合排序模型，
只有列联表、倾向得分匹配与聚类自助。术语：sulfenylation = 次磺酰化（-SOH），
sulfinylation = 亚磺酰化（-SO2H）。

================================================================================
一 统一效应量（运行前写定）
================================================================================
每条论断都化成同一张 2×2 列联表：行 = 被主张的属性 A（属某 GO 类 / 侧翼带 K 或 R），
列 = 论文的阳性指示 y（被检出持硫化 / 属 cluster 1 / 四器官共有）。

    效应量 = log2 优势比 = log2[(n11 * n00) / (n10 * n01)]

零格用 Haldane–Anscombe 加 0.5 修正（**运行前写定**，不按表调整）。
作者自己的统计量（比例、倍数富集、条目计数）照原样另报一列，不用它做前后对比——
倍数富集在背景改变后不可比，优势比可比。

PERS-010 不是列联表论断（属性就是阳性本身），单独声明：

    效应量 = log2(观测共有数 / 独立模型期望共有数)

期望值由各器官检出率的独立乘积给出；基线用边际检出率，控制后用层内检出率。
读法：若控制后该比值落回 0，则"四器官共有的核心持硫化组"是检出深度的产物。

================================================================================
二 两种控制（都做，运行前写定）
================================================================================
可检出性协变量**直接 import** `run_cross_protease_detectability_probe` 的定义，
不重写：位点层用 VIS10（十个纯可见性特征），蛋白层用 PROT8
（`phase2_claim_cohorts.py` 声明的聚合口径）。倾向得分 = 以 y 为因变量、
协变量标准化后的**逻辑回归**（lbfgs, C=1.0, max_iter=1000）。
选逻辑回归而不是 HistGradientBoosting，是因为后者的默认值随 sklearn 版本变动，
本轮在沙箱 python3（sklearn 1.9.1）里跑，而项目 venv 是 1.6.1。

* **控制一 检出深度分层**：按倾向得分五分位分 5 层，层内算优势比，
  用 Mantel–Haenszel 合并。论文自报检出深度可得时（PERS-005/006 的 Num. Pept、
  PERS-010 的谱图计数），另按其四分位再做一版作为敏感性。
* **控制二 可检出性匹配**：在倾向得分的 logit 上 1:1 最近邻无放回匹配，
  卡尺 = 0.2 × SD(logit)（Rosenbaum–Rubin 的常规卡尺，运行前写定）。
  匹配后重算优势比。

**凡是做匹配都配一个同规模随机对照。** 这条是本项目 v3 深度匹配那一轮踩过的坑：
匹配队列与同规模随机对照都跨零时，结论是功效陈述，不是"混淆被控制变量带走了"。
随机对照 = 保留同一批被匹配上的阳性，从阴性里等量随机抽 20 次
（种子 20260915+i）；报 20 次的均值，区间由 5000 次自助给出，
每个自助复本同时重抽聚类与随机挑一次抽样，把两种随机性一起积进区间。

================================================================================
三 前置检查（施加控制前先跑，阈值取那个模块预先声明的，不改）
================================================================================
`scripts/audit_cleaning_and_grouping.py` 的两条：

* `cleaning_collinearity(y, keep=匹配掩码, covariate=倾向得分, strata=五分位)`
  ——匹配规则会不会反而抬高它要消除的那个协变量自己的 AUC（阈值
  COVARIATE_RISE_LIMIT=0.02）、减少有信息层数、或把超过 ORPHAN_LIMIT=5% 的阳性
  推进单类层。任一条成立即 BLOCKING，该口径**不能作主口径**。
* `permutation_degeneracy(y, groups)`——分组置换还能不能打乱标签
  （阈值 MOVED_RATIO_LIMIT=0.5）。
  **蛋白层论断的聚类单位就是蛋白，每组一个观测，组内置换按构造必然退化**；
  这不是失败而是该检查的正确判决，退化时改用全局置换作零假设核对。
  为省时间，组数 > 2000 时不向该模块传 groups（它会跑 5000 次聚类区间，
  27000 个单点组下不可行），审计里记 interval_skipped_for_cost。

================================================================================
四 区间
================================================================================
5000 次聚类自助，种子 **20260915**。聚类单位：蛋白层论断按蛋白，位点层论断按蛋白
整体重抽（同一蛋白的位点一起进出）。自助用多项式权重实现（对聚类抽多重度，
再算加权列联表），因此 5000 次在 27000 个单位上也能跑完。
区间**不覆盖**：GO 版本与注释来源的不确定性、作者原始质谱搜索的不确定性、
协变量定义的替代选择。

================================================================================
五 三类判定（阈值与读法运行前写定）
================================================================================
记基线效应 b，控制后效应 c（及其 95% 区间），同规模随机对照效应 r（及区间）。
主口径**声明为可检出性匹配**（它带随机对照，能把"不成立"和"不可判定"分开），
检出深度分层作次口径。

* **存活 survives**：sign(c)=sign(b)，c 区间不跨零，且 |c| ≥ 0.5·|b|
* **衰减存活 attenuated**：sign(c)=sign(b)，c 区间不跨零，但 |c| < 0.5·|b|
* **消失 vanishes**：c 区间跨零，**且** r 区间不跨零
  （同规模随机对照下效应还在，所以损失来自控制而不是样本量）
* **反转 reverses**：sign(c) ≠ sign(b)，且 c 区间不跨零
* **不可判定 undecidable**：c 区间跨零**且** r 区间也跨零（功效受限）；
  或前置检查判 BLOCKING；或基线本身复现不出

次口径与主口径判定不同类时记 caliber_disagreement=True，并把"存活"降为"衰减存活"、
把"消失"降为"不可判定"。**区间跨零一律不写成"不成立"**，按上面两条分开写。

================================================================================
六 不做的事
================================================================================
* SFI-004（Garrido Ruiz 2022，PDB 里氧化态 Cys 更暴露）**不施加本轮控制**：
  它的观测通道是晶体可解析性，不是质谱肽段可见性，VIS10/DIG25 按定义描述不了它；
  且位点级相对 SASA 既不在补充材料（bi2c00349_si_001.pdf 只有图）也无法在本环境
  重算（无 SASA 引擎、需取回 1124 个 PDB 结构）。判定记 out_of_instrument_scope。
* SFE-006 **不宣称复现基线**：作者原阳性集在 SOHSite 网站与 Yang 2014 补充表，
  两者本环境不可达。改做转移检验，见 cohort 函数的 docstring。
* 不改 results/v2_blocker_log.csv 与 results/v2_requirement_status.csv，
  不跑 build_artifact_manifest.py，不动 protocols/v3_qtrp_analysis_plan.json。

产物（文件名独立，不覆盖任何既有产物）：
  results/phase2_claim_retest.csv
  results/phase2_claim_retest_audit.json
  results/phase2_claim_retest_strata.csv

用法：
    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=scripts python3 \
        scripts/run_phase2_claims_under_detectability_control.py
"""
from __future__ import annotations

import collections
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np
from scipy.stats import fisher_exact
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

import phase2_claim_cohorts as cohorts_module
from audit_cleaning_and_grouping import (
    COVARIATE_RISE_LIMIT, MOVED_RATIO_LIMIT, ORPHAN_LIMIT,
    cleaning_collinearity, permutation_degeneracy,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
OUT_CSV = os.path.join(RESULTS, "phase2_claim_retest.csv")
OUT_AUDIT = os.path.join(RESULTS, "phase2_claim_retest_audit.json")
OUT_STRATA = os.path.join(RESULTS, "phase2_claim_retest_strata.csv")

SEED = 20260915
REPLICATES = 5000
N_STRATA = 5
CALIPER_SD = 0.2
RANDOM_DRAWS = 20
SURVIVE_RETENTION = 0.5
HALDANE = 0.5
GROUPS_FOR_PRECHECK_LIMIT = 2000

# 基线是否按论文自述口径复现，逐条写定（不靠规则推导——第一版用
# `claim_id != "SFE-006"` 推导，把 PERS-006 错标成已复现）。
# 判读：作者统计量的方向与显著性都能在我的装配上重现即 True；
# 本体/材料不可得而只能换口径的记 False，理由写在 reproduction_notes 与报告 §3。
BASELINE_REPRODUCED = {
    "PERS-005": True,    # 该类 54 个基因座 36 个被检出 对 作者 57 个里 30 个，方向与量级一致
    "PERS-011": True,    # cluster 1 1401 对 1405；方向与显著性一致，倍数因 DAVID 分母口径偏低
    "PERS-012": True,    # 头条条目参考侧 2.18% 对 2.3%；"55 个条目"未再检验
    "PERS-010": True,    # 共有集 207/207 与 Supplementary Data 5 声明逐个对上；分母 1854 未复现
    "SNO-021": True,     # 位点集 982 对 983，差 1 个无 UniProt 对应的基因号
    "PERS-006": False,   # MapMan 本体不可得，GO 等价集下 11.4% 对 作者 58.8%，且基线不显著
    "SFE-006": False,    # 作者阳性集不可达，改为转移检验
    "SFI-004": False,    # 位点级 SASA 不在补充材料
}

SFI004_EXCLUSION = {
    "claim_id": "SFI-004",
    "verdict": "out_of_instrument_scope",
    "reason_zh": "观测通道不匹配：该论断来自 PDB 结构普查，混淆通道是晶体可解析性，"
                 "不是质谱肽段可见性，预先声明的 VIS10/DIG25 按定义描述不了它；"
                 "且位点级相对 SASA 不在补充材料（bi2c00349_si_001.pdf 仅含图 S1–S4），"
                 "本环境也无 SASA 引擎可重算 1124 个结构。",
    "baseline_reproduced": False,
    "counts_from_paper": {"structures": 1124, "CSO": 1171, "CSD": 469, "OCS": 382, "CYS": 7103},
}


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# ------------------------------------------------------------------ 统计工具

def log_odds_ratio(n11, n10, n01, n00):
    a, b, c, d = (n11 + HALDANE, n10 + HALDANE, n01 + HALDANE, n00 + HALDANE)
    return float(np.log2((a * d) / (b * c)))


def weighted_table(y, attribute, weights):
    n11 = float(np.dot(weights, (attribute == 1) & (y == 1)))
    n10 = float(np.dot(weights, (attribute == 1) & (y == 0)))
    n01 = float(np.dot(weights, (attribute == 0) & (y == 1)))
    n00 = float(np.dot(weights, (attribute == 0) & (y == 0)))
    return n11, n10, n01, n00


def make_weighted_table(y, attribute):
    """四个布尔掩码只算一次；自助的 5000 个复本里只换权重。"""
    masks = [((attribute == 1) & (y == 1)).astype(float),
             ((attribute == 1) & (y == 0)).astype(float),
             ((attribute == 0) & (y == 1)).astype(float),
             ((attribute == 0) & (y == 0)).astype(float)]
    stacked = np.vstack(masks)

    def table(weights):
        counts = stacked @ weights
        return float(counts[0]), float(counts[1]), float(counts[2]), float(counts[3])

    return table


def make_weighted_auc(y, scores):
    """加权 AUC 的向量化闭包。得分向量在自助里不变，所以排序与同值分块只做一次；
    每个复本只做 reduceat 与 cumsum，没有 Python 层循环。"""
    order = np.argsort(np.asarray(scores, dtype=float), kind="stable")
    y_sorted = np.asarray(y)[order]
    scores_sorted = np.asarray(scores, dtype=float)[order]
    starts = np.flatnonzero(np.concatenate([[True], scores_sorted[1:] != scores_sorted[:-1]]))
    positive_mask = (y_sorted == 1).astype(float)
    negative_mask = (y_sorted == 0).astype(float)

    def auc(weights):
        weights_sorted = np.asarray(weights, dtype=float)[order]
        block_positive = np.add.reduceat(weights_sorted * positive_mask, starts)
        block_negative = np.add.reduceat(weights_sorted * negative_mask, starts)
        total_positive, total_negative = block_positive.sum(), block_negative.sum()
        if total_positive <= 0 or total_negative <= 0:
            return float("nan")
        negative_before = np.concatenate([[0.0], np.cumsum(block_negative)[:-1]])
        concordant = float(np.dot(block_positive, negative_before + 0.5 * block_negative))
        return concordant / (total_positive * total_negative)

    return auc


def make_mh_log_odds_ratio(y, attribute, strata):
    """Mantel-Haenszel 的向量化闭包：层内四格用掩码矩阵一次算出。"""
    values = np.unique(strata)
    blocks = []
    for value in values:
        mask = strata == value
        blocks.append(np.vstack([
            (mask & (attribute == 1) & (y == 1)).astype(float),
            (mask & (attribute == 1) & (y == 0)).astype(float),
            (mask & (attribute == 0) & (y == 1)).astype(float),
            (mask & (attribute == 0) & (y == 0)).astype(float)]))
    stacked = np.vstack(blocks)

    def statistic(weights):
        counts = (stacked @ weights).reshape(len(values), 4)
        n11, n10, n01, n00 = counts[:, 0], counts[:, 1], counts[:, 2], counts[:, 3]
        total = n11 + n10 + n01 + n00
        usable = ((total > 0) & (n11 + n01 > 0) & (n10 + n00 > 0)
                  & (n11 + n10 > 0) & (n01 + n00 > 0))
        if not usable.any():
            return float("nan")
        safe = np.where(usable, total, 1.0)
        numerator = float(np.sum(np.where(usable, n11 * n00 / safe, 0.0)))
        denominator = float(np.sum(np.where(usable, n10 * n01 / safe, 0.0)))
        if numerator <= 0 or denominator <= 0:
            return float("nan")
        return float(np.log2(numerator / denominator))

    return statistic


def mh_log_odds_ratio(y, attribute, strata, weights):
    """Mantel–Haenszel 合并优势比，取 log2。只用两类都在的层。"""
    numerator = denominator = 0.0
    for value in np.unique(strata):
        mask = strata == value
        n11, n10, n01, n00 = weighted_table(y[mask], attribute[mask], weights[mask])
        total = n11 + n10 + n01 + n00
        if total <= 0:
            continue
        if (n11 + n01) <= 0 or (n10 + n00) <= 0 or (n11 + n10) <= 0 or (n01 + n00) <= 0:
            continue
        numerator += n11 * n00 / total
        denominator += n10 * n01 / total
    if numerator <= 0 or denominator <= 0:
        return float("nan")
    return float(np.log2(numerator / denominator))


def cluster_weight_matrix(groups, replicates, rng):
    """对聚类抽多重度，返回 (replicates, n_obs) 的权重；多项式等价于有放回重抽。"""
    unique, inverse = np.unique(groups, return_inverse=True)
    counts = rng.multinomial(len(unique), np.full(len(unique), 1.0 / len(unique)),
                             size=replicates)
    return counts[:, inverse]


def bootstrap_interval(statistic_fn, groups, replicates=REPLICATES, seed=SEED, chunk=250):
    rng = np.random.default_rng(seed)
    draws = []
    remaining = replicates
    while remaining > 0:
        take = min(chunk, remaining)
        weights = cluster_weight_matrix(groups, take, rng)
        for row in weights:
            value = statistic_fn(row.astype(float))
            if np.isfinite(value):
                draws.append(value)
        remaining -= take
    if not draws:
        return [float("nan"), float("nan")], 0
    return [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))], len(draws)


# ------------------------------------------------------------------ 倾向与匹配

def propensity_score(covariates, y):
    scaler = StandardScaler()
    matrix = scaler.fit_transform(covariates)
    model = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs")
    model.fit(matrix, y)
    score = model.predict_proba(matrix)[:, 1]
    score = np.clip(score, 1e-6, 1 - 1e-6)
    return score, np.log(score / (1 - score))


def quantile_strata(values, n_strata=N_STRATA):
    cuts = np.quantile(values, np.linspace(0, 1, n_strata + 1)[1:-1])
    return np.clip(np.searchsorted(cuts, values, side="right"), 0, n_strata - 1)


def nearest_neighbour_match(logit, y, caliper):
    """1:1 最近邻无放回匹配，卡尺内向两侧展开，找到更近的就停。

    返回 (阳性索引, 与之配对的阴性索引)。阳性按倾向得分从高到低处理，
    先给最难匹配的那些找对子——这一处理顺序也写在模块 docstring 的口径里。
    """
    positives = np.flatnonzero(y == 1)
    negatives = np.flatnonzero(y == 0)
    order = np.argsort(logit[negatives], kind="stable")
    pool_index = negatives[order]
    pool_value = logit[pool_index]
    used = np.zeros(len(pool_index), dtype=bool)
    matched_positive, matched_negative = [], []
    for i in positives[np.argsort(-logit[positives], kind="stable")]:
        target = logit[i]
        j = int(np.searchsorted(pool_value, target))
        best, best_distance = -1, np.inf
        left, right = j - 1, j
        while True:
            distance_left = abs(pool_value[left] - target) if left >= 0 else np.inf
            distance_right = abs(pool_value[right] - target) if right < len(pool_value) else np.inf
            nearest = min(distance_left, distance_right)
            if nearest > caliper or nearest > best_distance:
                break
            if distance_left <= caliper and left >= 0:
                if not used[left] and distance_left < best_distance:
                    best, best_distance = left, distance_left
                left -= 1
            if distance_right <= caliper and right < len(pool_value):
                if not used[right] and distance_right < best_distance:
                    best, best_distance = right, distance_right
                right += 1
        if best >= 0:
            used[best] = True
            matched_positive.append(int(i))
            matched_negative.append(int(pool_index[best]))
    return np.asarray(matched_positive, dtype=int), np.asarray(matched_negative, dtype=int)


def random_control_effect(y, attribute, groups, positive_index, n_negatives,
                          replicates=REPLICATES, seed=SEED, draws=RANDOM_DRAWS):
    """同规模随机对照：同一批阳性 + 等量随机阴性，20 次抽样。"""
    negatives = np.flatnonzero(y == 0)
    if n_negatives > len(negatives):
        n_negatives = len(negatives)
    draw_index = []
    for i in range(draws):
        rng = np.random.default_rng(seed + i)
        picked = rng.choice(negatives, size=n_negatives, replace=False)
        draw_index.append(np.concatenate([positive_index, picked]))
    point = []
    for index in draw_index:
        n11, n10, n01, n00 = weighted_table(y[index], attribute[index],
                                            np.ones(len(index)))
        point.append(log_odds_ratio(n11, n10, n01, n00))
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(replicates):
        index = draw_index[int(rng.integers(0, draws))]
        sub_groups = groups[index]
        unique, inverse = np.unique(sub_groups, return_inverse=True)
        counts = rng.multinomial(len(unique), np.full(len(unique), 1.0 / len(unique)))
        weights = counts[inverse].astype(float)
        n11, n10, n01, n00 = weighted_table(y[index], attribute[index], weights)
        values.append(log_odds_ratio(n11, n10, n01, n00))
    return (float(np.mean(point)),
            [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))],
            int(n_negatives), [float(v) for v in point])


def global_permutation_null(y, attribute, seed=SEED, draws=200):
    """前置检查判"组内置换退化"时改用的全局置换零假设核对。

    每次把 y 整体打乱再重算同一个效应量；若估计量本身无结构性偏倚，
    零分布应以 0 为中心。这不是功效检验，只是核对估计量与 Haldane 修正
    在该表规模下不会自己造出一个非零效应。
    """
    rng = np.random.default_rng(seed)
    observed = log_odds_ratio(*weighted_table(y, attribute, np.ones(len(y))))
    values = []
    labels = y.copy()
    for _ in range(draws):
        rng.shuffle(labels)
        values.append(log_odds_ratio(*weighted_table(labels, attribute, np.ones(len(labels)))))
    values = np.asarray(values, dtype=float)
    extreme = float(np.mean(np.abs(values) >= abs(observed)))
    return {"draws": int(draws), "seed": int(seed), "observed": float(observed),
            "null_mean": float(values.mean()),
            "null_interval": [float(np.quantile(values, 0.025)),
                              float(np.quantile(values, 0.975))],
            "permutation_p_two_sided": extreme,
            "null_centred_on_zero": bool(abs(values.mean()) < 0.1)}


# ------------------------------------------------------------------ 判定

def classify(baseline, controlled, controlled_interval, random_interval, blocked):
    if blocked:
        return "undecidable", "前置检查判 BLOCKING，该口径不作主口径"
    if not np.isfinite(controlled) or not np.isfinite(controlled_interval[0]):
        return "undecidable", "控制后效应不可估（层内单类或匹配后无可用列联表）"
    crosses_zero = controlled_interval[0] <= 0 <= controlled_interval[1]
    random_crosses = (not np.isfinite(random_interval[0])
                      or random_interval[0] <= 0 <= random_interval[1])
    same_sign = np.sign(controlled) == np.sign(baseline) or baseline == 0
    if not crosses_zero and not same_sign:
        return "reverses", "控制后符号反向且区间不跨零"
    if not crosses_zero and same_sign:
        if abs(controlled) >= SURVIVE_RETENTION * abs(baseline):
            return "survives", "符号一致、区间不跨零、保留基线效应的一半以上"
        return "attenuated", "符号一致、区间不跨零，但效应不足基线的一半"
    if crosses_zero and not random_crosses:
        return "vanishes", "控制后区间跨零，而同规模随机对照区间不跨零"
    return "undecidable", "控制后与同规模随机对照区间都跨零：在该样本量上不可判定"


DOWNGRADE = {"survives": "attenuated", "vanishes": "undecidable"}


# ------------------------------------------------------------------ 主流程

def run_contingency_claim(cohort):
    claim_id = cohort["claim_id"]
    y = cohort["y"].astype(int)
    attribute = cohort["attribute"].astype(int)
    groups = np.asarray(cohort["groups"])
    covariates = cohort["covariates"]
    record = {"claim_id": claim_id, "unit": cohort["unit"],
              "n_observations": int(len(y)), "n_positive": int(y.sum()),
              "n_attribute": int(attribute.sum()),
              "covariate_set": cohort["covariate_set"],
              "attribute_label": cohort["attribute_label"],
              "author_statistic": cohort["author_statistic"],
              "reproduction_notes": cohort["reproduction_notes"]}

    fast_table = make_weighted_table(y, attribute)
    n11, n10, n01, n00 = fast_table(np.ones(len(y)))
    record["baseline_table"] = {"attr_pos": n11, "attr_neg": n10,
                                "noattr_pos": n01, "noattr_neg": n00}
    baseline = log_odds_ratio(n11, n10, n01, n00)
    record["baseline_log2_or"] = baseline
    record["baseline_interval"], _ = bootstrap_interval(
        lambda w: log_odds_ratio(*fast_table(w)), groups)
    odds, p_value = fisher_exact([[int(n11), int(n10)], [int(n01), int(n00)]])
    record["baseline_fisher_p"] = float(p_value)
    # 作者框架里的可比数字
    with np.errstate(divide="ignore", invalid="ignore"):
        record["baseline_proportion_positive_in_attribute"] = float(n11 / (n11 + n10)) if (n11 + n10) else float("nan")
        record["baseline_fold_enrichment"] = float(
            (n11 / max(1.0, n11 + n01)) / max(1e-12, (n11 + n10) / len(y)))

    score, logit = propensity_score(covariates, y)
    strata = quantile_strata(score)
    fast_auc = make_weighted_auc(y, score)
    record["propensity_auc"] = fast_auc(np.ones(len(y)))
    record["propensity_auc_interval"], _ = bootstrap_interval(fast_auc, groups)

    # ---- 控制一：检出深度分层（倾向得分五分位）
    fast_mh = make_mh_log_odds_ratio(y, attribute, strata)
    mh = fast_mh(np.ones(len(y)))
    mh_interval, _ = bootstrap_interval(fast_mh, groups)
    record["stratified_log2_or"] = mh
    record["stratified_interval"] = mh_interval
    record["strata_detail"] = _strata_detail(y, attribute, strata, score)

    # ---- 论文自报深度的敏感性分层
    depth = cohort.get("observed_depth")
    if depth is not None and np.isfinite(depth).all() and len(np.unique(depth)) > 4:
        depth_strata = quantile_strata(depth, 4)
        fast_depth = make_mh_log_odds_ratio(y, attribute, depth_strata)
        record["observed_depth_log2_or"] = fast_depth(np.ones(len(y)))
        record["observed_depth_interval"], _ = bootstrap_interval(fast_depth, groups)
        record["observed_depth_strata"] = _strata_detail(y, attribute, depth_strata, depth)
    else:
        record["observed_depth_log2_or"] = None
        record["observed_depth_note"] = "论文未提供可用于分层的检出深度"

    # ---- 前置检查（施加匹配前）
    caliper = CALIPER_SD * float(np.std(logit))
    positive_index, negative_index = nearest_neighbour_match(logit, y, caliper)
    matched_index = np.concatenate([positive_index, negative_index]) if len(positive_index) else np.asarray([], dtype=int)
    keep = np.zeros(len(y), dtype=bool)
    keep[matched_index] = True
    n_groups = int(len(np.unique(groups)))
    precheck_groups = groups if n_groups <= GROUPS_FOR_PRECHECK_LIMIT else None
    collinearity = cleaning_collinearity(
        y, keep, covariate=score, strata=strata, groups=precheck_groups,
        label="%s 可检出性匹配规则 对 它要消除的倾向得分" % claim_id)
    collinearity["interval_skipped_for_cost"] = precheck_groups is None
    degeneracy = permutation_degeneracy(y, groups, label="%s 聚类单位" % claim_id)
    record["precheck_cleaning_collinearity"] = collinearity
    record["precheck_permutation_degeneracy"] = degeneracy
    if not degeneracy["informative"]:
        record["global_permutation_null"] = global_permutation_null(y, attribute)
        record["global_permutation_note"] = (
            "组内置换按构造退化（每个蛋白一个观测），按预先声明改用全局置换")
    blocked = bool(collinearity["blocking_reasons"])
    record["match_caliper"] = caliper
    record["n_matched_pairs"] = int(len(positive_index))
    record["matched_positive_retention"] = float(len(positive_index) / max(1, int(y.sum())))

    if len(positive_index) >= 10:
        m_y, m_a = y[matched_index], attribute[matched_index]
        m_groups = groups[matched_index]
        matched_table = make_weighted_table(m_y, m_a)
        matched = log_odds_ratio(*matched_table(np.ones(len(m_y))))
        matched_interval, _ = bootstrap_interval(
            lambda w: log_odds_ratio(*matched_table(w)), m_groups)
        record["matched_log2_or"] = matched
        record["matched_interval"] = matched_interval
        record["matched_propensity_auc"] = make_weighted_auc(
            m_y, score[matched_index])(np.ones(len(m_y)))
        random_point, random_interval, n_random, random_draws = random_control_effect(
            y, attribute, groups, positive_index, len(negative_index))
        record["random_control_log2_or"] = random_point
        record["random_control_interval"] = random_interval
        record["random_control_n_negatives"] = n_random
        record["random_control_draws"] = random_draws
    else:
        matched, matched_interval = float("nan"), [float("nan")] * 2
        random_interval = [float("nan")] * 2
        record["matched_log2_or"] = None
        record["matched_note"] = "匹配上的阳性不足 10 个，匹配口径不可用"
        record["random_control_log2_or"] = None

    primary, primary_reason = classify(baseline, matched, matched_interval,
                                       random_interval, blocked)
    secondary, secondary_reason = classify(baseline, mh, mh_interval,
                                           random_interval, blocked)
    disagreement = primary != secondary
    final = DOWNGRADE.get(primary, primary) if disagreement else primary
    record.update({
        "primary_caliber": "detectability_matching", "primary_verdict": primary,
        "primary_reason": primary_reason, "secondary_caliber": "depth_stratification",
        "secondary_verdict": secondary, "secondary_reason": secondary_reason,
        "caliber_disagreement": bool(disagreement), "verdict": final,
        "precheck_blocked": blocked,
    })
    return record


def _auc(y, scores):
    y = np.asarray(y)
    n1, n0 = int(y.sum()), int(len(y) - y.sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(np.asarray(scores, dtype=float), kind="stable")
    ranks = np.empty(len(y), dtype=float)
    sorted_scores = np.asarray(scores, dtype=float)[order]
    i = 0
    rank_values = np.empty(len(y), dtype=float)
    while i < len(sorted_scores):
        j = i
        while j + 1 < len(sorted_scores) and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        rank_values[i:j + 1] = (i + j) / 2.0 + 1.0
        i = j + 1
    ranks[order] = rank_values
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _weighted_auc(y, scores, weights):
    """加权 AUC = 加权 Mann–Whitney 统计；用排序后前缀和算，避免 O(n^2)。"""
    y = np.asarray(y)
    order = np.argsort(np.asarray(scores, dtype=float), kind="stable")
    ys, ws = y[order], np.asarray(weights, dtype=float)[order]
    scores_sorted = np.asarray(scores, dtype=float)[order]
    positive_weight = ws * (ys == 1)
    negative_weight = ws * (ys == 0)
    total_positive, total_negative = positive_weight.sum(), negative_weight.sum()
    if total_positive <= 0 or total_negative <= 0:
        return float("nan")
    negative_before = np.concatenate([[0.0], np.cumsum(negative_weight)[:-1]])
    # 同值组内按半权重计
    concordant = 0.0
    i = 0
    while i < len(ys):
        j = i
        while j + 1 < len(ys) and scores_sorted[j + 1] == scores_sorted[i]:
            j += 1
        block_positive = positive_weight[i:j + 1].sum()
        block_negative = negative_weight[i:j + 1].sum()
        concordant += block_positive * (negative_before[i] + 0.5 * block_negative)
        i = j + 1
    return float(concordant / (total_positive * total_negative))


def _strata_detail(y, attribute, strata, values):
    rows = []
    for value in np.unique(strata):
        mask = strata == value
        n11, n10, n01, n00 = weighted_table(y[mask], attribute[mask], np.ones(int(mask.sum())))
        rows.append({"stratum": int(value), "n": int(mask.sum()),
                     "positives": int(y[mask].sum()),
                     "attr_pos": int(n11), "attr_neg": int(n10),
                     "noattr_pos": int(n01), "noattr_neg": int(n00),
                     "informative": bool(n11 + n01 > 0 and n10 + n00 > 0
                                         and n11 + n10 > 0 and n01 + n00 > 0),
                     "log2_or": log_odds_ratio(n11, n10, n01, n00),
                     "covariate_mean": float(np.mean(values[mask]))})
    return rows


def run_pers010(cohort):
    """PERS-010：观测共有数 对 独立模型期望共有数，效应量取二者比值的 log2。"""
    keys = cohort["keys"]
    y = cohort["y"].astype(int)
    groups = np.asarray(cohort["groups"])
    n_organs = cohort["n_organs"]
    covariates = cohort["covariates"]
    record = {"claim_id": "PERS-010", "unit": "protein", "n_observations": int(len(y)),
              "n_positive": int(y.sum()), "n_attribute": None,
              "covariate_set": cohort["covariate_set"],
              "attribute_label": cohort["attribute_label"],
              "author_statistic": cohort["author_statistic"],
              "reproduction_notes": cohort["reproduction_notes"],
              "per_organ": cohort["per_organ"], "union_size": cohort["union_size"],
              "shared_size": cohort["shared_size"],
              "declared_shared": cohort["declared_shared"],
              "shared_matches_declared": cohort["shared_matches_declared"]}

    rate = np.asarray([cohort["per_organ"][o] / cohort["union_size"]
                       for o in ("liver", "kidney", "muscle", "brain")], dtype=float)

    def effect(weights, strata=None):
        total = weights.sum()
        observed = float(np.dot(weights, y == 1))
        if strata is None:
            expected = total * float(np.prod(rate))
        else:
            expected = 0.0
            for value in np.unique(strata):
                mask = strata == value
                block = weights[mask].sum()
                if block <= 0:
                    continue
                share = n_organs[mask]
                # 层内各器官检出率用该层平均检出器官数/4 的独立近似
                per_organ = float(np.dot(weights[mask], share) / (4.0 * block))
                expected += block * per_organ ** 4
        if expected <= 0:
            return float("nan")
        return float(np.log2(max(observed, HALDANE) / expected))

    baseline = effect(np.ones(len(y)))
    record["baseline_log2_or"] = baseline
    record["baseline_interval"], _ = bootstrap_interval(lambda w: effect(w), groups)
    record["baseline_shared_fraction"] = float(y.mean())
    record["baseline_expected_shared"] = float(len(y) * np.prod(rate))

    score, logit = propensity_score(covariates, y)
    strata = quantile_strata(score)
    fast_auc = make_weighted_auc(y, score)
    record["propensity_auc"] = fast_auc(np.ones(len(y)))
    record["propensity_auc_interval"], _ = bootstrap_interval(fast_auc, groups)
    record["stratified_log2_or"] = effect(np.ones(len(y)), strata)
    record["stratified_interval"], _ = bootstrap_interval(lambda w: effect(w, strata), groups)
    record["strata_detail"] = [
        {"stratum": int(v), "n": int((strata == v).sum()),
         "shared": int(y[strata == v].sum()),
         "mean_n_organs": float(np.mean(n_organs[strata == v])),
         "propensity_mean": float(np.mean(score[strata == v]))}
        for v in np.unique(strata)]

    depth = cohort["observed_depth"]
    depth_strata = quantile_strata(depth, 4)
    record["observed_depth_log2_or"] = effect(np.ones(len(y)), depth_strata)
    record["observed_depth_interval"], _ = bootstrap_interval(
        lambda w: effect(w, depth_strata), groups)

    caliper = CALIPER_SD * float(np.std(logit))
    positive_index, negative_index = nearest_neighbour_match(logit, y, caliper)
    matched_index = np.concatenate([positive_index, negative_index])
    keep = np.zeros(len(y), dtype=bool)
    keep[matched_index] = True
    n_groups = int(len(np.unique(groups)))
    collinearity = cleaning_collinearity(
        y, keep, covariate=score, strata=strata,
        groups=groups if n_groups <= GROUPS_FOR_PRECHECK_LIMIT else None,
        label="PERS-010 可检出性匹配规则 对 它要消除的倾向得分")
    collinearity["interval_skipped_for_cost"] = n_groups > GROUPS_FOR_PRECHECK_LIMIT
    record["precheck_cleaning_collinearity"] = collinearity
    record["precheck_permutation_degeneracy"] = permutation_degeneracy(
        y, groups, label="PERS-010 聚类单位")
    blocked = bool(collinearity["blocking_reasons"])
    record["precheck_blocked"] = blocked
    record["n_matched_pairs"] = int(len(positive_index))
    record["match_caliper"] = caliper
    record["matched_positive_retention"] = float(len(positive_index) / max(1, int(y.sum())))

    m_y = y[matched_index]
    m_organs = n_organs[matched_index]
    m_groups = groups[matched_index]

    def matched_effect(weights):
        total = weights.sum()
        observed = float(np.dot(weights, m_y == 1))
        per_organ = float(np.dot(weights, m_organs) / (4.0 * total))
        expected = total * per_organ ** 4
        return float(np.log2(max(observed, HALDANE) / expected)) if expected > 0 else float("nan")

    matched = matched_effect(np.ones(len(m_y)))
    matched_interval, _ = bootstrap_interval(matched_effect, m_groups)
    record["matched_log2_or"] = matched
    record["matched_interval"] = matched_interval
    record["matched_propensity_auc"] = make_weighted_auc(
        m_y, score[matched_index])(np.ones(len(m_y)))

    negatives = np.flatnonzero(y == 0)
    draw_index = [np.concatenate([positive_index,
                                  np.random.default_rng(SEED + i).choice(
                                      negatives, size=len(negative_index), replace=False)])
                  for i in range(RANDOM_DRAWS)]
    point = []
    for index in draw_index:
        w = np.ones(len(index))
        observed = float(np.dot(w, y[index] == 1))
        per_organ = float(np.dot(w, n_organs[index]) / (4.0 * w.sum()))
        expected = w.sum() * per_organ ** 4
        point.append(float(np.log2(max(observed, HALDANE) / expected)))
    rng = np.random.default_rng(SEED)
    values = []
    for _ in range(REPLICATES):
        index = draw_index[int(rng.integers(0, RANDOM_DRAWS))]
        sub = groups[index]
        unique, inverse = np.unique(sub, return_inverse=True)
        counts = rng.multinomial(len(unique), np.full(len(unique), 1.0 / len(unique)))
        w = counts[inverse].astype(float)
        observed = float(np.dot(w, y[index] == 1))
        per_organ = float(np.dot(w, n_organs[index]) / (4.0 * w.sum()))
        expected = w.sum() * per_organ ** 4
        if expected > 0:
            values.append(float(np.log2(max(observed, HALDANE) / expected)))
    random_interval = [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))]
    record["random_control_log2_or"] = float(np.mean(point))
    record["random_control_interval"] = random_interval
    record["random_control_n_negatives"] = int(len(negative_index))
    record["random_control_draws"] = [float(v) for v in point]

    primary, primary_reason = classify(baseline, matched, matched_interval,
                                       random_interval, blocked)
    secondary, secondary_reason = classify(baseline, record["stratified_log2_or"],
                                           record["stratified_interval"],
                                           random_interval, blocked)
    disagreement = primary != secondary
    record.update({
        "primary_caliber": "detectability_matching", "primary_verdict": primary,
        "primary_reason": primary_reason, "secondary_caliber": "depth_stratification",
        "secondary_verdict": secondary, "secondary_reason": secondary_reason,
        "caliber_disagreement": bool(disagreement),
        "verdict": DOWNGRADE.get(primary, primary) if disagreement else primary,
    })
    return record


def main():
    started = time.time()
    records, strata_rows = [], []
    for claim_id, builder in cohorts_module.COHORTS.items():
        print("[build]", claim_id, flush=True)
        cohort = builder()
        record = run_pers010(cohort) if claim_id == "PERS-010" else run_contingency_claim(cohort)
        for extra in ("n_positive_set", "n_cluster", "n_cluster_agi", "n_over", "n_under",
                      "declared_sites", "mapped_sites", "unmapped_genes",
                      "position_mismatch", "n_proteins", "transfer_cohort",
                      "n_positive_transfer", "n_observed_unmodified_transfer"):
            if extra in cohort:
                record[extra] = cohort[extra]
        records.append(record)
        for row in record.get("strata_detail", []):
            strata_rows.append(dict(claim_id=claim_id, caliber="propensity_quintile", **row))
        for row in record.get("observed_depth_strata", []) or []:
            strata_rows.append(dict(claim_id=claim_id, caliber="observed_depth_quartile", **row))
        print("   ", claim_id, record["verdict"],
              "baseline=%.4f" % record["baseline_log2_or"],
              "matched=%s" % record.get("matched_log2_or"), flush=True)

    records.append(dict(SFI004_EXCLUSION))

    fieldnames = ["claim_id", "verdict", "primary_verdict", "secondary_verdict",
                  "caliber_disagreement", "precheck_blocked", "unit", "n_observations",
                  "n_positive", "n_attribute", "attribute_label", "covariate_set",
                  "author_statistic_text", "baseline_reproduced",
                  "baseline_log2_or", "baseline_ci_low", "baseline_ci_high",
                  "baseline_fisher_p", "baseline_fold_enrichment",
                  "propensity_auc", "propensity_auc_ci_low", "propensity_auc_ci_high",
                  "stratified_log2_or", "stratified_ci_low", "stratified_ci_high",
                  "observed_depth_log2_or", "observed_depth_ci_low", "observed_depth_ci_high",
                  "matched_log2_or", "matched_ci_low", "matched_ci_high",
                  "n_matched_pairs", "matched_positive_retention", "matched_propensity_auc",
                  "random_control_log2_or", "random_control_ci_low", "random_control_ci_high",
                  "primary_reason", "secondary_reason"]
    table = []
    for record in records:
        author = record.get("author_statistic") or {}
        row = {
            "claim_id": record["claim_id"], "verdict": record["verdict"],
            "primary_verdict": record.get("primary_verdict", ""),
            "secondary_verdict": record.get("secondary_verdict", ""),
            "caliber_disagreement": record.get("caliber_disagreement", ""),
            "precheck_blocked": record.get("precheck_blocked", ""),
            "unit": record.get("unit", ""), "n_observations": record.get("n_observations", ""),
            "n_positive": record.get("n_positive", ""), "n_attribute": record.get("n_attribute", ""),
            "attribute_label": record.get("attribute_label", ""),
            "covariate_set": record.get("covariate_set", ""),
            "author_statistic_text": author.get("text", ""),
            "baseline_reproduced": BASELINE_REPRODUCED[record["claim_id"]],
            "baseline_log2_or": _round(record.get("baseline_log2_or")),
            "baseline_ci_low": _round((record.get("baseline_interval") or [None, None])[0]),
            "baseline_ci_high": _round((record.get("baseline_interval") or [None, None])[1]),
            "baseline_fisher_p": record.get("baseline_fisher_p", ""),
            "baseline_fold_enrichment": _round(record.get("baseline_fold_enrichment")),
            "propensity_auc": _round(record.get("propensity_auc")),
            "propensity_auc_ci_low": _round((record.get("propensity_auc_interval") or [None, None])[0]),
            "propensity_auc_ci_high": _round((record.get("propensity_auc_interval") or [None, None])[1]),
            "stratified_log2_or": _round(record.get("stratified_log2_or")),
            "stratified_ci_low": _round((record.get("stratified_interval") or [None, None])[0]),
            "stratified_ci_high": _round((record.get("stratified_interval") or [None, None])[1]),
            "observed_depth_log2_or": _round(record.get("observed_depth_log2_or")),
            "observed_depth_ci_low": _round((record.get("observed_depth_interval") or [None, None])[0]),
            "observed_depth_ci_high": _round((record.get("observed_depth_interval") or [None, None])[1]),
            "matched_log2_or": _round(record.get("matched_log2_or")),
            "matched_ci_low": _round((record.get("matched_interval") or [None, None])[0]),
            "matched_ci_high": _round((record.get("matched_interval") or [None, None])[1]),
            "n_matched_pairs": record.get("n_matched_pairs", ""),
            "matched_positive_retention": _round(record.get("matched_positive_retention")),
            "matched_propensity_auc": _round(record.get("matched_propensity_auc")),
            "random_control_log2_or": _round(record.get("random_control_log2_or")),
            "random_control_ci_low": _round((record.get("random_control_interval") or [None, None])[0]),
            "random_control_ci_high": _round((record.get("random_control_interval") or [None, None])[1]),
            "primary_reason": record.get("primary_reason", record.get("reason_zh", "")),
            "secondary_reason": record.get("secondary_reason", ""),
        }
        table.append(row)

    import csv as _csv
    with open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(table)
    if strata_rows:
        keys = sorted({k for row in strata_rows for k in row})
        with open(OUT_STRATA, "w", encoding="utf-8-sig", newline="") as handle:
            writer = _csv.DictWriter(handle, fieldnames=keys)
            writer.writeheader()
            writer.writerows(strata_rows)

    verdict_counts = collections.Counter(r["verdict"] for r in records)
    audit = {
        "script": "scripts/run_phase2_claims_under_detectability_control.py",
        "script_sha256": sha256_of(os.path.abspath(__file__)),
        "cohorts_module_sha256": sha256_of(os.path.join(ROOT, "scripts", "phase2_claim_cohorts.py")),
        "precheck_module_sha256": sha256_of(os.path.join(ROOT, "scripts", "audit_cleaning_and_grouping.py")),
        "probe_module_sha256": sha256_of(os.path.join(ROOT, "scripts", "run_cross_protease_detectability_probe.py")),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "python": sys.version.split()[0], "platform": platform.platform(),
        "numpy": np.__version__,
        "predeclared": {
            "estimand": "log2 odds ratio of the 2x2 table (claimed attribute x author's positive label); "
                        "Haldane-Anscombe 0.5 correction for zero cells",
            "pers010_estimand": "log2(observed shared / expected shared under a detection-independence model)",
            "covariates": {"site_level": "VIS10 imported from run_cross_protease_detectability_probe",
                           "protein_level": "PROT8 declared in phase2_claim_cohorts"},
            "propensity_model": "LogisticRegression(C=1.0, max_iter=1000, lbfgs) on standardised covariates",
            "control_1": "detection-depth stratification: %d propensity quintiles, Mantel-Haenszel pooling" % N_STRATA,
            "control_2": "detectability matching: 1:1 nearest neighbour without replacement on the "
                         "propensity logit, caliper %.2f x SD" % CALIPER_SD,
            "random_control": "same-size random control, %d draws, seeds %d..%d; interval integrates "
                              "cluster resampling and draw choice" % (RANDOM_DRAWS, SEED, SEED + RANDOM_DRAWS - 1),
            "bootstrap": {"replicates": REPLICATES, "seed": SEED,
                          "unit": "protein for protein-level claims; protein for site-level claims "
                                  "(sites resampled together)",
                          "does_not_cover": ["GO release and annotation source",
                                             "the authors' own MS search uncertainty",
                                             "alternative covariate definitions"]},
            "verdict_rules": {
                "survives": "sign preserved, CI excludes 0, |controlled| >= %.1f x |baseline|" % SURVIVE_RETENTION,
                "attenuated": "sign preserved, CI excludes 0, |controlled| < %.1f x |baseline|" % SURVIVE_RETENTION,
                "vanishes": "controlled CI crosses 0 AND same-size random control CI excludes 0",
                "reverses": "sign flipped and CI excludes 0",
                "undecidable": "controlled CI crosses 0 AND random control CI crosses 0 (power-limited); "
                               "or a precheck is BLOCKING; or the baseline was not reproducible",
                "primary_caliber": "detectability_matching",
                "secondary_caliber": "depth_stratification",
                "disagreement_rule": "survives->attenuated, vanishes->undecidable",
            },
            "precheck_thresholds": {"ORPHAN_LIMIT": ORPHAN_LIMIT,
                                    "COVARIATE_RISE_LIMIT": COVARIATE_RISE_LIMIT,
                                    "MOVED_RATIO_LIMIT": MOVED_RATIO_LIMIT,
                                    "groups_for_precheck_limit": GROUPS_FOR_PRECHECK_LIMIT},
        },
        "verdict_counts": dict(verdict_counts),
        "baseline_reproduced": BASELINE_REPRODUCED,
        "claims": records,
        "excluded": [SFI004_EXCLUSION],
        "terminology": {"sulfenylation": "次磺酰化 (-SOH)", "sulfinylation": "亚磺酰化 (-SO2H)"},
    }
    with open(OUT_AUDIT, "w", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2, default=_default)
    print("wrote", os.path.relpath(OUT_CSV, ROOT), os.path.relpath(OUT_AUDIT, ROOT))
    print("verdicts", dict(verdict_counts))


def _round(value, digits=4):
    if value is None or isinstance(value, str):
        return "" if value is None else value
    try:
        if not np.isfinite(value):
            return ""
    except TypeError:
        return value
    return round(float(value), digits)


def _default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        value = float(obj)
        return value if np.isfinite(value) else None
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    if isinstance(obj, set):
        return sorted(obj)
    raise TypeError(repr(obj)[:120])


if __name__ == "__main__":
    main()
