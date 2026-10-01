#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第二阶段补完：B 层十四条已发表论断在可检出性控制下重做。

第二阶段（`scripts/run_phase2_claims_under_detectability_control.py`）做了 A 层八条，
得存活 5 / 消失 1 / 反转 0 / 不可判定 1 / 口径不适用 1。本脚本补完
`results/ptm_site_preference_claims.csv` 里 `retestable_now=True` 的另外 **14 条 B 层**。

**沿用第二阶段已预声明的口径，本脚本不另立任何估计量、协变量或判定阈值。**
统计机制（log2 优势比与 Haldane–Anscombe 0.5 校正、倾向得分逻辑回归、
五分位 Mantel–Haenszel 分层、1:1 卡尺匹配、同规模随机对照、5000 次聚类自助
种子 20260915、两条前置检查）全部**直接 import** 那个脚本的函数，
一行都不重写；本文件只做三件事：调度、B 层特有的两处预声明补充、产物落盘。

本文件里没有任何模型被称为 persulfidation site predictor；这里不拟合排序模型。
术语：sulfenylation = 次磺酰化（-SOH），sulfinylation = 亚磺酰化（-SO2H）。

================================================================================
一 B 层特有的两处预声明补充（运行前写定）
================================================================================
第二阶段的八条论断全是 `claim_direction == positive_preference`，且都有作者自报的
效应量。B 层不是，因此需要补两条读法，**补的是读法，不是估计量**：

**补充 1 阴性论断的判定映射。** B 层有一条 `null_no_preference`（PERS-002，
Longen 2016 说 CXXC/CXC 不是持硫化的优先靶点）。作者主张的是"没有偏好"，
所以第二阶段那套"符号一致且保留一半以上"的规则按定义套不上去。映射写定为：

  * 基线区间**不**跨零  -> `undecidable`：作者的零结果在本装配里就不成立，
    控制后的结果没有解释力（与第二阶段 PERS-006 因基线不显著而判不可判定同型）。
  * 基线区间跨零、匹配后区间也跨零 -> `null_survives`：零结果在控制后仍是零结果。
    **阴性论断只能被确认到某个分辨率**，因此这一栏必须同时报匹配后区间的半宽
    （`null_resolution_log2_or`），不得只写"存活"。
  * 基线区间跨零、匹配后区间不跨零 -> `null_broken_by_control`：控制反而显出了
    作者说不存在的偏好。
  合并计数时 `null_survives` 计入"存活"一类，`null_broken_by_control` 计入"反转"一类，
  并在报告里单列，不混同于正向论断。

**补充 2 检出深度分层的层定义。** PERS-002/PERS-003 的 qPerS-SID 有四个供体臂，
论文自报深度只有 1..4 四个取值，四分位在四个取值上退化，因此这两条的
"论文自报深度分层"按**臂计数原值分 4 层**，不取四分位。其余论断的队列没有论文
自报深度，该列留空。

**补充 3 基线方向与论断方向的一致性检查。** 第二阶段的 `classify()` 比的是
"控制后的符号是否与**基线**符号一致"，它隐含假设基线的方向就是论断主张的方向——
A 层八条恰好都满足，所以那一轮不需要这一条。B 层不满足：转移队列里可能出现
基线符号与论断主张**相反**且区间不跨零的情形，此时 `classify()` 会给出 "survives"，
但它说的是"基线那个效应在控制后还在"，**不是"论断成立"**。因此补一条判读：

  * `claim_direction == positive_preference` 而基线 log2 优势比 < 0 且基线区间不跨零
    -> 论断的判定记 `baseline_contradicts_claim`，**单列一类，不计入存活/消失/反转/
    不可判定**；`classify()` 原本给出的那个判定另存一列
    （`verdict_by_phase2_rule`），因为它仍然回答"基线效应能否经受控制"这个问题。
  * 基线区间跨零（根本没有基线效应）时不触发本条，沿用第二阶段判定。

这一条是**第一次跑完之后发现的缺陷修正**，不是按结果调口径：判读只看"基线符号与
论断主张的方向是否一致"，与是哪条论断、结果好不好无关，换任何一批论断都同样适用。
更正前后的判定都留档（报告 §7）。

================================================================================
二 每个匹配分析都配同规模随机对照（一条不省）
================================================================================
第二阶段的规则原样执行：随机对照 = 保留同一批被匹配上的阳性、从阴性里等量随机抽
20 次（种子 20260915+i），区间由 5000 次自助给出。**这是判决"消失"与"功效不足"
的唯一依据**，因此凡是跑了匹配的论断都有这一列；没有这一列的论断只有两类：
口径不适用与材料不可达（下面逐条列名）。

声明的次属性（PERS-003 的 ±10 窗口、PERS-009 的 {K,R} 子集、SFE-007 的
"紧邻缺正电"那一半、SFE-011 的糖酵解、SNO-016 的 ±10 窗口）用**同一个**倾向得分与
同一组匹配对（倾向模型是 y ~ 协变量，与属性无关，所以匹配对不变），
只换 2x2 的行定义，因此次属性的结果与主属性逐位可比。

================================================================================
三 判 out_of_instrument_scope 的五条，逐条说明（不硬套）
================================================================================
任务书要求：若某条论断的结构使既有估计量不适用，逐条说明并判
`out_of_instrument_scope`。下面五条逐条给的是**结构性**理由——这些理由与材料
是否取得无关，取到材料也不会改变判定。详见 OUT_OF_SCOPE。

================================================================================
四 材料不可达的两条
================================================================================
PERS-007 与 PERS-008 都需要 Li 2024 RSC Chem Biol 的 ESI（锌指型与富集条目），
本轮**所有外网出口不可用**（`scripts/fetch_phase2b_claim_supplements.py` 的日志里
七次请求全部 SSL 握手失败，含 pypi.org，即不是域名放行问题），
且本树没有可替代的锌指持硫化队列可作转移检验。记 `blocked_material_unreachable`，
**既不算"不成立"也不算"不可判定"**，单独一类。

================================================================================
五 不做的事
================================================================================
* 不覆盖第二阶段的任何产物（本轮所有文件名带 phase2b 前缀）。
* 不改 results/v2_blocker_log.csv 与 results/v2_requirement_status.csv，
  不跑 build_artifact_manifest.py，不动 protocols/v3_qtrp_analysis_plan.json，
  不删除任何盘上文件。
* 付费墙材料不取、不缓存、不分发（Fu 2020 ARS 的在手补充表本脚本不读取内容用于分发）。

产物（文件名独立）：
  results/phase2b_claim_retest.csv              B 层 14 条
  results/phase2b_claim_retest_combined.csv     与第二阶段 8 条合并共 22 条，标注轮次
  results/phase2b_claim_retest_audit.json       脚本 sha256、预声明读法、前置检查判决
  results/phase2b_claim_retest_strata.csv       分层明细

用法：
    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=scripts \\
        /Users/lyuguohao/opt/persulfidation_v2/venv/bin/python \\
        scripts/run_phase2b_claims_backfill.py
"""
from __future__ import annotations

import collections
import csv
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np

import phase2b_claim_cohorts as cohorts_module
from audit_cleaning_and_grouping import (
    COVARIATE_RISE_LIMIT, MOVED_RATIO_LIMIT, ORPHAN_LIMIT,
)
from run_phase2_claims_under_detectability_control import (
    CALIPER_SD, HALDANE, N_STRATA, RANDOM_DRAWS, REPLICATES, SEED,
    SURVIVE_RETENTION, _round, _default, _strata_detail, bootstrap_interval,
    log_odds_ratio, make_mh_log_odds_ratio, make_weighted_table,
    run_contingency_claim, sha256_of, weighted_table,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
OUT_CSV = os.path.join(RESULTS, "phase2b_claim_retest.csv")
OUT_COMBINED = os.path.join(RESULTS, "phase2b_claim_retest_combined.csv")
OUT_AUDIT = os.path.join(RESULTS, "phase2b_claim_retest_audit.json")
OUT_STRATA = os.path.join(RESULTS, "phase2b_claim_retest_strata.csv")
PHASE2_CSV = os.path.join(RESULTS, "phase2_claim_retest.csv")
PHASE2_AUDIT = os.path.join(RESULTS, "phase2_claim_retest_audit.json")
FETCH_LOG = os.path.join(ROOT, "external", "intake", "phase2b_claim_supp",
                         "phase2b_supp_fetch_log.json")

# 基线复现的三类判读，运行前写定。
# B 层与 A 层最大的结构差别在这里：A 层八条每条都有作者自报的效应量（比例、倍数、
# 条目计数），B 层十四条里有七条 stat_type 就是 none_reported——论文根本没给这条论断的
# 效应量。因此"复现基线"必须分成三档，不能只报 True/False：
#   reproduced   作者自报的效应量被重算出来且方向与量级一致
#   scale_only   作者没给这条论断的效应量（none_reported / count_only），
#                只能核对队列规模与作者自报的位点/蛋白数，且已核对一致
#   not_reproduced  材料不可达（转移检验），或重算值与作者值不一致
BASELINE_REPRODUCTION = {
    "PERS-002": "scale_only",      # 作者未给检验；队列 724 阳性 对 论文 742 肽 / 540 蛋白
    "PERS-003": "scale_only",      # 同上，同一队列
    "PERS-009": "scale_only",      # 作者只给定性保守性图谱，论断行 n 留空
    "SNO-016": "not_reproduced",   # 作者数据本轮不可达，转移检验
    "SFE-007": "not_reproduced",   # 作者阳性集不可达（与 SFE-006 同源），转移检验
    "SFE-008": "not_reproduced",   # 同上，且谷胱甘肽化那一路无队列
    "SFE-011": "not_reproduced",   # 作者位点表与油菜蛋白组不可达，转移检验
    "PERS-001": "not_reproduced",  # 付费墙摘要级，基序内容不在可得材料里
    "SFI-005": "not_reproduced",   # 位点级氢键伙伴不在补充材料
    "SFE-009": "not_reproduced",   # 材料本轮不可达；且论断无可复现的效应量
    "SFE-010": "not_reproduced",   # 同上
    "SNO-017": "not_reproduced",   # 同上
    "PERS-007": "not_reproduced",  # ESI 不可达
    "PERS-008": "not_reproduced",  # ESI 不可达
}

# 结构性不适用的五条。理由与材料是否取得无关。
OUT_OF_SCOPE = [
    {
        "claim_id": "PERS-001", "verdict": "out_of_instrument_scope",
        "reason_zh": "被主张的属性在任何可得材料里都不存在：Fu 2020 是付费墙论文，"
                     "本树只读到摘要级信息，摘要只写 \"unique consensus motifs\"，"
                     "没有给出基序内容、频率或位置。预声明的估计量要求一个二值属性，"
                     "而这条论断的属性无法判读。位点表在手（付费墙，不分发）也不解决这一点——"
                     "从在手位点表自己推出一个基序再去检验它，检的是我们推的基序，"
                     "不是作者的论断，且属性会由数据定出来，违反'属性运行前写定'。",
        "what_would_unlock_it": "作者方法节或补充材料里写明的基序（需要全文访问权）",
        "counts_from_paper": {"sites": 1547, "proteins": 994},
    },
    {
        "claim_id": "SFI-005", "verdict": "out_of_instrument_scope",
        "reason_zh": "观测通道不匹配，与第二阶段 SFI-004 判定同源同理："
                     "该论断来自 PDB 结构普查，混淆通道是晶体可解析性，"
                     "不是质谱肽段可见性，预声明的 VIS10/PROT8 按定义描述不了它。"
                     "此外属性本身（氧化态 Cys 的氢键伙伴是 Thr 还是 Ser、Arg 还是 Lys）"
                     "需要逐个结构的氢键几何，补充材料 bi2c00349_si_001.pdf 只有图 S1–S4，"
                     "本环境既无氢键判定引擎也未取回 1124 个结构。",
        "what_would_unlock_it": "另造一把尺（按分辨率、B 因子、占有率匹配）+ 结构与氢键几何引擎",
        "counts_from_paper": {"CSD": 469, "CYS": 7103},
    },
    {
        "claim_id": "SFE-009", "verdict": "out_of_instrument_scope",
        "reason_zh": "论断的结构是'一组特征对判别有内在作用'（预测器特征有效性），"
                     "不是'某个具名属性在阳性集里富集'。预声明的估计量是一张 2x2 表的 "
                     "log2 优势比，需要一个二值属性；这条论断没有点名任何残基或位置，"
                     "十四个 AAIndex 理化属性 + PSAAP 是一个连续特征块。"
                     "把它硬塞进 2x2 只能靠自己造阈值，那是另立口径。"
                     "（另记：Xu 2016 的补充材料本轮也取不到，但这不是判定理由。）",
        "what_would_unlock_it": "一个判别力口径的估计量（匹配前后 ΔAUC），并为它单独写预声明",
        "counts_from_paper": {"positive_peptides": 900, "negative_peptides": 6856},
    },
    {
        "claim_id": "SFE-010", "verdict": "out_of_instrument_scope",
        "reason_zh": "与 SFE-009 同型：论断是'邻域理化指数可预测次磺酰化位点'，"
                     "即 2203 个特征的整体判别力，没有点名属性，"
                     "论断表里 n_positive 本身就因为正文只给特征数而留空。"
                     "预声明的 2x2 优势比对它不适用。"
                     "（另记：Al-Barakati 2018 的补充信息本轮也取不到，但这不是判定理由。）",
        "what_would_unlock_it": "同 SFE-009，需要判别力口径的估计量",
        "counts_from_paper": {"features": 2203, "physicochemical_properties": 14},
    },
    {
        "claim_id": "SNO-017", "verdict": "out_of_instrument_scope",
        "reason_zh": "被主张的属性是连续量且作者没给阈值：'−10 至 −2 区段相邻残基的 SASA "
                     "偏好高于非 SNO 位点'，正文未给数值与检验，也未给判'高'的切点。"
                     "要进预声明的 2x2 就得自己造一个 SASA 切点——SFE-001 与 SNO-002 那类"
                     "论断的作者自带切点（RSA>25%、RSA>10%），这一条没有。"
                     "此外位点级 SASA 需要结构与 SASA 引擎，本环境不具备"
                     "（与第二阶段 SFI-004 同一限制）。",
        "what_would_unlock_it": "连续属性口径的估计量（匹配后均值差）+ SASA 引擎",
        "counts_from_paper": {},
    },
]

# 材料不可达的两条。
BLOCKED_MATERIAL = [
    {
        "claim_id": "PERS-008", "verdict": "blocked_material_unreachable",
        "reason_zh": "需要 Li 2024 RSC Chem Biol 的 ESI Table S3（持硫化锌指蛋白名单）"
                     "才能装配阳性集，本轮所有外网出口不可用"
                     "（fetch 日志里七次请求全部 SSL 握手失败，含 pypi.org，"
                     "因此不是域名放行问题，也未尝试镜像站或伪装 User-Agent）。"
                     "属性侧还需要 GO:0016567/GO:0006511 一类条目的后代集，"
                     "而 QuickGO 同样不可达，只用直接注释会重犯第二阶段已更正的那个错误。"
                     "本树没有可替代的锌指持硫化队列，因此也没有转移检验的路。",
        "what_would_unlock_it": "外网出口恢复后重跑 scripts/fetch_phase2b_claim_supplements.py",
        "counts_from_paper": {},
    },
    {
        "claim_id": "PERS-007", "verdict": "blocked_material_unreachable",
        "reason_zh": "同一篇的 ESI（十个已发表持硫化数据集的锌指型汇编）不可达，"
                     "且属性（CCCC / CCCH / CCHC 配位型）需要逐蛋白的锌指型判读，"
                     "该判读只在作者的 ESI 里，本环境无法重建。"
                     "论断本身还是跨物种对比，即使拿到 ESI 也需要按物种分别装配。",
        "what_would_unlock_it": "同 PERS-008；拿到 ESI 后按物种内 2x2 逐个装配",
        "counts_from_paper": {"datasets_compiled": 10},
    },
]

NULL_VERDICT_BUCKET = {"null_survives": "survives", "null_broken_by_control": "reverses"}


def crosses_zero(interval):
    if interval is None or not np.isfinite(interval[0]) or not np.isfinite(interval[1]):
        return True
    return interval[0] <= 0 <= interval[1]


def secondary_attribute_effect(cohort, record, attribute):
    """次属性：换 2x2 的行定义，倾向得分与匹配对完全不变（倾向模型与属性无关）。"""
    y = cohort["y"].astype(int)
    groups = np.asarray(cohort["groups"])
    table = make_weighted_table(y, attribute)
    baseline = log_odds_ratio(*table(np.ones(len(y))))
    baseline_interval, _ = bootstrap_interval(lambda w: log_odds_ratio(*table(w)), groups)
    out = {"attribute_label": cohort.get("attribute_secondary_label", ""),
           "n_attribute": int(attribute.sum()),
           "baseline_log2_or": baseline, "baseline_interval": baseline_interval}
    index = record.get("_matched_index")
    if index is not None and len(index):
        matched_table = make_weighted_table(y[index], attribute[index])
        matched = log_odds_ratio(*matched_table(np.ones(len(index))))
        matched_interval, _ = bootstrap_interval(
            lambda w: log_odds_ratio(*matched_table(w)), groups[index])
        out["matched_log2_or"] = matched
        out["matched_interval"] = matched_interval
        out["retention_ratio"] = (abs(matched) / abs(baseline)) if baseline else None
    return out


def declared_depth_stratification(cohort, record):
    """补充 2：按供体臂计数原值分层的 Mantel–Haenszel 合并。"""
    depth = cohort.get("observed_depth_declared")
    if depth is None:
        return None
    y = cohort["y"].astype(int)
    attribute = cohort["attribute"].astype(int)
    groups = np.asarray(cohort["groups"])
    strata = np.asarray(depth, dtype=int)
    statistic = make_mh_log_odds_ratio(y, attribute, strata)
    point = statistic(np.ones(len(y)))
    interval, _ = bootstrap_interval(statistic, groups)
    return {"label": cohort.get("observed_depth_label", ""),
            "strata_definition": "供体臂计数原值（1..4），不取四分位",
            "log2_or": point, "interval": interval,
            "detail": _strata_detail(y, attribute, strata, np.asarray(depth, dtype=float))}


def classify_null(baseline_interval, matched_interval):
    if not crosses_zero(baseline_interval):
        return ("undecidable",
                "阴性论断读法：基线区间不跨零，作者的零结果在本装配里就不成立，"
                "控制后的结果没有解释力")
    if crosses_zero(matched_interval):
        return ("null_survives",
                "阴性论断读法：基线与匹配后区间都跨零，零结果在控制后仍是零结果"
                "（只能确认到该区间半宽的分辨率）")
    return ("null_broken_by_control",
            "阴性论断读法：基线区间跨零但匹配后不跨零，控制显出了作者说不存在的偏好")


def run_one(claim_id, builder):
    cohort = builder()
    record = run_contingency_claim(cohort)
    record["is_transfer"] = bool(cohort.get("is_transfer"))
    record["claim_direction"] = cohort.get("claim_direction", "positive_preference")
    record["transfer_cohort"] = cohort.get("transfer_cohort", "")
    for extra in ("n_positive_assembled", "n_observed_unmodified", "n_proteins",
                  "off_residue_rows", "missing_accessions", "n_positive_transfer",
                  "n_observed_unmodified_transfer", "n_sulfenyl_only", "n_sno_only",
                  "n_overlap_dropped", "n_positive_loci"):
        if extra in cohort:
            record[extra] = cohort[extra]

    # 次属性用同一批匹配对；重跑匹配索引（run_contingency_claim 没有把它带出来，
    # 但匹配只依赖倾向得分与卡尺，两者都是确定性的，所以重算给出同一组对子）
    from run_phase2_claims_under_detectability_control import (
        nearest_neighbour_match, propensity_score,
    )
    score, logit = propensity_score(cohort["covariates"], cohort["y"].astype(int))
    caliper = CALIPER_SD * float(np.std(logit))
    positive_index, negative_index = nearest_neighbour_match(
        logit, cohort["y"].astype(int), caliper)
    assert int(len(positive_index)) == int(record["n_matched_pairs"]), (
        "匹配对数与 run_contingency_claim 不一致：%s" % claim_id)
    record["_matched_index"] = (np.concatenate([positive_index, negative_index])
                               if len(positive_index) else np.asarray([], dtype=int))

    if cohort.get("attribute_secondary") is not None:
        record["secondary_attribute"] = secondary_attribute_effect(
            cohort, record, np.asarray(cohort["attribute_secondary"]).astype(int))

    declared_depth = declared_depth_stratification(cohort, record)
    if declared_depth is not None:
        record["declared_depth_stratification"] = declared_depth
        record["observed_depth_log2_or"] = declared_depth["log2_or"]
        record["observed_depth_interval"] = declared_depth["interval"]
        record.pop("observed_depth_note", None)

    baseline = float(record["baseline_log2_or"])
    baseline_interval = record.get("baseline_interval")
    if crosses_zero(baseline_interval):
        record["baseline_direction_agrees"] = "no_baseline_effect"
    elif record["claim_direction"] == "null_no_preference":
        record["baseline_direction_agrees"] = "no"
    else:
        record["baseline_direction_agrees"] = "yes" if baseline > 0 else "no"

    if record["claim_direction"] == "null_no_preference":
        verdict, reason = classify_null(baseline_interval, record.get("matched_interval"))
        record["verdict_by_phase2_rule"] = record["verdict"]
        record["verdict"] = verdict
        record["primary_reason"] = reason
        interval = record.get("matched_interval") or [float("nan")] * 2
        record["null_resolution_log2_or"] = (
            float((interval[1] - interval[0]) / 2.0)
            if np.isfinite(interval[0]) and np.isfinite(interval[1]) else None)
    elif record["baseline_direction_agrees"] == "no":
        # 补充 3：基线符号与论断主张的方向相反且区间不跨零
        record["verdict_by_phase2_rule"] = record["verdict"]
        record["verdict"] = "baseline_contradicts_claim"
        record["primary_reason"] = (
            "补充 3：论断主张正向偏好，而本队列的基线 log2 优势比为 %.4f、"
            "区间 [%.4f, %.4f] 不跨零，方向相反——第二阶段规则给出的 '%s' 说的是"
            "'基线那个（相反方向的）效应在控制后还在'，不是论断成立，因此单列一类"
            % (baseline, baseline_interval[0], baseline_interval[1],
               record["verdict_by_phase2_rule"]))

    record["baseline_reproduction"] = BASELINE_REPRODUCTION[claim_id]
    record["_cohort_notes"] = cohort["reproduction_notes"]
    record.pop("_matched_index", None)
    return record, cohort


FIELDNAMES = ["claim_id", "round", "verdict", "verdict_bucket", "claim_direction",
              "baseline_direction_agrees", "verdict_by_phase2_rule",
              "is_transfer", "transfer_cohort", "primary_verdict", "secondary_verdict",
              "caliber_disagreement", "precheck_blocked", "unit", "n_observations",
              "n_positive", "n_attribute", "attribute_label", "covariate_set",
              "author_statistic_text", "baseline_reproduction",
              "baseline_log2_or", "baseline_ci_low", "baseline_ci_high",
              "baseline_fisher_p", "baseline_fold_enrichment",
              "propensity_auc", "propensity_auc_ci_low", "propensity_auc_ci_high",
              "stratified_log2_or", "stratified_ci_low", "stratified_ci_high",
              "observed_depth_log2_or", "observed_depth_ci_low", "observed_depth_ci_high",
              "matched_log2_or", "matched_ci_low", "matched_ci_high",
              "n_matched_pairs", "matched_positive_retention", "matched_propensity_auc",
              "random_control_log2_or", "random_control_ci_low", "random_control_ci_high",
              "retention_ratio", "null_resolution_log2_or",
              "secondary_attribute_label", "secondary_attribute_matched_log2_or",
              "primary_reason", "secondary_reason"]


def bucket_of(verdict):
    if verdict in NULL_VERDICT_BUCKET:
        return NULL_VERDICT_BUCKET[verdict]
    return verdict


def row_of(record, round_tag):
    author = record.get("author_statistic") or {}
    baseline = record.get("baseline_log2_or")
    matched = record.get("matched_log2_or")
    retention = None
    if baseline not in (None, 0) and matched is not None:
        try:
            retention = abs(float(matched)) / abs(float(baseline))
        except (TypeError, ZeroDivisionError):
            retention = None
    secondary = record.get("secondary_attribute") or {}
    return {
        "claim_id": record["claim_id"], "round": round_tag,
        "verdict": record["verdict"], "verdict_bucket": bucket_of(record["verdict"]),
        "claim_direction": record.get("claim_direction", ""),
        "baseline_direction_agrees": record.get("baseline_direction_agrees", ""),
        "verdict_by_phase2_rule": record.get("verdict_by_phase2_rule", ""),
        "is_transfer": record.get("is_transfer", ""),
        "transfer_cohort": record.get("transfer_cohort", ""),
        "primary_verdict": record.get("primary_verdict", ""),
        "secondary_verdict": record.get("secondary_verdict", ""),
        "caliber_disagreement": record.get("caliber_disagreement", ""),
        "precheck_blocked": record.get("precheck_blocked", ""),
        "unit": record.get("unit", ""), "n_observations": record.get("n_observations", ""),
        "n_positive": record.get("n_positive", ""),
        "n_attribute": record.get("n_attribute", ""),
        "attribute_label": record.get("attribute_label", ""),
        "covariate_set": record.get("covariate_set", ""),
        "author_statistic_text": author.get("text", ""),
        "baseline_reproduction": record.get("baseline_reproduction", ""),
        "baseline_log2_or": _round(baseline),
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
        "matched_log2_or": _round(matched),
        "matched_ci_low": _round((record.get("matched_interval") or [None, None])[0]),
        "matched_ci_high": _round((record.get("matched_interval") or [None, None])[1]),
        "n_matched_pairs": record.get("n_matched_pairs", ""),
        "matched_positive_retention": _round(record.get("matched_positive_retention")),
        "matched_propensity_auc": _round(record.get("matched_propensity_auc")),
        "random_control_log2_or": _round(record.get("random_control_log2_or")),
        "random_control_ci_low": _round((record.get("random_control_interval") or [None, None])[0]),
        "random_control_ci_high": _round((record.get("random_control_interval") or [None, None])[1]),
        "retention_ratio": _round(retention, 3),
        "null_resolution_log2_or": _round(record.get("null_resolution_log2_or")),
        "secondary_attribute_label": secondary.get("attribute_label", ""),
        "secondary_attribute_matched_log2_or": _round(secondary.get("matched_log2_or")),
        "primary_reason": record.get("primary_reason", record.get("reason_zh", "")),
        "secondary_reason": record.get("secondary_reason", ""),
    }


def main():
    started = time.time()
    records, strata_rows = [], []
    for claim_id, builder in cohorts_module.COHORTS.items():
        print("[build]", claim_id, flush=True)
        record, cohort = run_one(claim_id, builder)
        records.append(record)
        for row in record.get("strata_detail", []):
            strata_rows.append(dict(claim_id=claim_id, caliber="propensity_quintile", **row))
        depth_block = record.get("declared_depth_stratification")
        if depth_block:
            for row in depth_block["detail"]:
                strata_rows.append(dict(claim_id=claim_id, caliber="donor_arm_count", **row))
        print("   ", claim_id, record["verdict"],
              "baseline=%.4f" % record["baseline_log2_or"],
              "matched=%s" % record.get("matched_log2_or"),
              "random=%s" % record.get("random_control_log2_or"), flush=True)

    for entry in OUT_OF_SCOPE + BLOCKED_MATERIAL:
        item = dict(entry)
        item["baseline_reproduction"] = BASELINE_REPRODUCTION[item["claim_id"]]
        item["primary_reason"] = item["reason_zh"]
        records.append(item)

    table = [row_of(record, "phase2b") for record in records]
    assert len(table) == 14, "B 层必须恰好 14 行，实得 %d" % len(table)

    with open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(table)

    # 合并表：第二阶段 8 行按原样带过来，缺的列留空，round 标 phase2
    with open(PHASE2_CSV, encoding="utf-8-sig", newline="") as handle:
        phase2_rows = list(csv.DictReader(handle))
    phase2_audit = json.load(open(PHASE2_AUDIT, encoding="utf-8"))
    phase2_direction = {"PERS-005": "positive_preference", "PERS-006": "positive_preference",
                        "PERS-010": "positive_preference", "PERS-011": "positive_preference",
                        "PERS-012": "positive_preference", "SNO-021": "positive_preference",
                        "SFE-006": "positive_preference", "SFI-004": "positive_preference"}
    combined = []
    for source in phase2_rows:
        row = {key: "" for key in FIELDNAMES}
        for key, value in source.items():
            if key in row:
                row[key] = value
        row["claim_id"] = source["claim_id"]
        row["round"] = "phase2"
        row["verdict_bucket"] = bucket_of(source["verdict"])
        row["claim_direction"] = phase2_direction[source["claim_id"]]
        row["is_transfer"] = source["claim_id"] == "SFE-006"
        row["transfer_cohort"] = ("fps2020_ath_sulfenyl"
                                  if source["claim_id"] == "SFE-006" else "")
        # 第二阶段的 baseline_reproduced 是 True/False 布尔；映射进本轮的三档，
        # A 层每条都有作者自报效应量，所以 True -> reproduced，False -> not_reproduced
        row["baseline_reproduction"] = ("reproduced"
                                        if source["baseline_reproduced"] == "True"
                                        else "not_reproduced")
        # 补充 3 同样施加到第二阶段八行（从盘上的区间重算，不假设它们都一致）
        try:
            point = float(source["baseline_log2_or"])
            low, high = float(source["baseline_ci_low"]), float(source["baseline_ci_high"])
            row["baseline_direction_agrees"] = (
                "no_baseline_effect" if low <= 0 <= high else ("yes" if point > 0 else "no"))
        except ValueError:
            row["baseline_direction_agrees"] = ""
        try:
            if source["matched_log2_or"] and source["baseline_log2_or"]:
                row["retention_ratio"] = _round(
                    abs(float(source["matched_log2_or"]))
                    / abs(float(source["baseline_log2_or"])), 3)
        except ValueError:
            pass
        combined.append(row)
    combined.extend(table)
    assert len(combined) == 22, "合并表必须恰好 22 行，实得 %d" % len(combined)
    with open(OUT_COMBINED, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(combined)

    if strata_rows:
        keys = sorted({k for row in strata_rows for k in row})
        with open(OUT_STRATA, "w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=keys)
            writer.writeheader()
            writer.writerows(strata_rows)

    counts_b = collections.Counter(r["verdict"] for r in records)
    counts_bucket_b = collections.Counter(bucket_of(r["verdict"]) for r in records)
    counts_all = collections.Counter(row["verdict"] for row in combined)
    counts_bucket_all = collections.Counter(row["verdict_bucket"] for row in combined)
    reproduction_b = collections.Counter(BASELINE_REPRODUCTION.values())
    reproduction_all = collections.Counter(row["baseline_reproduction"] for row in combined)

    fetch_log = json.load(open(FETCH_LOG, encoding="utf-8")) if os.path.exists(FETCH_LOG) else {}
    audit = {
        "script": "scripts/run_phase2b_claims_backfill.py",
        "script_sha256": sha256_of(os.path.abspath(__file__)),
        "cohorts_module_sha256": sha256_of(os.path.join(ROOT, "scripts", "phase2b_claim_cohorts.py")),
        "phase2_driver_sha256": sha256_of(os.path.join(
            ROOT, "scripts", "run_phase2_claims_under_detectability_control.py")),
        "phase2_cohorts_sha256": sha256_of(os.path.join(ROOT, "scripts", "phase2_claim_cohorts.py")),
        "precheck_module_sha256": sha256_of(os.path.join(ROOT, "scripts", "audit_cleaning_and_grouping.py")),
        "probe_module_sha256": sha256_of(os.path.join(
            ROOT, "scripts", "run_cross_protease_detectability_probe.py")),
        "fetch_module_sha256": sha256_of(os.path.join(
            ROOT, "scripts", "fetch_phase2b_claim_supplements.py")),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "python": sys.version.split()[0], "platform": platform.platform(),
        "numpy": np.__version__,
        "inherits_predeclared_from": {
            "script": "scripts/run_phase2_claims_under_detectability_control.py",
            "sections": "一 统一效应量 / 二 两种控制 / 三 前置检查 / 四 区间 / 五 三类判定",
            "audit": "results/phase2_claim_retest_audit.json",
            "phase2_predeclared": phase2_audit.get("predeclared"),
        },
        "phase2b_additional_predeclared": {
            "null_claim_mapping": {
                "undecidable": "baseline CI excludes 0 (the author's null does not reproduce here)",
                "null_survives": "baseline CI crosses 0 AND matched CI crosses 0; "
                                 "must be quoted together with null_resolution_log2_or",
                "null_broken_by_control": "baseline CI crosses 0 but matched CI excludes 0",
                "bucketing": NULL_VERDICT_BUCKET,
            },
            "declared_depth_strata": "PERS-002/PERS-003: donor-arm count 1..4 as raw strata, "
                                     "not quartiles (quartiles degenerate on four values)",
            "baseline_direction_consistency": {
                "rule": "a positive_preference claim whose baseline log2 OR is negative with an "
                        "interval excluding zero is recorded as baseline_contradicts_claim and is "
                        "NOT counted as survives/vanishes/reverses/undecidable; the Phase 2 rule's "
                        "own verdict is kept in verdict_by_phase2_rule",
                "why": "Phase 2's classify() compares the controlled sign against the BASELINE "
                       "sign, which silently assumes the baseline direction is the claimed one; "
                       "all eight A-tier baselines happened to satisfy that, B-tier transfers do not",
                "found": "after the first B-tier run; the rule is direction-only and claim-agnostic",
                "also_applied_to_phase2_rows": "recomputed from the on-disk Phase 2 intervals, "
                                               "not assumed",
            },
            "flank_attribute_convention": {
                "window": cohorts_module.FLANK_WINDOW,
                "fraction": cohorts_module.FLANK_FRACTION,
                "min_hits": cohorts_module.FLANK_MIN_HITS,
                "window_sensitivity": cohorts_module.FLANK_WINDOW_SENS,
                "min_hits_sensitivity": cohorts_module.FLANK_MIN_HITS_SENS,
                "hydrophobic_set": sorted(cohorts_module.HYDROPHOBIC),
                "applies_to": ["PERS-003", "PERS-009", "SNO-016"],
                "note": "claims whose authors named offsets (PERS-002, SFE-007, SFE-008) "
                        "use the named offsets instead",
            },
            "baseline_reproduction_categories": {
                "reproduced": "the author's own statistic was recomputed and agrees in direction and magnitude",
                "scale_only": "the paper reports no effect size for this claim "
                              "(stat_type none_reported / count_only); only the cohort scale was checked",
                "not_reproduced": "material unreachable (transfer test) or the recomputed value disagrees",
            },
        },
        "precheck_thresholds": {"ORPHAN_LIMIT": ORPHAN_LIMIT,
                                "COVARIATE_RISE_LIMIT": COVARIATE_RISE_LIMIT,
                                "MOVED_RATIO_LIMIT": MOVED_RATIO_LIMIT},
        "bootstrap": {"replicates": REPLICATES, "seed": SEED, "n_strata": N_STRATA,
                      "caliper_sd": CALIPER_SD, "random_draws": RANDOM_DRAWS,
                      "haldane": HALDANE, "survive_retention": SURVIVE_RETENTION},
        "verdict_counts_phase2b": dict(counts_b),
        "verdict_bucket_counts_phase2b": dict(counts_bucket_b),
        "verdict_counts_combined_22": dict(counts_all),
        "verdict_bucket_counts_combined_22": dict(counts_bucket_all),
        "baseline_reproduction_phase2b": dict(reproduction_b),
        "baseline_reproduction_combined_22": dict(reproduction_all),
        "baseline_reproduction": BASELINE_REPRODUCTION,
        "network_status_this_round": {
            "verdict": "no outbound HTTPS egress",
            "evidence": "every request in the fetch log failed at the TLS handshake "
                        "(SSL UNEXPECTED_EOF_WHILE_READING), including pypi.org, "
                        "so this is not a domain-allowlist decision",
            "fetch_log": os.path.relpath(FETCH_LOG, ROOT),
            "supplements_attempted": [e.get("key") for e in fetch_log.get("supplements", [])],
            "proteomes_attempted": [e.get("key") for e in fetch_log.get("proteomes", [])],
            "no_user_agent_spoofing": True, "no_mirrors_or_caches": True,
        },
        "claims": records,
        "out_of_instrument_scope": OUT_OF_SCOPE,
        "blocked_material_unreachable": BLOCKED_MATERIAL,
        "not_touched": ["results/v2_blocker_log.csv", "results/v2_requirement_status.csv",
                        "protocols/v3_qtrp_analysis_plan.json",
                        "scripts/build_artifact_manifest.py",
                        "results/phase2_claim_retest.csv (read only)"],
        "terminology": {"sulfenylation": "次磺酰化 (-SOH)", "sulfinylation": "亚磺酰化 (-SO2H)"},
    }
    with open(OUT_AUDIT, "w", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2, default=_default)

    print("wrote", os.path.relpath(OUT_CSV, ROOT), os.path.relpath(OUT_COMBINED, ROOT),
          os.path.relpath(OUT_AUDIT, ROOT))
    print("phase2b verdicts", dict(counts_b))
    print("combined buckets", dict(counts_bucket_all))
    print("baseline reproduction (22)", dict(reproduction_all))


if __name__ == "__main__":
    main()
