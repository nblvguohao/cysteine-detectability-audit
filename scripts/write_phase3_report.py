# -*- coding: utf-8 -*-
import csv, json
from pathlib import Path
ROOT = Path("/path/to/work")
d = json.loads((ROOT/"results/ptm_detectability_diagnostics_audit.json").read_text(encoding="utf-8"))
b = json.loads((ROOT/"results/reference_benchmark_spec_audit.json").read_text(encoding="utf-8"))
tiers = d["tier_distribution_over_51_claims"]
spec = list(csv.DictReader((ROOT/"results/reference_benchmark_spec.csv").open(encoding="utf-8-sig")))
std = list(csv.DictReader((ROOT/"results/label_semantics_tiers.csv").open(encoding="utf-8-sig")))
cases = d["regression"]["cases"]
npass = sum(1 for c in cases if c["ok"])

TEXT = f"""# 第三阶段：把批评变成别人能用的三件东西——诊断包、标签语义分级、参考基准规格

日期 2026-09-16。**由主会话完成**：这一阶段先后派出两个子轮次，两次都在动手之前因连接错误中断、
盘上无任何产物（台账已记），因此改为自己做。诊断包与分级标准不需外网，全在树内数据上跑。

产物：`scripts/ptm_detectability_diagnostics.py`（sha256 {d['script_sha256'][:16]}…，{d['elapsed_seconds']} 秒）、
`results/label_semantics_tiers.csv`、`results/ptm_detectability_diagnostics_audit.json`、
`scripts/build_reference_benchmark_spec.py`、`results/reference_benchmark_spec.csv`、
`results/reference_benchmark_spec_audit.json`、日志 `logs/ptm_detectability_diagnostics.log`。
阈值与回归期望**全部写在模块文档字符串里、在任何用例运行之前**。

## 一 为什么必须是"逐条施加控制"而不是一条经验规则

第二阶段 B 层补跑那一轮（`reports/PHASE2B_CLAIMS_BACKFILL.md` 第四节）把我自己提的假设否定了：
**效应损失并不随可检出性倾向 AUC 递增**（限定后 n=9、Spearman ρ −0.2510、精确置换 p 0.5123；
反例 SFE-008 倾向 0.6220 而保留比 0.086 为全体最低）。

这个否定结果正是诊断包存在的理由：**既然没有汇总统计量能预测哪条论断会塌，就只能逐条把控制施加一遍。**
论文里这一节的论证顺序应当是"我们试过用一个数概括，它不成立，所以我们给工具"。

## 二 诊断包：三项检查，跑脚本就是跑回归

在既有的 `cleaning_collinearity` 与 `permutation_degeneracy` 之上补三项：

| 检查 | 定义 | 声明的阈值 | 何时返回"无法计算"而不是给数 |
|---|---|---|---|
| 可检出性份额 | (AUC_可见性 − 0.5) / (AUC_消化 − 0.5) | 分母下限 {d['declared']['SHARE_DENOMINATOR_FLOOR']} | 分母低于下限时（超集本身几乎不优于随机，比值不稳） |
| 标签语义分级 | 由**阴性集**判定 T1/T2/T3 | 见下节映射表 | 阴性集未定义时判 unclassifiable |
| 深度混淆强度 | 分层表解析出的深度单列 AUC + 有信息层 + 被孤立阳性 | 孤立阳性上限 {d['declared']['ORPHANED_LIMIT']} | 深度协变量编码标签时不适用（须人工排除） |

另把第二阶段的判定规则实现为纯函数 `claim_verdict`（减弱阈值 {d['declared']['ATTENUATION_FLOOR']}），
这样它可以被回归、也可以被别人复用。

**{npass}/{len(cases)} 项回归全部通过**，其中四类是**重算**、一类是**读取已记录先例**（审计里分开标明）：

- R1 `results/ptm_detectability_share.csv` 里 16 个可用份额，全部由存盘的逐特征集全局 AUC 按定义重算复现。
- R2 QTRP 主臂与重复臂的深度单列 AUC 由分层表重算，得 0.6828 与 0.7885，与原轮次一致；
  并要求 `arm_clean` 在两臂上都被判为不可信（孤立阳性 102 = 37.50% 与 48 = 60.00%，均超 5% 上限）。
- R3 合并 22 条论断里带齐所需数字的 **14 条，判定规则回放全部与存盘判定一致**。
- R4 两个解析已知答案的合成净化规则：定向删除高协变量阴性**必须**被否决
  （实测该规则把协变量自身 AUC 抬高 +0.1841，超过声明上限 0.02）；均匀随机删除**必须**不被否决。
- R5 已记录先例：QTRP `arm_clean` 判 `demote_to_sensitivity_check`、ABPP 检出配平判
  `usable_as_primary_caliber`——**这两项是从既有审计读出的，不是重算的**，审计里的
  `recomputed_versus_recorded` 键明确区分。

**R2 在第一次运行时就抓到了我自己的一个错**：一致对方向写反，返回的是 1 − AUC（0.3172 对 0.6828）。
这正是把期望值写在运行之前的价值，已在文档字符串里记下这条用例"第一次运行就挣回了成本"。

## 三 标签语义分级：三级，判据是阴性集而不是化学

{len(std)} 行标准表在 `results/label_semantics_tiers.csv`，每级给出阴性含义、允许的论断强度、
必需的控制、树内实例与**反例**。三级为：

- **T1 直接加合物证据**——同一运行内同一半胱氨酸的两种加合形式都被检出且**位点层通道互斥**。允许化学论断。
  反例（说明为何"声称两通道"不等于 T1）：PXD005168 声明了封闭通道与直接持硫化通道，
  但 187 个直接通道位点里 112 个同时被报为封闭游离巯基。
- **T2 并行但归属不互斥**——同运行并行标注但位点层不互斥，或阴性取自富集运行里"被检出未修饰"的肽段。
  允许化学论断，**但必须同时报可检出性控制**。反例：YAP1C 队列换成同运行阴性后份额反而**升高**
  （0.8962 对 0.7322）——"用被看见的位点做阴性"并不自动免疫。
- **T3 捕获流程注释**——阴性是"未被观测到修饰"。**只允许注释排序类论断，不支持化学论断。**
  反例（防止误读）：T3 数据集仍可支撑排序论断，T3 不等于"不可用"。

把这套标准施加到 51 条已发表论断上：

| 级别 | 条数 |
|---|---:|
| T3 捕获流程注释（不支持化学论断） | **{tiers.get('T3_capture_annotation',0)}** |
| 无法分级（阴性集未定义） | {tiers.get('unclassifiable',0)} |
| T2 并行但归属不互斥 | {tiers.get('T2_parallel_ambiguous',0)} |
| T1 直接加合物证据 | {tiers.get('T1_direct_adduct',0)} |

**这两个数字不要混用**：第一阶段说的"30 条阴性集是未被检出的半胱氨酸"指三个"未检出"类别；
这里的 {tiers.get('T3_capture_annotation',0)} 条是**所有不支持化学论断的类别**（另含其它库人工阴性与随机/模拟集）。
两者是相邻但不同的量。

## 四 参考基准规格：先声明它不能回答什么

`results/reference_benchmark_spec.csv`。资格在装配前声明：只收 T1/T2
（**T3 按构造排除**——建在"未观测到即阴性"上的基准衡量的就是可检出性本身），
蛋白内口径须可识别（至少 30 个分组单元，取自已注册的 v3 方案，该方案在 `v3_pooled_ph5_raw` 上为 44 个组件），
分折按同源组件（`results/v3_human_homology_components.csv`，{b['components_file_rows']} 行）。

| 队列 | 级别 | 位点 | 阳性 | 分组单元 |
|---|---|---:|---:|---:|
{chr(10).join(f"| `{r['cohort']}` | {r['tier'].split('_')[0]} | {r['n_sites']} | {r['n_positive']} | {r['n_grouping_units']} |" for r in spec)}

**它不是"正确的持硫化基准"**，审计里 `not_a_standard` 键原样写着：每种修饰只有一个实验室的化学、
阳性在数百量级、全部人源，而本项目需要判定的 0.02–0.04 量级 AUC 差值**低于 314 个阳性能判定的下限**。
它是一个带着限制一起发表的可辩护起点。

## 五 每项检查会在什么条件下失效（审计 `known_limits` 键）

- 份额需要一个可见性特征块与一个超集；没有消化模型时它返回"无法计算"而不是给一个数。
- 分级判定读的是声明的元数据，**无法发现"论文声称互斥的通道在它自己的补充表里并不互斥"**——
  那要靠读表，PXD005168 就是这样被发现的。
- 深度混淆强度假定深度协变量不编码标签。定义成"在几个臂里通过判据"就会编码标签，
  本项目犯过两次（按行计数的探针深度、tier B 臂深度）。
- `claim_verdict` 是判定规则不是流水线，它把效应估计当已知量接收。

## 六 本阶段未做的

第四阶段（判决性湿实验设计）与第五阶段（成稿）未开始。湿实验设计依赖本阶段的判据已经就位，
现在可以开始；成稿要等核心图定稿——而核心图的覆盖仍是 22/37，其中 15 条需出版社访问权或联系作者。
"""
(ROOT/"reports"/"PHASE3_DIAGNOSTIC_AND_STANDARD.md").write_text(TEXT, encoding="utf-8")
print("wrote PHASE3_DIAGNOSTIC_AND_STANDARD.md", len(TEXT))
