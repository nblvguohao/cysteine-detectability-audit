# -*- coding: utf-8 -*-
import csv, json
from pathlib import Path
ROOT = Path("/Users/lyuguohao/Documents/Codex/https-chatgpt-com-s-t-6aa28ddd9e60819199c57fc1065a53b8/outputs/2026-09-12_persulfidation_upgrade")
rows = list(csv.DictReader((ROOT/"results/phase2e_claim_retest.csv").open(encoding="utf-8-sig")))
a = json.loads((ROOT/"results/phase2e_claim_retest_audit.json").read_text(encoding="utf-8"))
fa = json.loads((ROOT/"results/structural_site_features_audit.json").read_text(encoding="utf-8"))
cov = json.loads((ROOT/"results/phase2e_coverage.json").read_text(encoding="utf-8"))
R = {r["claim_id"] + "|" + r["caliber"]: r for r in rows}
def n(k, f, d=4):
    v = R[k][f]
    try: return f"{float(v):.{d}f}"
    except (TypeError, ValueError): return str(v)
def pc(k, f):
    return f"{float(R[k][f]):.1%}"
def iv(k, p):
    return f"[{n(k, p+'_ci_low')}, {n(k, p+'_ci_high')}]"
_mauc = [float(r["matched_propensity_auc"]) for r in rows]
_mret = [float(r["matched_positive_retention"]) for r in rows]
_deg = [r["claim_id"] + "|" + r["caliber"] for r in rows
        if "DEGENER" in r["precheck_permutation_verdict"].upper()]

TEXT = f"""# 第二阶段 e 轮：五条结构论断，覆盖 {cov['coverage_before_2e']}/37 → {cov['tested_in_retestable']}/37

日期 2026-09-16。两个脚本：`scripts/build_structural_site_features.py`
（sha256 {fa['script_sha256'][:16]}…，{fa['elapsed_minutes']} 分钟，沙箱 python）算结构特征，
`scripts/run_phase2e_claims.py`（sha256 {a['script_sha256'][:16]}…，{a['elapsed_minutes']} 分钟，项目 venv）跑再检验。
驱动逻辑**从 d 轮导入而非复制**（`run_phase2d_claims.run_one`，sha256 {a['driver_sha256'][:16]}…），
估计量一字未改。产物 `results/structural_site_features.csv`（{fa['sites_total']} 行）、
`results/phase2e_claim_retest.csv`（10 行）、两份审计、`results/phase2e_coverage.json`。

## 一 本轮的发现：一条论断能被再检验到什么程度，取决于它的**属性定义写得多精确**

这不是事先设想的结论。本轮五条论断的属性都来自同一张结构表、同一套仪器、同一控制，
差别只在**作者当初把属性定义写到了什么程度**。结果分成清楚的三档：

| 论断 | 作者写明了什么 | 组成核对（阳性 / 阴性 对作者报的） | 判定 |
|---|---|---|---|
| `SFE-001` | **量与阈值都写明**（相对可及面积 > 25%） | **{pc('SFE-001|primary','attribute_rate_positive')} / {pc('SFE-001|primary','attribute_rate_negative')}** 对 >60% / <30% | **存活** |
| `SNO-012` | 阈值写明但是**绝对面积**（硫原子 ≤1.0 Å²） | {pc('SNO-012|primary','attribute_rate_positive')} 埋藏 对 ~35% | 控制显出效应（口径相关） |
| `SNO-002` | **阈值没写** | {pc('SNO-002|primary','attribute_rate_positive')} / {pc('SNO-002|primary','attribute_rate_negative')} 对 29% / 23% | 不可判定 |
| `SNO-001` | 类别写明但工具（DSSP）不可复现 | {pc('SNO-001|primary','attribute_rate_positive')} / {pc('SNO-001|primary','attribute_rate_negative')} 对 40% / 29% | 不可判定 |
| `SNO-009` | 窗口写明、统计量是逐位频率 | {pc('SNO-009|primary','attribute_rate_positive')} / {pc('SNO-009|primary','attribute_rate_negative')} 对"升高 10 个百分点" | 不可判定 |

**这是给文章的第三根轴**：第一阶段说阴性集的构造决定你测到什么，d 轮说富集类论断的约束是背景选择，
这一轮说**属性定义的精确度决定论断能否被再检验**。三根轴都是可量化的，都不靠断言。

## 二 `SFE-001`：本项目最扎实的一条

作者用 NetSurfP 序列预测器，本轮用 AlphaFold 模型加 Shrake-Rupley 实算——**完全不同的方法**，
而组成落在 {pc('SFE-001|primary','attribute_rate_positive')} 对 {pc('SFE-001|primary','attribute_rate_negative')}，
作者报的是 >60% 对 <30%。阈值沿用作者明写的 25%，没有调过。

| 口径 | 基线 | 区间 | 匹配后 | 区间 | 同规模随机对照 | 保留比 |
|---|---:|---|---:|---|---:|---:|
| 主口径（{R['SFE-001|primary']['n_observations']} 观测 / {R['SFE-001|primary']['n_positive']} 阳性） | {n('SFE-001|primary','baseline_log2_or')} | {iv('SFE-001|primary','baseline')} | {n('SFE-001|primary','matched_log2_or')} | {iv('SFE-001|primary','matched')} | {n('SFE-001|primary','random_control_log2_or')} | **{n('SFE-001|primary','retention_ratio',3)}** |
| pLDDT≥70 敏感性 | {n('SFE-001|plddt70','baseline_log2_or')} | {iv('SFE-001|plddt70','baseline')} | {n('SFE-001|plddt70','matched_log2_or')} | {iv('SFE-001|plddt70','matched')} | {n('SFE-001|plddt70','random_control_log2_or')} | {n('SFE-001|plddt70','retention_ratio',3)} |

**组成用异方法复现、效应在可检出性控制下存活、敏感性口径同向**——三件都成立的只有这一条。
可以在文中作为"论断被正确建立起来时长什么样"的正面例子，与 `SFE-006` 的减半、
`PERS-008` 的背景假象形成对照。

## 三 `SNO-012`：控制把作者的"略有暴露富集"放大了，但判定标签随口径变

作者说约 35% 的 NO-Cys 硫原子埋藏、整体略有暴露富集。本轮 {pc('SNO-012|primary','attribute_rate_positive')} 埋藏
（高于作者的 35%，与声明过的 Shrake-Rupley 绝对面积差异一致），
阴性侧 {pc('SNO-012|primary','attribute_rate_negative')}——**NO-Cys 比未修饰 Cys 更不埋藏**。

- 主口径：基线 {n('SNO-012|primary','baseline_log2_or')} {iv('SNO-012|primary','baseline')} 跨零（作者的"略有"复现为零结果），
  匹配后 {n('SNO-012|primary','matched_log2_or')} {iv('SNO-012|primary','matched')} **不跨零**，
  同规模随机对照 {n('SNO-012|primary','random_control_log2_or')} {iv('SNO-012|primary','random_control')} 跨零 → `null_broken_by_control`
- pLDDT≥70：基线 {n('SNO-012|plddt70','baseline_log2_or')} {iv('SNO-012|plddt70','baseline')} **本身就不跨零** → 判 `undecidable`（作者的零结果在该口径下不成立）

**两个口径的判定标签不同，实质方向完全一致**：NO-Cys 更暴露，控制后更明显。
标签之差只因基线是否跨零，必须这样写，不能只报主口径那个更漂亮的标签。

## 四 Doulias 那三条：不可判定，且两条的组成对不上

`SNO-002`（暴露）：**作者没写埋藏/暴露的阈值**，借用 Yang 的 25% 后组成是
{pc('SNO-002|primary','attribute_rate_positive')} 对 {pc('SNO-002|primary','attribute_rate_negative')}，
而作者报 29% 对 23%——绝对水平只有一半，差值也从 6 个百分点缩到
{abs(float(R['SNO-002|primary']['attribute_rate_positive'])-float(R['SNO-002|primary']['attribute_rate_negative']))*100:.1f} 个。
按运行前写定的读法，这说明**该阈值不可复原**，结论只能是不可判定，不能说论断不成立。

`SNO-001`（α 螺旋过表达）：{pc('SNO-001|primary','attribute_rate_positive')} 对 {pc('SNO-001|primary','attribute_rate_negative')}，
方向与作者一致（40% 对 29%）但幅度小得多，基线 {n('SNO-001|primary','baseline_log2_or')} {iv('SNO-001|primary','baseline')} 跨零。
它的次属性（coil 欠表达）**方向与论断相反**：作者称 coil 欠表达，本轮得基线
{n('SNO-001|primary','secondary_baseline_log2_or')}（阳性侧 coil 略多），控制后 {n('SNO-001|primary','secondary_matched_log2_or')}。
pLDDT≥70 口径下次属性同向（基线 {n('SNO-001|plddt70','secondary_baseline_log2_or')}，
控制后 {n('SNO-001|plddt70','secondary_matched_log2_or')}）。
P-SEA 与 DSSP 的指派规则不同，这一条的替代成本最高。

`SNO-009`（0 至 +3 窗内 coil 升高）：{pc('SNO-009|primary','attribute_rate_positive')} 对 {pc('SNO-009|primary','attribute_rate_negative')}，
即**下降 {abs(float(R['SNO-009|primary']['attribute_rate_positive'])-float(R['SNO-009|primary']['attribute_rate_negative']))*100:.1f} 个百分点**，
而作者报升高 10 个百分点。基线 {n('SNO-009|primary','baseline_log2_or')} {iv('SNO-009|primary','baseline')}。
这一条在 P-SEA 下完全不复现。

## 五 全部十行

| 论断 | 口径 | 判定 | 观测 | 阳性 | 基线 | 基线 95% 区间 | 匹配后 | 匹配后 95% 区间 | 随机对照 | 保留比 | 倾向 AUC | 阳性率(+/−) |
|---|---|---|---:|---:|---:|---|---:|---|---:|---:|---:|---|
"""
ORDER = [f"{c}|{k}" for c in ("SFE-001","SNO-002","SNO-012","SNO-001","SNO-009")
         for k in ("primary","plddt70")]
for k in ORDER:
    r = R[k]
    TEXT += (f"| `{r['claim_id']}` | {r['caliber']} | {r['verdict']} | {r['n_observations']} | "
             f"{r['n_positive']} | {n(k,'baseline_log2_or')} | {iv(k,'baseline')} | "
             f"{n(k,'matched_log2_or')} | {iv(k,'matched')} | "
             f"{n(k,'random_control_log2_or')} | {n(k,'retention_ratio',3)} | "
             f"{n(k,'propensity_auc',4)} | {pc(k,'attribute_rate_positive')} / {pc(k,'attribute_rate_negative')} |\n")

TEXT += f"""
匹配质量（从全表十行算出）：匹配后倾向 AUC {min(_mauc):.4f}–{max(_mauc):.4f}，
阳性保留 {min(_mret):.3f}–{max(_mret):.3f}；匹配前置检查十行全部 `usable_as_primary_caliber`。
组内置换退化 {len(_deg)} 行{("：" + "、".join("`"+d+"`" for d in _deg)) if _deg else ""}。

## 六 结构特征这一步的口径与代价

{fa['counts']['proteins']} 个蛋白取 AlphaFold 模型（{fa['source']}），
**{fa['counts']['model_missing']} 个没有模型、{fa['counts']['length_mismatch']} 个序列长度与 UniProt 不符**，
这些蛋白的全部 Cys（含阳性）被剔除。得 {fa['sites_total']} 个 Cys 位点，
其中 pLDDT≥70 的 {fa['sites_plddt_ge_70']} 个。
`SFE-001` 主口径剔除 {R['SFE-001|primary']['n_sites_dropped_no_structure']} 个位点
（阳性 {R['SFE-001|primary']['n_positives_lost']} 个），Doulias 队列一个没丢。

SASA 用 **biotite 的 Shrake-Rupley**（`{fa['sasa_method'][:60]}…`）——freesasa 只有源码包、
本沙箱无法编译；两种算法的绝对面积略有差异，对 `SNO-012` 的绝对阈值有影响，已随该条声明。
二级结构用 **biotite 的 P-SEA**（本沙箱没有 DSSP 二进制），三态输出与 Doulias 报的一致但指派规则不同。
相对可及面积除以 {fa['thresholds']['cys_max_asa_A2']} Å²（{fa['thresholds']['cys_max_asa_source']}）。

pLDDT 的处理方式是运行前定的：**主口径不按 pLDDT 过滤**，因为按置信度过滤等于对有序性做条件，
而有序性与被测的可及性相关；pLDDT≥70 单列为声明过的敏感性口径。

## 七 覆盖与仍然做不了的

覆盖 **{cov['coverage_before_2e']}/{cov['retestable']} → {cov['tested_in_retestable']}/{cov['retestable']}**
（本轮 {cov['added_this_round']} 条，其中 {cov['added_in_retestable']} 条在可再检验集内；
`SNO-012` 在该集之外，是额外收获）。

仍未检验 {len(cov['still_untested_retestable'])} 条：

| 原因 | 条数 | 论断 |
|---|---:|---|
| pKa 无预测器 | {len(cov['still_untested_by_reason']['pKa_no_predictor'])} | {', '.join('`'+c+'`' for c in cov['still_untested_by_reason']['pKa_no_predictor'])} |
| 固有反应性无预测器 | {len(cov['still_untested_by_reason']['intrinsic_reactivity_no_predictor'])} | {', '.join('`'+c+'`' for c in cov['still_untested_by_reason']['intrinsic_reactivity_no_predictor'])} |
| d 轮判为不可再检验 | {len(cov['still_untested_by_reason']['judged_not_testable_2d'])} | {', '.join('`'+c+'`' for c in cov['still_untested_by_reason']['judged_not_testable_2d'])} |
| 阴性集未定义或不在仪器范围 | {len(cov['still_untested_by_reason']['placeholder_only_negative_set_undefined_or_out_of_scope'])} | {', '.join('`'+c+'`' for c in cov['still_untested_by_reason']['placeholder_only_negative_set_undefined_or_out_of_scope'])} |

另有 Marino 的 pKa 类论断 `SNO-010`、`SNO-011`、`SNO-015` 本就在可再检验集之外，同样卡在 pKa。
**pKa 一共卡住 5 条**（集内 2 条、集外 3 条），本轮按纪律没有用替代量顶替——
用别的量会改变检验的是什么。

## 八 本轮未做

未算 pKa、保守度、固有反应性。未重跑任何既有论断。
AlphaFold 模型存在会话临时目录、**未进本树**；实际用到的每个模型的版本号与 sha256 记在
`results/structural_site_features_audit.json` 的 `models_used` 里（{len(fa['models_used'])} 条）。
"""
p = ROOT/"reports"/"PHASE2E_STRUCTURAL_CLAIMS.md"
p.write_text(TEXT, encoding="utf-8")
print("wrote", p.name, len(TEXT))
