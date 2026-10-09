#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第二阶段主图：已发表论断在可检出性控制前后的效应量对照。

只读 `results/phase2_claim_retest.csv` 与 `results/phase2_claim_retest_audit.json`，
不重算任何统计量——图里的每个数都能在那两个文件里逐个找到。

左图 每条论断一行：空心点 = 基线 log2 优势比（论文自己的口径），
实心点 = 可检出性匹配后的 log2 优势比，两点间连线表示控制带来的移动；
灰色叉 = 同规模随机对照（用来把"消失"和"不可判定"分开）；横线 = 95% 区间。
颜色按三类判定，零线是语义零点（优势比 1）。

右图 倾向得分 AUC：只用预先声明的可检出性协变量（位点层 VIS10、蛋白层 PROT8）
把论文的阳性集与其阴性/背景集分开的能力。0.5 是无信息线。这一栏是混淆的大小，
不是模型性能，**不得读成任何 predictor 的表现**。

用法：
    PYTHONDONTWRITEBYTECODE=1 python3 scripts/plot_phase2_claim_retest.py
产物：results/phase2_claim_retest_fig.png / .pdf
"""
from __future__ import annotations

import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
IN_CSV = os.path.join(RESULTS, "phase2_claim_retest.csv")
IN_AUDIT = os.path.join(RESULTS, "phase2_claim_retest_audit.json")
OUT_PNG = os.path.join(RESULTS, "phase2_claim_retest_fig.png")
OUT_PDF = os.path.join(RESULTS, "phase2_claim_retest_fig.pdf")

VERDICT_ZH = {"survives": "存活", "attenuated": "衰减存活", "vanishes": "消失",
              "reverses": "反转", "undecidable": "不可判定",
              "out_of_instrument_scope": "口径不适用"}
VERDICT_COLOR = {"survives": "#1b6ca8", "attenuated": "#5aa9d6", "vanishes": "#b8860b",
                 "reverses": "#8b2f8b", "undecidable": "#8a8f98",
                 "out_of_instrument_scope": "#c7cbd1"}
ORDER = ["survives", "attenuated", "vanishes", "reverses", "undecidable",
         "out_of_instrument_scope"]
SHORT = {
    "PERS-005": "PERS-005 糖酵解通路过半被检出\n(Aroca 2017, 拟南芥)",
    "PERS-006": "PERS-006 过表达者多属初级代谢\n(Aroca 2017, 拟南芥)",
    "PERS-010": "PERS-010 四器官共有核心持硫化组\n(Bithi 2021, 小鼠)",
    "PERS-011": "PERS-011 富集于类囊体\n(Jurado-Flores 2023, 拟南芥)",
    "PERS-012": "PERS-012 富集于羧酸代谢\n(Matamoros 2024, 菜豆)",
    "SNO-021": "SNO-021 基序以赖氨酸为核心\n(Wang 2023, 弓形虫)",
    "SFE-006": "SFE-006 侧翼富集 K/R（转移检验）\n(Bui 2016 → 拟南芥次磺酰化队列)",
    "SFI-004": "SFI-004 氧化态 Cys 更暴露\n(Garrido Ruiz 2022, PDB)",
}


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def load():
    with open(IN_CSV, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    audit = json.load(open(IN_AUDIT, encoding="utf-8"))
    rows.sort(key=lambda r: (ORDER.index(r["verdict"]) if r["verdict"] in ORDER else 99,
                             r["claim_id"]))
    return rows, audit


def main():
    rows, audit = load()
    try:
        from kernel import apply_figure_style  # noqa: F401
    except Exception:  # noqa: BLE001
        pass
    for family in ("PingFang SC", "Heiti SC", "Arial Unicode MS", "Songti SC"):
        if family in {f.name for f in matplotlib.font_manager.fontManager.ttflist}:
            plt.rcParams["font.family"] = family
            break
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                         "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "legend.fontsize": 7, "axes.unicode_minus": False})

    n = len(rows)
    figure = plt.figure(figsize=(9.6, 0.86 * n + 2.6))
    grid = figure.add_gridspec(1, 2, width_ratios=[2.45, 1.0], wspace=0.06)
    left = figure.add_subplot(grid[0, 0])
    right = figure.add_subplot(grid[0, 1])

    positions = np.arange(n)[::-1]
    for position, row in zip(positions, rows):
        verdict = row["verdict"]
        colour = VERDICT_COLOR.get(verdict, "#8a8f98")
        baseline = _float(row["baseline_log2_or"])
        matched = _float(row["matched_log2_or"])
        b_low, b_high = _float(row["baseline_ci_low"]), _float(row["baseline_ci_high"])
        m_low, m_high = _float(row["matched_ci_low"]), _float(row["matched_ci_high"])
        random_point = _float(row["random_control_log2_or"])

        if np.isfinite(baseline) and np.isfinite(matched):
            left.annotate("", xy=(matched, position), xytext=(baseline, position),
                          arrowprops=dict(arrowstyle="-|>", color=colour, lw=1.0,
                                          alpha=0.55, shrinkA=3.2, shrinkB=3.2))
        if np.isfinite(b_low):
            left.plot([b_low, b_high], [position + 0.16] * 2, color="#5b6068", lw=1.0,
                      solid_capstyle="butt", zorder=2)
        if np.isfinite(baseline):
            left.plot([baseline], [position + 0.16], marker="o", mfc="white",
                      mec="#2b2f36", mew=1.1, ms=5.4, zorder=3)
        if np.isfinite(m_low):
            left.plot([m_low, m_high], [position - 0.16] * 2, color=colour, lw=1.6,
                      solid_capstyle="butt", zorder=2)
        if np.isfinite(matched):
            left.plot([matched], [position - 0.16], marker="o", color=colour, ms=6.0,
                      zorder=4)
        if np.isfinite(random_point):
            left.plot([random_point], [position - 0.16], marker="x", color="#5b6068",
                      ms=5.0, mew=1.1, zorder=5)
        if not np.isfinite(baseline):
            left.text(0.02, position, "基线不可复现 / 口径不适用", ha="left", va="center",
                      fontsize=7, color="#5b6068", style="italic")

    left.axvline(0.0, color="#2b2f36", lw=0.9, zorder=1)
    left.set_yticks(positions)
    left.set_yticklabels([SHORT.get(r["claim_id"], r["claim_id"]) for r in rows])
    left.set_xlabel("效应量 log2 优势比（被主张的属性 × 论文阳性标签）")
    left.set_title("已发表论断的效应量：论文口径（空心）→ 可检出性匹配后（实心）",
                   loc="left")
    left.margins(y=0.05)
    left.set_ylim(-0.75, n - 0.25)

    for position, row in zip(positions, rows):
        auc = _float(row["propensity_auc"])
        low, high = _float(row["propensity_auc_ci_low"]), _float(row["propensity_auc_ci_high"])
        colour = VERDICT_COLOR.get(row["verdict"], "#8a8f98")
        if np.isfinite(low):
            right.plot([low, high], [position] * 2, color=colour, lw=1.6,
                       solid_capstyle="butt")
        if np.isfinite(auc):
            right.plot([auc], [position], marker="o", color=colour, ms=6.0)
            right.annotate("%.3f" % auc, xy=(auc, position), xytext=(0, 6.5),
                           textcoords="offset points", ha="center", fontsize=6.5,
                           color="#2b2f36")
        else:
            right.text(0.505, position, "不适用", ha="left", va="center", fontsize=7,
                       color="#5b6068", style="italic")
    right.axvline(0.5, color="#2b2f36", lw=0.9)
    right.set_yticks(positions)
    right.set_yticklabels([])
    right.set_xlabel("可检出性协变量单独的 AUC")
    right.set_title("混淆有多大：只用可检出性\n就能把阳性集排到前面", loc="left")
    right.set_ylim(-0.75, n - 0.25)
    right.set_xlim(0.45, 1.0)

    handles = [plt.Line2D([], [], marker="o", color=VERDICT_COLOR[v], ls="none", ms=6,
                          label="%s（%s）" % (VERDICT_ZH[v], v))
               for v in ORDER if any(r["verdict"] == v for r in rows)]
    handles += [
        plt.Line2D([], [], marker="o", mfc="white", mec="#2b2f36", mew=1.1, ls="none",
                   ms=5.4, label="基线：论文自述口径"),
        plt.Line2D([], [], marker="x", color="#5b6068", ls="none", ms=5, mew=1.1,
                   label="同规模随机对照（功效参照）"),
        plt.Line2D([], [], color="#5b6068", lw=1.3, label="95% 区间（聚类自助 5000 次）"),
    ]
    # 图例放在整幅图下方、坐标轴标题之外，避免与 x 轴标签叠字
    figure.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.045),
                  ncol=3, frameon=False, handletextpad=0.5, columnspacing=1.6)

    figure.subplots_adjust(left=0.29, right=0.985, top=0.92, bottom=0.135)
    figure.savefig(OUT_PNG, dpi=300, bbox_inches="tight", pad_inches=0.12)
    figure.savefig(OUT_PDF, bbox_inches="tight", pad_inches=0.12)
    print("wrote", os.path.relpath(OUT_PNG, ROOT), os.path.relpath(OUT_PDF, ROOT))
    print("rows", n, "verdicts", {r["claim_id"]: r["verdict"] for r in rows})


if __name__ == "__main__":
    main()
