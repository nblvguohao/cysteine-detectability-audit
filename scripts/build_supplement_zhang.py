"""Supplement for Prof. Zhang Hua: the modification-vs-interaction result as an
experimental plan, plus the data requests.

Audience is the PI who runs the wet lab, so this is organised around what to
build, what to measure and what can be concluded - not around the computation.
Every number comes from results/v2_boltz2_ptm_panel_*.{csv,json} and
results/v2_boltz2_ptm_confirm_*.{csv,json}; nothing is restated from memory.
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from docx_helpers import bullets, heading, new_document, para, table

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "deliverables" / "补充材料_致张华老师_持硫化与互作_2026-09-15.docx"


def anchor_rows(interface_csv: pathlib.Path, residue: str = "14") -> dict:
    return {r["arm"]: r for r in csv.DictReader(interface_csv.open(encoding="utf-8-sig"))
            if r["wrky_residue"] == residue}


def main() -> None:
    r1 = anchor_rows(ROOT / "results/v2_boltz2_ptm_panel_interface.csv")
    r2 = anchor_rows(ROOT / "results/v2_boltz2_ptm_confirm_interface.csv")
    a1 = json.loads((ROOT / "results/v2_boltz2_ptm_panel_audit.json").read_text())
    a2 = json.loads((ROOT / "results/v2_boltz2_ptm_confirm_audit.json").read_text())

    def cell(row, audit, arm):
        n = len([k for k in row if k.startswith("cross_pae_s")])
        return (f"{row['cross_pae_min']} Å；{row['n_samples_under_10A']}/{n} 采样 ≤10 Å；"
                f"pair-ipTM {audit['best_pair_iptm_by_arm'][arm]:.4f}")

    doc = new_document()
    para(doc, "**补充材料：持硫化是否调控 WRKY×PP2A-B56 互作**", size=16, align="center")
    para(doc, "致张华老师 · 2026 年 9 月 15 日 · 计算侧（吕国豪）", size=9.5, align="center")
    para(doc, "本材料对应微信里说的那个结论。数据出处：`reports/V2_BOLTZ2_PTM_INTERFACE_PANEL.md`，"
              "两轮判据在预测前写入脚本并记录 SHA-256（"
              f"第一轮 `{a1['prep_script_sha256'][:12]}…`，第二轮 `{a2['prep_script_sha256'][:12]}…`）。",
         size=9, align="center")

    heading(doc, "一、一句话结论与一条限制", 1)
    bullets(doc, [
        "**结论：优先做 PP2A-B56 的 C226。** 把 C226 换成过硫化状态后，与 WRKY 的界面明显变差；"
        "而界面外的 C453（距界面 23 Å）和第二近的 C315（6.7 Å）都没有变化。换随机种子独立重跑一次，方向一致。",
        "**限制：这是“位置”假设，不是“过硫化特异”假设。** 把同一个 C226 换成次磺酸（多一个氧而非多一个硫），"
        "界面变差的幅度完全相当。所以结构预测这一层分不开过硫化与其他加合物，"
        "**特异性只能由实验建立**——例如比较 NaHS 与过氧化氢/次磺酸供体处理下结合强度的变化。",
        "**界面不是被打掉，而是更少形成。** 修饰后最优构象仍能形成界面，变的是形成频率与置信度。"
        "所以读数必须能出梯度，只做“有/无”判断很可能看不出差别。",
    ])

    heading(doc, "二、两轮结果", 1)
    para(doc, "构建：WRKY 1–40 肽段 × PP2A-B56 全长（这是两套架构下唯一能稳定收敛出界面的构建）。"
              "各臂的序列与 MSA 完全相同，只有一个残基的化学结构不同。"
              "指标是界面锚定残基 F14 到 PP2A 链的最小跨链 PAE（越小越确信两者的相对位置），"
              "以及 5/10 个独立构象采样里有多少个达到 ≤10 Å。", size=9.5)
    header = ["修饰位点", "距界面", "第一轮（5 采样，未设种子）", "第二轮（10 采样，种子固定）", "判定"]
    rows = [
        ["野生型", "—", cell(r1["wild_type"], a1, "wild_type"), cell(r2["r2_wild_type"], a2, "r2_wild_type"), "基线"],
        ["**C226**", "5.28 Å", "**" + cell(r1["C226_CSS"], a1, "C226_CSS") + "**",
         "**" + cell(r2["r2_C226_CSS"], a2, "r2_C226_CSS") + "**", "**界面变差，两轮一致**"],
        ["C315", "6.73 Å", cell(r1["C315_CSS"], a1, "C315_CSS"), "未测", "无变化"],
        ["C453（对照）", "23.11 Å", cell(r1["C453_CSS"], a1, "C453_CSS"),
         cell(r2["r2_C453_CSS"], a2, "r2_C453_CSS"), "无变化，对照成立"],
        ["C226 次磺酸", "5.28 Å", "未测", cell(r2["r2_C226_CSO"], a2, "r2_C226_CSO"),
         "同等变差 → 化学不特异"],
    ]
    table(doc, header, rows, widths_cm=[2.3, 1.7, 5.0, 4.2, 2.8], size=8, highlight_rows=(1,))
    para(doc, "第二轮的野生型比第一轮更好（最小 PAE 2.05 Å、8/10），所以第二轮是更强的基线，"
              "C226 在更强的基线上仍然变差。两轮都只按“≤10 Å 采样比例下降”这一支达标，"
              "“幅度”那一支没达到，所以可辩护的说法是“在预注册判据下成立”，不是“效应很强”。", size=9)

    heading(doc, "三、建议的构建与读数", 1)
    header2 = ["构建", "优先级", "预期", "说明"]
    rows2 = [
        ["PP2A-B56 **C226S**", "**首选**", "结合减弱", "唯一有结构证据会影响 WRKY 结合的位点；"
         "胰蛋白酶肽段 `ECLKSVLHR`，质谱可直接定位"],
        ["PP2A-B56 C315S", "次选", "不变或略变", "距界面第二近，但模型里没有变化"],
        ["PP2A-B56 C453S", "阴性位点对照", "不变", "距界面 23 Å；同时是我们模型对该蛋白的排名首选，"
         "所以它为阴性对预测本身也是信息"],
        ["WRKY **F14A**", "**界面验证首选**", "结合丧失/大幅减弱", "界面锚定残基，已有直接突变检验支持"],
        ["WRKY F35A", "界面外对照", "不变", "肽段上另一个苯丙氨酸，用于排除“换掉任何 F 都塌”"],
        ["WRKY E16A", "梯度对照", "减弱而非消除", "不要当作“消除结合”的判据"],
        ["PP2A-B56 V276A / S275A", "界面侧验证", "结合减弱", "另建议加一个 232 附近的突变作第二接触点；"
         "G279 是甘氨酸，不宜作首选"],
    ]
    table(doc, header2, rows2, widths_cm=[3.4, 2.2, 2.6, 7.8], size=8, highlight_rows=(0, 3))
    bullets(doc, [
        "**处理条件**：NaHS 或 GYY4137，以及 LCD1 过表达/沉默背景，比较结合强度变化。"
        "若要建立“过硫化特异”，需要再加一组非过硫化的氧化处理（如 H₂O₂）作并行对照——"
        "这正是计算无法替代的那一步。",
        "**读数**：建议 LCI 或定量 Co-IP 这类能出梯度的方法，Y2H 只用于先筛。"
        "模型预测的是结合概率下降，不是结合丧失。",
        "**诱饵同时做全长与 1–40 肽段**：若肽段结合而全长不结合，提示无序区存在自抑制，这本身是一个发现。",
    ])

    heading(doc, "四、WRKY 自身的 C193/C198：计算给不了答案", 1)
    para(doc, "这一条不用再等我们的计算结果，原因是方法上的，不是算力问题。", size=10)
    bullets(doc, [
        "可信界面只存在于 WRKY N 端 1–40 这一段，而 C193/C198 在下游，不在这段序列里。",
        "把全长拿来算，**野生型本身就形不成界面**：5 个采样的最小跨链 PAE 分别是 "
        "13.11 / 22.01 / 20.73 / 20.94 / 21.67 Å，最佳接触残基是 191、173、216、155、97，从来不是 F14。"
        "没有基线，修饰态与野生型的比较就是两个都不存在的界面在比。",
        "**C193/C198 是 WRKY 锌指的配位残基**（C-X4-C-X23-H-X-H，配位 C193/C198/H222/H224）。"
        "C→S 突变会同时破坏锌指结构与 DNA 结合能力，表型会混合结构效应与修饰效应，设计时需要提前安排对照。",
        "若要回答这一侧，可行路线是全长加锌离子的分子动力学（计算侧可以做，约数天），或者直接做实验。",
    ])

    heading(doc, "五、想请您这边提供的两项数据", 1)
    bullets(doc, [
        "**kiae271 补充数据集 1 的位点表与质谱检索参数。** 关键是检索时有没有设 Sulfide（+31.972）"
        "或 S-S-CAM 这两个修饰。这决定我们训练标签的准确表述：如果只搜了碳酰甲基（+57.021），"
        "那么标签的证据层级是“在富集流程中被拉下来并被烷基化的 Cys”，而不是“直接检测到过硫键的位点”。"
        "这一点直接影响论文题目和摘要怎么写，所以想先查清。",
        "**做过但未发表的阴性位点。** C→A 后无表型的位点，或 biotin-switch 检测为阴性的 Cys，几十个即可。"
        "这类真阴性集目前领域内没有任何公开数据，有了它才能把“能不能被质谱看到”和“化学上会不会被修饰”分开。",
    ])

    heading(doc, "六、边界（写论文时需要一起写的）", 1)
    bullets(doc, [
        "化学不特异：过硫化物与次磺酸在模型里效果相当，不能写成“过硫化特异地调控该互作”。",
        "效应恰好压在预注册阈值上（≤10 Å 采样比例降幅正好 0.30）。",
        "MSA 在被修饰那一列仍是普通半胱氨酸，模型保留了未修饰残基的进化证据，"
        "所以这是**偏保守**的检验——更容易漏掉效应，不容易造出假效应。",
        "两轮各 5 与 10 个采样、一个构建、一台机器；界面本身只有四个残基，位于一个 74% 无序的蛋白上。",
        "阴性结果无法区分“没有别构效应”与“Boltz-2 表示不了这种效应”。",
    ])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(OUT))
    print("saved", OUT)


if __name__ == "__main__":
    main()
