"""Cohorts for the five claims whose attribute is structural and has an explicit author threshold.

Attributes come from `results/structural_site_features.csv`
(`scripts/build_structural_site_features.py`), which computes per-cysteine solvent accessibility
and secondary structure from AlphaFold models. Everything else - negatives, grouping, covariates,
estimand - is unchanged from phase 2/2b/2d.

**Predeclared, fixed before the run.**

* Negatives are every other cysteine of the proteins carrying a positive, which is what Yang 2014,
  Doulias 2010 and Marino & Gladyshev's second control all state. Grouping is by protein, as in
  every earlier site-level cohort here.
* A site enters the cohort only if the structural table has a value for it. A protein with no
  usable AlphaFold model loses ALL its cysteines including its positives; the count is reported and
  is the honest coverage cost of this round, not something to be hidden in a footnote.
* Thresholds are the authors' own where stated:
  - `SFE-001` relative accessibility > 0.25, stated as ">25% RSA" by Yang 2014.
  - `SNO-012` sulfur-atom area <= 1.0 A^2, stated as the buried criterion by Marino & Gladyshev.
  - `SNO-002` Doulias 2010 does NOT state its buried/exposed cutoff. The same 0.25 RSA threshold is
    borrowed from Yang and declared as a substitution. The authors' reported composition (29%
    exposed among SNO sites, 23% among unmodified) is the check: if the borrowed threshold gives
    composition far from that, the threshold is not recoverable and the verdict must say so.
  - `SNO-001` helix is P-SEA class `a`, coil is class `c`; the authors report 40% vs 29% helix and
    32% vs 39% coil.
  - `SNO-009` coil fraction over positions 0..+3 >= 0.5; the authors report a 10-percentage-point
    rise in coil frequency across that window.
* `SNO-001` runs coil as the SECONDARY attribute, because the claim has two halves (helix
  over-represented, coil under-represented) and the second is a depletion, so its sign is expected
  to be negative.
* The pLDDT sensitivity is a separate caliber, not a filter on the primary: restricting to
  confident residues conditions on order, and order correlates with the accessibility being
  measured. Primary uses every site with a structural value.
"""
from __future__ import annotations

import csv
import os

import numpy as np

from phase2_claim_cohorts import site_feature_matrix
from phase2d_claim_cohorts import (
    _doulias_sites_and_sequences, _rows, _yang_sites_and_sequences,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
FEATURES = os.path.join(RESULTS, "structural_site_features.csv")
MARINO_SITES = os.path.join(RESULTS, "marino2010_nocys_sites.csv")

YANG_RSA = 0.25
MARINO_SG_BURIED = 1.0
COIL_WINDOW_FRACTION = 0.5

_FEATURE_CACHE = {}


def _features():
    if not _FEATURE_CACHE:
        for r in _rows(FEATURES):
            _FEATURE_CACHE[(r["accession"], int(r["position"]))] = r
    return _FEATURE_CACHE


def _marino_sites_and_sequences():
    import json
    cache = json.load(open(
        os.path.join(ROOT, "external", "intake", "manual_fetch_2026-09-16",
                     ".phase2d_uniprot_seq_cache.json"), encoding="utf-8"))
    sites, sequences = set(), {}
    for r in _rows(MARINO_SITES):
        if r["status"] != "confirmed_cys":
            continue
        acc, pos = r["accession"], int(r["position"])
        seq = cache.get(acc) or cache.get(acc.split("-")[0])
        if not seq or pos > len(seq) or seq[pos - 1] != "C":
            continue
        sequences[acc] = seq
        sites.add((acc, pos))
    return sites, sequences


def _structural_site_cohort(sites, sequences, claim_id, attribute_of, attribute_label,
                            author_statistic, notes, direction,
                            secondary_of=None, secondary_label=None, plddt_only=False):
    """Like phase2d._site_cohort, but the attribute is read from the structural table."""
    feats = _features()
    proteins = sorted({a for a, _ in sites})
    keys, y, attr, sec, dropped, positives_lost = [], [], [], [], 0, 0
    for acc in proteins:
        seq = sequences[acc]
        for i, ch in enumerate(seq):
            if ch != "C":
                continue
            pos = i + 1
            row = feats.get((acc, pos))
            value = None if row is None else attribute_of(row)
            if row is None or value is None or (plddt_only and row["plddt_ge_70"] != "1"):
                dropped += 1
                positives_lost += 1 if (acc, pos) in sites else 0
                continue
            keys.append((acc, pos))
            y.append(1 if (acc, pos) in sites else 0)
            attr.append(int(value))
            if secondary_of is not None:
                sec.append(int(secondary_of(row)))
    y = np.asarray(y, dtype=int)
    used_proteins = sorted({a for a, _ in keys})
    cohort = {
        "claim_id": claim_id, "unit": "site", "keys": keys, "y": y,
        "attribute": np.asarray(attr, dtype=int),
        "groups": np.asarray([a for a, _ in keys]),
        "covariates": site_feature_matrix(sequences, keys, "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "claim_direction": direction, "attribute_label": attribute_label,
        "author_statistic": author_statistic,
        "reproduction_notes": notes + [
            f"结构表覆盖：{len(keys)} 个 Cys 进入队列，{dropped} 个因无可用模型或该量缺失被剔除"
            f"（其中阳性 {positives_lost} 个）；蛋白 {len(used_proteins)}/{len(proteins)}"
            + ("；本行是 pLDDT>=70 的敏感性口径" if plddt_only else "")],
        "n_proteins": len(used_proteins), "sequences": sequences,
        "n_sites_dropped_no_structure": dropped, "n_positives_lost": positives_lost,
    }
    if secondary_of is not None:
        cohort["attribute_secondary"] = np.asarray(sec, dtype=int)
        cohort["attribute_secondary_label"] = secondary_label
    return cohort


def _rsa_exposed(row):
    return None if row["rsa"] == "" else float(row["rsa"]) > YANG_RSA


def _sg_buried(row):
    return None if row["sg_sasa"] == "" else float(row["sg_sasa"]) <= MARINO_SG_BURIED


def _is_helix(row):
    return None if row["sse"] == "" else row["sse"] == "a"


def _is_coil(row):
    return None if row["sse"] == "" else row["sse"] == "c"


def _coil_window(row):
    v = row["sse_coil_frac_0to3"]
    return None if v == "" else float(v) >= COIL_WINDOW_FRACTION


def cohort_sfe001(plddt_only=False):
    """SFE-001 Yang 2014: sulfenylation prefers solvent-exposed cysteines (>25% RSA)."""
    sites, sequences = _yang_sites_and_sequences()
    return _structural_site_cohort(
        sites, sequences, "SFE-001", _rsa_exposed,
        f"相对可及面积 > {YANG_RSA:.0%}（作者明写的阈值）",
        {"kind": "composition", "value": ">60% vs <30%", "n_positive": 1000, "n_negative": None,
         "text": ">60% 的次磺酰化 Cys 相对 RSA>25%（NetSurfP 预测），同批蛋白内未修饰 Cys <30%"},
        ["阴性用同蛋白其余 Cys，与作者口径一致",
         "属性替代：作者用 NetSurfP 序列预测器，本轮用 AlphaFold 模型加 Shrake-Rupley 实算，"
         "阈值沿用作者明写的 25% RSA",
         "作者报的组成（阳性 >60% 暴露、阴性 <30%）是本轮阈值与方法选择的核对点"],
        "positive_preference", plddt_only=plddt_only)


def cohort_sno002(plddt_only=False):
    """SNO-002 Doulias 2010: SNO sites slightly more exposed than unmodified cysteines."""
    sites, sequences = _doulias_sites_and_sequences()
    return _structural_site_cohort(
        sites, sequences, "SNO-002", _rsa_exposed,
        f"相对可及面积 > {YANG_RSA:.0%}（作者未给阈值，借用 Yang 的）",
        {"kind": "composition", "value": "29% vs 23% exposed", "n_positive": 139, "n_negative": 732,
         "text": "SNO 71% 埋藏（n=99）/29% 暴露（n=40）；未修饰 77% 埋藏（n=561）/23% 暴露（n=171）"},
        ["阴性用同蛋白其余 Cys，与作者口径一致",
         "**作者没写埋藏/暴露的阈值**，此处借用 Yang 明写的 25% RSA，已声明为替代",
         "作者报的组成（阳性 29% 暴露、阴性 23%）是阈值可否复原的核对点：若相差很远，"
         "结论只能说该阈值不可复原",
         "属性替代：作者用有实验结构的那部分位点，本轮统一用 AlphaFold 模型"],
        "positive_preference", plddt_only=plddt_only)


def cohort_sno012(plddt_only=False):
    """SNO-012 Marino 2010, a NULL-ish claim: about 35% of NO-Cys sulfur atoms are buried."""
    sites, sequences = _marino_sites_and_sequences()
    return _structural_site_cohort(
        sites, sequences, "SNO-012", _sg_buried,
        f"硫原子可及面积 ≤ {MARINO_SG_BURIED} Å²（作者明写的埋藏判据）",
        {"kind": "composition", "value": "~35% buried", "n_positive": 70, "n_negative": None,
         "text": "约 35% 的 NO-Cys 硫原子暴露值 ≤1.0 Å²（判为埋藏），尽管整体略有暴露富集"},
        ["作者主参照集（1000 个随机真核 PDB 蛋白的 Cys）无法重建，改用作者自己写明的第二对照"
         "——数据集内其余未修饰 Cys",
         "属性替代：作者用同源建模，本轮用 AlphaFold；**阈值是绝对面积 1.0 Å²，"
         "而 Shrake-Rupley 与 Lee-Richards 的绝对值略有差异**，这一点必须随数字引用",
         "这是一条组成陈述而非偏好陈述，判定走 classify_null 分支"],
        "null_no_preference", plddt_only=plddt_only)


def cohort_sno001(plddt_only=False):
    """SNO-001 Doulias 2010: SNO sites over-represented in helix, under-represented in coil."""
    sites, sequences = _doulias_sites_and_sequences()
    return _structural_site_cohort(
        sites, sequences, "SNO-001", _is_helix, "位点落在 α 螺旋（P-SEA 类 a）",
        {"kind": "composition_pvalue", "value": "helix 40% vs 29% (P<0.02)",
         "n_positive": 139, "n_negative": 561,
         "text": "螺旋 40% 对未修饰 29%（P<0.02）；coil 32% 对 39%（P<0.01）；β 折叠 28% 对 32%"},
        ["阴性用同蛋白其余 Cys，与作者口径一致",
         "属性替代：作者用 DSSP 类别，本轮用 biotite 的 P-SEA 指派（本沙箱没有 DSSP 二进制），"
         "两者都给螺旋/折叠/coil 三态但指派规则不同",
         "论断有两半：螺旋过表达（主属性）与 coil 欠表达（次属性，符号应为负）"],
        "positive_preference", secondary_of=_is_coil,
        secondary_label="位点落在 coil（论断的另一半，作者称欠表达，符号应为负）",
        plddt_only=plddt_only)


def cohort_sno009(plddt_only=False):
    """SNO-009 Doulias 2010: coil frequency rises ~10 points across positions 0 to +3."""
    sites, sequences = _doulias_sites_and_sequences()
    return _structural_site_cohort(
        sites, sequences, "SNO-009", _coil_window,
        f"0 至 +3 位里 coil 占比 ≥ {COIL_WINDOW_FRACTION:.0%}（作者点名的窗口）",
        {"kind": "composition_pvalue", "value": "+10 percentage points (P<0.001)",
         "n_positive": 139, "n_negative": 561,
         "text": "coil 频率在 0 至 +3 位升高约 10 个百分点（P<0.001），同时 β 折叠频率下降（P<0.001）"},
        ["阴性用同蛋白其余 Cys，与作者口径一致",
         "作者报的是窗口内 coil 频率的逐位升高，本仪器是 2×2，"
         "故属性预先声明为「窗口内 coil 占比 ≥50%」，方向一致但不是同一统计量",
         "属性替代：P-SEA 指派而非 DSSP"],
        "positive_preference", plddt_only=plddt_only)


COHORTS = {
    "SFE-001": cohort_sfe001,
    "SNO-002": cohort_sno002,
    "SNO-012": cohort_sno012,
    "SNO-001": cohort_sno001,
    "SNO-009": cohort_sno009,
}
