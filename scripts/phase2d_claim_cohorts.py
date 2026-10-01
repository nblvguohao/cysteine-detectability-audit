"""Cohorts for the claims unblocked by the 2026-09-16 material drop.

`reports/PHASE2C_SFE006_REPRODUCTION.md` closed gap G1. This module builds the rest of the
claims whose ATTRIBUTE is computable from sequence or from reference-proteome GO annotation.
Twelve further claims need structure (SASA, pKa, secondary structure, conservation) and are NOT
built here; they are listed in the runner's NEEDS_STRUCTURE table with what each requires.

**Predeclared, fixed before the run.**

* Negative sets follow each paper's own rule, and where a paper's rule is ambiguous or
  unreproducible the substitution is named in the cohort's `reproduction_notes`:
  - Yang 2014's pLogo background is not stated unambiguously in the text (proteome background
    versus unmodified cysteines of the same proteins). We use the latter, which is also the rule
    Bui 2016 states for the same data, and declare it.
  - Marino & Gladyshev's reference set is cysteines of 1000 randomly chosen eukaryotic PDB
    proteins, which cannot be reconstructed. The paper names a SECOND control - all known
    non-modified cysteines of the dataset proteins - and that is what we use.
  - Li 2024's enrichment background is DAVID's whole-genome default. We report that (the mouse
    reference proteome) as primary AND the 472 total zinc-finger proteins as a better-matched
    secondary background, because a whole-genome background is the construction this project's
    phase 1 identified as the one that inflates detectability share.
* Every site-level cohort keeps the authors' attribute offsets verbatim; nothing is re-tuned.
* Covariates are the registered sets: VIS10 for site-level, PROT8 for protein-level.
* Grouping for the cluster bootstrap is the protein for site-level cohorts, and the homology
  component for protein-level cohorts.
* GO attributes are term-with-descendants via QuickGO, so a protein annotated only to a child
  term still counts - matching how the enrichment tools the papers used behave.

Two claims are recorded as NOT TESTABLE rather than forced into a contingency table, with the
reason in the runner's NOT_TESTABLE table:
  * SFE-003 is a qualitative cross-modification motif comparison with no author statistic and no
    negative class of its own.
  * PERS-007's attribute (zinc-finger coordination type) is published ONLY for the positives -
    the authors' total-ZF table carries no type column - so the negative side of the 2x2 does not
    exist in the source. Deriving types ourselves would test our classifier, not their claim.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import os
import re
import zipfile

import numpy as np

from phase2_claim_cohorts import (
    go_with_descendants, protein_feature_matrix, read_proteome_tsv, site_feature_matrix,
    flank_residue_flag, homology_components,
)
from phase2b_claim_cohorts import HYDROPHOBIC, flank_count_flag, read_fasta_gz

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
PROTEOMES = os.path.join(ROOT, "external", "proteomes")
DROP = os.path.join(ROOT, "external", "intake", "manual_fetch_2026-09-16")
SEQ_CACHE = os.path.join(DROP, ".phase2d_uniprot_seq_cache.json")

YANG_SITES = os.path.join(RESULTS, "phase2c_yang2014_sulfenyl_sites.csv")
DOULIAS_SITES = os.path.join(RESULTS, "doulias2010_sno_sites.csv")
MARINO_SITES = os.path.join(RESULTS, "marino2010_nocys_sites.csv")
AKTER = os.path.join(DROP, "akter2018", "41589_2018_116_MOESM33_ESM.xlsx")
LI_ZIP = os.path.join(DROP, "li2024", "d3cb00106g2_suppl.zip")
LI_S2 = "Table S2 - Total ZF list and sorting.xlsx"

# --- author-named offsets, verbatim ---------------------------------------------------------
SFE002_GLU_OFFSETS = [-4, -3, 1, 3, 4, 5]
SFE002_LYS_OFFSETS = [-6, -5, -2, 6]
SNO005_GLY_OFFSET = [-1]
SNO006_NEG_PLUS3 = [3]
SNO006_NEG_MINUS1 = [-1]
SNO014_ASP_OFFSETS = [-1, 1]
SNO014_GLU_OFFSETS = [-3]
NEGATIVE_RESIDUES = {"D", "E"}

# --- GO terms, predeclared one per claim -----------------------------------------------------
GO_OXIDOREDUCTASE = "GO:0016491"      # SNO-008 primary: "largely composed of oxidoreductases"
GO_TRANSFERASE = "GO:0016740"         # SNO-008 secondary: "and transferases"
GO_EXOSOME = "GO:0070062"             # SFI-001: "highly enriched in exosomes"
GO_GLYCOLYSIS = "GO:0006096"          # SFI-002 secondary, author-named
GO_UBIQUITIN_CATABOLIC = "GO:0006511"  # PERS-008: "ubiquitin-mediated proteolysis"


def _cache():
    return json.load(open(SEQ_CACHE, encoding="utf-8")) if os.path.exists(SEQ_CACHE) else {}


def _rows(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def _sheet(path_or_bytes, sheet):
    import openpyxl
    src = io.BytesIO(path_or_bytes) if isinstance(path_or_bytes, bytes) else path_or_bytes
    wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
    ws = wb[sheet]
    header, out = None, []
    for row in ws.iter_rows(values_only=True):
        if header is None:
            header = [str(c).strip() if c is not None else "" for c in row]
            continue
        out.append(dict(zip(header, row)))
    wb.close()
    return out


def _li_sheet(sheet):
    with zipfile.ZipFile(LI_ZIP) as z:
        return _sheet(z.read(LI_S2), sheet)


def _site_cohort(sites, sequences, claim_id, attribute, attribute_label,
                 author_statistic, notes, direction, secondary=None, secondary_label=None):
    """All cysteines of the proteins carrying a positive; positives are `sites`."""
    proteins = sorted({a for a, _ in sites})
    keys, y = [], []
    for acc in proteins:
        seq = sequences[acc]
        for i, ch in enumerate(seq):
            if ch == "C":
                keys.append((acc, i + 1))
                y.append(1 if (acc, i + 1) in sites else 0)
    y = np.asarray(y, dtype=int)
    cohort = {
        "claim_id": claim_id, "unit": "site", "keys": keys, "y": y,
        "attribute": attribute(sequences, keys),
        "groups": np.asarray([a for a, _ in keys]),
        "covariates": site_feature_matrix(sequences, keys, "VIS10"),
        "covariate_set": "VIS10", "observed_depth": None,
        "claim_direction": direction,
        "attribute_label": attribute_label,
        "author_statistic": author_statistic,
        "reproduction_notes": notes,
        "n_proteins": len(proteins), "sequences": sequences,
    }
    if secondary is not None:
        cohort["attribute_secondary"] = secondary(sequences, keys)
        cohort["attribute_secondary_label"] = secondary_label
    return cohort


def _protein_cohort(positives, proteome, claim_id, go_terms, attribute_label,
                    author_statistic, notes, direction, background=None,
                    secondary_terms=None, secondary_label=None):
    """Protein-level GO attribute. `background` restricts the negative class when given."""
    expanded = go_with_descendants(go_terms)
    wanted = set()
    for t in go_terms:
        wanted |= expanded[t]
    universe = sorted(background) if background else sorted(proteome)
    universe = [a for a in universe if a in proteome and proteome[a]["seq"]]
    positives = {a for a in positives if a in proteome}
    keys = [a for a in universe if a in proteome]
    for a in sorted(positives):
        if a not in proteome:
            continue
        if a not in keys:
            keys.append(a)
    y = np.asarray([1 if a in positives else 0 for a in keys], dtype=int)
    attribute = np.asarray([1 if (proteome[a]["go"] & wanted) else 0 for a in keys], dtype=int)
    sequences = {a: proteome[a]["seq"] for a in keys}
    components = homology_components(sequences, keys)
    cohort = {
        "claim_id": claim_id, "unit": "protein", "keys": keys, "y": y,
        "attribute": attribute,
        "groups": np.asarray([components[a] for a in keys]),
        "covariates": protein_feature_matrix(sequences, keys),
        "covariate_set": "PROT8", "observed_depth": None,
        "claim_direction": direction,
        "attribute_label": attribute_label,
        "author_statistic": author_statistic,
        "reproduction_notes": notes,
        "n_proteins": len(keys), "sequences": sequences,
        "go_terms": go_terms,
    }
    if secondary_terms:
        sec_expanded = go_with_descendants(secondary_terms)
        sec_wanted = set()
        for t in secondary_terms:
            sec_wanted |= sec_expanded[t]
        cohort["attribute_secondary"] = np.asarray(
            [1 if (proteome[a]["go"] & sec_wanted) else 0 for a in keys], dtype=int)
        cohort["attribute_secondary_label"] = secondary_label
        cohort["go_terms_secondary"] = secondary_terms
    return cohort


# ---------------------------------------------------------------------------- site-level
def _yang_sites_and_sequences():
    rows = _rows(YANG_SITES)
    sites = {(r["accession"], int(r["position"])) for r in rows}
    sequences = read_fasta_gz(os.path.join(PROTEOMES, "hsa.fasta.gz"))
    sites = {(a, p) for a, p in sites if a in sequences and p <= len(sequences[a])}
    return sites, sequences


def cohort_sfe002():
    """SFE-002 Yang 2014: Glu over-represented at -4/-3/+1/+3/+4/+5, Lys at -6/-5/-2/+6."""
    sites, sequences = _yang_sites_and_sequences()
    return _site_cohort(
        sites, sequences, "SFE-002",
        lambda s, k: flank_residue_flag(s, k, SFE002_GLU_OFFSETS, {"E"}),
        "作者点名偏移 −4/−3/+1/+3/+4/+5 上出现 Glu",
        {"kind": "plogo", "value": None, "n_positive": 1000, "n_negative": None,
         "text": "pLogo：Glu 在 −4/−3/+1/+3/+4/+5 显著过表达，Lys 在 −6/−5/−2/+6 显著过表达（p<0.05）"},
        ["阳性取自 Yang 2014 Supplementary Dataset 2，与 SFE-006 复现同一队列",
         "阴性用同蛋白其余 Cys：作者 pLogo 的背景在正文里没写清是蛋白组背景还是同蛋白未修饰 Cys，"
         "此处取后者，与 Bui 2016 对同一数据陈述的规则一致，差异已声明",
         "两半（Glu 与 Lys）分别作主属性与次属性，主判定取 Glu 那半（偏移数更多）"],
        "positive_preference",
        secondary=lambda s, k: flank_residue_flag(s, k, SFE002_LYS_OFFSETS, {"K"}),
        secondary_label="作者点名偏移 −6/−5/−2/+6 上出现 Lys（论断的另一半）")


def _doulias_sites_and_sequences():
    rows = [r for r in _rows(DOULIAS_SITES) if r["canonical_position"]]
    cache = _cache()
    sequences = {}
    sites = set()
    for r in rows:
        acc, pos = r["accession"], int(r["canonical_position"])
        seq = cache.get(acc)
        if not seq or pos > len(seq) or seq[pos - 1] != "C":
            continue
        sequences[acc] = seq
        sites.add((acc, pos))
    return sites, sequences


def cohort_sno004():
    """SNO-004 Doulias 2010, a NULL claim: hydrophobicity of SNO sites equals unmodified Cys."""
    sites, sequences = _doulias_sites_and_sequences()
    return _site_cohort(
        sites, sequences, "SNO-004",
        lambda s, k: flank_count_flag(s, k, HYDROPHOBIC),
        "±5 窗内疏水残基占比 ≥30%（作者的疏水指数用同一方向）",
        {"kind": "mean_sd", "value": "0.03±0.69 vs 0.1±0.77", "n_positive": 309, "n_negative": 382,
         "text": "疏水指数 0.03±0.69（n=309）对 0.1±0.77（n=382）；309 个中 139 疏水区、170 亲水区"},
        ["阴性用同蛋白其余 Cys，与作者口径一致",
         "作者报的是连续疏水指数，本仪器是 2×2，故属性取本项目既有的疏水窗口标志（窗 ±5、阈 30%），"
         "方向一致但不是同一统计量，差异已声明",
         "这是一条阴性论断（作者自称无差别），判定走 classify_null 分支"],
        "null_no_preference")


def cohort_sno005():
    """SNO-005 Doulias 2010: the top-scoring motif has glycine exclusively at -1."""
    sites, sequences = _doulias_sites_and_sequences()
    return _site_cohort(
        sites, sequences, "SNO-005",
        lambda s, k: flank_residue_flag(s, k, SNO005_GLY_OFFSET, {"G"}),
        "−1 位为 Gly（作者最强基序的专一位）",
        {"kind": "motif_pvalue", "value": "n=37 (12%)", "n_positive": 37, "n_negative": None,
         "text": "最强基序 n=37（12%）、−1 位专一 Gly、P≤0.001"},
        ["阴性用同蛋白其余 Cys，与作者口径一致",
         "作者是基序富集（motif-x 类），本仪器是单偏移 2×2：检的是该基序的**专一位**本身",
         "作者的 n=37 是命中该基序的阳性数，不是全部阳性，故基线规模不可比，比较的是方向与控制后变化"],
        "positive_preference")


def cohort_sno006():
    """SNO-006 Doulias 2010: 2nd and 3rd motifs carry a negative residue at +3 and at -1."""
    sites, sequences = _doulias_sites_and_sequences()
    return _site_cohort(
        sites, sequences, "SNO-006",
        lambda s, k: flank_residue_flag(s, k, SNO006_NEG_PLUS3, NEGATIVE_RESIDUES),
        "+3 位为 Asp 或 Glu（作者第二强基序的专一位）",
        {"kind": "motif_pvalue", "value": "n=31 (10%) / n=25 (8%)", "n_positive": 56,
         "n_negative": None,
         "text": "第二基序 n=31（10%）+3 位负电；第三基序 n=25（8%）−1 位负电；均 P≤0.001"},
        ["阴性用同蛋白其余 Cys，与作者口径一致",
         "两个基序分别作主属性（+3 负电）与次属性（−1 负电）",
         "与 SNO-005 同理：作者 n 是命中基序的阳性子集，比较方向与控制后变化而非规模"],
        "positive_preference",
        secondary=lambda s, k: flank_residue_flag(s, k, SNO006_NEG_MINUS1, NEGATIVE_RESIDUES),
        secondary_label="−1 位为 Asp 或 Glu（作者第三强基序的专一位）")


def cohort_sno014():
    """SNO-014 Marino & Gladyshev 2010: acidic residues over-represented around NO-Cys."""
    rows = [r for r in _rows(MARINO_SITES) if r["status"] == "confirmed_cys"]
    cache = _cache()
    sequences, sites = {}, set()
    for r in rows:
        acc, pos = r["accession"], int(r["position"])
        seq = cache.get(acc) or cache.get(acc.split("-")[0])
        if not seq or pos > len(seq) or seq[pos - 1] != "C":
            continue
        sequences[acc] = seq
        sites.add((acc, pos))
    return _site_cohort(
        sites, sequences, "SNO-014",
        lambda s, k: flank_residue_flag(s, k, SNO014_ASP_OFFSETS, {"D"}),
        "−1 或 +1 位为 Asp（作者点名位）",
        {"kind": "overrepresentation", "value": None, "n_positive": 70, "n_negative": None,
         "text": "Asp 在 −1 与 +1、Glu 在 −3 位过表达；作者参照集为 1000 个随机真核 PDB 蛋白的 Cys"},
        ["位点取自 mmc1.pdf 第 8–11 页 Table S1 的 56 条确认位点（订阅材料，不分发）",
         "作者主参照集（1000 个随机 PDB 蛋白的 Cys）无法重建；改用作者自己写明的第二对照——"
         "数据集内全部已知非修饰 Cys，差异已声明",
         "队列跨物种（人、鼠、大鼠、假单胞菌等），聚类单位仍是蛋白",
         "Glu 在 −3 位作次属性"],
        "positive_preference",
        secondary=lambda s, k: flank_residue_flag(s, k, SNO014_GLU_OFFSETS, {"E"}),
        secondary_label="−3 位为 Glu（作者点名位）")


# ------------------------------------------------------------------------- protein-level
def cohort_sno008():
    """SNO-008 Doulias 2010: SNO proteins are largely oxidoreductases and transferases."""
    rows = [r for r in _rows(DOULIAS_SITES) if r["canonical_position"]]
    positives = {r["accession"] for r in rows}
    proteome = read_proteome_tsv("mmu")
    return _protein_cohort(
        positives, proteome, "SNO-008", [GO_OXIDOREDUCTASE],
        "具氧化还原酶活性（GO:0016491 及其后代）",
        {"kind": "percentage", "value": "39%", "n_positive": 99, "n_negative": None,
         "text": "99 个具催化活性的 SNO 蛋白中氧化还原酶 39%、转移酶 17%"},
        ["阳性为 Doulias 位点表里带规范位置的蛋白；作者的 99 是其中具催化活性的子集",
         "背景为小鼠参考蛋白组——作者这条是组成比例陈述，没有给出同实验未修饰蛋白名单，"
         "所以只能用蛋白组背景，即 NEG_A 构造、标签语义 T3，这一点必须随数字引用",
         "转移酶（GO:0016740）作次属性"],
        "positive_preference",
        secondary_terms=[GO_TRANSFERASE],
        secondary_label="具转移酶活性（GO:0016740 及其后代）——作者点名的第二类")


def cohort_sno008_secondary_terms():
    return [GO_TRANSFERASE]


def _akter_positives():
    accs = set()
    for sheet in ("A549", "HeLa"):
        for r in _sheet(AKTER, sheet):
            a = str(r.get("Uniprot Accession #") or "").strip().split("-")[0]
            if re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]", a):
                accs.add(a)
    return accs


def cohort_sfi001():
    """SFI-001 Akter 2018: nearly half of sulfinylated proteins are enriched in exosomes."""
    positives = _akter_positives()
    proteome = read_proteome_tsv("hsa")
    return _protein_cohort(
        positives, proteome, "SFI-001", [GO_EXOSOME],
        "注释到胞外外泌体（GO:0070062 及其后代）",
        {"kind": "percentage", "value": "47.3%", "n_positive": None, "n_negative": None,
         "text": "47.3% 的 S-亚磺酰化蛋白高度富集于外泌体（正文另一处写 ~50%）"},
        ["阳性为 MOESM33 的 A549 与 HeLa 两张位点表里全部唯一登录号（同位点表，蛋白级聚合）",
         "背景为人源参考蛋白组，与作者的富集背景一致（人类全蛋白组注释）",
         "作者未给 p 值，只给百分比；本轮给的是优势比与区间"],
        "positive_preference")


def cohort_sfi002():
    """SFI-002 Akter 2018: sulfinylated proteins are enriched in redox and glycolytic processes."""
    positives = _akter_positives()
    proteome = read_proteome_tsv("hsa")
    return _protein_cohort(
        positives, proteome, "SFI-002", [GO_OXIDOREDUCTASE],
        "具氧化还原酶活性（GO:0016491 及其后代）——作者点名的第一类",
        {"kind": "category_list", "value": None, "n_positive": None, "n_negative": None,
         "text": "富集类别含 oxidation-reduction、cell-cell adhesion、RNA processing、glycolysis、"
                 "核输入、脂肪酸 β 氧化；正文未给 p 值"},
        ["阳性与背景同 SFI-001",
         "作者点名六类而未给 p 值，运行前预先声明只判**第一类**（氧化还原），"
         "糖酵解（GO:0006096）作次属性；其余四类不判，避免多重比较后挑结果",
         "与 SFI-001 共用同一阳性集，两条不是独立检验"],
        "positive_preference",
        secondary_terms=[GO_GLYCOLYSIS],
        secondary_label="注释到糖酵解（GO:0006096 及其后代）——作者点名的另一类")


def cohort_sfi002_secondary_terms():
    return [GO_GLYCOLYSIS]


def cohort_pers008():
    """PERS-008 Li 2024: persulfidated zinc fingers are enriched in ubiquitin-mediated proteolysis."""
    ssh = {str(r.get("SSH ZFs") or "").strip() for r in _li_sheet("MEF SSH ZFs total")}
    total = {str(r.get("Entry") or "").strip() for r in _li_sheet("MEF Total ZFs")}
    ssh = {a for a in ssh if re.fullmatch(r"[A-Z0-9]{6,10}", a)}
    total = {a for a in total if re.fullmatch(r"[A-Z0-9]{6,10}", a)}
    proteome = read_proteome_tsv("mmu")
    return _protein_cohort(
        ssh, proteome, "PERS-008", [GO_UBIQUITIN_CATABOLIC],
        "注释到泛素依赖的蛋白分解（GO:0006511 及其后代）",
        {"kind": "pathway_enrichment", "value": None, "n_positive": len(ssh), "n_negative": None,
         "text": "GO 与 KEGG 均指向 ubiquitin-mediated proteolysis；阈值 0.05，正文未给具体 p 值"},
        ["阳性为 Table S2「MEF SSH ZFs total」里的登录号",
         "背景为小鼠参考蛋白组，与作者的 DAVID 默认全基因组背景一致（主口径）",
         "另跑一个以 472 个总锌指蛋白为背景的次口径——那是匹配得多的对照，"
         "而全基因组背景正是第一阶段判定会放大可检出性份额的那种构造",
         "作者用 KEGG 通路，本轮用 GO 的对应术语，差异已声明"],
        "positive_preference")


def pers008_restricted_background():
    total = {str(r.get("Entry") or "").strip() for r in _li_sheet("MEF Total ZFs")}
    return {a for a in total if re.fullmatch(r"[A-Z0-9]{6,10}", a)}


COHORTS = {
    "SFE-002": cohort_sfe002,
    "SNO-004": cohort_sno004,
    "SNO-005": cohort_sno005,
    "SNO-006": cohort_sno006,
    "SNO-014": cohort_sno014,
    "SNO-008": cohort_sno008,
    "SFI-001": cohort_sfi001,
    "SFI-002": cohort_sfi002,
    "PERS-008": cohort_pers008,
}
