#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""公开数据上的 ±7 位 K/R 位置剖面：背景口径决定这条论断能不能被问。

跑法（用项目 venv，需要 numpy）：

    cd "<本树根目录>"
    export PYTHONDONTWRITEBYTECODE=1
    /path/to/venv/bin/python \
        scripts/positional_kr_profile_public_2026-09-17.py

================================================================================
为什么做这一轮
================================================================================

以**前景自身的打乱**作背景（MoMo motif-x 的 `db_background: false`），
按构造分不开"切点几何"与"序列基序"。

按 `reports/MANUSCRIPT_POSITIONING_2026-09-17.md` 的主判据 **G1**，这条论断必须先有一条
**只用公开数据成立**的证明路径。本轮就是那条路径：在 census 里 7 个公开数据集
上，把同一张剖面算三遍，
**只换背景**。

================================================================================
一 判据（全部在运行前写定）
================================================================================

**窗口**：以每个半胱氨酸为中心，两侧各 7 位，共 15 位。中心位恒为 C。
偏移量 o ∈ {−7…−1, +1…+7}，共 14 个。

**读数**：`frac_KR(o)` = 该臂中偏移 o 处残基为 K 或 R 的比例。
**蛋白两端被截断的窗口，在缺失的那几个偏移上既不进分子也不进分母**（每个偏移单独报 n）。
若改成补一个哨兵字符，等于把"靠近末端"编码成"非 K/R"，会在两端人为压低比例。

**四个臂**（前景一个、背景三个）：

| 臂 | 定义 | 它回答的问题 |
|---|---|---|
| `FG` | census 中 `role == positive` 的位点 | 被报告为修饰的位点长什么样 |
| `BG_proteome` | 该物种蛋白组**全部**半胱氨酸 | 与"所有半胱氨酸"比，有没有偏好（多数论文的隐含背景） |
| `BG_observed` | census 中 `role == observed_unmodified` 的位点 | 与"同一实验里被检出但未报修饰"的半胱氨酸比，还剩多少偏好 |
| `BG_shuffled` | 把每个 FG 窗口的 14 个侧翼残基就地打乱 | MoMo 默认背景（厂商口径）会给出什么 |

`BG_shuffled` 的期望值在每个偏移上都等于该臂侧翼总均值——**这正是它的问题**：
它度量的是"偏离自身均值"，不是"偏离蛋白组"。把它和 `BG_proteome` 并排放，
就能看出 motif-x 的默认答的是另一个问题。打乱种子 20260915。

**统计量**：`ratio(o) = frac_KR_FG(o) / frac_KR_BG(o)`。区间不含 1.0 记为有效应。

**区间**：蛋白聚类自助 **5,000 次、种子 20260915**（房规 2、3，按蛋白整体重抽样，不按位点）。
- `BG_proteome` **不重抽样**：它是全蛋白组的完全枚举，不是样本。只重抽 FG。
- `BG_observed`、`BG_shuffled` 与 FG **联合重抽样**：同一个蛋白同时贡献两臂的窗口，
  按蛋白整体抽，保留二者的相关性。
主区间 95%；另出一条 **Bonferroni(14 个偏移) = 99.643%** 的敏感性区间，同表并列。
不做"只报校正后"的取舍——预写的判定看的是**两个背景在同一批偏移上是否一致**，
不是任何单个偏移是否显著。

================================================================================
二 预写分支（每个数据集独立判定，不跨数据集合并）
================================================================================

令 S_prot = {o : FG/BG_proteome 的 95% 区间不含 1.0}，
   S_obs  = {o : FG/BG_observed  的 95% 区间不含 1.0}。

- **B1 信号存活**：S_prot 非空，且 S_prot ⊆ S_obs 时同向 → 位置偏好不能用"被检出"解释。
- **B2 信号是检出属性**：S_prot 非空，而 S_prot ∩ S_obs = ∅ → 与蛋白组比有偏好，
  与"同一实验里检出但未修饰的半胱氨酸"比没有。**这是稿件主张 3 的预期结果。**
- **B3 不可判定**：S_prot 为空 → 这个样本量上问不了（同时报区间宽度）。
- **B4 反向**：S_obs 里存在与 S_prot 同一偏移但方向相反者 → 原样报，标为反转。
- **B5 部分**：以上都不是（S_prot 与 S_obs 部分重叠）→ 报重叠与不重叠的偏移各是哪些。

**B_shuffle（对所有数据集另判一次）**：存在偏移 o 使 FG/BG_shuffled 与 FG/BG_proteome 的
95% 区间**不相交** → 两个背景给出不相容的答案，只报打乱背景（MoMo 默认）会误述效应量。

**`natcomm2023_ath_sno` 没有 `observed_unmodified` 行**（census 实测 0 条），
它的 `BG_observed` 臂记 `unavailable`，判定只走 B3/B_shuffle，**不得**因此算作 B2 的支持。

================================================================================
三 锚点与阳性对照（房规 11、12）
================================================================================

- **A1 蛋白组半胱氨酸总数**必须复现已存盘的数（`results/invisible_cysteine_map_audit.json`）：
  hsa **262,194**、ath **308,561**。这两个数取自既有产物，本轮不重新定义。
  任一不符 → **REFUSED**，不写产物。
- **A2 中心位必须是 C**：每个抽出的窗口中心必须是 C。逐数据集报 `frac_site_is_cys`；
  任一数据集低于 **0.90** → **REFUSED**（说明位点编号与序列不同源）。
- **A3 打乱背景的构造检验**：`BG_shuffled` 的 14 个偏移的 `frac_KR` 极差应 < 0.02
  （按构造它们同分布）。超出 → 打乱实现有误，**REFUSED**。

  **A3 更正（2026-09-17 运行中发现，原判据一字未改，上面那条作废）**
  第一次运行在 `fps2020_ath_sulfenyl` 上被 A3 拦下（极差 **0.0268** ≥ 0.02），
  **这是判据缺陷，不是实现缺陷**：0.02 是个绝对阈值，对 1,742 个窗口和 23,641 个窗口
  一视同仁。打乱之后 14 个偏移各是一个二项比例，p̄≈0.13、n=1,742 时单个偏移的标准差
  就有 0.008，14 个同分布比例的极差期望约 3.4σ ≈ 0.028——**按构造就该这么大**。
  原阈值会把小数据集一律判成实现有误。
  **改判据（从零分布推出，不看观测值调参）**：同时要求
  ① `shuffled_spread < FG_spread`（打乱必须把位置结构抹平到比前景自身更平）；
  ② `shuffled_spread <= 6 * sqrt(p̄(1−p̄)/n_min)`，n_min 取 14 个偏移中最小的有效窗口数。
  真正的实现错误（比如根本没打乱）会让 ① 直接失败，量级是几十个 σ，不会被 6σ 放过。
  两个数都写进审计，原始那次 REFUSE 记在审计的 `corrections` 里。

================================================================================
四 拟南芥的一处必须预先定死的规则
================================================================================

census 里拟南芥用的是 TAIR 位点号（`AT1G01050`），而蛋白组宇宙是 `external/proteomes/ath.fasta.gz`
（UP000006548，39,272 条，锚点 308,561 个 C）。映射表 `ath_uniprot_tair_map.tsv.gz` 里
**4,138 个 TAIR 号全部能映射，但其中 1,268 个对应多条 UniProt 记录**（同工型/冗余）。

**规则（前景与背景用同一条，不许对前景挑更合适的那条）**：每个 TAIR 位点取
**最长的那条序列**；长度相同则取字典序最小的 UniProt 登录号。选定之后，若该位置不是 C，
该位点**丢弃并计数**（`site_not_cys_in_chosen_sequence`），不去换一条同工型来凑。
这条规则的效果由 A2 直接检验。

**注意**：背景 `BG_proteome` 仍是整条 `ath.fasta.gz`（含全部 39,272 条）以对上锚点 308,561；
前景走上面的选择规则。两者的宇宙**不完全相同**，这一点写进审计的
`what_this_does_NOT_establish`，不藏起来。

================================================================================
五 这一轮不能建立什么
================================================================================

1. **不能建立 K/R 位置结构的成因**。它与胰酶切点几何一致，但本轮没有非胰酶酶解的实验对照
   （那仍是全稿缺的那个对照，见 9.17）。
2. **不能把 B2 读成"论断是错的"**。B2 说的是"在这个报告口径下问不出来"，
   与 `reports/MANUSCRIPT_POSITIONING_2026-09-17.md` 主张 1 一致。
3. **不能跨数据集合并计数当作独立证据**：census 的 7 个数据集来自 7 篇工作，
   同一篇内部的位点高度相关，本轮按数据集分别判定，**不给合并 p 值**。
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import sys
from collections import defaultdict

import warnings

import numpy as np

# 本机 venv（numpy 2.0.2 / macOS Accelerate）在大 matmul 上会抛 divide-by-zero /
# overflow / invalid 三条假警告——与 9.5 登记的那一类同源。本轮当场核过：同一批输入下
# `w @ m` 与 `np.einsum` 逐位相同（最大差 0.0），结果无 NaN、无 inf。核过之后才关掉。
warnings.filterwarnings("ignore", message=r".*encountered in matmul", category=RuntimeWarning)

csv.field_size_limit(10 ** 7)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CENSUS = os.path.join(ROOT, "results", "ptm_census_sites.csv")
ATH_FASTA = os.path.join(ROOT, "external", "proteomes", "ath.fasta.gz")
HSA_FASTA = os.path.join(ROOT, "external", "proteomes", "hsa.fasta.gz")
ATH_MAP = os.path.join(ROOT, "external", "proteomes", "ath_uniprot_tair_map.tsv.gz")

OUT_LONG = os.path.join(ROOT, "results", "positional_kr_profile_public_2026-09-17.csv")
OUT_SUM = os.path.join(ROOT, "results", "positional_kr_profile_public_summary_2026-09-17.csv")
OUT_AUDIT = os.path.join(ROOT, "results", "positional_kr_profile_public_audit.json")

W = 7                      # 每侧宽度
OFFSETS = [o for o in range(-W, W + 1) if o != 0]
SEED = 20260915
REPS = 5000
BONF = 14                  # 14 个偏移
KR = frozenset("KR")

# 锚点：取自 results/invisible_cysteine_map_audit.json，本轮不重新定义
ANCHOR_NCYS = {"hsa": 262194, "ath_tair": 308561}


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def refuse(msg):
    print("REFUSED: {}".format(msg), file=sys.stderr)
    sys.exit(3)


# --------------------------------------------------------------------------
# 蛋白组宇宙
# --------------------------------------------------------------------------

def read_fasta_gz(path):
    seqs = {}
    acc, buf = None, []
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc, buf = ln.split("|")[1], []
            else:
                buf.append(ln.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def ath_tair_choice():
    """TAIR -> 选定序列。规则：最长；同长取字典序最小的 Entry。"""
    fasta = read_fasta_gz(ATH_FASTA)
    best = {}
    with gzip.open(ATH_MAP, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            entry = r["Entry"]
            if entry not in fasta:
                continue
            seq = fasta[entry]
            for t in r["TAIR"].split(";"):
                t = t.strip()
                if not t:
                    continue
                cur = best.get(t)
                if (cur is None
                        or len(seq) > len(cur[1])
                        or (len(seq) == len(cur[1]) and entry < cur[0])):
                    best[t] = (entry, seq)
    return fasta, best


def windows_from_seq(seq, pos1):
    """1-based 位点 -> (窗口字符串, 有效掩码)；越界处标 None。"""
    i = pos1 - 1
    if i < 0 or i >= len(seq):
        return None, None
    out, mask = [], []
    for o in OFFSETS:
        j = i + o
        if 0 <= j < len(seq):
            out.append(seq[j]); mask.append(True)
        else:
            out.append(""); mask.append(False)
    return seq[i], (out, mask)


def proteome_background(seqs):
    """全部半胱氨酸的 ±7 窗口。返回 (kr[n,14] bool, valid[n,14] bool, n_cys)。"""
    kr_rows, ok_rows = [], []
    n_cys = 0
    for seq in seqs.values():
        L = len(seq)
        start = seq.find("C")
        while start != -1:
            n_cys += 1
            kr, ok = [], []
            for o in OFFSETS:
                j = start + o
                if 0 <= j < L:
                    ok.append(True); kr.append(seq[j] in KR)
                else:
                    ok.append(False); kr.append(False)
            kr_rows.append(kr); ok_rows.append(ok)
            start = seq.find("C", start + 1)
    return np.array(kr_rows, bool), np.array(ok_rows, bool), n_cys


# --------------------------------------------------------------------------
# 自助
# --------------------------------------------------------------------------

def cluster_matrix(cluster_ids, kr, ok):
    """按聚类把 (kr, ok) 汇总成 (n_clusters, 28) 的计数矩阵。"""
    uniq = {c: i for i, c in enumerate(sorted(set(cluster_ids)))}
    m = np.zeros((len(uniq), 2 * len(OFFSETS)))
    idx = np.array([uniq[c] for c in cluster_ids])
    for col in range(len(OFFSETS)):
        np.add.at(m[:, col], idx, (kr[:, col] & ok[:, col]).astype(float))
        np.add.at(m[:, len(OFFSETS) + col], idx, ok[:, col].astype(float))
    return m


def boot_counts(rng, m, reps):
    """聚类自助：抽 n_clusters 个聚类（有放回）等价于多项式权重。"""
    n = m.shape[0]
    w = rng.multinomial(n, np.full(n, 1.0 / n), size=reps).astype(float)
    return w @ m


def ci(arr, level):
    lo = (1.0 - level) / 2.0 * 100.0
    return float(np.percentile(arr, lo)), float(np.percentile(arr, 100.0 - lo))


def main():
    script_sha = sha256_of(os.path.abspath(__file__))
    rng = np.random.default_rng(SEED)

    print("载入蛋白组宇宙 ...", flush=True)
    hsa = read_fasta_gz(HSA_FASTA)
    ath_fasta, ath_choice = ath_tair_choice()
    print("  hsa {} 条 / ath {} 条（TAIR 选定 {} 个位点）".format(
        len(hsa), len(ath_fasta), len(ath_choice)), flush=True)

    bg_cache, anchors = {}, {}
    for key, seqs in (("hsa", hsa), ("ath_tair", ath_fasta)):
        kr, ok, n_cys = proteome_background(seqs)
        bg_cache[key] = (kr, ok)
        anchors[key] = n_cys
        want = ANCHOR_NCYS[key]
        status = "PASS" if n_cys == want else "FAIL"
        print("  A1 {:<9} {} 期望 {} 实测 {}".format(key, status, want, n_cys), flush=True)
        if n_cys != want:
            refuse("A1 蛋白组半胱氨酸总数不符：{} 期望 {} 实测 {}".format(key, want, n_cys))

    # ---------------- 读 census，抽窗口 ----------------
    rows = defaultdict(lambda: defaultdict(list))   # dataset -> role -> [(acc, site)]
    meta = {}
    with open(CENSUS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            d = r["dataset_id"]
            meta[d] = (r["chemistry_family"], r["species"], r["proteome_key"])
            rows[d][r["role"]].append((r["accession"], int(r["site"])))

    def seq_for(pkey, acc):
        if pkey == "hsa":
            return hsa.get(acc)
        got = ath_choice.get(acc)
        return got[1] if got else None

    def build(pkey, pairs):
        kr_rows, ok_rows, clusters = [], [], []
        n_total, n_is_cys, n_no_seq = 0, 0, 0
        for acc, pos in pairs:
            n_total += 1
            seq = seq_for(pkey, acc)
            if seq is None:
                n_no_seq += 1
                continue
            centre, win = windows_from_seq(seq, pos)
            if win is None or centre != "C":
                continue
            n_is_cys += 1
            res, mask = win
            kr_rows.append([(res[i] in KR) and mask[i] for i in range(len(OFFSETS))])
            ok_rows.append(list(mask))
            clusters.append(acc)
        if not kr_rows:
            return None
        return (np.array(kr_rows, bool), np.array(ok_rows, bool), clusters,
                n_total, n_is_cys, n_no_seq)

    long_rows, sum_rows = [], []
    a2, a3 = {}, {}

    for dataset in sorted(rows):
        family, species, pkey = meta[dataset]
        fg = build(pkey, rows[dataset]["positive"])
        if fg is None:
            continue
        fkr, fok, fclust, n_total, n_is_cys, n_no_seq = fg
        frac_cys = n_is_cys / n_total if n_total else 0.0
        a2[dataset] = {"n_positive_rows": n_total, "n_window_built": n_is_cys,
                       "n_accession_without_sequence": n_no_seq,
                       "frac_site_is_cys": round(frac_cys, 4)}
        print("  A2 {:<32} frac_site_is_cys {:.4f}  ({}/{})".format(
            dataset, frac_cys, n_is_cys, n_total), flush=True)
        if frac_cys < 0.90:
            refuse("A2 {} 的 frac_site_is_cys {:.4f} < 0.90".format(dataset, frac_cys))

        obs = build(pkey, rows[dataset].get("observed_unmodified", []))

        # --- BG_shuffled：就地打乱 14 个侧翼（只在有效位之间换） ---
        sh_kr = np.zeros_like(fkr)
        for i in range(fkr.shape[0]):
            valid = np.where(fok[i])[0]
            perm = rng.permutation(valid)
            sh_kr[i, valid] = fkr[i, perm]
        sh_ok = fok.copy()

        bkr, bok = bg_cache[pkey]
        bg_kr_n = (bkr & bok).sum(0).astype(float)
        bg_ok_n = bok.sum(0).astype(float)
        bg_frac = bg_kr_n / bg_ok_n

        fg_kr_n = (fkr & fok).sum(0).astype(float)
        fg_ok_n = fok.sum(0).astype(float)
        fg_frac = fg_kr_n / fg_ok_n

        sh_kr_n = (sh_kr & sh_ok).sum(0).astype(float)
        sh_frac = sh_kr_n / fg_ok_n
        spread = float(sh_frac.max() - sh_frac.min())
        fg_spread = float(fg_frac.max() - fg_frac.min())
        # A3 修订版：阈值从二项零分布推出，随样本量走（原绝对阈值 0.02 的作废理由见文档字符串）
        pbar = float(sh_kr_n.sum() / fg_ok_n.sum())
        n_min = float(fg_ok_n.min())
        thr = 6.0 * (pbar * (1.0 - pbar) / n_min) ** 0.5
        a3[dataset] = {"shuffled_spread": round(spread, 4),
                       "fg_spread": round(fg_spread, 4),
                       "threshold_6sigma": round(thr, 4),
                       "n_min_windows": int(n_min)}
        if not (spread < fg_spread and spread <= thr):
            refuse("A3 {} 打乱极差 {:.4f}（前景 {:.4f}，6σ 阈 {:.4f}）".format(
                dataset, spread, fg_spread, thr))

        nO = len(OFFSETS)

        # --- 自助 ---
        m_fg = cluster_matrix(fclust, fkr, fok)
        b_fg = boot_counts(rng, m_fg, REPS)
        fg_boot = b_fg[:, :nO] / np.maximum(b_fg[:, nO:], 1.0)

        m_sh = cluster_matrix(fclust, sh_kr, sh_ok)
        # FG 与 shuffled 同聚类同权重：用同一份权重才保留相关性
        n_cl = m_fg.shape[0]
        w = rng.multinomial(n_cl, np.full(n_cl, 1.0 / n_cl), size=REPS).astype(float)
        b_fg2, b_sh = w @ m_fg, w @ m_sh
        r_shuf = (b_fg2[:, :nO] / np.maximum(b_fg2[:, nO:], 1.0)) / \
                 np.maximum(b_sh[:, :nO] / np.maximum(b_sh[:, nO:], 1.0), 1e-12)

        if obs is not None:
            okr, ook, oclust = obs[0], obs[1], obs[2]
            all_cl = sorted(set(fclust) | set(oclust))
            cidx = {c: i for i, c in enumerate(all_cl)}
            mf = np.zeros((len(all_cl), 2 * nO))
            mo = np.zeros((len(all_cl), 2 * nO))
            fi = np.array([cidx[c] for c in fclust])
            oi = np.array([cidx[c] for c in oclust])
            for col in range(nO):
                np.add.at(mf[:, col], fi, (fkr[:, col] & fok[:, col]).astype(float))
                np.add.at(mf[:, nO + col], fi, fok[:, col].astype(float))
                np.add.at(mo[:, col], oi, (okr[:, col] & ook[:, col]).astype(float))
                np.add.at(mo[:, nO + col], oi, ook[:, col].astype(float))
            n_cl2 = len(all_cl)
            w2 = rng.multinomial(n_cl2, np.full(n_cl2, 1.0 / n_cl2), size=REPS).astype(float)
            bf, bo = w2 @ mf, w2 @ mo
            r_obs = (bf[:, :nO] / np.maximum(bf[:, nO:], 1.0)) / \
                    np.maximum(bo[:, :nO] / np.maximum(bo[:, nO:], 1.0), 1e-12)
            obs_frac = (okr & ook).sum(0) / np.maximum(ook.sum(0), 1)
            n_obs = okr.shape[0]
        else:
            r_obs, obs_frac, n_obs = None, None, 0

        r_prot = fg_boot / bg_frac[None, :]

        S_prot, S_obs, shuffle_conflict = [], [], []
        for k, o in enumerate(OFFSETS):
            rec = {
                "dataset_id": dataset, "chemistry_family": family, "species": species,
                "proteome_key": pkey,
                "offset": o,
                "n_fg_windows": int(fg_ok_n[k]), "fg_frac_KR": round(float(fg_frac[k]), 6),
                "n_bg_proteome": int(bg_ok_n[k]), "bg_proteome_frac_KR": round(float(bg_frac[k]), 6),
                "n_bg_observed": int(n_obs), "bg_observed_frac_KR":
                    (round(float(obs_frac[k]), 6) if obs_frac is not None else ""),
                "bg_shuffled_frac_KR": round(float(sh_frac[k]), 6),
            }
            for name, arr in (("proteome", r_prot), ("observed", r_obs), ("shuffled", r_shuf)):
                if arr is None:
                    rec["ratio_vs_" + name] = ""
                    rec["ci95_lo_vs_" + name] = rec["ci95_hi_vs_" + name] = ""
                    rec["cibonf_lo_vs_" + name] = rec["cibonf_hi_vs_" + name] = ""
                    continue
                point = float(fg_frac[k] / (bg_frac[k] if name == "proteome" else
                                            (obs_frac[k] if name == "observed" else sh_frac[k])))
                lo, hi = ci(arr[:, k], 0.95)
                blo, bhi = ci(arr[:, k], 1.0 - 0.05 / BONF)
                rec["ratio_vs_" + name] = round(point, 4)
                rec["ci95_lo_vs_" + name] = round(lo, 4)
                rec["ci95_hi_vs_" + name] = round(hi, 4)
                rec["cibonf_lo_vs_" + name] = round(blo, 4)
                rec["cibonf_hi_vs_" + name] = round(bhi, 4)
                if name == "proteome" and not (lo <= 1.0 <= hi):
                    S_prot.append((o, 1 if point > 1 else -1))
                if name == "observed" and not (lo <= 1.0 <= hi):
                    S_obs.append((o, 1 if point > 1 else -1))
            if r_obs is not None:
                pl, ph = ci(r_prot[:, k], 0.95)
                sl, sh_hi = ci(r_shuf[:, k], 0.95)
                if ph < sl or sh_hi < pl:
                    shuffle_conflict.append(o)
            else:
                pl, ph = ci(r_prot[:, k], 0.95)
                sl, sh_hi = ci(r_shuf[:, k], 0.95)
                if ph < sl or sh_hi < pl:
                    shuffle_conflict.append(o)
            long_rows.append(rec)

        sp = {o for o, _ in S_prot}
        so = {o for o, _ in S_obs}
        dirp = dict(S_prot)
        diro = dict(S_obs)
        if obs is None:
            branch = "B3_undecidable_no_observed_arm" if not sp else "B_observed_arm_unavailable"
        elif not sp:
            branch = "B3_undecidable"
        elif any(o in diro and diro[o] != dirp[o] for o in sp):
            branch = "B4_reversed"
        elif sp and not (sp & so):
            branch = "B2_detection_property"
        elif sp <= so:
            branch = "B1_survives"
        else:
            branch = "B5_partial"

        sum_rows.append({
            "dataset_id": dataset, "chemistry_family": family, "species": species,
            "n_fg_sites": int(fkr.shape[0]), "n_fg_proteins": len(set(fclust)),
            "n_observed_unmodified_sites": n_obs,
            "n_bg_proteome_cys": anchors[pkey],
            "branch": branch,
            "offsets_vs_proteome": ";".join("{}{:+d}".format("+" if d > 0 else "-", o)
                                            for o, d in sorted(S_prot)),
            "offsets_vs_observed": ("" if obs is None else
                                    ";".join("{}{:+d}".format("+" if d > 0 else "-", o)
                                             for o, d in sorted(S_obs))),
            "shuffle_conflict_offsets": ";".join(str(o) for o in shuffle_conflict),
            "B_shuffle": "backgrounds_disagree" if shuffle_conflict else "backgrounds_agree",
            "frac_site_is_cys": a2[dataset]["frac_site_is_cys"],
            "shuffled_frac_spread": a3[dataset]["shuffled_spread"],
            "fg_frac_spread": a3[dataset]["fg_spread"],
        })
        print("  {:<32} {:<28} vs蛋白组 {:<26} vs检出 {}".format(
            dataset, branch, sum_rows[-1]["offsets_vs_proteome"] or "(无)",
            sum_rows[-1]["offsets_vs_observed"] or "(无/不可用)"), flush=True)

    with open(OUT_LONG, "w", newline="", encoding="utf-8") as fh:
        wri = csv.DictWriter(fh, fieldnames=list(long_rows[0]))
        wri.writeheader(); wri.writerows(long_rows)
    with open(OUT_SUM, "w", newline="", encoding="utf-8") as fh:
        wri = csv.DictWriter(fh, fieldnames=list(sum_rows[0]))
        wri.writeheader(); wri.writerows(sum_rows)

    audit = {
        "script": os.path.basename(__file__),
        "script_sha256": script_sha,
        "interpreter": sys.version.split()[0],
        "numpy": np.__version__,
        "run_utc": __import__("datetime").datetime.utcnow().isoformat() + "Z",
        "seed": SEED, "bootstrap_reps": REPS, "window_half_width": W,
        "offsets": OFFSETS,
        "criteria_declared_before_run": True,
        "inputs": {
            "census": {"path": "results/ptm_census_sites.csv", "sha256": sha256_of(CENSUS)},
            "hsa_fasta": {"path": "external/proteomes/hsa.fasta.gz", "sha256": sha256_of(HSA_FASTA)},
            "ath_fasta": {"path": "external/proteomes/ath.fasta.gz", "sha256": sha256_of(ATH_FASTA)},
            "ath_tair_map": {"path": "external/proteomes/ath_uniprot_tair_map.tsv.gz",
                             "sha256": sha256_of(ATH_MAP)},
        },
        "anchors": {
            "A1_proteome_cysteine_counts": {"expected": ANCHOR_NCYS, "observed": anchors,
                                            "source_of_expected":
                                                "results/invisible_cysteine_map_audit.json",
                                            "status": "PASS"},
            "A2_centre_is_cysteine": a2,
            "A3_shuffled_frac_spread": a3,
        },
        "corrections": [
            {
                "date": "2026-09-17",
                "what": "A3 的阈值从绝对值 0.02 改为随样本量走的二项 6σ 阈，并加一条 shuffled_spread < fg_spread。",
                "why": ("第一次运行在 fps2020_ath_sulfenyl 上被拦下（打乱极差 0.0268 >= 0.02）。"
                        "复核后判定是判据缺陷而非实现缺陷：打乱后 14 个偏移各是一个二项比例，"
                        "p 约 0.13、n=1742 时单偏移标准差约 0.008，14 个同分布比例的极差期望约 3.4σ，"
                        "原绝对阈值会把所有小数据集判成实现有误。"),
                "how_the_new_threshold_was_chosen": "从二项零分布推出（6σ），没有看观测值调参。",
                "what_this_correction_does_NOT_change": [
                    "没有改动任何前景/背景的定义、窗口宽度、自助设置或分支判定规则。",
                    "没有改动 A1、A2 两道锚点。",
                    "第一次运行被 A3 拦下时没有写出任何产物，因此没有被覆盖的数字。",
                ],
            },
            {
                "date": "2026-09-17",
                "what": "关掉了 numpy matmul 的三条 RuntimeWarning。",
                "why": ("macOS Accelerate BLAS 的假警告，与 project notes 9.5 登记的同类。"
                        "关之前当场核过：`w @ m` 与 np.einsum 在同一批输入上逐位相同（最大差 0.0），"
                        "输出无 NaN、无 inf。"),
                "what_this_correction_does_NOT_change": ["没有改动任何数值。"],
            },
        ],
        "arms": {
            "FG": "census role == positive",
            "BG_proteome": "every cysteine in the species proteome universe (complete enumeration, not resampled)",
            "BG_observed": "census role == observed_unmodified (same experiment, detected, not reported modified)",
            "BG_shuffled": "the 14 flank residues of each FG window permuted in place, seed {}".format(SEED),
        },
        "what_this_does_NOT_establish": [
            "不能建立 K/R 位置结构的成因；本轮没有非胰酶酶解的实验对照。",
            "B2 不等于论断是错的；它说的是在这个报告口径下问不出来。",
            "不给跨数据集合并 p 值：7 个数据集来自 7 篇工作，同一篇内部高度相关。",
            "拟南芥前景走 TAIR 最长序列规则，背景是整条 ath.fasta.gz，两个宇宙不完全相同。",
        ],
    }
    with open(OUT_AUDIT, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, ensure_ascii=False, indent=2)

    print("\nwrote {} ({} rows)".format(os.path.relpath(OUT_LONG, ROOT), len(long_rows)))
    print("wrote {} ({} rows)".format(os.path.relpath(OUT_SUM, ROOT), len(sum_rows)))
    print("wrote {}".format(os.path.relpath(OUT_AUDIT, ROOT)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
