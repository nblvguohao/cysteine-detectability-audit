"""方向一（大样本）· 二硫键状态 + 秩次稳健版检验

两件事：
1. **二硫键**：从 AlphaFold 模型量 SG–SG 距离（<=2.5 A 判为成键）。
   这本身就是位点选择的一个候选解释——锁在二硫键里的 Cys 无法被持硫化/亚硝基化。
   同时它解释了 16 号脚本里 pKa 均值高达 30+ 的异常：PROPKA 对二硫键 Cys 给极端高值。
2. **秩次版**：连续特征（rel_sasa / pLDDT / pKa）改用**蛋白内归一化秩次**重做。
   秩次对离群值免疫，且正是功效分析里所说的「排名类特征」。

沿用 16 号脚本的蛋白内配对 + 胰蛋白酶可检出性匹配设计。
"""
import json, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT, AF = ROOT / "results", ROOT / "data" / "af_large"
SSDIST = 2.5
TYPES = (("sno", "S-亚硝基化"), ("so", "S-磺酰化"), ("ox", "可逆 Cys 氧化"))

import importlib.util
spec = importlib.util.spec_from_file_location("s16", Path(__file__).parent / "16_structural_matched_test.py")
s16 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s16)
detectable_positions, moment_test, bh = s16.detectable_positions, s16.moment_test, s16.bh


def sg_coords(acc):
    f = AF / f"AF-{acc}-F1-model_v6.pdb"
    if not f.exists(): return {}
    out = {}
    for ln in f.read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "SG" and ln[17:20].strip() == "CYS":
            out[int(ln[22:26])] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return out


def disulfide_flags(acc):
    sg = sg_coords(acc)
    if len(sg) < 2: return {p: 0 for p in sg}
    pos = sorted(sg); arr = np.array([sg[p] for p in pos])
    d = np.linalg.norm(arr[:, None, :] - arr[None, :, :], axis=-1)
    np.fill_diagonal(d, 1e9)
    return {p: int(d[i].min() <= SSDIST) for i, p in enumerate(pos)}


def norm_rank(vals):
    """蛋白内归一化秩次 (0..1)，并列取平均秩。"""
    n = len(vals)
    order = sorted(range(n), key=lambda i: vals[i])
    r = [0.0] * n; i = 0
    while i < n:
        j = i
        while j + 1 < n and vals[order[j + 1]] == vals[order[i]]: j += 1
        avg = (i + j) / 2.0
        for k in range(i, j + 1): r[order[k]] = avg
        i = j + 1
    return [x / (n - 1) if n > 1 else 0.5 for x in r]


def main():
    feats = json.load(open(OUT / "15_structural_features.json"))["proteins"]

    # 逐蛋白加上二硫键标记
    ndis = ntot = 0
    for loc, fr in feats.items():
        if fr.get("status") != "ok": continue
        fl = disulfide_flags(fr["acc"])
        for c in fr["cys"]:
            c["ssbond"] = fl.get(c["pos"], 0)
            ndis += c["ssbond"]; ntot += 1
    print(f"全部 Cys {ntot}，其中判为二硫键 {ndis} ({100*ndis/ntot:.1f}%)")

    res, flat = {}, []
    for typ, label in TYPES:
        rows = s16.build(typ, feats)
        per = {}

        # --- 1. 二硫键状态（原值检验） ---
        units = []
        for seq, fmap, mods in rows:
            keep = detectable_positions(seq)
            pos = sorted(p for p in fmap if p in keep)
            if len(pos) < 2: continue
            idx = {p: i for i, p in enumerate(pos)}
            mm = {idx[p] for p in mods if p in idx}
            if not mm: continue
            units.append(([float(fmap[p]["ssbond"]) for p in pos], mm))
        per["在二硫键中"] = moment_test(units)

        # --- 2. 秩次版连续特征 ---
        for name, key in (("rel_sasa 秩次", "rel_sasa"), ("pLDDT 秩次", "plddt"), ("pKa 秩次", "pka")):
            units = []
            for seq, fmap, mods in rows:
                keep = detectable_positions(seq)
                pos = sorted(p for p in fmap if p in keep and fmap[p].get(key) is not None)
                if len(pos) < 2: continue
                idx = {p: i for i, p in enumerate(pos)}
                mm = {idx[p] for p in mods if p in idx}
                if not mm or len(mm) == len(pos): continue
                units.append((norm_rank([float(fmap[p][key]) for p in pos]), mm))
            per[name] = moment_test(units)

        # --- 3. pKa 原值，剔除二硫键 Cys ---
        units = []
        for seq, fmap, mods in rows:
            keep = detectable_positions(seq)
            pos = sorted(p for p in fmap if p in keep and not fmap[p]["ssbond"]
                         and fmap[p].get("pka") is not None)
            if len(pos) < 2: continue
            idx = {p: i for i, p in enumerate(pos)}
            mm = {idx[p] for p in mods if p in idx}
            if not mm or len(mm) == len(pos): continue
            units.append(([float(fmap[p]["pka"]) for p in pos], mm))
        per["pKa（剔除二硫键）"] = moment_test(units)

        # --- 4. 描述性：修饰 vs 未修饰 Cys 的 pKa 分布（可检出、非二硫键） ---
        mv, uv = [], []
        for seq, fmap, mods in rows:
            keep = detectable_positions(seq)
            for p, c in fmap.items():
                if p not in keep or c["ssbond"] or c.get("pka") is None: continue
                (mv if p in mods else uv).append(c["pka"])
        desc = lambda v: {"n": len(v), "中位": round(float(np.median(v)), 2),
                          "均值": round(float(np.mean(v)), 2),
                          "P5-P95": [round(float(np.percentile(v, 5)), 2),
                                     round(float(np.percentile(v, 95)), 2)]} if v else None
        per["_pKa分布_已修饰"] = desc(mv); per["_pKa分布_未修饰"] = desc(uv)

        for k, v in per.items():
            if not k.startswith("_"): flat.append((f"{typ}|{k}", v["p"]))
        res[typ] = {"label": label, "features": per}

    qs = bh(flat)
    json.dump({"n_cys": ntot, "n_disulfide": ndis, "results": res, "_bh_q": qs},
              open(OUT / "17_disulfide_rank.json", "w"), ensure_ascii=False, indent=1)

    for typ, label in TYPES:
        print(f"\n{label}")
        for k, v in res[typ]["features"].items():
            if k.startswith("_"):
                print(f"  {k}: {v}"); continue
            print(f"  {k:20s} obs={v['obs']:>8} exp={v['exp']:>8} z={v['z']:>7} "
                  f"q={qs[f'{typ}|{k}']:.3g} n={v['n_modified']}")


if __name__ == "__main__":
    main()
