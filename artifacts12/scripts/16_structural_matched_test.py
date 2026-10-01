"""方向一（大样本）· 结构特征的蛋白内配对检验（先做胰蛋白酶可检出性匹配）

设计沿用 scripts/12：同一个蛋白内部，把被修饰的 Cys 与未被修饰的 Cys 相比。
零假设 = 该蛋白 n 个可检出 Cys 中随机抽 k 个被修饰。
统计量 T = 被修饰 Cys 的特征值之和；用超几何矩得到 E[T]、Var[T]，转 z。
这样物种、丰度、亚细胞定位、提取方法在蛋白内天然配对掉。

三块产出：
  A 结构特征（可及性 / 二级结构 / pLDDT / pKa / 局部电荷 / 8 Å 内 Cys 数）
    —— 每项都同时报「未匹配」与「可检出性匹配后」两个 z
  B 试点两条反直觉结论在千级样本上是否成立（0/12 在 loop、5/12 完全埋藏）
  C 四个稳健序列信号在按二级结构 / 埋藏分层后是否还在（限制条件一的排查）
所有 P 值报 Benjamini-Hochberg 校正后的 q。
"""
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
MINL, MAXL, MISSED, K = 7, 30, 2, 5
TYPES = (("sno", "S-亚硝基化"), ("so", "S-磺酰化"), ("ox", "可逆 Cys 氧化"))
SEQ_CLASSES = {"酸性 DE": set("DE"), "芳香 FWY": set("FWY"),
               "疏水 AVLIMFW": set("AVLIMFW"), "半胱氨酸 C": set("C"),
               "碱性 KRH": set("KRH")}


def phi(z): return 0.5 * (1 + math.erf(z / math.sqrt(2)))
def pval(z): return 2 * (1 - phi(abs(z)))


def bh(pairs):
    """pairs: [(key, p)] -> {key: q}（Benjamini-Hochberg）"""
    ps = sorted(pairs, key=lambda x: x[1])
    m, q, prev = len(ps), {}, 1.0
    for i in range(m - 1, -1, -1):
        k, p = ps[i]
        prev = min(prev, p * m / (i + 1))
        q[k] = round(min(prev, 1.0), 5)
    return q


def cut_sites(s):
    return [i + 1 for i, ch in enumerate(s)
            if ch in "KR" and not (i + 1 < len(s) and s[i + 1] == "P")]


def detectable_positions(s):
    cs = cut_sites(s)
    cuts = sorted(set([0] + cs + [len(s)]))
    ok = set()
    for i in range(len(cuts) - 1):
        for m in range(MISSED + 1):
            j = i + 1 + m
            if j >= len(cuts): break
            st, en = cuts[i], cuts[j]
            if MINL <= en - st <= MAXL:
                ok.update(range(st + 1, en + 1))
    return ok


def win_frac(s, pos, chars, k=K):
    lo, hi = max(0, pos - 1 - k), min(len(s), pos + k)
    w = s[lo:pos - 1] + s[pos:hi]
    return sum(1 for ch in w if ch in chars) / (len(w) or 1)


def moment_test(units):
    """units: [(values[list], modified_idx[set])] 每个蛋白一条。"""
    T = mu = var = 0.0; nmod = npro = 0
    for vals, mods in units:
        n, k = len(vals), len(mods)
        if k == 0 or n < 2 or k == n: continue
        m_ = sum(vals) / n
        s2 = sum((v - m_) ** 2 for v in vals) / n
        T += sum(vals[i] for i in mods)
        mu += k * m_
        var += k * (n - k) / (n - 1) * s2
        nmod += k; npro += 1
    sd = math.sqrt(var) if var > 0 else 0.0
    z = (T - mu) / sd if sd else 0.0
    return {"obs": round(T / nmod, 4) if nmod else None,
            "exp": round(mu / nmod, 4) if nmod else None,
            "diff": round((T - mu) / nmod, 4) if nmod else None,
            "z": round(z, 2), "p": pval(z) if sd else 1.0,
            "n_proteins": npro, "n_modified": nmod}


def build(typ, feats):
    """把修饰数据与结构特征合到一起。返回 [(seq, {pos: featdict}, modified_set)]"""
    data = json.load(open(OUT / f"09_mapped_{typ}.json"))["proteins"]
    rows = []
    for p in data:
        fr = feats.get(p["locus"])
        if not fr or fr.get("status") != "ok": continue
        fmap = {c["pos"]: c for c in fr["cys"]}
        mods = {m for m in p["modified"] if m in fmap}
        if not mods: continue
        rows.append((p["seq"], fmap, mods))
    return rows


def run_feature(rows, getter, restrict_detectable):
    units = []
    for seq, fmap, mods in rows:
        keep = detectable_positions(seq) if restrict_detectable else None
        pos = sorted(p for p in fmap if keep is None or p in keep)
        if len(pos) < 2: continue
        vals = []
        for p in pos:
            v = getter(fmap[p], seq, p)
            if v is None: vals = None; break
            vals.append(float(v))
        if vals is None: continue
        idx = {p: i for i, p in enumerate(pos)}
        mm = {idx[p] for p in mods if p in idx}
        if not mm: continue
        units.append((vals, mm))
    return moment_test(units)


def main():
    feats = json.load(open(OUT / "15_structural_features.json"))["proteins"]

    STRUCT = {
        "相对可及性 rel_sasa":      lambda c, s, p: c["rel_sasa"],
        "完全埋藏 (rel_sasa<=0.05)": lambda c, s, p: 1.0 if c["rel_sasa"] <= 0.05 else 0.0,
        "在 loop/coil 上":           lambda c, s, p: 1.0 if c["ss3"] == "-" else 0.0,
        "在 helix 上":               lambda c, s, p: 1.0 if c["ss3"] == "H" else 0.0,
        "在 sheet 上":               lambda c, s, p: 1.0 if c["ss3"] == "E" else 0.0,
        "pLDDT":                     lambda c, s, p: c["plddt"],
        "无序 (pLDDT<70)":           lambda c, s, p: 1.0 if c["plddt"] < 70 else 0.0,
        "pKa (PROPKA3)":             lambda c, s, p: c["pka"],
        "8A 内碱性残基数":            lambda c, s, p: c["nbasic"],
        "8A 内酸性残基数":            lambda c, s, p: c["nacid"],
        "8A 内净电荷":                lambda c, s, p: c["netq"],
        "8A 内其他 Cys 数":           lambda c, s, p: c["ncys8"],
    }

    result = {"A_structural": {}, "B_pilot_claims": {}, "C_stratified_sequence": {}}
    flat = []

    # ---- A 结构特征 ----
    for typ, label in TYPES:
        rows = build(typ, feats)
        per = {}
        for name, g in STRUCT.items():
            raw = run_feature(rows, g, False)
            mat = run_feature(rows, g, True)
            per[name] = {"raw": raw, "matched": mat}
            flat.append((f"A|{typ}|{name}", mat["p"]))
        result["A_structural"][typ] = {"label": label, "n_proteins_used": len(rows), "features": per}

    # ---- B 试点两条反直觉结论 ----
    for typ, label in TYPES:
        rows = build(typ, feats)
        keepf = detectable_positions
        n_mod = n_loop = n_bur = 0
        exp_loop = exp_bur = 0.0
        for seq, fmap, mods in rows:
            keep = keepf(seq)
            pos = [p for p in fmap if p in keep]
            mm = [p for p in mods if p in keep]
            if len(pos) < 2 or not mm: continue
            n_mod += len(mm)
            n_loop += sum(1 for p in mm if fmap[p]["ss3"] == "-")
            n_bur += sum(1 for p in mm if fmap[p]["rel_sasa"] <= 0.05)
            exp_loop += len(mm) * sum(1 for p in pos if fmap[p]["ss3"] == "-") / len(pos)
            exp_bur += len(mm) * sum(1 for p in pos if fmap[p]["rel_sasa"] <= 0.05) / len(pos)
        result["B_pilot_claims"][typ] = {
            "label": label, "n_modified": n_mod,
            "loop_obs_pct": round(100 * n_loop / n_mod, 1) if n_mod else None,
            "loop_exp_pct": round(100 * exp_loop / n_mod, 1) if n_mod else None,
            "buried_obs_pct": round(100 * n_bur / n_mod, 1) if n_mod else None,
            "buried_exp_pct": round(100 * exp_bur / n_mod, 1) if n_mod else None,
        }

    # ---- C 序列信号按结构分层 ----
    STRATA = {
        "全部": lambda c: True,
        "仅 helix": lambda c: c["ss3"] == "H",
        "仅 sheet": lambda c: c["ss3"] == "E",
        "仅 loop/coil": lambda c: c["ss3"] == "-",
        "仅埋藏 (<=0.25)": lambda c: c["rel_sasa"] <= 0.25,
        "仅暴露 (>0.25)": lambda c: c["rel_sasa"] > 0.25,
    }
    for typ, label in TYPES:
        rows = build(typ, feats)
        per = {}
        for sname, sfun in STRATA.items():
            per[sname] = {}
            for cname, chars in SEQ_CLASSES.items():
                units = []
                for seq, fmap, mods in rows:
                    keep = detectable_positions(seq)
                    pos = sorted(p for p in fmap if p in keep and sfun(fmap[p]))
                    if len(pos) < 2: continue
                    idx = {p: i for i, p in enumerate(pos)}
                    mm = {idx[p] for p in mods if p in idx}
                    if not mm or len(mm) == len(pos): continue
                    units.append(([win_frac(seq, p, chars) for p in pos], mm))
                r = moment_test(units)
                r["diff_pp"] = round(r["diff"] * 100, 2) if r["diff"] is not None else None
                per[sname][cname] = r
                flat.append((f"C|{typ}|{sname}|{cname}", r["p"]))
        result["C_stratified_sequence"][typ] = {"label": label, "strata": per}

    qs = bh(flat)
    result["_bh_q"] = qs
    result["_n_tests"] = len(flat)
    json.dump(result, open(OUT / "16_structural_matched.json", "w"), ensure_ascii=False, indent=1)
    print("检验总数", len(flat), "| 写出 16_structural_matched.json")

    for typ, label in TYPES:
        b = result["B_pilot_claims"][typ]
        print(f"\n{label}  n_mod={b['n_modified']}")
        print(f"  loop  实测 {b['loop_obs_pct']}%  蛋白内期望 {b['loop_exp_pct']}%")
        print(f"  埋藏  实测 {b['buried_obs_pct']}%  蛋白内期望 {b['buried_exp_pct']}%")
        for name in STRUCT:
            m = result["A_structural"][typ]["features"][name]["matched"]
            r = result["A_structural"][typ]["features"][name]["raw"]
            q = qs[f"A|{typ}|{name}"]
            print(f"  {name:26s} z_raw={r['z']:7.2f}  z_matched={m['z']:7.2f}  q={q:.4g}  n={m['n_modified']}")


if __name__ == "__main__":
    main()
