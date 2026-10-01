"""方向一（大样本）· 三种半胱氨酸修饰的序列组成检验 —— 解析法

蛋白内配对：每个蛋白从它的 n 个 Cys 中不放回地抽 k 个（k=实际修饰数）。
该零分布的均值与方差可精确求出：
    E[Σ] = k·μ_i        Var[Σ] = k·(n−k)/(n−1)·σ²_i
蛋白之间独立，求和后用中心极限定理得 z 与双侧 P。
另用置换法对其中一项做校验，确认近似可靠。
"""
import json, random, math
from pathlib import Path
OUT = Path(__file__).resolve().parents[1] / "results"
random.seed(20260910); K = 5
CLASSES = {"碱性 KRH": set("KRH"), "酸性 DE": set("DE"), "芳香 FWY": set("FWY"),
           "疏水 AVLIMFW": set("AVLIMFW"), "小 AGS": set("AGS"),
           "极性 STNQ": set("STNQ"), "脯氨酸 P": set("P"), "甘氨酸 G": set("G"),
           "半胱氨酸 C": set("C")}
def win(s, pos, k=K):
    lo, hi = max(0, pos-1-k), min(len(s), pos+k)
    return s[lo:pos-1] + s[pos:hi]
def phi(z): return 0.5*(1+math.erf(z/math.sqrt(2)))

def prep(typ):
    d = json.loads((OUT / f"09_mapped_{typ}.json").read_text(encoding="utf-8"))
    out = []
    for p in d["proteins"]:
        s = p["seq"]; mods = sorted(set(p["modified"]))
        comp = {}
        for c in p["cys"]:
            w = win(s, c); L = len(w) or 1
            comp[c] = {lab: sum(1 for ch in w if ch in aas)/L for lab, aas in CLASSES.items()}
        out.append((p["cys"], mods, comp))
    return d["label"], len(d["proteins"]), out

def analytic(prepd, lab):
    T = mu = var = 0.0; nmod = 0
    for cys, mods, comp in prepd:
        n, k = len(cys), len(mods)
        if k == 0 or n < 2: continue
        vals = [comp[c][lab] for c in cys]
        m = sum(vals)/n
        s2 = sum((v-m)**2 for v in vals)/n
        T += sum(comp[c][lab] for c in mods)
        mu += k*m
        var += k*(n-k)/(n-1)*s2 if n > 1 else 0.0
        nmod += k
    sd = math.sqrt(var) if var > 0 else 0.0
    z = (T-mu)/sd if sd else 0.0
    p = 2*(1-phi(abs(z))) if sd else 1.0
    return T/nmod, mu/nmod, z, p, nmod

def permute_check(prepd, lab, nperm=4000):
    obs = sum(comp[c][lab] for cys, mods, comp in prepd for c in mods)
    ge = le = 0
    for _ in range(nperm):
        t = 0.0
        for cys, mods, comp in prepd:
            if not mods: continue
            for c in random.sample(cys, len(mods)): t += comp[c][lab]
        if t >= obs: ge += 1
        if t <= obs: le += 1
    return min(1.0, 2*min(ge, le)/nperm)

allres = {}
for typ in ("sno", "so", "ox"):
    label, npro, prepd = prep(typ)
    res = {}
    for lab in CLASSES:
        o, e, z, p, nmod = analytic(prepd, lab)
        res[lab] = {"observed": round(o, 4), "expected": round(e, 4),
                    "diff_pp": round((o-e)*100, 2), "z": round(z, 2), "p_two_sided": p}
    allres[typ] = {"label": label, "n_proteins": npro, "n_modified": nmod, "classes": res}
    print("=" * 80)
    print(f"{label}   蛋白 {npro}，修饰位点 {nmod}")
    print("=" * 80)
    print(f"{'残基类别':<14}{'修饰位点':>9}{'随机期望':>10}{'差值':>10}{'z':>8}{'双侧 P':>12}")
    for lab, v in sorted(res.items(), key=lambda kv: kv[1]["p_two_sided"]):
        star = "  ***" if v["p_two_sided"] < 0.0056 else ("  *" if v["p_two_sided"] < 0.05 else "")
        pstr = f"{v['p_two_sided']:.2e}" if v["p_two_sided"] < 1e-4 else f"{v['p_two_sided']:.4f}"
        print(f"{lab:<14}{v['observed']*100:>8.1f}%{v['expected']*100:>9.1f}%"
              f"{v['diff_pp']:>+9.2f}pp{v['z']:>8.2f}{pstr:>12}{star}")
    print()

# 校验：对 sno 的极性一项用置换法复核
label, npro, prepd = prep("sno")
pp = permute_check(prepd, "极性 STNQ")
za = allres["sno"]["classes"]["极性 STNQ"]["p_two_sided"]
print(f"校验（sno · 极性 STNQ）：解析 P = {za:.2e}    置换 4000 次 P = {pp:.4f}"
      f"    {'一致（置换分辨率下限即 1/4000）' if pp <= 0.001 or abs(pp-za) < 0.02 else '不一致，需检查'}")
(OUT / "10_large_scale_motif.json").write_text(json.dumps(
    {"window": K, "method": "analytic (exact mean/variance of within-protein sampling) + normal approx",
     "permutation_check": {"type": "sno", "class": "极性 STNQ", "p_perm_4000": pp},
     "bonferroni_threshold": round(0.05/len(CLASSES), 5), "results": allres},
    indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n*** = 通过 Bonferroni（0.05/9 = {0.05/9:.4f}）")
print(f"[写出] {OUT/'10_large_scale_motif.json'}")
