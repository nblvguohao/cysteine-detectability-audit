"""按「可检出性」匹配后重跑 —— 检验碱性富集是否真的消失

做法：对每个蛋白做**计算机胰蛋白酶酶切**（K/R 之后切，P 之前不切），
允许 0–2 个漏切。只保留落在可检出肽段（长度 7–30 aa）内的 Cys，
再用同样的蛋白内配对检验比较。

若碱性富集在匹配后大幅缩小 → 证实它源于可检出性偏倚。
"""
import json, math, re
from pathlib import Path
OUT = Path(__file__).resolve().parents[1] / "results"
KRH = set("KRH"); K = 5
MINL, MAXL, MISSED = 7, 30, 2
CLASSES = {"碱性 KRH": set("KRH"), "酸性 DE": set("DE"), "芳香 FWY": set("FWY"),
           "疏水 AVLIMFW": set("AVLIMFW"), "极性 STNQ": set("STNQ"), "半胱氨酸 C": set("C")}
def phi(z): return 0.5*(1+math.erf(z/math.sqrt(2)))

def cut_sites(s):
    """胰蛋白酶切点（K/R 之后，后面不是 P）"""
    return [i+1 for i, ch in enumerate(s) if ch in "KR" and not (i+1 < len(s) and s[i+1] == "P")]

def detectable_positions(s):
    """返回落在可检出肽段内的残基位置集合（1-based）"""
    cuts = [0] + cut_sites(s) + ([len(s)] if (not cut_sites(s) or cut_sites(s)[-1] != len(s)) else [])
    cuts = sorted(set(cuts))
    ok = set()
    for i in range(len(cuts)-1):
        for m in range(MISSED+1):
            j = i+1+m
            if j >= len(cuts): break
            st, en = cuts[i], cuts[j]
            if MINL <= en-st <= MAXL:
                ok.update(range(st+1, en+1))
    return ok

def win(s, pos, k=K):
    lo, hi = max(0, pos-1-k), min(len(s), pos+k)
    return s[lo:pos-1] + s[pos:hi]

def analyse(prots, lab, restrict=None):
    T=mu=var=0.0; nmod=0; npro=0
    for p in prots:
        s=p["seq"]; cys=p["cys"]; mods=set(p["modified"])
        if restrict is not None:
            keep = restrict[p["locus"]]
            cys = [c for c in cys if c in keep]; mods = {m for m in mods if m in keep}
        n,k = len(cys), len(mods)
        if k==0 or n<2: continue
        vals=[ (lambda w: sum(1 for ch in w if ch in CLASSES[lab])/(len(w) or 1))(win(s,c)) for c in cys]
        m_=sum(vals)/n; s2=sum((v-m_)**2 for v in vals)/n
        idx={c:i for i,c in enumerate(cys)}
        T+=sum(vals[idx[c]] for c in mods); mu+=k*m_
        var+=k*(n-k)/(n-1)*s2 if n>1 else 0.0
        nmod+=k; npro+=1
    sd=math.sqrt(var) if var>0 else 0.0
    z=(T-mu)/sd if sd else 0.0
    return {"obs":T/nmod if nmod else 0,"exp":mu/nmod if nmod else 0,
            "diff_pp":round(((T-mu)/nmod)*100,2) if nmod else 0,
            "z":round(z,2),"p":2*(1-phi(abs(z))) if sd else 1.0,
            "n_proteins":npro,"n_modified":nmod}

for typ,label in (("sno","S-亚硝基化"),("so","S-磺酰化"),("ox","可逆 Cys 氧化")):
    d=json.loads((OUT/f"09_mapped_{typ}.json").read_text(encoding="utf-8"))
    prots=d["proteins"]
    keep={p["locus"]: detectable_positions(p["seq"]) for p in prots}
    # 覆盖率
    allc=sum(len(p["cys"]) for p in prots); okc=sum(sum(1 for c in p["cys"] if c in keep[p["locus"]]) for p in prots)
    allm=sum(len(p["modified"]) for p in prots); okm=sum(sum(1 for c in p["modified"] if c in keep[p["locus"]]) for p in prots)
    print("="*82); print(f"{label}")
    print(f"  全部 Cys {allc} 中，落在可检出肽段内 {okc}（{okc/allc*100:.0f}%）")
    print(f"  修饰 Cys {allm} 中，落在可检出肽段内 {okm}（{okm/allm*100:.0f}%）  ← 应接近 100%，可作自检")
    print("="*82)
    print(f"{'残基类别':<14}{'原分析 差值':>13}{'原 z':>8}{'匹配后 差值':>14}{'匹配后 z':>10}{'匹配后 P':>12}")
    res={}
    for lab in CLASSES:
        a=analyse(prots,lab); b=analyse(prots,lab,restrict=keep)
        pst=f"{b['p']:.1e}" if b['p']<1e-4 else f"{b['p']:.4f}"
        print(f"{lab:<14}{a['diff_pp']:>+12.2f}pp{a['z']:>8.1f}{b['diff_pp']:>+13.2f}pp{b['z']:>10.1f}{pst:>12}")
        res[lab]={"raw":a,"matched":b}
    print()
    (OUT/f"12_matched_{typ}.json").write_text(json.dumps(
        {"type":typ,"label":label,"tryptic":{"min_len":MINL,"max_len":MAXL,"missed":MISSED},
         "coverage":{"all_cys":allc,"detectable_cys":okc,"all_mod":allm,"detectable_mod":okm},
         "classes":res},indent=2,ensure_ascii=False),encoding="utf-8")
print("判读：若「碱性 KRH」的 z 在匹配后大幅下降，即证实原信号来自可检出性偏倚。")
