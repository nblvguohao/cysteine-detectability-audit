"""方向一（大样本）· 把 PlantPTMViewer 的位点映射到拟南芥蛋白序列

产出：每个蛋白的全部 Cys，标注哪些被该类修饰过。
三种修饰：S-亚硝基化(sno)、S-磺酰化(so)、可逆半胱氨酸氧化(ox)。
"""
import csv, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data" / "ptmviewer"; OUT = ROOT / "results"
PROT = Path("/private/tmp/claude-501/-Users-lyuguohao-Documents-Codex-2026-08-26-ban-zhuanghua/"
            "752a4588-c398-4a64-b83c-23230b6bc8f6/scratchpad/ptmv/ath_proteome.tsv")

def load_proteome():
    seqs = {}
    with open(PROT) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            oln = (r.get("Gene Names (ordered locus)") or "").strip()
            s = r.get("Sequence") or ""
            if not s: continue
            for g in oln.split():
                seqs.setdefault(g.upper(), []).append((int(r["Length"]), s))
    return seqs

def main():
    seqs = load_proteome()
    print(f"蛋白组载入 {len(seqs)} 个 locus")
    allres = {}
    for typ, label in (("sno", "S-亚硝基化"), ("so", "S-磺酰化"), ("ox", "可逆 Cys 氧化")):
        rows = [r for r in csv.DictReader(open(D / f"{typ}_allSPECIES.csv"))
                if r["species"] == "ath" and r["modified_aa"] == "C"]
        byprot = {}
        for r in rows:
            byprot.setdefault(r["locus"].upper(), {"len": int(r["length"]), "sites": set(),
                                                   "desc": r["description"]})["sites"].add(int(r["modification_pos"]))
        mapped, unmapped, lenmismatch, sitemismatch = [], 0, 0, 0
        for locus, info in byprot.items():
            cand = seqs.get(locus)
            if not cand: unmapped += 1; continue
            hit = next((s for L, s in cand if L == info["len"]), None)
            if hit is None: lenmismatch += 1; continue
            bad = [p for p in info["sites"] if not (0 < p <= len(hit) and hit[p-1] == "C")]
            if bad: sitemismatch += 1; continue
            cys = [i+1 for i, ch in enumerate(hit) if ch == "C"]
            if len(cys) < 2: continue
            mapped.append({"locus": locus, "length": len(hit), "seq": hit,
                           "modified": sorted(info["sites"]), "cys": cys,
                           "n_cys": len(cys), "desc": info["desc"][:60]})
        nsite = sum(len(m["modified"]) for m in mapped)
        print(f"\n{label} ({typ})")
        print(f"  原始位点 {len(rows)}，涉及蛋白 {len(byprot)}")
        print(f"  成功映射且 Cys 数>=2 的蛋白 {len(mapped)}，其中修饰位点 {nsite}，"
              f"未修饰 Cys {sum(m['n_cys'] for m in mapped)-nsite}")
        print(f"  丢弃：locus 找不到 {unmapped}，长度不符 {lenmismatch}，位点非 Cys {sitemismatch}")
        allres[typ] = {"label": label, "n_proteins": len(mapped), "n_modified": nsite,
                       "n_unmodified": sum(m["n_cys"] for m in mapped)-nsite}
        (OUT / f"09_mapped_{typ}.json").write_text(json.dumps(
            {"type": typ, "label": label, "proteins": mapped}, ensure_ascii=False), encoding="utf-8")
    (OUT / "09_mapping_summary.json").write_text(json.dumps(allres, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n[写出] {OUT}/09_mapped_*.json")

if __name__ == "__main__":
    main()
