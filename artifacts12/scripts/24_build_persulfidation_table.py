"""从 PRIDE PXD072089 建水稻持硫化位点表（PNAS 2026, 10.1073/pnas.2608150123）

输入：data/persulfidation/SS-all-peptides.tsv（PRIDE 该项目唯一的 SEARCH 文件，MaxQuant peptides 格式）
输出：
  24_mapped_osa_persulfidation.json   ← **主用**：只取单 Cys 肽段推出的位点（定位无歧义）
  24_mapped_osa_allpeptides.json      ← 对照：含 Cys 肽段里每个 C 都算位点
  24_localization_artifact.json       ← 两者在「8 Å 内其他 Cys 数」上的差异

★ 为什么必须只用单 Cys 肽段
  这份肽段表不带位点级修饰标注，位点只能由「肽段起点 + 肽段内 C 的偏移」推出。
  若一条肽段含 2 个以上 C，把每个 C 都当成修饰位点，就等于**人为把互相靠近的 Cys
  成对塞进修饰组**——同一胰蛋白酶肽段内的 Cys 相距 ≤30 aa，在三维上常落在 8 Å 内。
  实测后果：「8 Å 内其他 Cys 数」的 z 从 −2.61（仅单 Cys 肽段）翻成 +2.54（全部肽段），
  方向完全颠倒。这是与碱性残基那个胰蛋白酶假象同一族的**定位歧义假象**。
  多 Cys 肽段只占 11.1%，却足以翻转结论。
"""
import csv, json, collections
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
PEP = ROOT / "data" / "persulfidation" / "SS-all-peptides.tsv"


def sequences():
    seq = {}
    for r in csv.DictReader(open(ROOT / "data" / "uniprot" / "osa_proteome.tsv"), delimiter="\t"):
        if r["Sequence"]: seq[r["Entry"]] = r["Sequence"]
    extra = OUT / "24_extra_sequences.json"          # 参考蛋白组外、按 accession 单独取回的
    if extra.exists(): seq.update(json.load(open(extra)))
    return seq


def build(peptides, seq, single_cys_only):
    sites, stat = collections.defaultdict(set), collections.Counter()
    for r in peptides:
        n = int(r["C Count"])
        if n == 0: continue
        if single_cys_only and n != 1:
            stat["跳过_多Cys肽段"] += 1; continue
        acc = r["Leading razor protein"]; base = acc.split("-")[0]
        key = acc if acc in seq else (base if base in seq else None)
        if key is None: stat["无序列"] += 1; continue
        s = seq[key]; sp = r["Start position"]
        if not sp or sp.lower() in ("nan", ""): stat["无起点"] += 1; continue
        st = int(float(sp))
        if s[st - 1:st - 1 + len(r["Sequence"])] != r["Sequence"]:   # razor 蛋白坐标校验
            q = s.find(r["Sequence"])
            if q < 0: stat["肽段不在序列中"] += 1; continue
            st = q + 1
        for i, ch in enumerate(r["Sequence"]):
            if ch == "C": sites[key].add(st + i)
        stat["采用肽段"] += 1

    prots = []
    for key, ss in sites.items():
        s = seq[key]
        good = {p for p in ss if 0 < p <= len(s) and s[p - 1] == "C"}
        if not good: stat["位点非C"] += 1; continue
        cys = [i + 1 for i, ch in enumerate(s) if ch == "C"]
        if len(cys) < 2: stat["Cys<2"] += 1; continue
        prots.append({"locus": key, "acc": base_of(key), "length": len(s), "seq": s,
                      "modified": sorted(good), "cys": cys, "n_cys": len(cys)})
    return prots, dict(stat)


def base_of(k): return k.split("-")[0]


def main():
    seq = sequences()
    peptides = [r for r in csv.DictReader(open(PEP), delimiter="\t")
                if not r["Reverse"] and not r["Potential contaminant"]]
    cc = collections.Counter(int(r["C Count"]) for r in peptides if int(r["C Count"]) > 0)
    meta = {"含Cys肽段": sum(cc.values()), "单Cys肽段": cc[1],
            "多Cys肽段": sum(v for k, v in cc.items() if k > 1),
            "Cys数分布": dict(sorted(cc.items()))}

    for single, fn, tag in ((True,  "24_mapped_osa_persulfidation.json", "持硫化（水稻，仅单Cys肽段）"),
                            (False, "24_mapped_osa_allpeptides.json",    "持硫化（水稻，全部含Cys肽段）")):
        prots, stat = build(peptides, seq, single)
        json.dump({"osa_persulfidation": {
            "type": tag, "single_cys_peptides_only": single,
            "source": "PRIDE PXD072089 / PNAS 2026 (10.1073/pnas.2608150123)",
            "note": ("位点由单 Cys 肽段推出，定位无歧义" if single else
                     "含 Cys 肽段内每个 C 均计为位点——定位歧义，仅作对照，勿用于邻近 Cys 类特征"),
            "peptide_stats": meta, "build_stats": stat, "proteins": prots}},
            open(OUT / fn, "w"), ensure_ascii=False)
        print(f"{tag}: {len(prots)} 蛋白 / {sum(len(p['modified']) for p in prots)} 位点  {stat}")
    json.dump(meta, open(OUT / "24_peptide_stats.json", "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
