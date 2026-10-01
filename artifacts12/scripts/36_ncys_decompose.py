"""把「8 Å 内其他半胱氨酸更少」拆成两件事，settle 掉它对候选集定义的依赖

问题：该特征的符号会随候选集定义翻转——
  候选集 = 全部可检出 Cys        → 7 套里 6 套「更少」（z 最强 −11.1）
  候选集 = 仅单 Cys 肽段内的 Cys → 7 套里 4 套变成「更多」
这说明 ncys8 混了两种邻居，而单 Cys 肽段过滤正好只砍掉其中一种。

拆法：对每个 Cys，把 8 Å 内的其他 Cys 按**序列间隔**分成两类
  近邻 |Δpos| <= 30   —— 与该 Cys 可能落在同一胰蛋白酶肽段（序列邻近）
  远邻 |Δpos| >  30   —— 只能由三维折叠带来（结构邻近）
两类各自做蛋白内配对 + 可检出性匹配检验。若「更少」来自近邻、「更多」来自远邻，
则原结论应改写为两条，而不是一条。

距离用 SG–SG（缺 SG 时回退到 CB/CA），阈值 8 Å，与 15/21 号脚本一致。
"""
import csv, json, math, collections, importlib.util
from pathlib import Path
from multiprocessing import Pool

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
AFD = [ROOT / "data" / "af_large", ROOT / "data" / "af_cross"]
CUT, SEQGAP = 8.0, 30

def _load(n, f):
    s = importlib.util.spec_from_file_location(n, Path(__file__).parent / f)
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
s16 = _load("s16", "16_structural_matched_test.py")
s17 = _load("s17", "17_disulfide_and_rank.py")


def cys_atoms(acc):
    """acc -> {resseq: (x,y,z)}，优先 SG，回退 CB、CA"""
    for d in AFD:
        f = d / f"AF-{acc}-F1-model_v6.pdb"
        if f.exists(): break
    else:
        f = None
    if f is None:
        for d in AFD:
            g = d / f"{acc}.pdb"
            if g.exists(): f = g; break
    if f is None: return {}
    best = {}
    for ln in f.read_text().splitlines():
        if not ln.startswith("ATOM") or ln[17:20].strip() != "CYS": continue
        nm = ln[12:16].strip()
        if nm not in ("SG", "CB", "CA"): continue
        r = int(ln[22:26])
        rank = {"SG": 0, "CB": 1, "CA": 2}[nm]
        if r in best and best[r][0] <= rank: continue
        best[r] = (rank, (float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return {r: xyz for r, (_, xyz) in best.items()}


def split_counts(job):
    key, acc = job
    xyz = cys_atoms(acc)
    if len(xyz) < 2: return key, {}
    pos = sorted(xyz)
    out = {}
    for p in pos:
        near = far = 0
        x1, y1, z1 = xyz[p]
        for q in pos:
            if q == p: continue
            x2, y2, z2 = xyz[q]
            if math.dist((x1, y1, z1), (x2, y2, z2)) > CUT: continue
            if abs(q - p) <= SEQGAP: near += 1
            else: far += 1
        out[p] = (near, far)
    return key, out


def datasets():
    DS = []
    osa = json.load(open(OUT / "24_mapped_osa_persulfidation.json"))["osa_persulfidation"]["proteins"]
    DS.append(("持硫化（水稻）", "肽段推断", [(p["locus"], p["acc"], p["seq"], p["modified"]) for p in osa]))
    hp = json.load(open(OUT / "33_mapped_hsa_persulfidation.json"))["hsa_persulfidation"]["proteins"]
    DS.append(("持硫化（人）", "肽段推断", [(p["acc"], p["acc"], p["seq"], p["modified"]) for p in hp]))
    accmap = json.load(open(OUT / "14_accession_map.json"))
    loc2acc = {v: k for k, v in accmap.items()}
    for t, lab in (("sno", "S-亚硝基化（拟南芥）"), ("so", "S-磺酰化（拟南芥）"),
                   ("ox", "可逆 Cys 氧化（拟南芥）")):
        ps = json.load(open(OUT / f"09_mapped_{t}.json"))["proteins"]
        DS.append((lab, "数据库标注",
                   [(p["locus"], loc2acc.get(p["locus"]), p["seq"], p["modified"])
                    for p in ps if loc2acc.get(p["locus"])]))
    cre = json.load(open(OUT / "19_mapped_cre.json"))
    for t, lab in (("sno", "S-亚硝基化（莱茵衣藻）"), ("ox", "可逆 Cys 氧化（莱茵衣藻）")):
        DS.append((lab, "数据库标注",
                   [(p["acc"], p["acc"], p["seq"], p["modified"]) for p in cre[t]["proteins"]]))
    ani = json.load(open(OUT / "19_mapped_dbptm.json"))
    for t, lab in (("S-nitrosylation", "S-亚硝基化（人鼠）"), ("Glutathionylation", "谷胱甘肽化（人鼠）")):
        DS.append((lab, "数据库标注",
                   [(p["acc"], p["acc"], p["seq"], p["modified"]) for p in ani[t]["proteins"]]))
    pal = json.load(open(OUT / "27_mapped_palmitoylation.json"))["S-palmitoylation"]["proteins"]
    DS.append(("S-棕榈酰化（人鼠）", "数据库标注",
               [(p["acc"], p["acc"], p["seq"], p["modified"]) for p in pal]))
    return DS


def main():
    DS = datasets()
    jobs = {}
    for _, _, ps in DS:
        for key, acc, _, _ in ps:
            if acc and key not in jobs: jobs[key] = (key, acc)
    jobs = list(jobs.values())
    cache = OUT / "36_ncys_split.json"
    split = json.load(open(cache)) if cache.exists() else {}
    todo = [j for j in jobs if j[0] not in split]
    print(f"需算 {len(jobs)} 个蛋白（已缓存 {len(split)}，待算 {len(todo)}）", flush=True)
    if todo:
        done = 0
        with Pool(8) as pool:
            for key, d in pool.imap_unordered(split_counts, todo, chunksize=16):
                split[key] = {str(k): v for k, v in d.items()}
                done += 1
                if done % 1500 == 0:
                    tmp = cache.with_suffix(".tmp")
                    json.dump(split, open(tmp, "w")); tmp.replace(cache)
                    print(f"  {done}/{len(todo)}", flush=True)
        tmp = cache.with_suffix(".tmp"); json.dump(split, open(tmp, "w")); tmp.replace(cache)

    rows, flat = [], []
    for lab, src, ps in DS:
        for which, gi in (("近邻 |Δpos|≤30", 0), ("远邻 |Δpos|>30", 1), ("合计（原 ncys8）", None)):
            us_b, us_s = [], []
            for key, acc, seq, mods in ps:
                d = split.get(key)
                if not d: continue
                val = {int(k): (v[gi] if gi is not None else v[0] + v[1]) for k, v in d.items()}
                for keep, bucket in ((s16.detectable_positions(seq), us_b),
                                     (tryptic_single_cys(seq), us_s)):
                    pos = sorted(p for p in val if p in keep)
                    if len(pos) < 2: continue
                    idx = {p: i for i, p in enumerate(pos)}
                    mm = {idx[p] for p in mods if p in idx}
                    if not mm or len(mm) == len(pos): continue
                    bucket.append(([float(val[p]) for p in pos], mm))
            tb, ts = s16.moment_test(us_b), s16.moment_test(us_s)
            rows.append({"数据集": lab, "阳性来源": src, "邻居类型": which,
                         "z_全部可检出": tb["z"], "obs_全部": tb["obs"], "exp_全部": tb["exp"],
                         "n位点_全部": tb["n_modified"],
                         "z_单Cys肽段": ts["z"], "n位点_单Cys": ts["n_modified"]})
            if tb["p"] is not None: flat.append((f"{lab}|{which}|b", tb["p"]))
            if ts["p"] is not None: flat.append((f"{lab}|{which}|s", ts["p"]))

    q = s16.bh(flat)
    for r in rows:
        r["q_全部可检出"] = q.get(f"{r['数据集']}|{r['邻居类型']}|b")
        r["q_单Cys肽段"] = q.get(f"{r['数据集']}|{r['邻居类型']}|s")
    cols = ["数据集","阳性来源","邻居类型","obs_全部","exp_全部","z_全部可检出","q_全部可检出",
            "n位点_全部","z_单Cys肽段","q_单Cys肽段","n位点_单Cys"]
    with open(OUT / "36_ncys_decompose.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        for r in rows: w.writerow({c: r.get(c) for c in cols})
    json.dump({"n_tests": len(flat), "rows": rows},
              open(OUT / "36_ncys_decompose.json", "w"), ensure_ascii=False, indent=1)

    print(f"\n检验 {len(flat)} 项\n")
    print(f"{'数据集':22s}{'来源':10s}{'近邻 z':>9s}{'远邻 z':>9s}{'合计 z':>9s}")
    for lab, src, _ in DS:
        rs = {r["邻居类型"]: r for r in rows if r["数据集"] == lab}
        f_ = lambda k: (f"{rs[k]['z_全部可检出']:+.2f}" if rs.get(k) and rs[k]["z_全部可检出"] is not None else "—")
        print(f"{lab:22s}{src:10s}{f_('近邻 |Δpos|≤30'):>9s}{f_('远邻 |Δpos|>30'):>9s}{f_('合计（原 ncys8）'):>9s}")


def tryptic_single_cys(seq, minl=7, maxl=30):
    cuts = [0] + [i + 1 for i in range(len(seq) - 1)
                  if seq[i] in "KR" and seq[i + 1] != "P"] + [len(seq)]
    out = set()
    for a, b in zip(cuts, cuts[1:]):
        pep = seq[a:b]
        if not (minl <= len(pep) <= maxl): continue
        cs = [a + i + 1 for i, ch in enumerate(pep) if ch == "C"]
        if len(cs) == 1: out.add(cs[0])
    return out


if __name__ == "__main__":
    main()
