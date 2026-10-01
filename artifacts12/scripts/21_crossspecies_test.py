"""跨物种复制 · 把方向一的结构规律拿到别的物种上重跑

七套数据，三个谱系：
  拟南芥（陆生植物）  S-亚硝基化 / S-磺酰化 / 可逆 Cys 氧化   —— 已在 15 号脚本算过特征
  莱茵衣藻（绿藻）    S-亚硝基化 / 可逆 Cys 氧化              —— PlantPTMViewer，与拟南芥同库同流程
  人 + 小鼠（动物）   S-亚硝基化 / 谷胱甘肽化                  —— dbPTM 实验验证位点

检验设计与 16/17/18 完全一致：蛋白内配对 + 胰蛋白酶可检出性匹配，
连续特征同时报原值与蛋白内秩次，全部检验统一 BH 校正。
问题只有一个：**「修饰位点回避无序区」等规律，是不是只在拟南芥成立。**
"""
import json, math, warnings
from pathlib import Path
from multiprocessing import Pool
import importlib.util
import numpy as np

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"

def _load(n, f):
    s = importlib.util.spec_from_file_location(n, Path(__file__).parent / f)
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
s16 = _load("s16", "16_structural_matched_test.py")
s17 = _load("s17", "17_disulfide_and_rank.py")

MAXASA_CYS = 167.0
BASIC = {"ARG": ("NE","NH1","NH2"), "LYS": ("NZ",), "HIS": ("ND1","NE2")}
ACID  = {"ASP": ("OD1","OD2"), "GLU": ("OE1","OE2")}
R8, SSDIST = 8.0, 2.5


def feat_one(job):
    key, acc, seq, afdir = job
    f = Path(afdir) / f"AF-{acc}-F1-model_v6.pdb"
    if not f.exists(): return key, {"status": "no_model", "acc": acc}
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
        import pydssp
        st = PDBParser(QUIET=True).get_structure("m", str(f)); model = st[0]
        res = [r for r in list(model)[0] if r.id[0] == " "]
        if len(res) != len(seq): return key, {"status": "len_mismatch", "acc": acc}
        ShrakeRupley().compute(model, level="R")
        ss = pydssp.assign(pydssp.read_pdbtext(f.read_text()), out_type="c3")

        import propka.run, tempfile, shutil, os
        td = Path(tempfile.mkdtemp(prefix="pk_")); cwd = os.getcwd(); pka = {}
        try:
            shutil.copy(f, td / "m.pdb"); os.chdir(td)
            propka.run.single("m.pdb", optargs=("--quiet",))
            insum = False
            for line in (td / "m.pka").read_text(errors="ignore").splitlines():
                if line.startswith("SUMMARY OF THIS PREDICTION"): insum = True; continue
                if insum:
                    p = line.split()
                    if len(p) >= 4 and p[0] == "CYS":
                        try: pka[int(p[1])] = float(p[3])
                        except Exception: pass
                    elif line.startswith("-") or line.startswith("Free energy"): break
        except Exception: pass
        finally:
            os.chdir(cwd); shutil.rmtree(td, ignore_errors=True)

        bpos = [a.coord for r in res if r.get_resname() in BASIC for a in r if a.get_id() in BASIC[r.get_resname()]]
        apos = [a.coord for r in res if r.get_resname() in ACID  for a in r if a.get_id() in ACID[r.get_resname()]]
        bpos = np.asarray(bpos) if bpos else np.zeros((0,3))
        apos = np.asarray(apos) if apos else np.zeros((0,3))
        sg = {r.id[1]: r["SG"].coord for r in res if r.get_resname() == "CYS" and "SG" in r}
        ssb = {}
        if len(sg) >= 2:
            ps = sorted(sg); arr = np.array([sg[p] for p in ps])
            d = np.linalg.norm(arr[:,None,:] - arr[None,:,:], axis=-1); np.fill_diagonal(d, 1e9)
            ssb = {p: int(d[i].min() <= SSDIST) for i, p in enumerate(ps)}

        out = []
        for i, r in enumerate(res):
            if r.get_resname() != "CYS": continue
            pos = r.id[1]
            ref = sg.get(pos, r["CA"].coord if "CA" in r else None)
            if ref is None: continue
            out.append({"pos": pos, "rel_sasa": round(float(r.sasa)/MAXASA_CYS, 4),
                        "ss3": str(ss[i]), "plddt": round(float(np.mean([a.get_bfactor() for a in r])), 2),
                        "pka": pka.get(pos), "ssbond": ssb.get(pos, 0),
                        "nbasic": int((np.linalg.norm(bpos-ref,axis=1)<=R8).sum()) if len(bpos) else 0,
                        "nacid":  int((np.linalg.norm(apos-ref,axis=1)<=R8).sum()) if len(apos) else 0,
                        "ncys8": sum(1 for p2,c2 in sg.items() if p2!=pos and float(np.linalg.norm(c2-ref))<=R8)})
        return key, {"status": "ok", "acc": acc, "cys": out}
    except Exception as e:
        return key, {"status": "error", "acc": acc, "err": f"{type(e).__name__}: {e}"}


FEATS_BIN = {
    "落在无序区 (pLDDT<70)": lambda c: float(c["plddt"] < 70),
    "完全埋藏 (rel_sasa<=0.05)": lambda c: float(c["rel_sasa"] <= 0.05),
    "在 loop/coil 上": lambda c: float(c["ss3"] == "-"),
    "在 sheet 上": lambda c: float(c["ss3"] == "E"),
    "在预测二硫键中": lambda c: float(c["ssbond"]),
}
FEATS_CONT = {"相对可及性": "rel_sasa", "pLDDT": "plddt",
              "8A 内其他 Cys 数": "ncys8", "8A 内酸性残基数": "nacid",
              "pKa（剔除二硫键）": "pka"}


def units(prots, feats, valfun, rank=False, skip_ss=False):
    us = []
    for p in prots:
        fr = feats.get(p["key"])
        if not fr or fr["status"] != "ok": continue
        fmap = {c["pos"]: c for c in fr["cys"]}
        keep = s16.detectable_positions(p["seq"])
        pos = sorted(q for q in fmap if q in keep and valfun(fmap[q]) is not None
                     and not (skip_ss and fmap[q]["ssbond"]))
        if len(pos) < 2: continue
        idx = {q: i for i, q in enumerate(pos)}
        mm = {idx[q] for q in p["modified"] if q in idx}
        if not mm or len(mm) == len(pos): continue
        v = [float(valfun(fmap[q])) for q in pos]
        us.append((s17.norm_rank(v) if rank else v, mm))
    return us


def main():
    ath_feat = json.load(open(OUT / "15_structural_features.json"))["proteins"]
    cre = json.load(open(OUT / "19_mapped_cre.json"))
    ani = json.load(open(OUT / "19_mapped_dbptm.json"))
    osa = json.load(open(OUT / "24_mapped_osa_persulfidation.json"))

    DS = []
    # ★ 持硫化本身（水稻，PRIDE PXD072089 / PNAS 2026）——第一次不靠替代修饰
    DS.append(("水稻（持硫化）", "持硫化",
               [{"key": p["locus"], "acc": p["acc"], "seq": p["seq"], "modified": p["modified"]}
                for p in osa["osa_persulfidation"]["proteins"]], ROOT/"data"/"af_cross"))
    for t, lab in (("sno","S-亚硝基化"),("so","S-磺酰化"),("ox","可逆 Cys 氧化")):
        ps = [{"key": p["locus"], "acc": None, "seq": p["seq"], "modified": p["modified"]}
              for p in json.load(open(OUT / f"09_mapped_{t}.json"))["proteins"]]
        DS.append(("拟南芥（陆生植物）", lab, ps, None))
    for t, lab in (("sno","S-亚硝基化"),("ox","可逆 Cys 氧化")):
        ps = [{"key": p["acc"], "acc": p["acc"], "seq": p["seq"], "modified": p["modified"]}
              for p in cre[t]["proteins"]]
        DS.append(("莱茵衣藻（绿藻）", lab, ps, ROOT/"data"/"af_cross"))
    for t, lab in (("S-nitrosylation","S-亚硝基化"),("Glutathionylation","谷胱甘肽化")):
        ps = [{"key": p["acc"], "acc": p["acc"], "seq": p["seq"], "modified": p["modified"]}
              for p in ani[t]["proteins"]]
        DS.append(("人 + 小鼠（动物）", lab, ps, ROOT/"data"/"af_cross"))
    # ★ S-棕榈酰化（+238 Da 脂酰基）——与同库同物种的 S-亚硝基化(+29)、谷胱甘肽化(+305)
    #   构成「修饰基团体积」三点梯度，用来检验「大基团才需要暴露」这一假说
    palf = OUT / "27_mapped_palmitoylation.json"
    if palf.exists():
        ps = [{"key": p["acc"], "acc": p["acc"], "seq": p["seq"], "modified": p["modified"]}
              for p in json.load(open(palf))["S-palmitoylation"]["proteins"]]
        DS.append(("人 + 小鼠（动物）", "S-棕榈酰化", ps, ROOT/"data"/"af_cross"))
    # ★ 持硫化（人，PRIDE PXD044043 BTA 富集）——持硫化的第二个物种，用于跨物种复制
    hpf = OUT / "33_mapped_hsa_persulfidation.json"
    if hpf.exists():
        ps = [{"key": p["acc"], "acc": p["acc"], "seq": p["seq"], "modified": p["modified"]}
              for p in json.load(open(hpf))["hsa_persulfidation"]["proteins"]]
        DS.append(("人（持硫化）", "持硫化", ps, ROOT/"data"/"af_cross"))

    # --- 算新物种的结构特征（拟南芥直接复用 15 号结果，并补上二硫键标记）---
    cache = OUT / "21_features_cross.json"
    newf = json.load(open(cache)) if cache.exists() else {}
    jobs = {}
    for _, _, ps, afd in DS:
        if afd is None: continue
        for p in ps:
            if p["key"] not in newf: jobs[p["key"]] = (p["key"], p["acc"], p["seq"], str(afd))
    jobs = list(jobs.values())
    if jobs:
        print(f"待算新物种蛋白 {len(jobs)} 个（已缓存 {len(newf)}）", flush=True)
        done = 0
        with Pool(processes=8) as pool:
            for k, r in pool.imap_unordered(feat_one, jobs, chunksize=4):
                newf[k] = r; done += 1
                if done % 200 == 0:                      # 增量落盘，被中断可续跑
                    tmp = cache.with_suffix(".tmp")
                    json.dump(newf, open(tmp, "w"), ensure_ascii=False)
                    tmp.replace(cache)
                    print(f"  {done}/{len(jobs)}", flush=True)
        tmp = cache.with_suffix(".tmp")
        json.dump(newf, open(tmp, "w"), ensure_ascii=False); tmp.replace(cache)
    tally = {}
    for r in newf.values(): tally[r["status"]] = tally.get(r["status"], 0) + 1
    print("新物种特征状态:", tally, flush=True)

    for loc, fr in ath_feat.items():
        if fr.get("status") != "ok": continue
        if "ssbond" not in fr["cys"][0] if fr["cys"] else False: pass
        fl = s17.disulfide_flags(fr["acc"])
        for c in fr["cys"]: c["ssbond"] = fl.get(c["pos"], 0)

    rows, flat = [], []
    for lineage, mod, ps, afd in DS:
        feats = ath_feat if afd is None else newf
        n_ok = sum(1 for p in ps if feats.get(p["key"], {}).get("status") == "ok")
        for name, g in FEATS_BIN.items():
            m = s16.moment_test(units(ps, feats, g))
            rows.append({"谱系": lineage, "修饰": mod, "特征": name, "类型": "二分",
                         "obs": m["obs"], "exp": m["exp"], "差pp": round(m["diff"]*100, 2) if m["diff"] is not None else None,
                         "z_原值": m["z"], "z_秩次": None, "p_原值": m["p"],
                         "n_位点": m["n_modified"], "n_蛋白": m["n_proteins"], "n_有模型": n_ok})
            flat.append((f"{lineage}|{mod}|{name}|raw", m["p"]))
        for name, key in FEATS_CONT.items():
            skip = (key == "pka")
            g = lambda c, k=key: c.get(k)
            m = s16.moment_test(units(ps, feats, g, False, skip))
            r = s16.moment_test(units(ps, feats, g, True, skip))
            rows.append({"谱系": lineage, "修饰": mod, "特征": name, "类型": "连续",
                         "obs": m["obs"], "exp": m["exp"], "差pp": None,
                         "z_原值": m["z"], "z_秩次": r["z"], "p_原值": m["p"], "p_秩次": r["p"],
                         "n_位点": m["n_modified"], "n_蛋白": m["n_proteins"], "n_有模型": n_ok})
            flat += [(f"{lineage}|{mod}|{name}|raw", m["p"]), (f"{lineage}|{mod}|{name}|rank", r["p"])]

    q = s16.bh(flat)
    for r in rows:
        base = f"{r['谱系']}|{r['修饰']}|{r['特征']}"
        r["q_原值"] = q[f"{base}|raw"]; r["q_秩次"] = q.get(f"{base}|rank")
        if r["类型"] == "二分":
            r["判定"] = "显著" if r["q_原值"] < 0.05 else "阴性"
        else:
            ok = r["q_原值"] < 0.05 and r["q_秩次"] is not None and r["q_秩次"] < 0.05 \
                 and r["z_原值"] * r["z_秩次"] > 0
            r["判定"] = "稳健" if ok else ("脆弱（仅原值）" if r["q_原值"] < 0.05 else "阴性")

    import csv as _csv
    cols = ["谱系","修饰","特征","类型","obs","exp","差pp","z_原值","z_秩次","q_原值","q_秩次","n_位点","n_蛋白","n_有模型","判定"]
    with open(OUT / "21_crossspecies_summary.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = _csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow({c: r.get(c) for c in cols})
    json.dump({"n_tests": len(flat), "rows": rows}, open(OUT / "21_crossspecies_summary.json","w"),
              ensure_ascii=False, indent=1)

    print(f"\n检验 {len(flat)} 项\n")
    hdr = f"{'特征':24s}" + "".join(f"{l[:2]+'·'+m[:5]:>16s}" for l,m,_,_ in DS)
    print(hdr)
    for name in list(FEATS_BIN) + list(FEATS_CONT):
        cells = []
        for lineage, mod, _, _ in DS:
            r = next(x for x in rows if x["谱系"]==lineage and x["修饰"]==mod and x["特征"]==name)
            v = f"{r['差pp']:+.1f}pp" if r["差pp"] is not None else f"z{r['z_原值']:+.1f}"
            mark = "*" if r["q_原值"] < 0.05 else " "
            cells.append(f"{v}{mark}".rjust(16))
        print(f"{name:24s}" + "".join(cells))


if __name__ == "__main__":
    main()
