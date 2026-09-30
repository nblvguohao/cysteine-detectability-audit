# -*- coding: utf-8 -*-
import csv, hashlib, json, time
from pathlib import Path
ROOT = Path("/Users/lyuguohao/Documents/Codex/https-chatgpt-com-s-t-6aa28ddd9e60819199c57fc1065a53b8/outputs/2026-09-12_persulfidation_upgrade")
R = ROOT / "results"
DROP = "external/intake/manual_fetch_2026-09-16/"

ITEMS = [
 {"priority": "1 (关掉 G1，决定性)", "tag": "yang2014_ncomms5776",
  "paper": "Yang et al. 2014, Nat Commun 5:4776, Site-specific mapping and quantification of protein S-sulphenylation",
  "doi": "10.1038/ncomms5776", "url": "https://www.nature.com/articles/ncomms5776",
  "which_file": "Supplementary Data（位点级 S-次磺酰化表；页面下方 Supplementary Information 区，注意要的是 Data 不是 Supplementary Information PDF）",
  "unlocks": "SFE-001、SFE-002、SFE-003、SFE-004、SFE-005（5 条）；**并且是 SFE-006 的阳性来源**",
  "why_it_matters": "SFE-006 是唯一判消失的论断，现在只是转移检验。拿到这张表就能在原队列上复现基线，把全文头条从转移检验升级为复现——这是整份差距表里唯一的决定性项。",
  "i_tried": "PMC 只有作者稿 PDF 补充（NIHMS615802-supplement-1.pdf），不含数据表；Springer CDN 按 ncomms/41467 两套命名各试 15×6 组合均未命中；nature.com 重定向到未放行域。",
  "drop_to": DROP + "yang2014_ncomms5776/"},
 {"priority": "1 (关掉 G1，决定性)", "tag": "bui2016_sohsite",
  "paper": "Bui et al. 2016, BMC Genomics 17:9, SOHSite",
  "doi": "10.1186/s12864-015-2299-1", "url": "https://doi.org/10.1186/s12864-015-2299-1",
  "which_file": "Additional file 6 与 7（我按 MOESM1–8 试遍六种扩展名，只有 1/2/3/4/5/8 存在且都是图与性能表，**位点表应在 6 或 7**）；若仍没有，则取 SOHSite 网站的数据集下载",
  "unlocks": "SFE-006 的作者阳性集（与 Yang 2014 互为补充）",
  "why_it_matters": "同上：把唯一的消失结论从转移检验变成复现。两篇拿到一篇即可开工，两篇都有最稳。",
  "i_tried": "已取得 MOESM1/2/3/4/5/8（共 6 个 .docx，在 external/intake/phase2_g1_supp/），逐个解包核过——全是流程图、五折性能表、InterPro 结构域分布，无位点级表。",
  "drop_to": DROP + "bui2016_sohsite/"},
 {"priority": "2 (关掉 G2 覆盖，9 条)", "tag": "doulias2010_pnas",
  "paper": "Doulias et al. 2010, PNAS 107:16958, Structural profiling of endogenous S-nitrosocysteine residues",
  "doi": "10.1073/pnas.1008036107", "url": "https://www.pnas.org/doi/10.1073/pnas.1008036107",
  "which_file": "SI Appendix（含位点表）以及页面上单列的 Dataset/Table 文件（若有）",
  "unlocks": "SNO-001 至 SNO-009（9 条，本清单里单篇解锁最多）",
  "why_it_matters": "S-亚硝基化是唯一跨两个物种的家族，这 9 条是该家族论断的主体。缺了它，覆盖停在 22/37。",
  "i_tried": "pnas.org 已申请放行并获批，但出版社对自动客户端返回 HTTP 403——放行域名解决不了，必须浏览器手动下载。",
  "drop_to": DROP + "doulias2010_pnas/"},
 {"priority": "2 (关掉 G2 覆盖，6 条)", "tag": "marino2010_jmb",
  "paper": "Marino & Gladyshev 2010, J Mol Biol 395:844, Structural analysis of cysteine S-nitrosylation",
  "doi": "10.1016/j.jmb.2009.10.042", "url": "https://doi.org/10.1016/j.jmb.2009.10.042",
  "which_file": "Supplementary Table S1（Elsevier；PMC 只有作者稿、不含补充表，**需机构订阅**）",
  "unlocks": "SNO-010 至 SNO-015（6 条）",
  "why_it_matters": "这 6 条是'酸碱基序'那一支的原始论断，属于本领域被引最多的位点偏好陈述之一。",
  "i_tried": "PMC 作者稿无补充表；Elsevier 正文非开放获取。未尝试镜像站。",
  "drop_to": DROP + "marino2010_jmb/"},
 {"priority": "3 (补齐第三类化学)", "tag": "akter2018_nchembio",
  "paper": "Akter et al. 2018, Nat Chem Biol 14:995, DiaAlk chemical proteomics of cysteine sulfinic acid",
  "doi": "10.1038/s41589-018-0116-2", "url": "https://www.nature.com/articles/s41589-018-0116-2",
  "which_file": "Supplementary Dataset 1、2、3（.xlsx；PMC 里的确切文件名为 NIHMS1500332-supplement-Supplementary_Dataset_1/2/3.xlsx）",
  "unlocks": "SFI-001、SFI-002（2 条）",
  "why_it_matters": "亚磺酰化(-SO2H)一支目前只有 1 条对照被召回，结论只能写候选；补上这两条能让该家族至少有可再检验的实体。",
  "i_tried": "PMC 列出了三个 .xlsx 的确切文件名，但 /bin/ 路径返回 1.8 KB 下载拦截页；Springer CDN 按 41589 命名试遍未命中。",
  "drop_to": DROP + "akter2018_nchembio/"},
 {"priority": "3 (持硫化侧)", "tag": "li2024_rsc_chembiol",
  "paper": "Li et al. 2024, RSC Chem Biol, 锌指蛋白持硫化",
  "doi": "10.1039/d3cb00106g", "url": "https://doi.org/10.1039/d3cb00106g",
  "which_file": "ESI Table S3（持硫化锌指蛋白名单；金色开放获取，浏览器应可直接下载）",
  "unlocks": "PERS-007、PERS-008（2 条）",
  "why_it_matters": "这两条是持硫化家族里少数带明确属性的论断，且该 ESI 汇编自 10 个已发表持硫化数据集，对 G3 的数据集搜寻也有用。",
  "i_tried": "pubs.rsc.org 与 www.rsc.org 都已申请放行并获批：落地页对自动客户端 403，suppdata 侧 TLS 握手失败（SSL UNEXPECTED_EOF）。两条路都走不通。",
  "drop_to": DROP + "li2024_rsc_chembiol/"},
]

NOT_ACCESS = [
 {"item": "SFE-012 (Lu 2023, Microbiol Spectr)、SFI-006 (Castro 2024, Redox Biol)、SFI-003 (Wood 2003, Science)",
  "why": "**不是访问权问题**：这三篇根本没有位点级表（Wood 2003 与 Castro 2024 是结构/酶学论文，Lu 2023 的 hasSuppl=N）。不要为它们花时间；文中作为'连再检验都不可能'的证据使用。"},
 {"item": "SNO-018/019/020 (dbSNO 2.0)、SFE-005 (Yang 2014)、PERS-004 (Longen 2016)",
  "why": "**卡点是阴性集未定义，不是取不到**（PERS-004 的表本来就在手）。下载解决不了：要么联系作者问他们的阴性集，要么由我们自己定义——但那样检的就不是作者的论断了。建议在文中单列为'阴性集未声明'一类。"},
]

rows = []
for it in ITEMS:
    rows.append({k: it[k] for k in ("priority","tag","paper","doi","url","which_file","unlocks","why_it_matters","i_tried","drop_to")})
out = R / "manual_fetch_list.csv"
with out.open("w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
audit = {
 "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
 "new_computation": False,
 "items": len(rows),
 "claims_unlocked_if_all_fetched": 24,
 "network_attempts_this_round": {
   "egress_restored": "pypi.org 200, ebi.ac.uk 200 (previous round had no outbound HTTPS at all)",
   "domains_requested_and_granted": ["pubs.rsc.org", "www.rsc.org"],
   "still_unreachable_after_grant": {"pubs.rsc.org": "HTTP 403 to automated clients",
                                     "www.rsc.org": "TLS handshake failure (SSL UNEXPECTED_EOF)"},
   "pnas.org": "granted in an earlier round; publisher returns HTTP 403 to automated clients",
   "pmc_attachments": "PMC /bin/ paths return a 1.8 KB download interstitial for these records",
   "no_user_agent_spoofing_no_mirrors": True,
 },
 "partial_material_obtained": {
   "bui2016_sohsite": ["MOESM1","MOESM2","MOESM3","MOESM4","MOESM5","MOESM8"],
   "bui2016_contains_site_table": False,
   "note": "unpacked and checked: flowcharts, five-fold performance tables and an InterPro domain distribution; no site-level table",
 },
 "not_an_access_problem": NOT_ACCESS,
 "drop_directory": DROP,
 "outputs": {"results/manual_fetch_list.csv": sha(out)},
}
(R / "manual_fetch_list_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
print("items:", len(rows), "| claims unlocked if all fetched:", audit["claims_unlocked_if_all_fetched"])
for r in rows:
    print(f"  [{r['priority'][:14]:14s}] {r['tag']:22s} -> {r['unlocks'][:44]}")
