# -*- coding: utf-8 -*-
import csv, json
from pathlib import Path
ROOT = Path("/Users/lyuguohao/Documents/Codex/https-chatgpt-com-s-t-6aa28ddd9e60819199c57fc1065a53b8/outputs/2026-09-12_persulfidation_upgrade")
a = json.loads((ROOT/"results/manual_fetch_intake_audit.json").read_text(encoding="utf-8"))
rows = list(csv.DictReader((ROOT/"results/manual_fetch_intake.csv").open(encoding="utf-8-sig")))
by = {r["group"]: r for r in rows}

TEXT = f"""# 人工下载材料入库与结构核验（2026-09-16）

{a['files_ingested']} 个文件、{round(a['total_bytes']/1e6,1)} MB，已按原名存入 `{a['drop_directory']}` 下各自子目录。
**每个文件都被打开读过结构，不是按文件名推断**（方法见第五节）。
清单与逐文件哈希：`results/manual_fetch_intake.csv`、`results/manual_fetch_intake_files.csv`、
`results/manual_fetch_intake_audit.json`。本轮不含新计算、未做任何再检验。

## 一 结论：{a['n_unlocked']} 条论断现在可以开工，8 条仍差位点级材料

清单原本预期解锁 23 条。实到之后逐个核结构，**{a['n_unlocked']} 条确认可开工**，
差额全部来自同一篇（Doulias 2010），原因是材料本身的粒度，不是下载出了问题。

| 论文 | 文件 | 位点级？ | 位点表在哪 | 解锁 |
|---|---:|---|---|---|
| Yang 2014 | {by['yang2014']['n_files']} | **是** | `{by['yang2014']['site_table_location']}` | SFE-001～004 |
| Marino & Gladyshev 2010 | {by['marino2010']['n_files']} | **是** | `{by['marino2010']['site_table_location']}` | SNO-010～015 |
| Akter 2018 | {by['akter2018']['n_files']} | **是** | `{by['akter2018']['site_table_location']}` | SFI-001、SFI-002 |
| Li 2024 | {by['li2024']['n_files']} | 否（蛋白级，与论断单位一致） | `{by['li2024']['site_table_location']}` | PERS-007、PERS-008 |
| Doulias 2010 | {by['doulias2010']['n_files']} | **部分** | `{by['doulias2010']['site_table_location']}` | 仅 SNO-008 |
| Bui 2016 | {by['bui2016']['n_files']} | 否 | — | 无 |

## 二 G1 解锁了——头条可以从转移检验升级为原队列复现

**`41467_2014_BFncomms5776_MOESM785_ESM.xlsx`（Yang 2014 Dataset 2）就是那张位点表**：
`IDs from SulfenM` 953 行、`IDs from SulfenQ` 517 行、`Total` 1106 行，
列含 `ID / Gene Name / Protein Description / Modified Site / Modification / Peptide Sequence`。

这同时是 **SFE-006 的作者阳性来源**（SOHSite 的阳性取自 Yang 2014 与 RedoxDB）。
也就是说，当前唯一判"消失"的那条论断，现在可以在**作者自己的队列上**复现基线再施加控制——
这正是差距登记表里 G1 那一条，也是子刊审稿人会第一个要的东西。

顺带澄清 Bui 2016：八个附件我全部解包核过（S1 前 20 理化性质、S2 五折性能、S3 独立测试、
S4 GO、S5 KEGG、S6 InterPro，另两个是流程图），**确实没有位点表**——
之前取不到的 MOESM6/7 也只是 GO 与 KEGG 注释表。但这不影响 G1，因为 Yang 2014 已经够了。

## 三 Doulias 2010：粒度不够，需要你再花一分钟

四张表都是**蛋白级**（`Protein name / Uniprot Accession / MW / Unique Peptides / Peptide sequence`）。
肽段列里内嵌了 Cys 位置（例如 `C473VAYAESHDQALVGDK` 即 Cys473），
但本机可用的转换器只还原出 **23 个唯一 (登录号, 位置) 对**，而 SNO-001 作者报的阳性是 **139**。
`st03` 是金属配位 Cys（7 个蛋白），`st04` 是 Trx 敏感蛋白（78 行，无位点）。

因此：**SNO-008（功能富集，蛋白级）现在就能做**——合计 200 个唯一登录号足够；
**SNO-001～007 与 SNO-009 是位点级论断（二级结构、SASA、pKa、侧翼组成、基序），还不能做**。

最省事的补法：**用 Word 打开 `st01.doc` 与 `st02.doc`，另存为 `.docx` 或 `.xlsx`**，
表格结构就能完整读出。或者提供该文的 SI Appendix / Dataset——139 个位点应当在那里。

## 四 其余三篇的核验细节

- **Marino & Gladyshev 2010**：`mmc1.pdf` 第 **8–11 页**是 Table S1「NO-Cys 数据集」（Source / ID /
  Protein name / NO-Cys 位置），第 **12–15 页**是 Table S2「NO-Cys 周边暴露区」。
  你另外传的 `main.pdf` 是该文正文（Table 1/2/3）。需要从 PDF 解析成表，下一轮做。
  **订阅材料，按纪律不分发**，仅本地再检验用。
- **Akter 2018**：`MOESM33` 是主位点表（A549 1174 行、HeLa 1099 行，
  `Uniprot Accession # / Site # / Modified sequence / Mean Ratio`），
  另含 `A549_SOH_vs_SO2H` 103 行与 `HeLa_SOH_vs_SO2H` 89 行——
  这两张可直接支撑次磺酰化与亚磺酰化的对比。`MOESM32` 为 243/221 行，`MOESM34` 为 Srx WT/KO。
- **Li 2024**：`d3cb00106g2_suppl.zip` 里的 Table S1（Total proteome 5928 行、
  Persulfidated proteins 2436 行）与 Table S2（MEF Total ZFs 472、SSH ZFs 119，
  按 C4/C2H2/CCCH 分型 68/32/57）。是蛋白级，但 PERS-007/008 本就是蛋白级论断，**单位匹配**。

## 五 核验方法，以及本轮我自己纠正的一处

xlsx 用 openpyxl 读工作表名、行列数与表头；docx 解包读 `word/document.xml` 的表格行；
doc 用 textutil 转换后按 `\\x07` 单元格分隔符解析（txt/html/docx 三条路径都试过，
antiword / catdoc / libreoffice 本机均不存在）；pdf 用 pypdfium2 逐页抽文本并统计登录号与 Cys 记号；
zip 列出条目并解出其中的 xlsx 再同样处理。

**一处自纠**：首次统计 Cys 记号时，我用的模式要求数字后有词边界，
于是 `C473VAYAESHDQALVGDK` 这种"位置后紧跟序列字母"的写法被整批漏掉，
一度把 PNAS 那几张表判成完全没有位点信息。改用 `C(\\d{{1,4}})(?=[A-Z])` 重数后，
结论修正为**部分可得（23 对）**。这与本项目此前两次深度定义错误是同一类：
先看真实字符串长什么样，再写判据。
"""
p = ROOT/"reports"/"MANUAL_FETCH_INTAKE_2026-09-16.md"
p.write_text(TEXT, encoding="utf-8")
print("wrote", p.name, len(TEXT))
