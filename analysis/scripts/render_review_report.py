"""Render the curated review report (Chinese) from the curated table.

Inputs : results/tables/intratumoral_bacteria_bioinformatics_curated.csv
         data/processed/curated_authors.json
         analysis/out/doi_verification.json
Output : literature/intratumoral-bacteria-bioinformatics-review.md
"""

from __future__ import annotations

import csv
import json
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CUR = ROOT / "results" / "tables" / "intratumoral_bacteria_bioinformatics_curated.csv"
AUTH = json.loads((ROOT / "data" / "processed" / "curated_authors.json").read_text(encoding="utf-8"))
VER = json.loads((ROOT / "analysis" / "out" / "doi_verification.json").read_text(encoding="utf-8"))
OUT = ROOT / "literature" / "intratumoral-bacteria-bioinformatics-review.md"
OUT.parent.mkdir(parents=True, exist_ok=True)

rows = list(csv.DictReader(CUR.open(encoding="utf-8-sig")))


def surname(pmid: str) -> str:
    name = AUTH.get(pmid, "")
    return name.split()[0] if name else "—（无作者字段）"


def remark(r: dict) -> str:
    notes = []
    if r["doi"] and isinstance(VER.get(r["doi"]), dict) and VER[r["doi"]].get("retracted"):
        notes.append("**已撤稿**")
    if "bioRxiv" in r["journal"] or "medRxiv" in r["journal"]:
        notes.append("预印本")
    if r["article_type"] == "综述/评论":
        notes.append("综述/评论")
    return "；".join(notes)


def link(r: dict) -> str:
    t = r["title"].rstrip(".")
    if r["doi"]:
        return f"[{t}](https://doi.org/{r['doi']})"
    return t


def table(theme: str) -> str:
    head = ("| 年份 | 第一作者 | 期刊 | 标题 | 摘要中提及的公共数据 | 摘要含湿实验关键词 | 备注 | PMID |\n"
            "|---|---|---|---|---|---|---|---|\n")
    lines = []
    for r in rows:
        if r["theme"] != theme:
            continue
        wet = "是" if r["wetlab_mention_in_abstract"] == "有提及" else "否"
        lines.append(
            f"| {r['year']} | {surname(r['pmid'])} | {r['journal']} | {link(r)} | "
            f"{r['data_mentioned_in_abstract']} | {wet} | {remark(r)} | "
            f"[{r['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/) |")
    return head + "\n".join(lines) + "\n"


themes = list(OrderedDict.fromkeys(r["theme"] for r in rows))
n = len(rows)
retracted = [k for k, v in VER.items() if isinstance(v, dict) and v.get("retracted")]

doc = f"""# 瘤内菌（intratumoral bacteria / tumor microbiome）生物信息学分析文献清单

检索日期：2026-10-15 ｜ 数据库：PubMed（NCBI E-utilities）｜ 人工归类的候选文献 {n} 篇

---

## 一、先看这三条结论

**1. 这个方向的"纯生信 + 复用 TCGA/GEO"论文数量已经很大，但可以清楚分成两类。** 一类是**方法学与资源类**工作：从公共测序数据里把微生物读段"读"出来、去污染、建模、做成图谱或在线工具。另一类是**下游应用类**：拿现成的 TCGA 微生物丰度矩阵去做预后 signature、免疫浸润、免疫治疗反应预测。后者数量最多，2023 年起在 Frontiers、Scientific Reports、Microbiology Spectrum 这一类期刊上呈井喷式增长，典型套路是"TCGA 训练 + GEO 验证 + LASSO/Cox 建模 + CIBERSORT/ssGSEA + 列线图"。

**2. 引用前必须知道的坏消息：这个领域最重要的两篇 TCGA 微生物组论文已被撤稿。** Nature 2020 的 Poore 等（[Poore 2020](https://doi.org/10.1038/s41586-020-2095-1)）和 Nature Communications 2022 的 Hermida 等（[Hermida 2022](https://doi.org/10.1038/s41467-022-30512-3)）在 Crossref 上都标注为 RETRACTED。这两篇正是大量下游论文的数据来源与引用对象。

**3. 2023 年之后出现了一条独立且必须阅读的方法学质疑线。** Gihawi 等在 mBio 上以再分析证明"多数被报道的肿瘤相关微生物实际上根本不在样本里"，并指出 Poore 等及十余篇后续研究很可能无效（[Gihawi 2023](https://doi.org/10.1128/mbio.01607-23)）；同一团队随后把 TCGA 全部 5734 例全基因组测序样本重新分析，公开了细菌、病毒、古菌、真菌的读段计数供他人核查（[Gihawi 2025](https://doi.org/10.1126/scitranslmed.ads6335)）。对 Genomics England 8908 例 22 种癌症的分析只找到结直肠癌有可靠可区分的微生物群落（[Gihawi 2025](https://doi.org/10.1126/scitranslmed.ads6166)）；对英国 100,000 基因组计划 16369 例高深度肿瘤全基因组的分析同样发现，除口咽-消化道肿瘤（细菌、真菌、病毒、古菌乃至滴虫的多界群落）外，多数癌种的微生物信号与背景不可区分（[Dohlman 2026](https://doi.org/10.1016/j.cell.2026.04.015)）。Cancer Discovery 的新闻稿把这组结果概括为"微生物 DNA 并非多数肿瘤的普遍特征，结直肠癌是例外"（[Cancer Discov 2026](https://doi.org/10.1158/2159-8290.CD-NW2025-0104)）。

**这对你意味着什么**：如果你要做瘤内菌生信，肠癌/口腔癌/胃癌这类口咽-消化道肿瘤是高价值且可信度较高的方向；泛癌或低微生物量癌种（如胰腺、肺、乳腺的"signature 建模"）即使文章发得很多，方法学风险也显著更高，投稿或引用时最好主动交代去污染流程与阴性对照。相关的操作规范见本清单第 2 组。

---

## 二、方法学可靠性与争议（引用 TCGA 微生物结果前必读）

{table("2 方法学可靠性与争议")}

## 三、数据资源与分析工具

这一组是做瘤内菌生信最值得先读的，因为它们决定了"丰度矩阵从哪来、可不可信"。**注意：下表中 2020 年 Nature 与 2022 年 Nature Communications 那两篇已被撤稿（表中已标注），但它们在二手文献里被引用的频率极高，看到别人引用时要留个心。**

{table("1 数据资源与分析工具")}

## 四、泛癌 / 跨癌种生信分析

{table("3 泛癌/跨癌种生信分析")}

## 五、按癌种的应用类论文

### 5.1 肺癌与胸部肿瘤

{table("4 肺癌/胸部肿瘤")}

### 5.2 消化道肿瘤（食管、胃、结直肠、肝、胆胰）

{table("5 消化道肿瘤")}

### 5.3 头颈/口腔/鼻咽、神经、泌尿生殖与其他

{table("6 头颈/口腔/鼻咽/神经与泌尿生殖系统")}

## 六、综述与系统评价

{table("7 综述/系统评价/立场文章")}

## 七、做选题时的几个实操提示

**最常被复用的四类公共数据入口。** 一是 TCGA 经 cBioPortal 直接分发的微生物丰度矩阵（多篇论文用的是同一套"TCGA 微生物组"矩阵）；二是 Cancer Microbiome Atlas（TCMA），它专门处理污染问题，是判断"组织定植菌 vs 污染物"的对照资源（[Dohlman 2021](https://doi.org/10.1016/j.chom.2020.12.001)）；三是 GEO 里的肿瘤转录组/单细胞数据，用 Vulture、microbioTA 这类工具自己重新比对调用微生物（[Vulture 2024](https://doi.org/10.1093/gigascience/giad117)、[microbioTA 2023](https://doi.org/10.1093/nar/gkac851)）；四是 PRISM 这类新方法，同时产出宿主与微生物信号（[PRISM 2026](https://doi.org/10.1016/j.ccell.2026.01.007)）。

**三个最容易翻车的地方。** 宿主读段残留与试剂污染（HumanFilt 等宿主去除流程值得照抄，[HumanFilt 2026](https://doi.org/10.1128/msystems.00550-26)）；成分数据（compositional）变换会人为制造表型相关（有专门论文讨论这个争点，见第 2 组中 2025 年 mSystems 那篇）；以及下游关联检验的统计可靠性——2026 年 mSystems 的系统评估发现，TCGA 微生物谱与"细胞类型比例""患者生存"的关联在置换检验下多为统计假象，与基因表达的关联则中等可靠（[mSystems 2026](https://doi.org/10.1128/msystems.00180-26)）。这条结论直接决定了：**基于 TCGA 微生物丰度做的生存 signature，正是可靠性最低的那一类下游产物。**

**发表层面的现实。** 如果你想快速出题，第 5 组里"某癌种 + TCGA 微生物 + 免疫浸润 + 免疫治疗队列验证（IMvigor210、GDSC）"的套路在 2023–2026 年仍然可发表；但如果审稿人熟悉第 2 组文献，纯 signature 类工作会被追问去污染与阴性对照。稳妥组合是：口咽-消化道癌种 + 严格的宿主去除与污染对照 + 一小组实验验证（菌群原位检测或细胞/类器官干预），这正是你在需求里提到的"加了小部分湿实验"。

---

## 八、检索与筛选过程（可复现）

- 检索式共 7 条主题式查询，覆盖"intratumoral bacteria / tumor microbiome / cancer microbiome"等核心词与"TCGA / GEO / CPTAC / ICGC / bioinformatics / signature / deconvolution / single-cell"等数据与方法词，去重后得 2765 条 PubMed 记录。
- 机器初筛（必须同时命中"微生物主题词"和"生信方法词"）得到 1060 条候选，再按主题相关性/是否复用公共数据/是否综述打分排序，人工归类出上述 {n} 篇。
- 表中的"摘要中提及的公共数据"和"摘要含湿实验关键词"两列由脚本从 PubMed 摘要文本机械抽取，没有读全文。两列都可能不准：前者会漏（摘要未提 GEO 而正文用了 GEO），后者会高估（"we performed" 之类用于描述计算分析的措辞会被误判为湿实验）。判断某篇是否"纯生信"仍需看原文方法部分。
- 所有 DOI 均经 Crossref/OpenAlex 解析验证；其中两篇标注为撤稿（见第一节结论 2）。Urologiia 2026 那篇系统评价在 PubMed 中未登记 DOI，仅给出 PMID。

## 九、附带产物

- 精选表：`results/tables/intratumoral_bacteria_bioinformatics_curated.csv`
- 全量候选表：`results/tables/intratumoral_microbiome_bioinformatics_candidates.csv`
- 原始检索记录（2765 条，含摘要）：`data/processed/tumor_microbiome_pubmed_records.csv` / `.json`
- 脚本：`analysis/scripts/pubmed_search_tumor_microbiome_bioinformatics.py`、`screen_tumor_microbiome_candidates.py`、`fix_doi_extraction.py`、`merge_supplementary_records.py`、`build_curated_table.py`、`fetch_curated_authors.py`、`render_review_report.py`

> 表中"撤稿"提示基于 Crossref `update-to` 元数据；引用任一论文前建议再核对一次期刊页面与 PubMed 记录状态。
"""

OUT.write_text(doc, encoding="utf-8")
print("wrote", OUT, len(doc), "chars; themes:", themes)
print("retracted DOIs flagged:", retracted)
