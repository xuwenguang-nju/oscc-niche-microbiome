# 10-literature — 二手解读 → 一手文献的溯源与提取

## Purpose

用户提供的 17 条微信公众号"文献解读"链接是**二手材料**，不足以支撑"该加哪些分析"的决策。
本模块把这些解读回溯到 PubMed/PMC 原始记录，取得开放全文，按一手方法学重新评估。

## 为什么必须做这一步

复核中实际发生的事：

| 问题 | 后果 |
|---|---|
| 第 1 篇首轮误匹配成**胶质瘤**论文 | 若不校验，整篇"HNSC 可借鉴模块"的参考会错 |
| "13/13 篇有湿实验" | 关键词 `animal` 命中讨论里的 "animal models"；实际仅 5 篇（02/03/12/13/16） |
| TCMA 被我们称为 "RNA-seq-mined" | 与一篇已发表论文（PMID 35862975 称 WXS/WGS）**相反**；已改为中性表述并加待办 |

## Inputs

- 用户提供的微信链接（已人工去重至 17 篇）
- PubMed / PMC E-utilities（`eutils.ncbi.nlm.nih.gov`）
- Europe PMC REST（DOI 精确匹配用）

## Commands

```bash
# 1) 标题片段 → PMID/期刊/年/DOI/PMC
python3 resolve_papers.py

# 2) 修正误匹配 + 按标题片段逐条验证 + 取 PMC JATS XML
python3 fetch_fulltext.py

# 3) XML → Markdown（摘要+章节+图注+方法词表命中）
python3 to_text.py

# 4) 从方法正文提取队列规模/平台/软件/分析类别矩阵
python3 extract_analyses.py
```

全部在 `ssh:yinxu` 上跑（默认分析环境），仅需 Python 标准库 + 无第三方依赖
（`to_text.py` 用 `xml.etree`）。

## Outputs

| 路径 | 内容 |
|---|---|
| `literature/stimulus-papers/md/*.md` | 15 篇可读全文（13 篇有正文） |
| `literature/stimulus-papers/texts/*.xml` | PMC JATS 原始 XML（追溯用） |
| `literature/stimulus-papers/resolved_papers.json` | PMID/DOI/PMC/标题 |
| `literature/stimulus-papers/paper_analysis_matrix.json` | 队列/平台/软件/分析类别 |
| `results/reports/literature-based-expansion-plan.md` | 重估后的扩展方案 |

## 结果

15 篇解析成功；13 篇有正文；2 篇（06 HCC-MRS、11 breast-bile-acid）PMC 仅收录摘要。
期刊含 Hepatology、Mol Cancer、Gut Microbes、Genome Med、Microbiol Spectr 等。

## 结果改变参数

| 参数 | 默认 | 影响 |
|---|---|---|
| 检索式 | 标题片段 + `[Title]` | 放宽会引入误匹配（已发生一次） |
| 标题验证片段 | 每篇指定 | **必须保留**，这是防误匹配的唯一闸门 |
| `CATEGORIES` 关键词表 | 见脚本 | 决定矩阵覆盖面；过宽会像 `animal` 那样误判 |
| `min_align`/统计阈值 | 不适用（仅阅读） | — |

## Known limitations

- **只收 PMC 开放获取**；订阅期刊正文未取得，未尝试绕过付费墙。
- 关键词矩阵记录的是"**文中提及**"。不出现 ≠ 未做，故只能用于"这篇做了什么"，
  **不能**用于"这个领域没做什么"。
- 15 篇中仅精读 01、16 两篇（与本研究最直接相关）；其余依赖程序化提取，可能有遗漏。
- 17 条链接中 2 条触发微信验证墙，未绕过；1 篇在 PubMed 无法按标题验证。

## 待办

- [ ] 确认 TCMA 底层测序 assay（需 PMID 33382980 方法部分，非 OA）
- [ ] 若后续要引用这 15 篇中的具体统计量，须逐篇核对原文数值而非依赖本项目矩阵
