# 08-replication — 跨队列复现 / 种水平解析 / 隐藏类群（S8）

## Purpose

回答三个跟进问题，并明确哪些**无法**用现有数据回答。

1. S7 的「癌旁组织是瘤内菌最近的源」能否在独立队列复现？
2. S7 的「倒灌梯度」（癌表面比深部肿瘤更接近唾液）能否复现？
3. S7 中唯一肿瘤富集的属 Parvimonas，其信号是否来自 *P. micra*？未达属水平的 ASV 里
   是否还藏着别的肿瘤富集类群？

## Answerability（先看这一节）

| 问题 | 可回答性 |
|---|---|
| ① 源排序复现 | **可**（PRJNA1415359 n=31；PRJNA412445 仅 n=8，功效不足） |
| ② 倒灌梯度复现 | **不可**：项目内只有 PRJNA866676 同时采了深部肿瘤与癌表面 |
| ③ 种水平 | **可，但受限于精确匹配的可解析率（12–20%）** |
| ④ 隐藏类群 | **仅科水平可检验**（属/种水平对未注释 ASV 本就不可能） |

②的替代检验（病灶侧黏膜 vs 对侧黏膜何者更接近唾液）在报告中标注为**代理指标，非复现**。

## Inputs

| 输入 | 说明 |
|---|---|
| `data/genus/<study>/genus_table.csv` | 样本 × 属，由 `s8_build_tables.py` 从 ASV 表重建 |
| `data/genus/<study>/family_table.csv` | 样本 × 科（全部 ASV） |
| `data/genus/<study>/family_table_genus_unassigned.csv` | 样本 × 科（**仅**有科无属的 ASV） |
| `data/genus/<study>/species_table_Parvimonas.csv` | 样本 × 种标签（含 `Parvimonas_unresolved`） |
| `results/tables/oscc_niche_sample_manifest.csv` | run → patient_id / niche（UTF-8 BOM，须 `utf-8-sig`） |
| 新队列原始读段 | PRJNA412445：102 双端 MiSeq run，518 MB，ENA |

## Methods

### 新队列处理（PRJNA412445）
- 参数全部实测：2×251 bp、**无 5′ 引物残留**、Q≥25 至 248(R1)/245(R2)
  → `truncLen=245,240`（保留 99.6% 读段、95%+ 碱基），`maxEE=2`，`pool=FALSE`
- 结果：102/102 样本、**15,826 ASV、1.198 M 读段、合并率 87.6%、耗时 921 s（32 线程）**
- 宿主读段：0%（该队列 ASV 全部为细菌）

### 源排序统计
- 每患者按「到该患者肿瘤的距离」给各源排序（秩 1 = 最近），对均值秩做
  Wilcoxon signed-rank 检验（零假设：无偏好，均值秩 =(k+1)/2）
- Bray-Curtis 与 Aitchison 双口径并行报告
- 属过滤：流行率 ≥10%

### 种水平
- 用 `assignTaxonomy(minBoot=50)` 对 toSpecies 训练集取 Species 列（**精确匹配**）
- 不强行归属：无法解析者保留 `Parvimonas_unresolved`，并在所有结论中报告可解析比例

### 隐藏类群
- 只取「有科、无属」的 ASV，折叠到科
- 患者内 CLR + Wilcoxon + BH
- **两道过滤**：流行率 ≥20% **且** 中位相对丰度 ≥0.5%（后者必需——低丰度科在 CLR 下
  伪计数会主导，曾因此产生 >1 的荒谬占比）

## Commands

```bash
# 下载（通用，处理单/双端与 gzip 尾部损坏）
TCMA_ROOT=$PWD python3 scripts/fetch_manifest.py --study PRJNA412445 --threads 10
# 参数实测 + 引物残留检查
TCMA_ROOT=$PWD python3 scripts/profile_cohort.py --study PRJNA412445 --runs 6
# DADA2
Rscript scripts/run_dada2_paired.R --study PRJNA412445 --trunclen 245,240 --threads 32
# 注释（属 + 科 + 隐藏科）
Rscript scripts/assign_taxonomy_full.R --study PRJNA412445 --res-root $PWD \
        --ref-root ~/oscc-mb/reference --threads 24
# 建表
python3 scripts/s8_build_tables.py --study PRJNA412445 --asv ... --tax ... --asv-seqs ...
# 种水平
Rscript scripts/species_only.R --study-dir ... --genus Parvimonas --tax-file taxonomy_full.csv \
        --ref-species ~/oscc-mb/reference/silva/silva_nr99_v138.2_toSpecies_trainset.fa.gz
# 分析
python3 scripts/s8a_replicate_ranking.py --genus-dir data/genus --out results/tables
python3 scripts/s8b_species_micra.py     --genus-dir data/genus --out results/tables
python3 scripts/s8c_hidden_taxa.py       --genus-dir data/genus --out results/tables
```

## Result-changing parameters

| 参数 | 默认 | 说明 |
|---|---|---|
| `--trunclen` | 245,240 | 提高会丢样本；降低会缩短重叠（<12 bp 无法合并） |
| `minBoot`（种） | 50 | 提高会降低可解析率；此处是精确匹配，影响很小 |
| 流行率过滤 | 0.10 / 0.20 | 影响可检验类群数量，不改变方向 |
| `--min_rel`（S8c） | 0.005 | **必需**，见上 |

## Outputs

`results/tables/s8a_*`、`s8b_*`、`s8c_*`；图 `figures/main/F4-replication-and-species.png`
（由 `s8_plot.py` 出图，经 `14-figure-layout/scripts/figlib.render_at_width` 定到 180 mm，
不在本步直接调用）；报告 `results/reports/step8-replication-species-hidden.md`。

## Known limitations

- PRJNA412445 对源排序复现功效不足（n=8），阴性不可解读为反证。
- 倒灌梯度无第二套数据，未复现。
- 种水平可解析率仅 12–20%，结论限于可解析子集，不可外推到整个 Parvimonas 属。
- 隐藏类群只在科水平可检验。
- 412445 的「对侧正常黏膜」不等同于 866676 的「癌旁组织」。

## Environment

服务器 R 4.4.0 + DADA2 1.34.0（`~/oscc-mb/Rlib`）、Python 3.10.12 + pandas/scipy；
SILVA 138.2 toGenus 与 toSpecies 训练集在 `~/oscc-mb/reference/silva/`。
