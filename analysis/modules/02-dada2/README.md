# 02-dada2 — 逐队列 ASV 推断

> **2026-10-22 修订：宿主读段必须在属水平汇总前剔除。**
> `collapse_to_genus()` 原先保留全部 ASV，导致真核（宿主）扩增子进入分母。
> 实测各队列宿主读段占比差异极大：PRJNA1415359 0.0%、PRJNA597251 0.0%、
> PRJNA1415268 3.3%、**PRJNA666746 49.2%**。不剔除会让后者每个细菌属的相对丰度
> 被稀释约两倍，直接偏倚跨队列效应量。现规则：只保留 Kingdom ∈ {Bacteria, Archaea}，
> 细菌但无属级注释的读段保留为 `unassigned` 以便核对读段守恒。
> 四个队列的属表已用 `rederive_genus_tables.py` 以同一实现重算，
> 汇总见 `results/tables/genus_table_summary.csv`。

## PRJNA666746 被独立处理了两次（可重复性检验）

该队列先在本地 Windows 上跑完（DADA2 多线程不可用，1–2 核，**12.8 h**），
后在 Linux 计算节点上重跑（32 线程，**3.7 h**）。两次使用同一流程与同一下采样深度
（100,000 reads/样本），但各自抽取独立的随机子样本，因此不会完全相同。

| | ASVs | 非嵌合 reads | 平台 | 耗时 |
|---|---|---|---|---|
| 计算节点（权威，用于 S5） | 36,477 | 2,582,449 | Linux, 32 线程 | 3.7 h |
| 本地备份 | 36,428 | 2,572,088 | Windows, 1–2 核 | 12.8 h |

**属水平一致性（`check_run_concordance.py`）**：465 个共同属、100 个样本配对，
Bray-Curtis 相异度中位 **0.0144**（p95 0.026，最大 0.032），
Spearman ρ 中位 **0.993**（最低 0.960）。
产物：`results/tables/prjna666746_run_concordance.csv` 与 `..._summary.json`。

**目录约定**（重要，避免踩坑）：`results/dada2/PRJNA666746/` 是计算节点的权威结果；
本地重跑结果单独放在 `results/dada2/PRJNA666746_localrun/`。
两者同名文件曾互相覆盖，导致 `asv_table.csv` 与 `taxonomy_silva.csv` 来自不同运行而报错。
`rederive_genus_tables.py` 现已自动识别 `asv_table.csv` 的方向（ASV×样本 或 样本×ASV）。

## Purpose

把每个队列的原始读段独立处理成 ASV 表与（可选）分类注释。**六套队列绝不合并处理**——它们
在平台（MiSeq / HiSeq 2500 / Ion Torrent S5）、读长（251–460 bp）与深度（23k–968k reads/run）
上都不同，合并会把批次效应变成不可分离的混杂。跨队列比较在后续步骤用"分队列效应量 + 随机
效应 meta 分析"完成。

> **2026-10-27 修正：`cohort_params.yaml` 与实际运行不符，已订正。**
> Table 1 由该 yaml 生成，因此错报了两处，且两处都会误导复现：
>
> | 队列 | 曾记 | 实际（`run_status.json`） |
> |---|---|---|
> | PRJNA412445 | truncLen 236/236 | **245/240** |
> | PRJNA866676 | truncLen "无（变长）"、pool "pseudo" | **truncLen 350**、**pool FALSE** |
>
> 另外 `pool` 对**所有**队列都曾写作 `"pseudo"`/未同步：`run_dada2_paired.R` 与
> `run_dada2_single_end.R` 里 `pool = FALSE` 是**硬编码**的，yaml 的 `pool` 字段从未被读取。
> 现在 yaml 六队列一律 `pool: false`，并注明该字段仅为记录。
>
> `15-tables/build_main_tables.py` 也已修正：单端队列过去一律打印 "none (variable length)"，
> 即使实际用了固定 `truncLen`（866676 用 350）。Table 1 已重新生成。
> 教训：**以 `run_status.json` 为准**，yaml 是运行前的计划，运行后不会自动回写。

## Inputs

- `data/raw/16s/<study>/*.fastq.gz` — 上游模块 01-preprocess 的下载产物
- `analysis/modules/02-dada2/input/cohort_params.yaml` — 逐队列参数，全部来自 01-preprocess
  的实测（`trimLeft` 来自引物检测，`truncLen` 来自质量曲线，`single_end`/`ion_torrent`
  来自文件布局与长度分布）
- `data/reference/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz` — 139.0 MB，
  Zenodo 20955974，sha256 见 `data/reference/REFERENCE_SOURCES.md`
- `data/reference/ehomd/eHOMD_RefSeq_dada2_V15.22.fasta.gz` — 0.1 MB，Zenodo 7380292。
  口腔特异参考库；PRJNA412445 原文用的就是 HOMD，因此它是检验"SILVA 是否丢失口腔类群"的
  正确对照

## Methods

每个队列单独一次调用，流程固定为：

1. **filterAndTrim** — `maxN = 0`、`truncQ = 2`、`minLen = 150`、`maxLen = 500`、
   `rm.phix = FALSE`（Ion Torrent 无 PhiX）。`trimLeft` 与 `truncLen` 来自逐队列配置。
2. **learnErrors** — `nbases = 1e8`；双端时 R1/R2 分别学习。
3. **dada** — `pool = "pseudo"`。**Ion Torrent 队列额外传
   `HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32`**（DADA2 对 Ion Torrent 同聚物错误的
   必需设置；这两个参数经 `dada()` 的 `...` 传给底层比对器）。
4. **mergePairs** — 仅双端。
5. **removeBimeraDenovo** — `method = "consensus"`。
6. **assignTaxonomy**（独立脚本 `assign_taxonomy.R`，可对已有 ASV 表补做，无需重跑
   DADA2）— SILVA 138.2 为主，eHOMD V15.22 为对照。`taxLevels` 到 Genus。物种级不做：
   V3–V4 + Ion Torrent 上物种级不可靠，且跨队列不可比。
7. **collapse_to_genus** — 把 ASV 计数按 SILVA 属名聚合，得到 `genus_table_silva.csv`。
   **这是后续所有分析的分析单元。**

结果性参数（种子无关，但全部记录）：`trimLeft`、`truncLen`、`maxEE`、`pool`、
Ion Torrent 的两个比对参数。

## Software and data sources

- R 4.4.2（`C:\Program Files\R\R-4.4.2`，不改动系统库）
- DADA2 1.34.0（Bioconductor 3.20 Windows 二进制）安装于项目本地库 `analysis/r-lib`
- Biostrings 2.74.1、ShortRead 1.64.0、RcppParallel 5.1.11.2、IRanges 2.40.1、S4Vectors 0.44.0
- yaml 2.3.10、jsonlite 1.8.9
- 参考库：SILVA 138.2（可选，未下载时为 `unavailable`）
- 外部数据库：无

## Commands and scripts

```bash
# DADA2 推断（单队列；--max-files 限制 run 数用于飞行员）
Rscript analysis/modules/02-dada2/scripts/run_dada2.R --study PRJNA1415359 --no-multithread
Rscript analysis/modules/02-dada2/scripts/run_dada2.R --study PRJNA866676 --no-multithread

# 全部队列（多线程）
Rscript analysis/modules/02-dada2/scripts/run_dada2.R --all

# 分类注释 + 属水平聚合（可对已有 ASV 表补做，不重跑 DADA2）
Rscript analysis/modules/02-dada2/scripts/assign_taxonomy.R --study PRJNA1415359 --no-multithread
Rscript analysis/modules/02-dada2/scripts/assign_taxonomy.R --all

# Ion Torrent 的 maxEE 校准（默认 2；用 --out 写到临时目录避免覆盖正式结果）
Rscript analysis/modules/02-dada2/scripts/run_dada2.R --study PRJNA866676 \
    --max-files 2 --no-multithread --maxee 3 --out analysis/out/dada2_calib/me3
python analysis/modules/02-dada2/scripts/summarise_maxee_calibration.py
```

## Outputs

每队列一个目录 `results/dada2/<study>/`：

| 文件 | 内容 |
|---|---|
| `asv_table.rds` | 去嵌合后的序列表（DADA2 原生对象，列名为序列） |
| `asv_table.csv` | 同一份表的 CSV 版（列名为 ASV id，行名为样本） |
| `asv_seqs.fasta` | ASV 序列 |
| `taxonomy.rds` / `taxonomy.csv` | 仅当 SILVA 存在时 |
| `qc.csv` | 每样本：读入 → 过滤后 → 去噪后 → 非嵌合 → 占读入百分比 |
| `run_status.json` | 机器可读的成功/失败记录，含 trimLeft、truncLen、平台、ASV 数、耗时 |
| `taxonomy_silva.csv` / `.rds` | SILVA 138.2 注释到属 |
| `taxonomy_ehomd.csv` / `.rds` | eHOMD V15.22 注释到属（对照） |
| `taxonomy_compare.csv` | 每 ASV 的两个参考库属名是否一致，以及哪一方能注释 |
| `genus_table_silva.csv` | **属水平计数表（分析单元）**，行=样本，列=属 |
| `taxonomy_status.json` | 注释率（ASV 数 / 读段数两个口径）、属数、两库一致率 |

## 属水平注释的实测结果

对两个飞行员队列跑通了 SILVA 138.2 + eHOMD V15.22 双参考注释：

| 队列 | ASV | 属注释率（ASV） | 属注释率（读段） | 属数 | SILVA/eHOMD 一致率（ASV） | 一致率（读段） |
|---|---|---|---|---|---|---|
| PRJNA1415359 (MiSeq) | 318 | 85.8% | **100%** | 65 | 71.1% | **96.8%** |
| PRJNA866676 (Ion Torrent) | 935 | 67.8% | **100%** | 268 | 27.4% | 60.2% |

三个要点：

- **读段口径的注释率是 100%**，即每个读段都落到了某个属（含 `unassigned` 桶）。属水平
  分析不受注释缺失影响。
- **两库在 ASV 数上分歧大、在读段数上分歧小**（MiSeq：71.1% vs 96.8%）。说明分歧集中在
  低丰度 ASV，而占读段主体的高丰度类群两库一致。
- **eHOMD 能注释 SILVA 注释不了的 ASV**：MiSeq 队列中 301/318 (94.7%) 被 eHOMD 注释，
  SILVA 只有 273/318 (85.8%)；`taxonomy_compare.csv` 显示 31 个 ASV 是 "eHOMD only"、
  仅 1 个是 "SILVA only"。

生物学校验（MiSeq 队列 top 8 属）：*Streptococcus* 32,728、*Neisseria* 15,851、
*Rothia* 11,755、*Prevotella* 5,160、*Veillonella* 5,032、*Gemella* 4,638、
*Haemophilus* 4,043、*Granulicatella* 4,037 —— 全部是口腔常见属，符合预期。

## Ion Torrent 的 maxEE 校准（已完成）

`maxEE = 2` 是 Illumina 默认值，对 Ion Torrent 偏严。用同一批 2 个 run 比较三档
（`results/tables/dada2_maxee_calibration.csv`）：

| maxEE | 过滤保留率 | 非嵌合占读入 | ASV | ASV 相对 maxEE=2 的增幅 |
|---|---|---|---|---|
| 2 | 47.5% | 45.0% | 935 | — |
| 3 | 63.2% | 60.0% | 1002 | **+7.2%** |
| 4 | 72.4% | 68.7% | 1023 | +9.4% |

**结论：保留率接近翻倍，ASV 只增 9%。** 找回的读段是已检出序列的冗余拷贝，而不是新多样性。
**采纳 maxEE = 3**：拿到大部分保留率收益，同时把引入误差衍生 ASV 的风险保持在较低水平。
注意此校准只用了 2 个 run，全队列上应复核一次。

## 飞行员的 DADA2 结果（不是最终结果）

两条代码分支都跑通：

| 队列 | 分支 | trimLeft | truncLen | 过滤保留 | 合并率 | 非嵌合占读入 | ASV |
|---|---|---|---|---|---|---|---|
| PRJNA1415359 | MiSeq 双端 | 17 / 21 | 280 / 225 | 64.4% / 64.2% | 95% | **42.2% / 58.0%** | 318 |
| PRJNA1415359（错误参数对照） | MiSeq 双端 | 0 / 0 | 290 / 240 | 57.0% / 57.2% | 93% | 24.7% / 32.8% | — |
| PRJNA866676 | Ion Torrent 单端 | 0 | 0 | 26.9% / 47.5% | n/a | 24.4% / 45.0% | 935 |

两条值得记录的观察：

- **引物裁剪使保留率实质性提升**（非嵌合读段占读入 24.7%→42.2%）。这验证了 01-preprocess
  那一步的必要性。
- **嵌合率在两个平台上差异极大**：MiSeq 双端数据去掉 86% 的 ASV 序列（1923/2241），
  Ion Torrent 单端只去掉 7%（71/1006）。同一套嵌合算法、平台差异如此之大，说明嵌合率本身
  就是平台特征，**不能跨平台直接比较 ASV 数**。

## 一个已出现的生物学信号：Ion Torrent 队列带污染物特征

Ion Torrent 队列的属表里，`unassigned` 占 15.6%（9,045/57,927 读段），且第 4 位是
**Stenotrophomonas（4,641 读段）**——典型的试剂/环境污染属。这不是流程错误，而是低生物量
组织 16S 的常见现象，正好是 S2 等流行性去污染步骤要处理的对象。**该队列不能跳过 S2。**

## Limitations

- **ASV 数不可跨队列直接比较。** 平台、读长、深度、嵌合率都不同。跨队列分析必须在属水平
  并用分队列效应量。
- **Ion Torrent 属数（268）明显多于 MiSeq（65）**，但需注意两者样本数相同（各 2 个）而
  Ion 队列 ASV 更多、长尾更长。这不必然是噪声（>5 读段的属覆盖 99.8% 读段），但**必须经过
  S2 去污染和流行率过滤后才能纳入比较**。
- **maxEE 校准只用了 2 个 run**，全队列上应复核；`cohort_params.yaml` 里的 maxEE 尚未
  据此更新（当前靠 `--maxee` 命令行覆盖）。
- **PRJNA666746 的引物状态模糊**（R1 34%、R2 47%，低于 50% 阈值故未裁剪），两种设置应在
  全队列各跑一次再定案。
- **分类注释只在两个飞行员队列上跑过**，其余四套尚未运行（`--all` 需要先把全量读段下完）。
- **eHOMD 只做了注释对照，没有做属表聚合**；若最终采用 eHOMD 为主参考，需要扩展
  `collapse_to_genus` 的调用。
- **耗时**：`--no-multithread` 下 DADA2 约 250–315 s/样本（MiSeq）、约 480–620 s/样本
  （Ion Torrent）；注释约 130–170 s/队列。全量 898 runs 必须开多线程并按队列分块——本环境的
  `run_in_context` 对长时间原生进程不稳定，所有验收以"产物存在且非空"为准。
- 物种级鉴定未实现，且按设计不应做。
