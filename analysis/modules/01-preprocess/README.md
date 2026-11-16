# 01-preprocess — 原始读段获取与参数实测

## Purpose

把六套公共 OSCC 16S 队列的原始读段落到本地，并在**不做任何生物学假设的前提下**实测出
DADA2 所需的三个结果性参数：读段长度/质量曲线（决定 `truncLen`）、5′ 引物残留（决定
`trimLeft`）、以及平台类型（决定是否需要 Ion Torrent 专用分支）。

这个模块存在的理由：这三个参数是结果性的，套模板会静默降低数据质量。本项目就踩到了这一点
（见下文"实测发现"）。

## Inputs

- `results/tables/oscc_niche_sample_manifest.csv` — **898 个 run（档案全量）**的清单，含 ENA 的
  `fastq_ftp` 直连路径、生态位标签、患者编号、读段数、平台。由
  `analysis/scripts/build_niche_sample_manifest.py` 生成（S0 阶段）。

  > **档案 run ≠ 分析样本。** 进入 DADA2 的只有 **782 个（245 例患者）**：PRJNA1415359 只跑了
  > 125 个 run（32 肿瘤 + 32 癌旁 + 31 术前 + 30 术后唾液，来自 32 例患者），见
  > `analysis/modules/02-dada2/input/subset_PRJNA1415359.txt`。被排除的 116 个是 80 个只捐唾液、
  > 无组织的供者，加 18 例无肿瘤患者贡献的 36 个唾液 run；其余 5 个队列两者相同。
  > 论文与图件（Fig 1A、Fig 1F、Table 1）一律以 **782** 为准。

## Methods

**获取**（`fetch_reads.py`）。按清单从 ENA 的 HTTPS 镜像下载 FASTQ 到
`data/raw/16s/<study>/`。下载是幂等的：已存在且大小与远端一致的文件跳过；大小不一致的
删除重下。失败会写入 `data/processed/fetch_status.csv` 并以非零退出码返回，因此调用方可以
据此分支——本环境的 `run_in_context` 有时对失败仍返回 0，所以**验收标准是产物存在且非空，
不是退出码**。

**读段画像**（`profile_reads.py`）。每个队列抽 3 个文件、每文件前 20,000 条读段，统计长度
分布与逐循环平均 Phred。`truncLen` 取"平均质量首次跌破 Q20 的循环 −1"，再考虑配对重叠。

**引物残留**（`detect_primer_residue.py`）。用 IUPAC 感知的正则匹配 7 种常见 16S 引物
（27F / 338F / 341F / 515F / 806R / 926R / 1492R）及其反向互补，只检查读段**最前 26 bp**，
并区分"引物在 5′ 端"与"序列内部恰好含保守区"。超过 50% 的读段带引物就判定需要 trimLeft。

## 实测发现（三条，均改变了下游参数）

**1. `trimLeft` 必须逐队列实测，不能统一设为 0。**
两套队列的读段仍带引物：PRJNA1415268（R1 341F 81%、R2 806R 76%）与 PRJNA1415359
（R1 341F 85%、R2 806R 94%）。另四套已经切除（0%）。

**2. 引物残留会静默压低保留率。** 用 `trimLeft = 0` 试跑 PRJNA1415359 时，过滤只保留
57.0%/57.2% 的读段；改为实测的 `trimLeft = [17, 21]` 并相应把 `truncLen` 从 [290,240] 调为
[280,225] 后，过滤保留率升到 **64.4%/64.2%**，配对合并率 93%→95%，最终非嵌合读段占原始
读段比例从 24.7%/32.8% 升到 **42.2%/58.0%**。
关键点：trimLeft 会缩短可用读长，`truncLen` 必须同步下调，否则相当于从下游硬切。

**3. PRJNA866676 不是它对应论文的数据。** ENA 记录是**单端**（每 run 一个文件）、读长
**428 bp 可变（370–460）**，符合 Ion Torrent S5 XL；而它最常被引用的论文
（PMID 36046741）描述的是 **Illumina HiSeq 2500 双端、305 个样本、65 例患者**。
ENA 这套是 313 个 run、67 个患者编号。**结论：不能把该论文的方法学当作这套数据的依据**，
区段只能按实测扩增子长度（约 428 bp，与 V3–V4 一致）判断。

## 参考库获取（`fetch_references.py`）

两个参考库，用途不同：

- **SILVA 138.2 toGenus**（139.0 MB，Zenodo 20955974）——跨队列统一口径的默认参考库。
- **eHOMD V15.22**（0.1 MB，Zenodo 7380292）——Human Oral Microbiome Database，口腔特异。
  PRJNA412445 原文用的就是 HOMD，因此它是检验"SILVA 是否丢失口腔类群"的正确对照。

两个文件都记录了 URL、字节数与 sha256 到 `data/reference/REFERENCE_SOURCES.md`，可据此校验。
下载幂等。

**一件事没做成**：最初按 `silva_nr99_v138.2_toGenus_trainset.fa.gz` 去 Zenodo 4587955
（社区里最常被引用的记录号）取，返回 404——那个记录只到 138.1。实际可用的是 Zenodo 20955974，
它同时提供 138.2 的 toGenus 与 toSpecies 两个版本。

## Software and data sources

- Python 3.11.16
- 仅标准库（csv、gzip、json、re、urllib）——本模块不依赖第三方包
- 外部数据库：ENA portal API 与 `ftp.sra.ebi.ac.uk` HTTPS 镜像；访问日期 2026-10-20
- 下载规模（飞行员，6 runs/队列）：66 个文件，2.12 GB，全部成功
- 全量规模（实测 `fastq_bytes` 求和）：44.9 GB

## Commands and scripts

```bash
# 飞行员下载（每队列 6 个 run）
python analysis/modules/01-preprocess/scripts/fetch_reads.py --pilot --threads 6

# 全量下载（约 45 GB）
python analysis/modules/01-preprocess/scripts/fetch_reads.py --threads 6

# 只取某几套队列
python analysis/modules/01-preprocess/scripts/fetch_reads.py --studies PRJNA1415359 PRJNA666746

# 参数实测
python analysis/modules/01-preprocess/scripts/profile_reads.py
python analysis/modules/01-preprocess/scripts/detect_primer_residue.py

# DADA2 安装到项目本地 R 库（不改动系统 R）
Rscript analysis/modules/01-preprocess/scripts/install_dada2.R
```

## Outputs

| 路径 | 内容 |
|---|---|
| `data/raw/16s/<study>/*.fastq.gz` | 原始读段（只读，不修改） |
| `data/processed/fetch_status.csv` | 每个文件的下载状态、字节数、失败信息 |
| `results/tables/read_profile_by_cohort.csv` | 每队列 × 读向的长度分布与逐循环质量统计 |
| `results/tables/read_profile_summary.csv` | 每队列 × 读向的建议 `truncLen`、5′ 优势碱基、长度范围 |
| `results/tables/primer_residue.csv` | 每队列 × 读向的 5′ 引物命中率、引物名、建议 `trimLeft` |
| `results/tables/dada2_trimlen.json` | 机器可读的 trimLeft 配置（含支持证据） |
| `analysis/out/read_profile_report.txt` | 逐队列质量曲线与判定理由 |
| `analysis/out/primer_residue.log` | 引物检测输出 |
| `analysis/out/install_dada2.log` | DADA2 安装日志 |
| `data/reference/silva/*.fa.gz` | SILVA 138.2 toGenus（139.0 MB） |
| `data/reference/ehomd/*.fa.gz` | eHOMD V15.22（0.1 MB） |
| `data/reference/REFERENCE_SOURCES.md` | 两个参考库的 URL、字节数、sha256、许可与用途 |

## Limitations

- **DADA2 的 Windows 安装路径**：bioconda 没有 Windows 构建，pixi/conda 方案在此平台不可行
  （`pixi.toml` 保留作为该事实的记录，但它无法求解 `bioconductor-dada2`）。实际采用
  Bioconductor 官方 Windows 二进制（Bioc 3.20 → dada2 1.34.0），装进 `analysis/r-lib`。
- **PRJNA666746 的引物状态是模糊的**（R1 34%、R2 47%），低于 50% 阈值故未裁剪。可能是混合
  批次。两个设置需在全队列上各跑一次再定，不应凭阈值直接定案。
- **PRJNA866676 的 maxEE 需要放宽**：`maxEE = 2` 是 Illumina 的经验值，对 Ion Torrent 偏严。
  飞行员过滤只保留 26.9%/47.5%，其中相当一部分很可能是 maxEE 造成的而非真实低质量。需在
  全队列上比较 maxEE = 2/3/4 后再锁定，当前数值不应视为最终值。
- 质量画像只抽样（每队列 3 文件 × 2 万读段），足够决定参数但不足以刻画队列内部的样本间
  异质性。
- 引物检测只覆盖 7 种常见引物；若某队列用了非标准引物，会被误判为"已切除"。
