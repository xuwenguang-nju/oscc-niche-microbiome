# 06-host-coupling — 宿主侧耦合（S6）

> **2026-10-22 重要更正。** 本模块最初用 GSE227919（总 RNA-seq），原因是我把 TCMA 的
> 域名探错成 `tcma.galaxyproject.org`，据此判定"TCMA 不可得"。用户指出真实地址是
> **`thecancermicrobiomeatlas.org`**，本机与服务器均可访问。S6 已改用 TCMA/TCGA-HNSC 重做，
> GSE227919 降级为旁证。
>
> **当前生效版本：TCMA v2 / TCGA-HNSC**（177 aliquot：155 原发瘤 + 22 癌旁；185 属）
> - 方向检验：6/9 与 S5 一致（sign test p=0.51）
> - **配对检验（22 对）3 个属通过 FDR**：Actinomyces(q=0.011)、Veillonella(q=0.012)、
>   Lancefieldella(q=0.012，经旧名 Atopobium 近似匹配）
> - **Parvimonas 方向相反**（−0.33，p=0.97），S5 的肿瘤富集未在此平台复现
> - 详见 `results/reports/step6b-tcma-host-coupling.md`
>
> **S6-v3 按解剖部位分层（GDC primary_site）后的结论修正：**
> 分层**没有**救回 Parvimonas——在最干净的口腔组（78 瘤 / 17 癌旁 / 17 对）里
> log2FC = **−0.03（几乎正好为零，p=0.64）**。且它在口腔样本中丰度 1.55%、
> **检出率 90.5%**，所以既不是"部位混杂"也不是"灵敏度不足"。
> 相反，分层让另外三个属变强：Actinomyces q 0.011→**0.0003**、
> Lancefieldella 0.012→**0.0044**、Veillonella 0.012→0.0157，
> 因为喉部样本方向相反、会稀释口腔信号。
> 详见 `results/reports/step6c-site-stratified.md`。
> 限制：咽部 0 例癌旁、喉仅 2 例，二者基本不可检验。
>
> **S6-v4 宿主程序耦合（TCMA 相关面板，4 测定 × 50 特征 × 10 属）：**
> **Parvimonas 200 个检验单元 0 个通过 FDR（min q=0.083，max |R|=0.32）**，
> 与任何免疫/EMT/增殖/生长程序都不耦合；而 Prevotella(q=0.002)、Campylobacter(0.005)、
> Streptococcus(0.011)、Actinomyces(0.011)、Veillonella(0.018) 均有显著耦合。
> Campylobacter 的 20 个命中全在甲基化层面（全部负相关），提示微生物-宿主关联可能经
> 细胞组成（免疫浸润）间接产生，**是相关不是因果**。
> 详见 `results/reports/step6d-host-coupling-panels.md`。
>
> **S6-v5 舌部细分（61 瘤 / 16 癌旁 / 16 对）：**
> Parvimonas **方向翻正（+0.40）**，是 S6 中首次与 S5 同向，但 p=0.38 不显著；
> 按 d≈0.22 需约 170 对才有 80% 功效，故这是功效不足而非"无效应"的证据。
> 三个可复现属在舌部依旧成立（Actinomyces q=0.0005、Lancefieldella 0.0069、
> Veillonella 0.0228），强度与口腔组相当 → v3 的"口腔组变强"主因是**排除喉部**。
> 口底 0 例癌旁，无法再细分。详见 `results/reports/step6e-tongue-subsite.md`。
>
> 下面 GSE227919 部分保留为方法学记录（其微生物信号 47% 为环境/试剂背景，不足以做检验）。

## Purpose

回答"瘤内菌是否与宿主免疫/EMT 程序相关"。S5 的四个 16S 队列没有宿主转录组，
所以本节必须建立在**微生物与宿主来自同一份样本**的数据上。

## 数据可得性（已实测，见 `results/tables/s6_data_availability.json`）

| 方案 | 结果 |
|---|---|
| TCMA（TCGA RNA-seq 挖微生物） | 不可得：`tcma.galaxyproject.org` 域名无法解析，figshare 403，Zenodo 无记录 |
| 16S 队列 + 外部宿主组学 | 无同源样本，不做 |
| **GSE227919 总 RNA-seq** | ✅ 采用。一份文库同时含宿主转录组与原核读段 |

GSE227919：Control 18 / HkNR 17 / Dysplasia 22 / Cancer 9，共 66 例。

## Inputs

```
data/raw/gse227919/GSE227919_2023_20_03_expr_mat_microbe.txt.gz   948 taxa x 66
data/raw/gse227919/GSE227919_2023_20_03_expr_mat_raw.txt.gz       21510 genes x 66
data/raw/gse227919/GSE227919_2023_21_03_pml_annotations.csv.gz    66 x 31（含 Class）
results/tables/s5_meta_analysis.csv                                S5 的四队列结果
```

两个矩阵都是**空格分隔**（不是制表符），脚本用 `sniff_sep()` 自动判定。
微生物矩阵行 id 形如 `ti|1969841|org|Proteus_phage_...|accession|NC_...`，
用 `org\|([^|]+)\|` 抽取物种名再取属。

## Methods

1. **属水平汇总 + 背景剔除**：948 行 → 924 行为可解析的属级；
   `EXCLUDE` 剔除宿主/真菌/病毒行。**注意**：该正则最初含裸 `homo`，
   会误删 *Stenotrophomonas*（含 "homon"），连带丢 40% 读段；现为精确匹配。
   保留 916 taxa / 297 属，丢弃 32 行经核对为宿主/真菌/病毒。
2. **样本对齐**：微生物 id 用 `-`、宿主/注释用 `_`，`norm_sample()` 统一，66/66 全对齐。
3. **S6a 方向检验**：只比方向，不比效应量——16S 扩增子与 RNA-seq 挖出的微生物读段
   在深度、poly-A/宿主背景与分类偏倚上不可比。符号检验作整体判据。
4. **S6b 宿主耦合**：6 个**预先指定**的宿主程序（细胞毒/IFN-γ/检查点/p-EMT/增殖/上皮分化），
   各取一组标志基因做样本内 z-score 均值；与 S5 可复现属做 Spearman，Benjamini-Hochberg 校正。
   预指定而非事后挑选，是为了避免钓鱼。
5. **S6c 信号构成诊断**：把 297 属按口腔典型 / 环境试剂背景 / 其他分层，
   判断这条路线的信噪比是否足以做任何检验。

## Run

```bash
python3 discover_s6_sources.py      # 数据可得性探查
python3 probe_gse227919.py          # 系列元数据
python3 fetch_gse227919.py          # 下载 + gzip 完整性校验
python3 inspect_gse227919.py        # 结构探查
python3 s6_host_coupling.py         # S6a + S6b
python3 s6_signal_composition.py    # S6c 诊断
python3 plot_s6.py                  # 出图
```

环境：Python 3.10.12，pandas 2.3.3，scipy 1.15.3（`pip install --user`）。
DADA2 无关，本节不涉及。

## Outputs

```
results/tables/s6a_direction_check.csv
results/tables/s6a_group_abundance.csv
results/tables/s6b_genus_signature_corr.csv
results/tables/s6b_signal_composition.csv / .json
results/tables/s6_host_signature_scores.csv
results/tables/s6_summary.json
results/tables/s6_data_availability.json
results/tables/s6_gse227919_{probe,files,structure}.json
figures/s6-host-coupling-GSE227919.png
results/reports/step6-host-coupling.md
```

**附图 Fig. S2（`figures/supplementary/FigS2-host-coupling-site-platform.png`）**由
`s_combined_host_site_platform.py` 出图，经 `14-figure-layout` 的 `render_at_width` 定到 180 mm。
它原有 12 块 panel，已**去重为 7 块 A–G**：删掉的 5 块（旧 C/D/E/I/J）各自是与主图同一张表、
同一个统计量的重画——组成条（Fig. 6L 与 Fig. 7F，重复三次）、Parvimonas 按部位（Fig. 6J）、
按部位可检验性（Fig. 6K）、按队列的秩一致性（Fig. 7A），以及 Fig. 7C 前 12 名的嵌套子集（旧 J）。
保留的 7 块都只在本图出现。详见 `results/reports/figure-legends.md` 交叉核对第 4 条。

## 结论与限制

- **S6a**：5/9 方向一致，符号检验 p=1.00，与随机无异。*Parvimonas*（S5 中唯一四队列一致
  肿瘤升高的属）在此也呈升高方向（p=0.072）。
- **S6b**：54 个检验 **0 个**通过 FDR（最小 q=0.094）。
- **S6c（关键）**：RNA-seq 微生物信号中 **47% 为环境/试剂背景**，口腔典型菌仅 37%
  （榜首含 *Paeniglutamicibacter*、*Catellicoccus*、*Porphyrobacter*、*Alcanivorax* 等非口腔菌）。
- 因此 S6 的阴性**不构成对 16S 结论的反驳**，应表述为"在本项目条件下**不可判定**"
  （not evaluable），而非"已检验且不支持"。
- n=9 例 Cancer，任何组间检验功效都很低；本节全部数字只能作为方向性探针。

## 复现注意

- 脚本用 `project_root()` 向上搜索项目根（`results/` + `analysis/` 或 `results/` + `data/`），
  因为同一脚本要在项目内与计算节点 `~/oscc-mb/work/s6/` 两种布局下运行。
  早期用固定 `parents[N]` 时，项目里残留的 `analysis/results/` 空目录会被误判为根，已清除该残留。
- 服务器上传的脚本落在 run 的 `inputs/` 目录，需要 `cp` 到工作目录；
  用 `ls -t ... | head -1` 取最新版本，不要用 `| tail -1`（按 UUID 排序会取到旧版）。
