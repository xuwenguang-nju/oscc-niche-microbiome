# 13-function — 功能预测（S13，PICRUSt2）

## Purpose

在 PRJNA866676 的六生态位数据上做 KEGG/MetaCyc 功能预测，回答两个问题：

1. **这份预测可信吗？** PICRUSt2 的 NSTI 衡量每个 ASV 距参考基因组的距离。
   多数瘤内菌文章报 KEGG 结果却不报 NSTI，正好掩盖了预测在"外推"的情形。本模块把 NSTI 作为一手 QC 报告。
2. **功能谱能否复现 S7 的来源排序？** 若通路空间也显示"肿瘤最接近癌旁"，则来源归属结论
   不是 16S 分类学的产物，而是**同一批样本上的独立证据轴**。

## 环境

**PICRUSt2 2.6.3**，装在用户级 micromamba 前缀（无 root）：

```bash
# micromamba 静态二进制，无需 conda、无需 sudo
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
~/.local/bin/micromamba create -y -p ~/oscc-mb/envs/picrust2 \
    -c conda-forge -c bioconda -c defaults picrust2=2.6.3
```

3.8 GB，含 `epa-ng`、`gappa`、`hmmsearch`。**不需要 minimap2**。

> 曾误判"无 conda 就装不上 PICRUSt2"。实际上：(a) micromamba 用户级可装；
> (b) picrust2 在 bioconda 的 `noarch` 子目录，之前的检索只筛了 `linux-64` 所以显示为空。

## Inputs

| 文件 | 来源 |
|---|---|
| `PRJNA866676.prev2.fasta` | 由 `e1_build_inputs.py` 从 DADA2 ASV 表生成 |
| `PRJNA866676.prev2.biom` | 同上（BIOM v2.1，行=观测、列=样本） |
| `oscc_niche_sample_manifest.csv` | 项目内，生态位标签 |

### 输入过滤（必须声明）

| 规则 | 保留 ASV | 保留读段 |
|---|---|---|
| 全部 | 35,375 | 100% |
| **流行率 ≥ 2 样本（采用）** | **7,910（22.4%）** | **97.22%** |

丢弃 27,465 个单例 ASV。依据：单例在扩增子数据里主要是残余测序误差，且 EPA-NG 放置成本
随唯一序列数线性增长，而它们只携带 2.78% 读段。脚本内有 `<90% 读段则拒绝运行` 的守卫。

## Scripts

| 脚本 | 作用 | 运行位置 |
|---|---|---|
| `e1_build_inputs.py` | ASV 表 → fasta + BIOM（含过滤与守卫） | ssh:yinxu |
| *(picrust2_pipeline.py)* | 放置 + 预测 EC/KO + 推通路 | ssh:yinxu |
| `e2_pathways.py` | NSTI QC、PERMANOVA、初见源排序、配对差异 | ssh:yinxu |
| `e3_distance.py` | 距离层面的来源比较（排除癌表面） | ssh:yinxu |
| `e4_annotate.py` | 加 MetaCyc 名称 + 效应量汇总 | ssh:yinxu |
| `e5_plot.py` | F11 三面板 | local |

```bash
# 完整流程
python3 e1_build_inputs.py --study PRJNA866676 --min-prev 2
micromamba run -p ~/oscc-mb/envs/picrust2 picrust2_pipeline.py \
    -i PRJNA866676.prev2.biom -s PRJNA866676.prev2.fasta \
    -o out -p 32 --in_traits EC,KO --verbose
python3 e2_pathways.py && python3 e3_distance.py && python3 e4_annotate.py
```

## Outputs

`results/tables/`：`e2_qc_nsti.json`、`e2_permanova.json`、`e2_pathway_niche_ranking.csv`、
`e2_pathway_paired_tests.csv`、`e3_tumour_distances.csv`、`e3_distance_paired_tests.csv`、
`e4_pathway_tests_annotated.csv`、`e4_comparison_summary.csv`、`e4_top_pathways.md`、
`e_pathway_relabund.csv`

`figures/`：`FigS1-function-niches.png`（三面板：NSTI 质控、通路 PCoA、配对来源距离）

另见 `analysis/scripts/figS3_functional_shift.py` → `figures/supplementary/FigS3-functional-shift.png`：原发队列 PICRUSt2 通路火山图，以及校正 consortium score 后全部归零的配对图。数据源 `results/tables/fw2024_420/picrust2_kegg_pathway_tests.csv`，全表 `tableS4_primary_cohort_kegg_pathways.csv`。

## 关键结果

| 项 | 值 |
|---|---|
| 放置成功 | 7,608 / 7,910（**丢失 3.8%**） |
| 样本加权 NSTI 中位 | **0.020**（p95 0.041，最大 0.078） |
| ASV NSTI > 0.5 | 8 个（0.11%），**携带 0.003% 读段** |
| 通路 PERMANOVA | F=14.4，R²=0.190，p=0.001 |
| 癌旁为最近外部来源 | **20/37**（随机 1/4，binom **p=1.7e-4**） |
| 肿瘤到癌旁 < 到唾液 | 22/37（59.5%），Wilcoxon **p=0.023** |
| 倒灌梯度（癌表面更近唾液） | 24/36（66.7%），**p=0.062（未显著）** |
| 差异通路 \|中位差\| ≥0.5% | **0 / 544**（三个比较皆 0） |

## 结果改变参数

| 参数 | 取值 | 影响 |
|---|---|---|
| 流行率过滤 | **≥2 样本** | 不过滤则 35,375 ASV，EPA-NG 耗时剧增且单例放置不可靠 |
| NSTI 剔除阈值 | **PICRUSt2 默认 2.0** | 本数据无 ASV 超过，等于未剔除 |
| 参考谱（来源排序） | **留一法**（用**其他**患者的该生态位均值） | 不剔除自身会让患者与自己匹配，距离虚低 |
| 是否含 cancer_surface | **排除**（同病变） | 含入则"癌表面最近"会掩盖真正的外部来源 |
| 效应量阈值 | **0.5% 功能潜能** | 不设阈值则"277/544 显著"会被误读为实质差异 |
| 校正 | BH，按每个比较独立 | — |

## Limits

1. **基因组潜能 ≠ 活性**：无宏转录组/宏蛋白组。
2. 相对丰度，非绝对定量。
3. 仅 PRJNA866676；PRJNA412445 尚未做功能预测。
4. 通路为 unstratified 汇总，丢失物种贡献信息。
5. PICRUSt2 参考库偏肠道/体表，口腔厌氧菌覆盖有限；NSTI 好只说明被放置的 ASV 有近亲。
6. 癌表面 vs 深部配对仅 36 对，倒灌检验功效低。

## 已知踩坑（都会静默出错或直接崩）

1. **BIOM v2.1 的 rows/columns 必须带 `"metadata": null`**，否则 `Table.from_json` 抛 `KeyError: 'metadata'`
2. biom 约定**行=观测、列=样本**；反向会得到 0 行匹配
3. DADA2 的 CSV 行索引是**序列**不是 ASV id，必须按位置映射而不是按 id join
4. 必须用 `micromamba run`（或把环境 bin 加入 PATH），否则 `place_seqs.py` 找不到
5. `weighted_nsti.tsv.gz` 在 `KO_metagenome_out/` 子目录，不在输出根目录
6. 这个版本**必须同时提供 `-i`（BIOM）和 `-s`（fasta）**，只给 fasta 直接打印 usage 退出
