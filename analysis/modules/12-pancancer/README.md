# 12-pancancer — 泛癌上下文（S12）

## Purpose

把 HNSC 放进其他 TCGA 肿瘤里，检验"OSCC 的瘤内群落是口腔来源的"这一论断在项目层面是否成立。

**核心逻辑**：若瘤内菌来自口腔，则上呼吸消化道肿瘤（HNSC/ESCA/STAD）应以口腔菌为主，
下消化道肿瘤（COAD/READ）应以肠道菌为主，且方向相反。TCMA 用**同一 pipeline**处理这五个
项目（同一参考库、同一去污染），所以项目间差异不是处理差异造成的——这是用 TCMA 而非
另找 16S 队列的关键理由。

## Inputs

| 文件 | 来源 | 说明 |
|---|---|---|
| `data/raw/tcma/pancancer.genus.relabund.txt` | TCMA 网页导出（v2 2022） | 185 属 × 620 样本 |
| `data/raw/tcma/pancancer.sample_metadata.txt` | 同上 | 620 × 83 字段（project / AnatomicSite / sample_type） |
| `data/raw/tcma/pancancer.phylum.relabund.txt` | 同上 | 门水平（备用，未分析） |

**样本对齐**：612→616 可用（丢弃 4 个总丰度为 0 的样本：ESCA 2、READ 1、STAD 1）。
矩阵列与元数据行**完全一致**（620/620 交集），已核验。

**范围限制**：网页版只含 5 个 TCGA 项目。TCMA 完整库（1772 例、12 癌种）在 Duke Data
Repository，本环境无公开 API，**未获取**。故本节实为"上呼吸消化道 vs 胃肠道"。

## Methods

| 步骤 | 方法 | 参数 / 依据 |
|---|---|---|
| 组成 | 各项目主要属平均相对丰度 | top-12 属 |
| 排序 | PCoA on Bray-Curtis（属水平相对丰度） | 自实现（Gower 中心化 + 特征分解） |
| 组间差异 | PERMANOVA，999 次置换 | Anderson 2001；报告 F、p、R² |
| 口腔/肠道轴 | 先验指标集求和后相减 | **指标集不从未数据学习**，避免循环 |
| 组间检验 | Kruskal-Wallis | 轴值、*Parvimonas* 丰度 |
| 灵敏度 | 只取 Primary Tumor 重算 | 排除正常组织稀释 |

### 口腔 / 肠道指示属（先验声明）

```
口腔 (11): Streptococcus, Veillonella, Fusobacterium, Prevotella, Neisseria,
           Rothia, Actinomyces, Haemophilus, Porphyromonas, Capnocytophaga, Parvimonas
肠道 (10): Bacteroides, Faecalibacterium, Roseburia, Blautia, Eubacterium,
           Ruminococcus, Alistipes, Akkermansia, Barnesiella, Odoribacter
```

所选 21 个属在两个数据集中**全部命中**（11/11、10/10），无须替换。

## Scripts

```bash
python3 s12_pancancer.py   # 全部统计 → s12_*.csv/json
python3 s12_plot.py        # F10 三面板
```

## Outputs

| 文件 | 内容 |
|---|---|
| `s12_project_composition.csv` | top-12 属 × 5 项目 |
| `s12_pcoa_coordinates.csv` | 620 样本 PCoA1-3 + 项目/部位/样本类型 |
| `s12_oral_gut_axis.csv` | 项目级口腔/肠道轴（均值、比值、检出属数中位） |
| `s12_oral_gut_per_sample.csv` | 逐样本轴值（作图用） |
| `s12_focus_genera_by_project.csv` | 8 个关注属的均值与检出率 |
| `s12_sample_type_sensitivity.csv` | 肿瘤 vs 正常分层重算 |
| `s12_pancancer_summary.json` | 全部检验统计量 |
| `figures/f10-pancancer-context.png` | a PCoA · b 口腔−肠道轴 · c *Parvimonas* |

## 结果改变参数

| 参数 | 取值 | 影响 |
|---|---|---|
| 零丰度样本 | 丢弃 4 个 | 保留会引入 NaN 列 |
| PCoA 距离 | Bray-Curtis | 换 Aitchison 会改变轴解释度但结论方向应保持 |
| 指标属清单 | 21 个先验属 | 换清单会改变轴数值；结论（符号翻转）对合理的口腔/肠道属清单应稳健 |
| 置换次数 | 999 | 影响 p 的分辨率，不影响 F |

## 关键结果

- PERMANOVA（项目）：**F=37.9，p=0.001，R²=0.199**；PCoA1+2 = 37.6%
- 口腔−肠道轴：HNSC **+74.9**、ESCA +78.8、STAD +59.6、COAD **−47.9**、READ **−44.7**；
  Kruskal-Wallis **H=370.5，p=6.4e-79**
- *Parvimonas*：HNSC **1.55%（检出 80%）** vs 其余 0.31–0.42%（36–56%）；
  **H=93.6，p=2.2e-19**
- 只取原发瘤重算：结论完全不变
- **HNSC 癌旁组织口腔菌占比 74.2% ≈ 肿瘤 75.3%，*Parvimonas* 1.89% > 1.51%**
  → 与 S7 的"瘤内菌最接近癌旁"**独立互证**

## Limits

1. 仅 5 个 TCGA 项目（缺 LUAD/LUSC/BRCA/PAAD 等真正的泛癌谱）
2. 指标集是先验的；"口腔 vs 肠道"并非二元（STAD 居中即为证据）
3. 属级相对丰度，无法区分 *P. micra* 等种
4. PCoA 前两轴仅 37.6%
5. 正常样本的 HNSC n=22、READ n=4，分层的正常组结论功效有限
