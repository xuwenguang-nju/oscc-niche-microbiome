# 07-six-niche — 六生态位来源归属（S7，PRJNA866676）

## Purpose

项目里唯一的**同患者多生态位**队列：肿瘤 / 癌旁 / 癌表面 / 对侧正常黏膜 / 唾液 / 舌苔，
共 313 样本、67 例患者，其中 **26 例六生态位齐全**。用来回答三件事：

1. 瘤内菌群能否与口腔各生态位区分？
2. 它最像哪个生态位（可能的来源）？
3. 是否存在**倒灌梯度**——暴露的癌表面是否比深部肿瘤更像唾液？

## 数据与流程

```
PRJNA866676  313 runs, Ion Torrent 单端, 8.2 GB
  ↓  fetch_prjna866676.py        （URL 前缀补全 + gzip 校验 + 容错恢复）
  ↓  profile_866676.py           实测: 中位 426 bp, Q>=25 至 349 位
  ↓  run_dada2_single_end.R      truncLen=350, maxEE=3, BAND_SIZE=32, pool=FALSE
  ↓  assign_taxonomy.R           SILVA 138.2 -> 1,550 属, 读段层面 100% 有属级
  ↓  s7_niche_analysis.py        Q1/Q2/Q3
  ↓  s7b_decontam_sensitivity.py 去污染敏感性
  ↓  s7c_s5_genera.py            S5 十属在本队列的复现
  ↓  s7_plot.py                  两张图
```

**产出**：313 样本 × 35,375 ASV，17,293,002 条非嵌合读段，中位保留 **94.84%**；
宿主污染仅 **0.16%** 真核读段。

## 分析口径（重要）

- **只报距离与秩，不报混合比例**。S3 的模拟已证明混合比例在这些数据里不可辨识，
  任何"瘤内菌 X% 来自唾液"的表述都是错的。
- **所有关键比较都做患者内配对**（六生态位来自同一人），PERMANOVA 用
  patient-blocked 置换，否则零假设是错的、p 值会虚低。
- 只用出现率 ≥20% 的属（129 个），避免稀有噪声属造成"生态位可分"的假象。

## 去污染敏感性（必读）

各生态位试剂/环境背景菌占比极不均衡：**癌旁 47.6%、肿瘤 13.3%、癌表面 8.8%、
对侧 8.0%、唾液 0.5%、舌苔 0.8%**（其中癌旁的 Stenotrophomonas 就占 42.6%）。
该研究**无阴性/空白对照**，无法用 decontam 频率法。因此剔除 13 个已知试剂污染属
（129→116 属）后重算，结果：

| 结论 | 前 | 后 | 翻转？ |
|---|---|---|---|
| 癌旁最近 | −0.109, p=0.0035 | **−0.158, p<0.0001** | 否，更强 |
| 癌表面较近 | −0.041, p=0.080 | **−0.049, p=0.033** | 否，转显著 |
| 对侧较远 | +0.065, p=0.0022 | **+0.104, p<0.0001** | 否，更强 |
| 舌苔最远 | +0.080, p<0.0001 | **+0.096, p<0.0001** | 否，更强 |
| 倒灌梯度 | p=2.8e-05 | **p=4.9e-06** | 否，更强 |

**全部加强，无一翻转** → 结论不是污染伪影。

## 结果摘要

- PERMANOVA（患者内配对）：pseudo-F=16.9，p=0.0005，R²=0.275
- 源排序（1=最近）：癌旁 **1.78** > 癌表面 2.25 > 唾液 2.86 > 舌苔 3.78 > 对侧 3.93
- 倒灌梯度：到唾液的中位距离，深部肿瘤 **0.706** vs 癌表面 **0.463**；
  **81% 患者癌表面更近，p=4.9e-06**
- **S5 十属在本队列 10/10 方向一致**（肿瘤 vs 癌旁）
- **Parvimonas** 生态位梯度：肿瘤 **2.00%** > 癌旁 0.78% > 癌表面 0.68% >
  唾液 0.36% > 舌苔 0.19% > 对侧正常 **0.12%**（比任何口腔生态位高 6–17 倍）
- 其余九属在肿瘤中低于所有口腔生态位 → 多为"肿瘤内口腔菌整体变少"，
  而非主动清除

## Outputs

```
results/tables/s7_niche_composition.csv         6 生态位 × 129 属
results/tables/s7_permanova.csv                 配对/非配对 PERMANOVA
results/tables/s7_pairwise_niche_distances.csv  15 对生态位距离
results/tables/s7_source_ranking.csv           逐患者源排序
results/tables/s7_source_ranking_tests.csv     逐源检验
results/tables/s7_backflow_gradient.csv        倒灌梯度
results/tables/s7_discriminative_genera.csv    生态位标志属
results/tables/s7b_contaminant_burden.csv      污染占比
results/tables/s7b_before_after.csv            去污染前后
results/tables/s7c_s5_genera_by_niche.csv      S5 十属 × 六生态位
results/tables/s7c_s5_direction_check.csv      方向一致性
results/tables/s7_summary.json / s7b_summary.json / s7c_summary.json
figures/s7-six-niche-source-attribution.png
figures/s7-s5-genera-niche-profile.png
results/reports/step7-six-niche.md
```

## 不做的分析（明确排除）

- **不做定量来源分解**（如 SourceTracker/FEAST 的百分比）：S3 已证明不可辨识。
  若将来有绝对定量或更强的先验，再考虑。
- 不做跨平台合并（本队列是 Ion Torrent 单端，与四个 MiSeq 队列不合并建模，
  只在方向层面比较）。
- 不做时间序列/纵向推断（横断面设计）。

## 限制

1. 单一队列（中国人群、Ion Torrent），无外部验证。
2. 无阴性/空白对照，污染只能用"剔除已知污染属"处理。
3. 横断面 → 倒灌梯度与"表面接触"一致，但不能定因果。
4. 癌旁最接近肿瘤也可能只是解剖邻近的必然结果，不等于"来源"。
5. 26 例六生态位齐全；部分患者缺位，配对检验的 n 在 27–37 之间。
