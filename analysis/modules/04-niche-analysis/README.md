# 04-niche-analysis — 去污染、残差菌属、生态位归属排序

## Purpose

按 S3 模拟的结论，实施两个**通过验证的**读数，并完成没有阴性对照时的去污染：

- **S2 去污染**：多生态位等流行性（equiprevalence）。六套公共队列**没有任何一套提交了
  阴性对照**（已逐一扫描确认为 0），因此 `decontam` 的对照组法不可用。替代方案：在同一研究
  的多个生态位之间，污染物会等比例出现，真正定植菌的检出率会随生态位变化。
- **S4a 残差菌属（H2）**：输出**一个菌属列表**，不是分数。S3 已证明任何"肿瘤特异比例"在此
  数据条件下不可辨识。
- **S4b 生态位归属排序（H1）**：输出**排序**，不是贡献百分比。同样由 S3 决定。

## Inputs

- `results/dada2/<study>/genus_table_silva.csv` — 02-dada2 的属水平计数表（行=样本，列=属）
- `results/tables/oscc_niche_sample_manifest.csv` — run → 生态位 / 患者编号的映射

## Methods

**样本对齐。** 属表的行名形如 `SRR37008826_1_F_filt.fastq.gz`，取其 run 号与清单匹配；
若同一 run 的两条 mate 都在表中则保留第一条。

**深度过滤。** 属注释读段总数 < 1000 的样本被剔除（`DEPTH_MIN`）。

**`unassigned` 处理。** 保留在测序深度里（它代表未注释读段，是真实数据的一部分），但
**不参与**任何生物学检验。

**流行率过滤。** 只在 ≥10% 样本中检出的属进入检验（`PREVALENCE_MIN`），避免把稀有长尾
当作检验对象——S3 已表明这类类群在过离散抽样下会整类消失。

**变换。** 计数 → CLR（centred log-ratio，伪计数 0.5）。16S 属表是成分数据，直接对原始比例
做相关会被测序深度驱动。距离同时报告 Aitchison（CLR 上的欧氏距离）。

### S2：等流行性去污染

对每个属，构造「生态位 × 是否检出」的列联表，做卡方检验，BH 校正。某属若**无法拒绝**
"检出率与生态位无关"（q ≥ 0.05），标记为 `equiprevalent`——候选污染物。
同时输出各生态位的检出率，便于人工核对。

注意这是**记录而非静默剔除**：结果表保留全部属及其判定，是否剔除由下游显式决定。

### S4a：双阈值残差菌属（H2）

对每个患者（同时具备 tumour / adjacent_normal_tissue / preop_saliva 三者）：

1. 计算两个配对 log-ratio：`log(tumour) − log(adjacent normal)` 与
   `log(tumour) − log(preop saliva)`，按属逐一计算。
2. 跨患者对每个属的 log-ratio 做**单样本 t 检验**（自由度 = 患者数 − 1）。
3. 两个配对各自 BH 校正。
4. **残差 = 在两个配对中都显著（q < 0.05）且方向相同**。上调记为 `tumour_enriched`，
   下调记为 `tumour_depleted`。

为什么必须跨患者聚合：S3 的实测显示患者内估计的 SD 为 0.06–0.11，聚合到队列层面降到
0.03–0.07（方差缩减因子 1.2–4.2）。逐患者下结论不可靠。

### S4b：生态位归属排序（H1）

对每个患者的肿瘤样本，计算到各生态位质心的 Aitchison 距离，两种参照：

- `vs_cohort_centroid`：其他患者该生态位的质心（**剔除本人**，避免自我参照泄漏）
- `vs_own_patient`：该患者自己的正常黏膜 / 术前唾液（临床上更有意义的对照）

跨患者用 **Friedman 检验**（配对、非参数）判断各生态位距离是否有差异，并按**中位距离排序**
输出"最近 → 最远"。**不输出任何贡献百分比。**

## Software and data sources

- Python 3.11.16
- numpy 2.4.6、pandas 2.3.3、scipy 1.17.1
- 上游：DADA2 1.34.0 + SILVA 138.2（见 02-dada2）
- 外部数据库：无

## Commands and scripts

```bash
python analysis/modules/04-niche-analysis/scripts/niche_attribution.py \
    --study PRJNA1415359
```

## Outputs

| 路径 | 内容 |
|---|---|
| `results/tables/s2_equiprevalence_<study>.csv` | 每属的生态位检出率、卡方 p、BH q、是否等流行 |
| `results/tables/s4a_residual_taxa_<study>.csv` | 每属的两个 log-ratio 均值、q 值、方向、是否残差 |
| `results/tables/s4b_distance_to_niche_<study>.csv` | 每患者肿瘤到各生态位质心的距离 |
| `results/tables/s4b_distance_own_patient_<study>.csv` | 每患者肿瘤到自身各生态位参照的距离 |
| `results/tables/sample_niche_map_<study>.csv` | 样本 → 生态位 / 患者的对齐结果（供核查） |
| `analysis/out/s4b_niche_ranking_<study>.json` | 排序结果与 Friedman 检验统计量 |

## Limitations

- **等流行性判定依赖生态位数量与样本量。** PRJNA1415359 的 4 个生态位、约 241 个样本足够；
  但单个生态位样本数少时卡方检验功效低，会把污染物误判为"随生态位变化"。
- **S2 只标记不剔除。** 本模块输出判定，剔除策略需在下游显式声明并做敏感性分析（去污染前 /
  后两版结果）。
- **CLR 的伪计数选择影响稀有属。** 伪计数 0.5 对 65–270 个属的表格是常规选择，但对
  极端稀疏的属仍不稳定；S3 的流行率过滤是主要防线。
- **S4a 的 t 检验假设患者间独立**，未校正年龄/性别/吸烟/HPV 等协变量——这些信息在公共队列
  的样本级元数据里大多缺失（已核实）。
- **S4b 的质心法假设生态位内部足够同质。** 若某生态位内部异质性很大（例如唾液跨个体差异），
  质心距离会低估真实离散度。
- **未实现**：跨队列随机效应 meta 分析（S5）、宿主侧耦合（S6）。本模块只做单队列 S2/S4。
