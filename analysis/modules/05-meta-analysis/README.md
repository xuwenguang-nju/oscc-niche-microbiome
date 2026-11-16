# 05-meta-analysis — 跨队列复现（S5）

> **2026-10-22 更新：四队列结果（宿主已剔除）。**
> 队列：PRJNA1415359（首尔，MiSeq）、PRJNA1415268（首尔，MiSeq）、
> PRJNA597251（青岛，HiSeq）、PRJNA666746（印度，HiSeq；下采样 100k 后端到端在
> Linux 计算节点完成，DADA2 1.34.0 与本地一致）。
>
> **主结论未变且更强**：|d_z| ≥ 0.8 的菌属中，**没有任何一个在 ≥2 个队列复现**
> （各队列大效应数为 0、0、1、2，彼此不重叠，最大 |d_z| 分别 0.655/0.790/0.877/0.825）。
> 在四个队列均可检验的 80 个属中，10 个方向一致且合并 q<0.05。
> **Parvimonas 是唯一在全部四队列均为肿瘤富集的属**（+0.561/+0.559/+0.310/+0.109，
> 合并 +0.756，q=1.4e-04，I²=0），但效应量都停留在中等水平。
> 与之相对，Staphylococcus、Alloprevotella、Veillonella、Schaalia、Actinomyces 等
> 在四队列一致地在肿瘤中**降低**（合并 -0.65 ~ -1.06）。
>
> 注意：`Collinsella` 仅在 2 个队列可检验，其 -0.99 的合并效应不应与四队列结论等同看待。
> 生产命令：`cross_cohort_meta.py`（默认已含四队列）。

## Purpose

回答一个问题：**PRJNA1415359 里"效应量小、无残差菌属"这个阴性结果，是真实的，还是那个队列
特殊？**

单队列阴性无法区分两种情形：
- 确实没有肿瘤特异信号
- 这个队列特殊（人群、采样、平台）

因此用三套独立配对队列做复现检验（PRJNA1415268 首尔 47 对 / PRJNA597251 青岛 24 对 /
PRJNA666746 印度 50 对）。这不是"再找几个显著菌属"，而是**验证阴性本身可复现**——两者的
结论含义完全不同。

## Inputs

- `results/dada2/<study>/genus_table_silva.csv` — 各队列的属水平表（上游 02-dada2）
- `results/tables/oscc_niche_sample_manifest.csv` — run → 生态位 / 患者编号

## Methods

**逐队列效应量。** 每个属，对每例患者计算配对 log-ratio：

```
log((tumour + 0.5)/sum) − log((adjacent normal + 0.5)/sum)
```

跨患者取均值即为效应估计，标准误 = SD/√n。这与 04-niche-analysis 的 S4a 是同一估计量，
但**只保留 tumour vs adjacent-normal 这一个对比**——因为四套队列只有它共有（只有
PRJNA1415359 有唾液）。

同时计算 Cohen's `d_z` = mean/SD（配对设计下的标准化效应量），这是判断"效应大小"而非
"是否显著"的关键。

**跨队列合并。** DerSimonian–Laird 随机效应模型，报 Cochran's Q 与 I²。

为什么用随机效应而非固定效应：四套队列在国家、平台（MiSeq / HiSeq 2500）、读长
（251–301 bp）、深度（2.3 万–97 万 reads/run）上都不同。固定效应模型会断言"存在一个共同
真效应"，而这个设计不支持该断言。

**复现判据（预先设定）。**

| 判据 | 含义 |
|---|---|
| `max|d_z|` per cohort | 该队列中最大效应有多大 |
| `|d_z| ≥ 0.8 的属数` per cohort | 是否存在大效应菌属 |
| `large_effect_replicated` | 在 **≥2 个队列**中都达 |d_z| ≥ 0.8 的属（真正值得追的候选） |
| `direction_consistent` | 所有队列效应同向 |
| `median I²` | 跨队列异质性 |

## Software and data sources

- Python 3.11.16、numpy 2.4.6、pandas 2.3.3、scipy 1.17.1
- 上游：DADA2 1.34.0 + SILVA 138.2（02-dada2）
- 统计方法：DerSimonian & Laird (1986) 随机效应；Cohen's d_z 配对标准化

## Commands and scripts

```bash
# 单队列冒烟测试（会优雅提示"少于两个队列"）
python analysis/modules/05-meta-analysis/scripts/cross_cohort_meta.py \
    --studies PRJNA1415359

# 完整复现
python analysis/modules/05-meta-analysis/scripts/cross_cohort_meta.py \
    --studies PRJNA1415359 PRJNA1415268 PRJNA597251 PRJNA666746
```

前置：各队列必须已跑完 02-dada2 的 `run_dada2.R` 与 `assign_taxonomy.R`。

## Outputs

| 路径 | 内容 |
|---|---|
| `results/tables/s5_effect_<study>.csv` | 单队列每属的效应、SE、Cohen's d_z、q、流行率 |
| `results/tables/s5_meta_analysis.csv` | 跨队列合并：pooled effect、I²、Q、方向一致性、大效应复现数 |
| `results/tables/s5_effect_by_cohort.csv` | 长表（genus × study），供森林图使用 |
| `analysis/out/s5_replication_summary.json` | 复现判据的汇总结论 |

## Limitations

- **只有 tumour vs adjacent-normal 一个对比能跨队列合并。** 唾液对比只有 PRJNA1415359 有；
  PRJNA866676 有癌表面与舌苔但没有配对正常组织。因此"倒灌梯度"检验无法在此模块完成。
- **平台差异会引入队列间异质性**，I² 会偏高。这是真实的、不可消除的差异，不是模型缺陷；
  随机效应模型已把 τ² 计入。
- **"adjacent normal" 的定义跨研究不一致**（同侧黏膜 / 对侧黏膜 / 距肿瘤距离未记录）。这是
  公共数据的固有限制，已核实在样本级元数据中确实缺失。
- **未纳入协变量**（年龄、性别、吸烟、槟榔、HPV/p16）——这些信息在公共队列的样本级元数据
  中大多缺失，无法校正。
- **复现判据依赖 |d_z| ≥ 0.8 这一阈值。** 0.8 是 Cohen 的"大效应"惯例值，但临床意义阈值
  应结合具体菌属的效应与测量噪声判断；本模块同时输出连续分布以便自行判断。
- **未做批次校正**。设计上选择"分队列估计 + 随机效应合并"而非合并原始数据，因此不需要批次
  校正，但也意味着不能做跨队列的样本级联合建模。
