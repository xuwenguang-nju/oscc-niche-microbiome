# 03-simulation — 检测下限与可辨识性模拟

## Purpose

在做任何真实数据分析**之前**，用已知真值回答两个问题：在六套公共 OSCC 队列
实际的测序深度和过离散水平下，

1. 能否估计出"瘤内菌群中不被已知源生态位解释的质量**比例**"（定量拆分）？
2. 能否识别出构成这部分未被解释质量的**菌属集合**（残差菌属）？

结论（本模块的主要交付物）：

- **读出一（定量比例）不可辨识，必须从研究设计中移除。** 三个目标函数（L1 / L2 / chi²）都无法还原真值，且失败是结构性的而非调参问题。L2 在构造上对"无法解释的质量"失明（`p@S` 在新增菌属处恒为 0，残差对 `p` 无导数，`sum(p)` 不受向下压力）；chi² 走向另一极端（稀有类群主导，估计饱和到约 0.9）；L1 对未被解释的质量敏感，但它自身的零分布（把一个**源样本**拿去对着同一套源参考打分）就已达到 0.25–0.71 量级，与信号同阶，偏差无法与信号分离。用真实口径（源参考跨患者合并、被打分患者留一剔除）后，`moderate_c80` 场景：肿瘤样本 0.297，源样本零分布 0.382——**零分布高于信号**。
- **读出二（残差菌属集合）成立。** p_novel = 0.10 时，全部五个深度下 sensitivity = 1.00、specificity = 1.00；加入更贴近真实的渗漏情形（肿瘤富集菌在邻近黏膜中也可检出 25%）后仍为 1.00 / 1.00（968k 队列为 0.923 / 1.00）。
- **真正的限制是"可检出性"，不是统计功效。** 能通过流行率过滤的肿瘤特异菌属数量随深度单调上升（p_novel = 0.10 时：23k → 9 个、43k → 10 个、105k → 11 个、130k → 11 个、968k → 14 个）。**深度只影响能检出多少菌，不影响检出后的判断准确性。**

因此研究设计改为：**残差菌属集合 + 定性归属排序**（肿瘤菌群最像哪个生态位），不再报告任何混合百分比。

## Inputs

本模块不读取外部数据；全部输入由脚本按固定种子生成。运行参数记录在
`analysis/modules/03-simulation/input/s3_simulation_settings.json`（含种子、
属数、过离散参数、深度表、浓度表、p_novel 网格、求解器说明、包版本）。

深度取值来自 ENA 实测的每 run 读段数中位数（见
`analysis/out/niche_technical_details.txt`）：

| 队列 | median reads/run |
|---|---|
| PRJNA412445 | 23,428 |
| PRJNA597251 | 42,602 |
| PRJNA866676 | 104,655 |
| PRJNA1415359 | 129,505 |
| PRJNA666746 | 967,658 |

## Methods

**真值构造。** 150 个属，其中 135 个为"源可解释"（G_SRC）、15 个为"源中完全不存在"
（G_NOVEL），因此残差菌属的判定真值无歧义。G_NOVEL 从群落丰度**中段**（秩 10–70）
抽取，而非极端稀有尾部——真实 OSCC 富集菌（*Fusobacterium*、*Capnocytophaga*、
*Treponema*、*Prevotella*）位于 0.1–5% 量级；真正稀有的类群在过离散抽样下会整类
消失，任何方法都无法评分，故通过流行率过滤（肿瘤样本中读段 ≥10 且检出率 ≥50%）
从灵敏度指标中剔除，并显式说明。

**生态位相似度。** 每个源生态位为患者特异的 Dirichlet(c·base) 抽样，c ∈ {20, 80, 400}。
报告的是**实际达成的**唾液 vs 相邻黏膜 Bray–Curtis 距离（0.63 / 0.42 / 0.23），而非
假定的相似度。

**汇（肿瘤）构造。** 唾液 0.40 + 相邻黏膜 0.40 + 龈下菌斑 0.20 的加权和，占 (1 − p_novel)，
其余 p_novel 来自 G_NOVEL 组成；p_novel ∈ {0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30}。

**计数抽样。** Dirichlet-multinomial，θ = 50（属水平 16S 的过离散），按上表深度抽样。

**队列结构。** 一个队列 = 20 例配对患者。源参考谱为**跨患者合并**并**剔除被打分患者本身**
（留一共参考），以对应真实分析中必须避免的自我参照泄漏。

**求解器。** 迭代重加权最小二乘（IRLS）+ 半单纯形投影，纯 numpy 实现。不使用 scipy
的 SLSQP：它在某些沙箱环境中会原生崩溃（0xC06D007F），且该环境中这类崩溃不会反映为
非零退出码。在相同输入上，本求解器与 SLSQP 的精确解**一致到四位小数**（完全可解释的
汇与含 20% 未知质量的汇两种情形均已核对）。

**目标函数比较。** 三种权重均实现并输出（`s3_objective_comparison.csv`）；
L1 为本模块采用的权重。Bray–Curtis 目标被明确**否决**：它要求两个组成各自归一到 1，
而本模型的已解释混合在构造上小于 1，分母 (sum y + sum m) 会把目标变成奖励放大 sum(p)
的比值，从而悄悄改变被估计的量。

**残差菌属判定。** 对每例患者计算两个配对 log-ratio（肿瘤 vs 相邻黏膜、肿瘤 vs 唾液），
再对跨患者的 log-ratio 做单样本 t 检验，BH 校正。一个菌属必须**在两个配对中同时显著
且同向**才判定为残差。

**固定参数。** α = 0.05；θ = 50；流行率过滤 10 reads / 50% 检出；每格 20 次重复；
队列规模 20。全部记录在 settings JSON 中。

## Software and data sources

- Python interpreter: 3.11.16
- Python package `numpy`: 2.4.6
- Python package `pandas`: 2.3.3
- Python package `scipy`: 1.17.1
- Python package `matplotlib`: 3.11.2
- 外部数据库：无。深度参数取自 ENA portal API 于 2026-10-20 的实测查询
- Wisp/model runtime metadata: unavailable

## Commands and scripts

```bash
# 由于沙箱中长时运行会原生中止，扫描按场景分块执行（每块约 3.5 分钟）
python analysis/modules/03-simulation/scripts/simulate_source_attribution.py \
    --reps 20 --patients 20 --scenario distinct_c20
python analysis/modules/03-simulation/scripts/simulate_source_attribution.py \
    --reps 20 --patients 20 --scenario moderate_c80
python analysis/modules/03-simulation/scripts/simulate_source_attribution.py \
    --reps 20 --patients 20 --scenario near_identical_c400
# 更贴近真实的渗漏情形（肿瘤富集菌在邻近黏膜可检出 25%）
python analysis/modules/03-simulation/scripts/simulate_source_attribution.py \
    --reps 20 --patients 20 --scenario moderate_c80 --leak 0.25

# 汇总各分块为最终表
python analysis/modules/03-simulation/scripts/simulate_source_attribution.py --assemble-only

# 图
python analysis/modules/03-simulation/scripts/plot_simulation_calibration.py
```

烟雾测试：`--quick`（4 次重复 × 8 例患者，约 50 秒）。

## Outputs

| 路径 | 内容 |
|---|---|
| `results/tables/s3_objective_comparison.csv` | L1 / L2 / chi² 的偏差对比（支撑"读出一"被否决） |
| `results/tables/s3_replicates_<scenario>.csv` | 每个场景的分块重复结果（`--assemble-only` 的输入） |
| `results/tables/s3_simulation_replicates.csv` | 2800 行合并重复（4 场景 × 5 深度 × 7 p_novel × 20 重复） |
| `results/tables/s3_simulation_calibration.csv` | 每格汇总：估计值、零分布、超额、灵敏度、特异度、可检出菌属数 |
| `results/tables/s3_simulation_mdc.csv` | 每场景 × 深度的最小可检出比例（要求 ≥50% 的目标菌属可检出） |
| `results/tables/s3_perpatient_variability.csv` | 患者内 SD vs 队列 SD，量化"聚合的必要性" |
| `figures/s3-simulation-calibration.png` | 两联图：(a) 比例不可辨识；(b) 真实限制是可检出性，且随深度变化 |

## Limitations

- **真实数据的源–汇相似度未知。** 本模块扫描了 0.23–0.63 的实际 Bray–Curtis 范围，但真实
  OSCC 的唾液/黏膜/肿瘤两两相似度需由 S1–S2 处理后的数据实测后再回填本图。
- **θ = 50 是假定值。** 真实过离散需从数据估计；θ 更小会使稀有类群更不稳定，可能进一步
  缩小可检出集合。
- **漏检渗漏。** `--leak` 只模拟了均匀渗漏（每个 G_NOVEL 属在三个源里各获得 p_novel/3 的
  质量）。真实情形的渗漏可能不均等，且可能只出现在相邻黏膜而非唾液。
- **G_NOVEL 的丰度分布是构造的。** 中段秩 + 对数正态抖动，并非真实 OSCC 富集菌的实测
  分布；这会影响"可检出菌属数"的绝对值，但不影响"灵敏度=1.00 一旦可检出"这一结论。
- **只模拟了属水平。** 物种水平解析在本设计中被排除（V3-V4 + Ion Torrent 不可靠），
  故未评估。
- **模拟不能替代实测的批次效应。** 六套队列的平台/区段/深度差异远大于本模块覆盖的范围，
  跨队列比较的可行性须由步骤 3 的真实数据检验。
