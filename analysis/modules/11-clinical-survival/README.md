# 11-clinical-survival — 临床关联、预后、新颖性、平台一致性（S11）

## Purpose

回答四个此前未做的问题，全部基于 TCGA-HNSC：

| 子分析 | 问题 | 结论 |
|---|---|---|
| A2 | 菌群指标与分期/分级/HPV/吸烟/部位有关吗 | 18 项检验 BH 后全阴（最低 q=0.218） |
| A3 | 菌群多样性/丰度能预测 OS 吗 | 单因素弱关联，调整分期后消失（LRT p>0.05） |
| A4 | 我们复现的属在文献里新吗 | 13 个中 5 个无 OSCC 报道（含 2 个改名属） |
| D | TCMA 与 16S 在属级别可比吗 | ρ=0.70–0.76，可比；偏差是类群特异的 |

**D 是本章最重要的一项**：它把 S6 的"平台不可比"改写为"平台可比、偏差有限且可定位"。

## Inputs

| 文件 | 来源 | 说明 |
|---|---|---|
| GDC `/cases` API | api.gdc.cancer.gov | 528 例 TCGA-HNSC：分期、分级、部位、吸烟、饮酒、生存 |
| cBioPortal `hnsc_tcga` | cbioportal.org/api（`clinical_data`） | HPV_STATUS_ISH / HPV_STATUS_P16（GDC 无此字段） |
| `s6v2_tcma_genus.csv` | 本项目 S6（TCMA 网页导出） | 177 aliquot × 185 属 |
| `oscc_niche_sample_manifest.csv` | 本项目 S1 | 生态位标签 |
| `data/genus/*/genus_table.csv` | 本项目 S2/D | 三套 16S 属表 |
| `ormckb_oscc_genus_votes.csv` | 项目既有 | ORMCKB OSCC 属级文献计数 |
| `s5_meta_analysis.csv` | 本项目 S5 | 13 个可复现属 |

## Methods（含决策依据）

### A1 主表
- 177 aliquot → **155 患者**；OS 可用 **156/177 样本** → 去重后 **134 患者、73 事件**
- 分期正则 `Stage\s*(IV|III|II|I)`，无尾部 `\b`（否则 "IVA" 中 IV 后匹配失败，曾漏掉 79 例）
- `sex_at_birth`（GDC 新字段）替代已弃用的 `gender`
- `primary_site` 是 **case 顶层字段**，不在 `diagnoses` 下

### A2 临床关联
Kruskal-Wallis（3+ 组）/ Mann-Whitney U（2 组），BH 校正跨 18 项检验。

### A3 生存
- 单因素 Cox（连续，每 1 SD）、KM 中位分割 + log-rank、多因素 Cox（分期序数+年龄+性别）、
  **嵌套模型 LRT**（临床模型 vs +菌群指标）、ΔAIC
- **分期必须序数编码**：`factor(stage)` 时 Stage I 的 7 例零死亡 → 完全分离，
  HR 4×10⁷、CI Inf。已用序数 + early/late 二分做敏感性
- **必须按 `submitter_id` 去重**，否则同一患者多 aliquot 造成伪重复

### A4 新颖性
S5 的 13 个方向一致属 → 查 ORMCKB 票数 → 判是否有既往 OSCC 报道、方向是否一致。
**改名属（Schaalia←Actinomyces、Lancefieldella←Atopobium）需人工判读**，"无报道"应写成
"no prior report under this name"。

### D 平台一致性
- 属级相对丰度，两边各自按样本归一后取均值
- Spearman ρ，**并做流行率过滤（两侧均 ≥20%）**——零膨胀稀有菌会扭曲秩相关
  （Genome Medicine 2023 的方法学教训）
- 另报 top-20 丰度重叠（Jaccard）与秩差最大的属，使偏差模式可见
- **非配对设计**：TCMA 与 16S 无共同样本，故混合了平台效应与队列效应

## Scripts

```bash
# 全部在 ssh:yinxu 的 ~/oscc-mb/work/clin/ 下运行
python3 probe_clinical.py          # 探查 GDC 字段填充率（一次性）
python3 s11_build_master.py        # A1 主表 → master.csv
python3 s11_clinical_assoc.py      # A2 → s11_clinical_assoc.csv
Rscript s11_survival.R             # A3 → s11_survival_*.csv（需 survival 包）
python3 s11_genus_concordance.py   # D  → s11_genus_concordance*.csv
# 本地
python3 s11_novelty.py             # A4（纯项目文件）
python3 s11_plot.py                # F6 + F9
```

## Outputs

`results/tables/`：`master.csv`、`s11_clinical_assoc.csv`、`s11_survival_{univariable,km,multivariable,lrt,stage_sensitivity}.csv`、
`s11_novelty_assessment.csv`、`s11_genus_concordance{,_disagreement}.csv`

`figures/`：`f6-clinical-prognosis-TCGA-HNSC.png`、`f9-platform-concordance.png`

## Result-changing parameters

| 参数 | 取值 | 影响 |
|---|---|---|
| 分期编码 | **序数 1–4** | 因子编码会完全分离，HR 不可解释 |
| 生存去重键 | **submitter_id** | 不去重则 n=177 伪重复，KM 分组数错 |
| 流行率阈值（D） | **≥20%** | 不过滤 ρ 偏低（0.58 vs 0.70），受零膨胀驱动 |
| KM 分割 | 中位数 | 不用最优 cutpoint，避免选择偏倚 |
| 多重校正 | BH，按分析族分别校正 | 18 项（A2）/ 3 项（A3）|

## Limits

1. HPV 仅 28 例 → HPV 结论不可下
2. TCMA 无癌旁 → 全部为跨患者横断面比较
3. D 的 ρ 非配对，平台与队列效应混杂
4. TCGA-HNSC 以咽/喉为主，真口腔癌比例低
5. C-index 0.50–0.56，菌群指标几乎无区分度
