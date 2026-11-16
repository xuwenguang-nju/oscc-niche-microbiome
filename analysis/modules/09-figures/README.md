# 09-figures — 论文图件

## Purpose

生成稿件正文引用的图件（稿件索引见 `results/reports/manuscript-index.md`）。每张图一个脚本，一个函数对应一个面板，
便于单独重渲染某个面板时不必重跑全部。

## Inputs

| 文件 | 用途 |
|---|---|
| `results/tables/oscc_niche_sample_manifest.csv` | run → patient_id / niche（UTF-8 BOM，须 `utf-8-sig`） |
| `results/tables/s3_simulation_mdc.csv` | F2 需在此基础上作图 |
| `results/tables/s7_*.json` `s7_source_ranking_tests.csv` `s7b_before_after.csv` | F3 |
| `results/tables/s5_meta_analysis.csv` | F5 |
| `results/tables/s8*` `species_Parvimonas_*.csv` | F6 |
| `results/tables/s6v4_coupling_long.csv` | F7 |

## Commands

```bash
python analysis/modules/09-figures/scripts/fig1_study_design.py
```

在**本地**持久运行时执行（服务器 Python 无 matplotlib）。用 `run_in_context` 跑绘图脚本会触发
native crash 0xC06D007F 且返回 exit 0 —— 这一点在 S3 阶段已确认，本模块所有绘图都在
`python`（local）里跑。

## Outputs

| 图 | 脚本 | 面板 | 状态 |
|---|---|---|---|
| `figures/fig1-study-design.png` | `fig1_study_design.py` | a 队列×生态位覆盖 · b 分析流程 · c 证据地图 | 已完成，已肉眼复核 |

其余（F2–F7）目前直接复用 `analysis/modules/0x-*/scripts/plot_*.py` 产出的 PNG，
尚未迁入本模块统一出图。

## Figure 1 的设计取舍（值得记录）

- **panel a 用"覆盖矩阵"而非"队列地图"**：读者真正需要知道的不是地理，而是
  "六个队列各自采了哪些生态位、缺口在哪"。地理信息压进行标签副行。
- **panel b 的流程框标注了每步由哪个队列支撑**，因为本文的核心约束是"有些问题只有一个队列
  能回答"；把它放在方法图里，比放在 Discussion 里更早建立预期。
- **panel c 是自我限制的显式化**：每条结论画出"可检验队列数"与"确认队列数"，n=1 的条目
  直接标注 `single cohort — cannot be replicated`。这比在正文用文字交代更难被读者忽略。

## 已知未做

- F1 未画地理分布或患者配对结构（配对结构在 S7/S8 的正文表格里已给出，图上冗余）。
- 未做灰度/色盲友好版本；当前配色（蓝阶 + 红色强调）在 deuteranopia 下需复核。
- 未输出矢量格式；投稿前应改出 PDF/EPS（matplotlib `savefig(..., format="pdf")`）。

## Environment

本地 Python 3.10 + matplotlib（`D:\AI-Wisp-Science\OSCC_Crotonylation\.pixi\envs\default`）。
样式来自 `figure-style` skill 的 `apply_figure_style()`；`panel_crops()` 产出的检查裁切图属
临时垃圾，已删除，不应出现在 `figures/` 旁。
