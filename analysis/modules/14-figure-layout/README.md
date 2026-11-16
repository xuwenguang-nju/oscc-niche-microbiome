# 14-figure-layout — 图件定版与物理宽度控制

## Purpose

把项目里所有图件统一到期刊可用的物理宽度，并固定交付目录：

| 目录 | 内容 | 宽度 |
|---|---|---|
| `figures/main/` | **8 张主图** F1–F8 | 180 mm 双栏（S1 例外 120 mm） |
| `figures/supplementary/` | **7 张附图** S1–S7 | 180 mm；S1 为 120 mm 单栏 |
| `figures/legacy/` | 16 张重排前的原始输出 | 超宽，仅存档，**不要用于投稿** |

清单见 `results/tables/figure_manifest.csv`（含每张图的标题与图注草稿）。

---

## 为什么需要这个模块

项目里所有绘图脚本都写成 `figsize=(7.09, ...)`（=180 mm）然后
`savefig(..., bbox_inches="tight")`。**`tight` 会按内容重新计算画布**，导致交付的 PNG
普遍比标称宽度多出 15–90 mm：

| 实测 | 张数 |
|---|---|
| ≤182 mm（合格） | **3 / 16** |
| 超宽 190–269 mm | **13 / 16** |

最宽的是 `s8-replication-and-species`（**269 mm**）。

**关键认识：只改 `figsize` 没用**，因为 `tight` 会覆盖它。必须走
"渲染 → 量实际 tight bbox → 迭代修正" 的闭环。

---

## 机制（`figlib.render_at_width`）

不改任何绘图脚本的排版代码，而是：

1. 用 `runpy` 执行源脚本，同时 monkey-patch `plt.figure` 把 `figsize` 乘一个 `scale`
2. 拦下脚本的 `savefig`，拿到它本来要保存的 Figure
3. 强制换成 Agg canvas，量 `fig.get_tightbbox()`
4. 按目标宽度反推 `scale`，迭代至误差 ≤0.8 mm
5. 用**精确 bbox** 保存到指定路径

**字号不缩放**——pt 值保持作者设定，所以画布缩小后文字相对变大，这正是栏宽所需的效果。
也正因如此宽度对 `scale` 不是线性函数，必须迭代。

```python
from figlib import render_at_width
render_at_width("plot_forest.py", "figures/main/F5-cross-cohort-meta.png", mm=180)
```

多图脚本用 `which=N` 选第几张（`s7_plot.py` 存两张，`s11_plot.py` 存两张）。

### 收敛失败是特性，不是 bug

若"缩小画布反而让 tight bbox 变宽"，说明**有固定宽度的元素**（通常是超长单行标题）。
此时函数抛错而不是硬压，因为继续缩只会把字挤碎：

```
RuntimeError: width does not converge: shrinking the canvas makes the tight bbox
wider (231 -> 231 mm). This figure's density needs a layout change, not a scale change.
```

实际遇到的两次都是长 suptitle 钉住宽度，把标题换行即解决（F2、F5、S2、S5）。

---

## 布局改动（A 为主、B 为例外）

多数图只需缩放；**三个横排面板在 180 mm 下会被挤碎**，因此改了版式：

| 图 | 原布局 | 新布局 | 理由 |
|---|---|---|---|
| **F4** | 1×3（10.6 in 宽） | **3 列 × 4 行刚性网格**（AB / CD / EFG / HIJ，B、C 各跨 2 列） | panel a 含两个队列的块状点阵，压到 1/3 宽后块头标题撞进 b 的 y 轴；后续又发现"每块 panel 宽度不等、字母对不齐"，改用等宽列网格（见 `output/f4-grid-alignment.md`） |
| **S5** | 1×3 | **2×2** | panel b 是配对柱状图，需宽度 |
| **S1** | 1×3 | **3×1 纵向，120 mm** | 三面板回答彼此独立的问题，压扁会看不清；改单栏后字号反而更大 |

其余图保持原布局，仅缩放。

---

## 交付前的目视复核（每张都做过）

自动碰撞检测（`figure-style` 的 `non-tick text collisions`）不足以发现所有问题，
以下都是**看图**才发现的：

| 图 | 发现的问题 | 修法 |
|---|---|---|
| F1 | panel a 两行行标签逐行重叠 | 增高 panel a、字号 5.8→5.1、`linespacing=0.95` |
| F3 | 三个 panel 标题横向相撞；"median" 图例压数据；panel c 注释压 Saliva 标签 | 标题压成单行；删图例（红线自明）；注释移到轴下 |
| F4 | 三列布局崩坏（标题互叠、y 轴标签跨 panel） | 改 2×2 |
| F5 | panel c 标题超出右边界；两处图例挤在同一带 | 标题缩短；图例分层下移 |
| F8 | 三标题相撞；panel c x 轴标签被右边界截断 | 标题缩短；标签去掉 "mean relative"→"mean rel." |
| S5 | 三列挤压 | 改 2×2 |
| S1 | panel b 图例压 x 轴标签；panel c 注释贴轴 | 图例与注释再下移，纵向间距 0.62→0.85 |

---

## 一个顺带发现的数据口径矛盾

F5 的 panel a 原本用 `n_cohorts == 4` 筛选（10 属），而 S5 报告与 S11 新颖性分析用
`n_cohorts >= 3`（13 属）。**图会与全文其他文档不一致。**

已统一为 **k≥3**，并在图注写明 "3+ of 4 cohorts required"。
这是图标题自动显示 "10 genera" 时暴露的——**图与报告的交叉核对能抓出文字审阅漏掉的问题**。

---

## Scripts

```bash
python3 relayout_main.py     # 一次性重排 8 张主图 -> figures/main/
```

附图用同样的 `render_at_width` 逐张调用（见本模块 README 下方记录）。

## Outputs

| 路径 | 内容 |
|---|---|
| `figures/main/F1..F8-*.png` | 8 张主图，330 dpi |
| `figures/supplementary/S1..S7-*.png` | 7 张附图，330 dpi |
| `figures/legacy/*.png` | 16 张旧输出（超宽，存档） |
| `results/tables/figure_manifest.csv` | 图件清单 + 标题 + 图注草稿 |

## 已知限制

1. **仅有 PNG**（330 dpi）。投稿通常要 PDF/EPS，需再加一步矢量导出——
   做法：同样的 `render_at_width` 把后缀改 `.pdf`，`savefig` 自动走矢量后端。
2. **图注是草稿**，写清了"这张图证明什么"，但未按目标期刊格式（未编号引用、未含统计量符号规范）。
3. **重跑绘图脚本仍会往 `figures/` 根目录写临时件**——脚本默认输出到 `figures/`，而
   `figlib` 刻意不吞掉那次 `savefig`。根目录已于本轮清空（5 个文件归档到 `figures/legacy/`），
   但只要重跑某个脚本，根目录就会再出现同名文件。处理办法二选一：改各脚本的输出路径，
   或每次重跑后手动归档。**归档时不得覆盖既有存档**：`figures/legacy/` 里已有
   `fig1-study-design.png`（2026-10-23 的横版旧档），与本轮那份竖版中间件同名但内容不同，
   所以本轮那份存为 `fig1-study-design_2026-10-26.png`。今后遇到同名一律加日期后缀，
   不要覆盖。
4. S1 是唯一的单栏纵向图，若期刊要求全部同宽，需重新排成 180 mm 横排。
5. **`_run_capturing()` 会覆盖正式交付文件。** 它刻意让脚本自己的 `savefig` 正常执行
   （某些脚本保存后会 `out.stat()` 打印文件大小，吞掉写入会让它们报错）。对大多数脚
   本无害——它们写到 `figures/` 根目录——但 **S2 的脚本直接写 `figures/supplementary/`**，
   所以任何用 `_run_capturing` 做的几何核查都会把已定宽的文件替换成未缩放版
   （实测：176.0 mm → 187.7 mm）。若用 `_run_capturing` 验证某张图，**之后必须再跑一次
   `render_at_width` 恢复交付文件**，并核对该图尺寸。

6. **本环境新建 python 进程无法渲染 matplotlib（2026-10-27 发现）。** 任何 `savefig`
   （Agg/PDF 皆然）在**新建进程**中都会以 `Windows fatal exception: code 0xc06d007f`
   （= `ERROR_PROC_NOT_FOUND`，DLL 延迟加载找不到入口点）原生崩溃；连
   `plt.subplots()` + `plot()` + `savefig()` 这样的最小脚本也崩。常驻的 python REPL 不受影响。
   因此本模块的出图**只能走 REPL**（`wisp_python`），不要用 `run_in_context` 或 shell 调
   `python script.py` 来渲染。修复方式未知（疑为 pixi 环境二进制与 Python ABI 不同步）。
   另注：在 REPL 里**不要** `importlib.reload(pathlib)`——matplotlib 各模块在导入时已捕获
   `pathlib.Path` 类对象，reload 后旧类仍被引用，会报
   `type object 'Path' has no attribute '_flavour'`。
