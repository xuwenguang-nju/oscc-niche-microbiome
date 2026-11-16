#!/usr/bin/env python3
"""F5 — cross-cohort meta-analysis: what replicates, and what does not.

Eight panels. The conclusion is deliberately two-sided: a set of genera replicates its
DIRECTION across cohorts with moderate effects, but no genus reaches a large effect in
more than one cohort, so single-cohort "large effect" claims do not survive.

A  volcano over all 146 testable genera; the 13 that replicate are highlighted and the
   number of cohorts testing each genus is the marker size
B  forest: per-cohort Cohen's dz and the random-effects pooled effect, 13 genera
C  why a 1-cohort study is not evidence: max |dz| per cohort vs the large-effect threshold
D  how novel these genera are against prior OSCC reports
E  effect matrix: every genus x cohort cell
F  heterogeneity (I^2) of the pooled estimates
G  the |dz| distribution by number of cohorts, showing large values are single-cohort
H  single-cohort |dz| vs pooled effect: the inflation that pooling removes

Two panels of the previous version were dropped as redundant: a bar chart counting
cohorts up vs down (already readable from B and E) and a three-point curve of how the
replicating set shrinks with the k-cohort threshold (already stated in B's title).

Reads:  results/tables/s5_meta_analysis.csv, s5_effect_by_cohort.csv,
        s11_novelty_assessment.csv
Writes: figures/main/F5-cross-cohort-meta.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "14-figure-layout" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))   # volcano_panel
from figlib import panel_letter_fixed, project_root  # noqa: E402

ROOT = project_root(Path(__file__).resolve().parent)
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures"
# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
exec(compile(open(_FIGSTYLE / "runtime.py", encoding="utf-8").read(),
              "runtime.py", "exec"))

apply_figure_style(sizes=(7.0, 6.2, 5.4))

meta = pd.read_csv(TAB / "s5_meta_analysis.csv")
eff = pd.read_csv(TAB / "s5_effect_by_cohort_full.csv")   # carries per-cohort se and n
DZ = pd.read_csv(TAB / "s5_pooled_dz_and_logratio.csv")   # pooled d_z, same scale as panel H
nov = pd.read_csv(TAB / "s11_novelty_assessment.csv")

K = 3
REP = meta[(meta.n_cohorts >= K) & (meta.pooled_q_bh < 0.05) &
           meta.direction_consistent].sort_values("pooled_effect")
STUDIES = ["PRJNA1415359", "PRJNA1415268", "PRJNA597251", "PRJNA666746"]
SLAB = {"PRJNA1415359": "1415359", "PRJNA1415268": "1415268",
        "PRJNA597251": "597251", "PRJNA666746": "666746"}
SCOL = {"PRJNA1415359": "#1f4e79", "PRJNA1415268": "#3182bd",
        "PRJNA597251": "#a6611a", "PRJNA666746": "#7b3294"}
BLUE, RED, ORANGE = "#2166ac", "#b2182b", "#a6611a"
LARGE = 0.8                     # |d_z| threshold for a "large" effect (panels D, G, H)
# Panel C plots a log2 ratio, so the d_z threshold is converted with the median per-cohort SD of
# the log2 ratio (1.15 ln-equivalent = 1.66 log2), giving 0.8 x 1.66 = 1.33. The shaded band is a
# reference for "large" only; no genus in this panel is claimed to fall inside or outside it.
LARGE_LOG = 1.33
L2 = float(np.log(2))

fig = plt.figure(figsize=(7.09, 7.05))
gs = fig.add_gridspec(3, 3, left=0.134, right=0.985, top=0.968, bottom=0.058,
                      wspace=0.78, hspace=0.66,
                      height_ratios=[1.42, 1.0, 1.0])
# The volcano needs the width of two cells; every other panel keeps one.
ax_volcano = fig.add_subplot(gs[0, 0:2])
axes = np.array([
    [ax_volcano, None, fig.add_subplot(gs[0, 2])],
    [fig.add_subplot(gs[1, j]) for j in range(3)],
    [fig.add_subplot(gs[2, j]) for j in range(3)],
], dtype=object)


def lbl(ax, letter, text):
    ax.set_title(text, fontsize=6.2, multialignment="center", pad=5)
    panel_letter_fixed(ax, letter, case="upper")


# ------------------------------------------------------------------ A volcano
# Drawn first because it is the selection step: 146 testable genera, the 13 that
# replicate highlighted, cohort count encoded as marker size.
from volcano_panel import draw_volcano  # noqa: E402

draw_volcano(axes[0, 0], letter="A", panel_letter=panel_letter, tab=TAB)

# ------------------------------------------------------------------ B forest
ax = axes[1, 0]
gen = REP.genus.tolist()
ys = np.arange(len(gen))
# Panel C uses ONE scale throughout: per-cohort mean log-ratio and the pooled log-ratio. An
# earlier version plotted per-cohort Cohen's d_z against a pooled log-ratio, mixing a
# standardised and an unstandardised effect on the same axis.
for i, g in enumerate(gen):
    d = eff[(eff.genus == g) & (eff.study.isin(STUDIES))]
    # per-cohort and pooled effects are natural-log ratios in the source table; the panel is in
    # log2 to match Table S1, Fig. 4L and the volcano in panel A
    for _, r in d.iterrows():
        e_, s_ = r["effect"] / L2, r["se"] / L2
        ax.plot([e_ - 1.96 * s_, e_ + 1.96 * s_], [i, i],
                color=SCOL[r["study"]], lw=0.8, alpha=0.9, zorder=2)
        ax.scatter(e_, i, s=11, color=SCOL[r["study"]], edgecolor="0.25",
                   linewidth=0.25, zorder=3)
    row = REP[REP.genus == g].iloc[0]
    pe_, ps_ = row["pooled_effect"] / L2, row["pooled_se"] / L2
    ax.scatter(pe_, i, marker="D", s=26, color="0.12", zorder=5)
    ax.plot([pe_ - 1.96 * ps_, pe_ + 1.96 * ps_],
            [i, i], color="0.12", lw=1.1, zorder=4)
ax.axvspan(-LARGE_LOG, LARGE_LOG, color="0.93", zorder=0)
ax.axvline(0, color="0.35", lw=0.8, zorder=2)
ax.set_yticks(ys)
ax.set_yticklabels(gen, fontsize=4.9, style="italic")
ax.set_xlabel("mean log$_2$ ratio, tumour vs adjacent\n(± 95% CI; diamonds pooled)",
              fontsize=5.4)
ax.set_ylim(-2.9, len(gen) - 0.2)
# the plot no longer spans the old ln range, so the panel is padded to its own data and the
# reference band is drawn wide enough that it visibly brackets the largest pooled effect
_x = np.concatenate([eff[eff.genus.isin(gen)]["effect"].values / L2,
                     REP[REP.genus.isin(gen)]["pooled_effect"].values / L2])
_pad = 0.12 * (np.nanmax(_x) - np.nanmin(_x))
ax.set_xlim(np.nanmin(_x) - _pad, np.nanmax(_x) + _pad)
handles = [plt.Line2D([], [], marker="o", ls="", color=SCOL[k], ms=3.0,
                      label=SLAB[k]) for k in SCOL]
handles.append(plt.Line2D([], [], marker="D", ls="", color="0.12", ms=3.4, label="pooled"))
ax.legend(handles=handles, frameon=True, framealpha=1.0, edgecolor="none",
          facecolor="white", fontsize=3.8, loc="lower center", ncol=2,
          handletextpad=0.15, columnspacing=0.6, labelspacing=0.10,
          borderaxespad=0.10)
lbl(ax, "C", f"{len(gen)} genera, {K}+ cohorts")

# ------------------------------------------------------------------ B max effect per cohort
ax = axes[1, 1]
rows = []
for st in STUDIES:
    d = eff[eff.study == st]
    rows.append((st, d["cohens_dz"].abs().max(), int((d["cohens_dz"].abs() >= LARGE).sum())))
xs = np.arange(len(rows))
ax.bar(xs, [r[1] for r in rows], color=[SCOL[r[0]] for r in rows], edgecolor="0.3",
       lw=0.4, width=0.62, zorder=3)
for x, (_, m, n) in zip(xs, rows):
    ax.text(x, m + 0.035, f"{m:.2f}\nn={n}", ha="center", fontsize=4.6, color="0.25")
ax.axhline(LARGE, color=RED, lw=1.0, ls=(0, (4, 1.6)), zorder=4)
ax.set_xticks(xs)
ax.set_xticklabels([SLAB[r[0]] for r in rows], rotation=42, ha="right", fontsize=5.0)
ax.set_ylabel("max |Cohen's $d_z$| in cohort", fontsize=5.8)
ax.set_ylim(0, max(r[1] for r in rows) * 1.30)
lbl(ax, "D", "No large effect replicates")

# ------------------------------------------------------------------ C novelty
ax = axes[0, 2]
n = nov.copy()
n["prior"] = np.where(n["oscc_papers"].fillna(0) > 0, "reported", "none")
cnt = n.groupby("prior")["genus"].size().reindex(["none", "reported"]).fillna(0)
ax.bar([0, 1], cnt.values, color=["#bdbdbd", BLUE], edgecolor="0.3", lw=0.4,
       width=0.55, zorder=3)
for x, v in zip([0, 1], cnt.values):
    ax.text(x, v + 0.18, f"{int(v)}", ha="center", fontsize=5.6, color="0.25")
ax.set_xticks([0, 1])
ax.set_xticklabels(["no prior\nOSCC report", "previously\nreported"], fontsize=5.0)
ax.set_ylabel("genera", fontsize=5.8)
ax.set_ylim(0, max(cnt.values) * 1.28)
ax.set_xlabel("prior OSCC literature", fontsize=5.6)
lbl(ax, "B", f"{int(cnt.loc['none'])}/{int(cnt.sum())} never reported")

# ------------------------------------------------------------------ D effect matrix
ax = axes[1, 2]
M = np.full((len(gen), len(STUDIES)), np.nan)
for i, g in enumerate(gen):
    for j, st in enumerate(STUDIES):
        v = eff[(eff.genus == g) & (eff.study == st)]["cohens_dz"]
        if len(v):
            M[i, j] = v.iloc[0]
im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-1.6, vmax=1.6)
ax.set_xticks(np.arange(len(STUDIES)))
ax.set_xticklabels([SLAB[s] for s in STUDIES], rotation=42, ha="right", fontsize=4.8)
ax.set_yticks(np.arange(len(gen)))
ax.set_yticklabels(gen, fontsize=4.4, style="italic")
for i in range(len(gen)):
    for j in range(len(STUDIES)):
        if np.isfinite(M[i, j]):
            ax.text(j, i, f"{M[i,j]:.1f}", ha="center", va="center", fontsize=3.8,
                    color="white" if abs(M[i, j]) > 0.95 else "0.15")
cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.03)
cb.ax.tick_params(labelsize=4.2)
cb.set_label("$d_z$", fontsize=4.6)
ax.set_xlabel("cohort", fontsize=5.6)
lbl(ax, "E", "Direction, not size")

# ------------------------------------------------------------------ E heterogeneity
ax = axes[2, 0]
ys = np.arange(len(gen))
cols = [RED if v > 50 else ORANGE if v > 25 else BLUE for v in REP.I2_pct]
ax.barh(ys, REP.I2_pct.values, color=cols, edgecolor="0.3", lw=0.35, height=0.6, zorder=3)
ax.set_yticks(ys)
ax.set_yticklabels(gen, fontsize=4.4, style="italic")
ax.invert_yaxis()
ax.set_xlabel("I² heterogeneity (%)", fontsize=5.8)
ax.set_xlim(0, 100)
ax.text(0.97, 0.03, "low I² = same direction and size", transform=ax.transAxes,
        fontsize=4.3, ha="right", va="bottom", color="0.35", style="italic")
lbl(ax, "F", "Most are homogeneous")

# ------------------------------------------------------------------ F |dz| by cohort count
ax = axes[2, 1]
data, labels = [], []
for k in [1, 2, 3, 4]:
    d = meta[meta.n_cohorts == k]["max_abs_dz"].dropna()
    if len(d):
        data.append(d.values); labels.append(f"k={k}\n(n={len(d)})")
bp = ax.boxplot(data, positions=np.arange(len(data)), widths=0.5, patch_artist=True,
                medianprops=dict(color="0.15", lw=1.0), boxprops=dict(lw=0.5),
                whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5),
                flierprops=dict(ms=1.4, mfc="0.6", mec="none"))
for patch in bp["boxes"]:
    patch.set_facecolor("#cfe0ee"); patch.set_edgecolor("0.35")
ax.axhline(LARGE, color=RED, lw=1.0, ls=(0, (4, 1.6)), zorder=4)
ax.set_xticks(np.arange(len(data)))
ax.set_xticklabels(labels, fontsize=4.6)
ax.set_xlabel("cohorts testing the genus", fontsize=5.6)
ax.set_ylabel("max |$d_z$| across cohorts", fontsize=5.8)
lbl(ax, "G", "Big effects are single-cohort")

# ------------------------------------------------------------------ H inflation
ax = axes[2, 2]
# Panel H compares two standardised quantities: the largest per-cohort |d_z| and the pooled |d_z|.
# It previously compared max |d_z| against the pooled LOG-RATIO, which are different scales.
ks3 = DZ[DZ.k >= K].copy()   # DZ already carries max_abs_dz and pooled_dz
ax.scatter(ks3.max_abs_dz.abs(), ks3.pooled_dz.abs(), s=13, color=BLUE,
           alpha=0.55, edgecolor="none", zorder=3)
lim = float(max(ks3.max_abs_dz.abs().max(), ks3.pooled_dz.abs().max())) * 1.08
ax.plot([0, lim], [0, lim], color="0.55", lw=0.9, ls=(0, (4, 1.6)), zorder=2)
ax.axhline(LARGE, color=RED, lw=0.9, ls=":", zorder=2)
ax.axvline(LARGE, color=RED, lw=0.9, ls=":", zorder=2)
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel("max single-cohort |$d_z$|", fontsize=5.8)
ax.set_ylabel("pooled |$d_z$|", fontsize=5.8)
ax.text(0.03, 0.97, "below the diagonal:\npooling shrinks it", transform=ax.transAxes,
        fontsize=4.4, ha="left", va="top", color="0.35", style="italic")
lbl(ax, "H", "Pooling removes inflation")

out = FIG / "manuscript-v6" / "supplementary" / "FigS2-cross-cohort-meta.png"
fig.savefig(out, dpi=330, facecolor="white")
print(f"wrote {out.relative_to(ROOT)}")
try:
    panel_crops(fig)
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
