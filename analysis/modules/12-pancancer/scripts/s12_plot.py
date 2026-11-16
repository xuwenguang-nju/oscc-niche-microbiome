#!/usr/bin/env python3
"""F8 — pan-cancer context: is the oral signal specific to the aerodigestive tract?

Nine panels on the TCMA v2 five-project release (HNSC, ESCA, STAD, COAD, READ). The point
is that the HNSC microbial profile is not "a tumour microbiome"; it is an oral community,
and it separates cleanly from the colorectal one.

A  PCoA of genus profiles, coloured by project
B  the oral-vs-gut indicator axis per project
C  Parvimonas, the oral marker, by project
D  genus x project mean abundance (the 12 most variable genera)
E  per-sample distribution of the oral-minus-gut axis
F  sample-type sensitivity: tumour vs normal within each project
G  number of oral genera detected per sample
H  prevalence vs mean abundance for the oral markers, by project
I  project composition as a stacked bar

Reads:  results/tables/s12_pcoa_coordinates.csv, s12_oral_gut_axis.csv,
        s12_project_composition.csv, s12_focus_genera_by_project.csv,
        s12_sample_type_sensitivity.csv, s12_oral_gut_per_sample.csv
Writes: figures/manuscript-v6/main/Fig7-pan-cancer.png
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
from figlib import panel_letter_fixed, project_root  # noqa: E402

ROOT = project_root(Path(__file__).resolve().parent)
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures" / "manuscript-v6" / "main"
# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
exec(compile(open(_FIGSTYLE / "runtime.py", encoding="utf-8").read(),
              "runtime.py", "exec"))

apply_figure_style(sizes=(7.0, 6.2, 5.4))

pcoa = pd.read_csv(TAB / "s12_pcoa_coordinates.csv", index_col=0)
axis = pd.read_csv(TAB / "s12_oral_gut_axis.csv")
comp = pd.read_csv(TAB / "s12_project_composition.csv")
focus = pd.read_csv(TAB / "s12_focus_genera_by_project.csv")
sens = pd.read_csv(TAB / "s12_sample_type_sensitivity.csv")
per = pd.read_csv(TAB / "s12_oral_gut_per_sample.csv")

PROJ = ["HNSC", "ESCA", "STAD", "COAD", "READ"]
PCOL = {"HNSC": "#b2182b", "ESCA": "#e66101", "STAD": "#7b3294",
        "COAD": "#2166ac", "READ": "#1b7837"}
AERO = ["HNSC", "ESCA", "STAD"]
GUT = ["COAD", "READ"]
BLUE, RED, ORANGE, GREEN = "#2166ac", "#b2182b", "#a6611a", "#1b7837"

fig = plt.figure(figsize=(7.09, 7.05))
gs = fig.add_gridspec(3, 3, left=0.134, right=0.985, top=0.968, bottom=0.058,
                      wspace=0.62, hspace=0.68)
axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(3)] for i in range(3)])


def lbl(ax, letter, text):
    ax.set_title(text, fontsize=6.2, multialignment="center", pad=5)
    panel_letter_fixed(ax, letter, case="upper")


# ------------------------------------------------------------------ A PCoA
ax = axes[0, 0]
for p in PROJ:
    d = pcoa[pcoa.project == p]
    ax.scatter(d["PCoA1"], d["PCoA2"], s=4.0, color=PCOL[p], alpha=0.65,
               linewidths=0, zorder=3, label=p)
ax.set_xlabel("PCoA1", fontsize=5.8)
ax.set_ylabel("PCoA2", fontsize=5.8)
ax.legend(frameon=False, fontsize=4.5, loc="upper right", handletextpad=0.1,
          labelspacing=0.14, borderaxespad=0.2, markerscale=2.4)
lbl(ax, "A", "Projects separate")

# ------------------------------------------------------------------ B oral-gut axis
ax = axes[0, 1]
d = axis.set_index("project").reindex(PROJ)
xs = np.arange(len(PROJ))
cols = [PCOL[p] for p in PROJ]
ax.bar(xs, d["oral_minus_gut_pct"], color=cols, edgecolor="0.3", lw=0.4, width=0.62, zorder=3)
ax.axhline(0, color="0.3", lw=0.8, zorder=4)
for x, (v, n) in enumerate(zip(d["oral_minus_gut_pct"], d["n"])):
    ax.text(x, v + (3.5 if v >= 0 else -3.5), f"{v:+.0f}", ha="center",
            va="bottom" if v >= 0 else "top", fontsize=4.8, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels(PROJ, rotation=42, ha="right", fontsize=5.0)
ax.set_ylabel("oral − gut indicator genera (%)", fontsize=5.8)
ax.set_ylim(-72, 108)
lbl(ax, "B", "Signature flips by site")

# ------------------------------------------------------------------ C Parvimonas
ax = axes[0, 2]
f = focus[focus.genus == "Parvimonas"].set_index("project").reindex(PROJ)
ax.bar(xs, f["mean_pct"], color=cols, edgecolor="0.3", lw=0.4, width=0.62, zorder=3)
for x, v in zip(xs, f["mean_pct"]):
    ax.text(x, v + 0.028, f"{v:.2f}", ha="center", fontsize=4.8, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels(PROJ, rotation=42, ha="right", fontsize=5.0)
ax.set_ylabel("$\\it{Parvimonas}$ mean %", fontsize=5.8)
ax.set_ylim(0, f["mean_pct"].max() * 1.30)
lbl(ax, "C", "Parvimonas peaks in HNSC")

# ------------------------------------------------------------------ D genus x project
ax = axes[1, 0]
M = comp.set_index("genus")[PROJ]
M = M.reindex(M.max(axis=1).sort_values(ascending=False).index).head(12)
Z = M.div(M.max(axis=1), axis=0)
im = ax.imshow(Z.values, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
ax.set_xticks(np.arange(len(PROJ)))
ax.set_xticklabels(PROJ, rotation=42, ha="right", fontsize=4.7)
ax.set_yticks(np.arange(len(Z)))
ax.set_yticklabels(Z.index, fontsize=4.3, style="italic")
ax.set_xlabel("project", fontsize=5.6)
cb = fig.colorbar(im, ax=ax, fraction=0.024, pad=0.045)
cb.ax.tick_params(labelsize=4.2)
cb.set_label("relative to genus max", fontsize=4.5)
lbl(ax, "D", "Top genera, by project")

# ------------------------------------------------------------------ E per-sample axis
ax = axes[1, 1]
data, labs = [], []
for p in PROJ:
    v = per[per.project == p]["oral_minus_gut_pct"].dropna().to_numpy()
    if len(v):
        data.append(v); labs.append(p)
bp = ax.boxplot(data, positions=np.arange(len(data)), widths=0.55, patch_artist=True,
                medianprops=dict(color="0.15", lw=1.0), boxprops=dict(lw=0.5),
                whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5),
                flierprops=dict(ms=1.2, mfc="0.6", mec="none"))
for patch, p in zip(bp["boxes"], labs):
    patch.set_facecolor(PCOL[p]); patch.set_alpha(0.55); patch.set_edgecolor("0.35")
ax.axhline(0, color="0.3", lw=0.8, zorder=4)
ax.set_xticks(np.arange(len(labs)))
ax.set_xticklabels(labs, rotation=42, ha="right", fontsize=4.9)
ax.set_ylabel("oral − gut axis, per sample", fontsize=5.8)
lbl(ax, "E", "Projects differ, with overlap")

# ------------------------------------------------------------------ F sample-type sensitivity
ax = axes[1, 2]
t = sens[sens.sample_type == "Primary Tumor"].set_index("project").reindex(PROJ)
n = sens[sens.sample_type == "Solid Tissue Normal"].set_index("project").reindex(PROJ)
xs2 = np.arange(len(PROJ))
ax.bar(xs2 - 0.18, t["oral_minus_gut_pct"], width=0.34, color=[PCOL[p] for p in PROJ],
       edgecolor="0.3", lw=0.35, zorder=3, label="tumour")
ax.bar(xs2 + 0.18, n["oral_minus_gut_pct"], width=0.34, color=[PCOL[p] for p in PROJ],
       alpha=0.42, edgecolor="0.3", lw=0.35, zorder=3, label="normal")
ax.axhline(0, color="0.3", lw=0.8, zorder=4)
ax.set_xticks(xs2)
ax.set_xticklabels(PROJ, rotation=42, ha="right", fontsize=4.9)
ax.set_ylabel("oral − gut axis (%)", fontsize=5.8)
ax.legend(frameon=False, fontsize=4.4, loc="upper right", handlelength=1.0,
          labelspacing=0.2, borderaxespad=0.2)
ax.text(0.02, 0.06, "normal tissue shows it too", transform=ax.transAxes,
        fontsize=4.3, color="0.35", style="italic")
lbl(ax, "F", "Not a tumour-only effect")

# ------------------------------------------------------------------ G oral genera count
ax = axes[2, 0]
d = axis.set_index("project").reindex(PROJ)
ax.bar(xs, d["median_n_oral_genera"], color=cols, edgecolor="0.3", lw=0.4,
       width=0.62, zorder=3)
for x, v in zip(xs, d["median_n_oral_genera"]):
    ax.text(x, v + 0.16, f"{v:.0f}", ha="center", fontsize=4.9, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels(PROJ, rotation=42, ha="right", fontsize=5.0)
ax.set_ylabel("oral genera detected\nper sample (median)", fontsize=5.6)
ax.set_ylim(0, d["median_n_oral_genera"].max() * 1.26)
lbl(ax, "G", "Oral genera everywhere")

# ------------------------------------------------------------------ H prevalence vs abundance
ax = axes[2, 1]
for p in PROJ:
    d = focus[focus.project == p]
    ax.scatter(d["mean_pct"], d["prevalence_pct"], s=13, color=PCOL[p], alpha=0.75,
               edgecolor="none", zorder=3, label=p)
ax.set_xscale("symlog", linthresh=0.01)
ax.set_xlabel("mean % abundance", fontsize=5.8)
ax.set_ylabel("% samples detected", fontsize=5.8)
ax.legend(frameon=False, fontsize=4.2, loc="upper left", handletextpad=0.1,
          labelspacing=0.12, borderaxespad=0.2, markerscale=1.4, ncol=1,
          columnspacing=0.6)
lbl(ax, "H", "Detected vs abundant")

# ------------------------------------------------------------------ I stacked composition
ax = axes[2, 2]
top = comp.set_index("genus")[PROJ]
top = top.reindex(top.mean(axis=1).sort_values(ascending=False).index).head(8)
mat = top.div(top.sum(axis=0), axis=1) * 100
bottoms = np.zeros(len(PROJ))
CMAP = plt.get_cmap("tab20")
for i, (g, row) in enumerate(mat.iterrows()):
    ax.bar(np.arange(len(PROJ)), row.values, bottom=bottoms, color=CMAP(i % 20),
           edgecolor="white", lw=0.3, width=0.66, zorder=3, label=g[:14])
    bottoms += row.values
ax.set_xticks(np.arange(len(PROJ)))
ax.set_xticklabels(PROJ, rotation=42, ha="right", fontsize=4.9)
ax.set_ylabel("% of top-8 genus signal", fontsize=5.6)
ax.set_ylim(0, 100)
ax.legend(frameon=False, fontsize=3.6, loc="center left", bbox_to_anchor=(1.01, 0.5),
          handletextpad=0.2, labelspacing=0.12, borderaxespad=0.0)
lbl(ax, "I", "Different composition")

out = FIG / "Fig7-pan-cancer.png"
fig.savefig(out, dpi=330, facecolor="white")
print(f"wrote {out.relative_to(ROOT)}")
try:
    panel_crops(fig)
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
