#!/usr/bin/env python3
"""Figure for S6-v4: host-programme coupling from the TCMA correlation panels.

Panel A  Heatmap of -log10(min q) per genus x assay. Parvimonas is the only genus in
         the panel with no significant coupling anywhere.
Panel B  Parvimonas's own strongest couplings vs the 0.05 line, to make the null
         explicit rather than a footnote.

Run: python3 s6v4_plot.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
SKILL = _FIGSTYLE
sys.path.insert(0, str(SKILL))
from runtime import apply_figure_style, panel_crops, panel_letter

ROOT = Path(__file__).resolve().parents[4]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

per = pd.read_csv(TAB / "s6v4_genus_summary.csv")
pvt = pd.read_csv(TAB / "s6v4_parvimonas_top.csv")
summary = json.loads((TAB / "s6v4_summary.json").read_text(encoding="utf-8"))

ASSAYS = ["Methylation27K", "RNAseq", "RPPA", "miRNA"]
mat = (per.pivot(index="genus", columns="assay", values="min_q")
          .reindex(columns=[a for a in ASSAYS if a in per["assay"].unique()]))
mat = mat.loc[mat.min(axis=1).sort_values().index]          # best-coupled genus on top
nlq = -np.log10(mat.clip(lower=1e-4))

apply_figure_style(sizes=(9.4, 7.0, 6.4))
fig = plt.figure(figsize=(9.4, 3.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.2], wspace=0.42)

# ---------------- panel a ----------------
ax = fig.add_subplot(gs[0, 0])
vmax = float(np.nanmax(nlq.values))
im = ax.imshow(nlq.values, cmap="viridis", vmin=0, vmax=vmax, aspect="auto")
ax.set_xticks(np.arange(nlq.shape[1]))
ax.set_xticklabels([c.replace("Methylation27K", "Methyl-\nation") for c in nlq.columns],
                   fontsize=6.2)
ax.set_yticks(np.arange(nlq.shape[0]))
ax.set_yticklabels(nlq.index, fontsize=6.4)
# highlight Parvimonas row
for i, g in enumerate(nlq.index):
    if g == "Parvimonas":
        ax.add_patch(plt.Rectangle((-0.5, i - 0.5), nlq.shape[1], 1, fill=False,
                                   edgecolor="#b2182b", lw=1.6, zorder=5))
# annotate hits
for i in range(nlq.shape[0]):
    for j in range(nlq.shape[1]):
        q = mat.values[i, j]
        txt = "" if pd.isna(q) else ("ns" if q >= 0.05 else f"{q:.3f}")
        if txt:
            ax.text(j, i, txt, ha="center", va="center", fontsize=4.9,
                    color="white" if nlq.values[i, j] > vmax * 0.55 else "0.15")
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cb.set_label("$-\\log_{10}$(min q)", fontsize=6.4)
cb.ax.tick_params(labelsize=5.8)
ax.set_title("Best host coupling per genus x assay\n"
             "(cell = min q over that panel's 50 features)",
             multialignment="center")
panel_letter(ax, "A", case="upper")

# ---------------- panel b ----------------
ax = fig.add_subplot(gs[0, 1])
p = pvt.head(12).copy()
p["label"] = p["assay"].str.replace("Methylation27K", "Meth") + "  " + p["feature"]
p = p.sort_values("q", ascending=False)
yy = np.arange(len(p))
cols = ["#b2182b" if q < 0.05 else "0.72" for q in p["q"]]
ax.barh(yy, -np.log10(p["q"]), color=cols, height=0.64, zorder=3)
ax.axvline(-np.log10(0.05), color="#b2182b", lw=1.0, ls=(0, (4, 2)), zorder=4)
ax.text(-np.log10(0.05), len(p) - 0.4, " FDR 0.05", fontsize=5.9, color="#b2182b",
        va="top", ha="left")
ax.set_yticks(yy)
ax.set_yticklabels(p["label"], fontsize=5.5)
ax.set_xlabel("$-\\log_{10}$ q  (best 12 couplings of Parvimonas)", fontsize=6.8)
ax.set_title(f"Parvimonas couples with nothing\n"
             f"(0/200 cells; min q = {summary['parvimonas']['min_q']:.3f})", multialignment="center")
panel_letter(ax, "B", case="upper")
ax.margins(y=0.06)

out = FIG / "s6-host-coupling-panels-TCGA-HNSC.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print(f"file: {out.name} {out.stat().st_size:,} bytes")
try:
    panel_crops(fig)
    print("crops: .cache/figure-style/")
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
