#!/usr/bin/env python3
"""Figure for S6-v2 (TCMA / TCGA-HNSC host-side coupling).

Panel A  Direction test: the S5-replicating genera in TCMA (155 tumour vs 22 adjacent
         normal), log2 fold change, coloured by agreement with the 16S direction.
         Approximate matches (genus renamed/split since the TCMA reference) are marked.
Panel B  Paired per-patient test (22 patients with both aliquots): -log10 q, with the
         FDR line. Three genera pass.
Panel C  Sample-type composition, so the 155-vs-22 imbalance is visible rather than
         hidden.

Run: TCMA_ROOT=<dir> python3 s6v2_plot.py
"""
from __future__ import annotations

import json
import os
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


def locate_root() -> Path:
    env = os.environ.get("TCMA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results" / "tables").is_dir():
            return cand
    return here.parents[4]


ROOT = locate_root()
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

d = pd.read_csv(TAB / "s6v2_direction_check.csv")
summary = json.loads((TAB / "s6v2_summary.json").read_text(encoding="utf-8"))
AGREE, DISAGREE, MISSING = "#b2182b", "#2166ac", "0.72"

apply_figure_style(sizes=(9.6, 7.2, 6.4))
fig = plt.figure(figsize=(4.7, 6.4))
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.55], hspace=0.75)

# ---------------- panel a ----------------
ax = fig.add_subplot(gs[0, 0])
det = d[d["detected"].astype(bool)].copy().sort_values("log2fc")
y = np.arange(len(det))
cols = [AGREE if a else DISAGREE for a in det["direction_agrees"]]
ax.barh(y, det["log2fc"], color=cols, height=0.62, zorder=3)
ax.axvline(0, color="0.4", lw=0.9, zorder=4)
labels = [g + (" *" if ap else "") for g, ap in zip(det["genus"], det["approximate_match"])]
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=6.6)
ax.set_xlabel("log$_2$ (tumour / adjacent normal)\nTCMA, TCGA-HNSC", fontsize=6.8)
na, nt = int(det["direction_agrees"].sum()), len(det)
ax.set_title(f"{na}/{nt} agree with the 16S direction\n"
             f"(sign test p = {summary['s6a']['sign_test_p']:.2f}; "
             f"* = genus renamed since TCMA)", multialignment="center", loc="center")
panel_letter(ax, "A", case="upper")
ax.margins(y=0.05)
ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 0.9)

# ---------------- panel b ----------------
ax = fig.add_subplot(gs[1, 0])
p = det[det["q_paired"].notna()].copy()
p["nlq"] = -np.log10(p["q_paired"])
p = p.sort_values("nlq")
yy = np.arange(len(p))
bc = [AGREE if q < 0.05 else MISSING for q in p["q_paired"]]
ax.barh(yy, p["nlq"], color=bc, height=0.62, zorder=3)
ax.axvline(-np.log10(0.05), color="#b2182b", lw=1.0, ls=(0, (4, 2)), zorder=4)
ax.text(-np.log10(0.05), len(p) - 0.2, " FDR 0.05", fontsize=5.9, color="#b2182b",
        va="top", ha="left")
ax.set_yticks(yy)
ax.set_yticklabels(p["genus"], fontsize=6.2)
ax.set_xlabel("$-\\log_{10}$ q (paired Wilcoxon,\nn = 22 patients)", fontsize=6.8)
ax.set_title(f"{int((p['q_paired'] < 0.05).sum())} genera pass FDR\n"
             f"in the paired subset", multialignment="center", loc="center")
panel_letter(ax, "B", case="upper")
ax.margins(y=0.08)

# ---------------- panel c ----------------
ax = fig.add_subplot(gs[2, 0])
counts = summary["sample_type_counts"]
order = [("Primary Tumor", "#7b3294"), ("Solid Tissue Normal", "#a6dba0")]
left = 0.0
for key, col in order:
    v = counts.get(key, 0)
    ax.barh([0], [v], left=left, color=col, height=0.45, zorder=3)
    ax.text(left + v / 2, 0, str(v), ha="center", va="center", fontsize=7,
            color="white", fontweight="bold")
    left += v
ax.set_yticks([])
ax.set_xlabel("samples", fontsize=6.8)
ax.set_title("Sample types\n(177 HNSC aliquots)", multialignment="center", loc="center")
panel_letter(ax, "C", case="upper")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in order]
ax.legend(handles, ["Primary tumour", "Adjacent normal"], frameon=False, fontsize=5.6,
          loc="upper right", bbox_to_anchor=(1.0, 1.28), ncol=1,
          handlelength=1.0, labelspacing=0.25)
ax.set_ylim(-0.4, 0.4)

out = FIG / "s6-host-coupling-TCGA-HNSC.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print(f"file: {out.name} {out.stat().st_size:,} bytes")
try:
    panel_crops(fig)
    print("crops: .cache/figure-style/")
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
