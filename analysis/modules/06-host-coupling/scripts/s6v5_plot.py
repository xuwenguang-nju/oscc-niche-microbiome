#!/usr/bin/env python3
"""Figure for S6-v5: tongue-only and finer oral subsites.

Panel A  Parvimonas log2FC across every cohort grouping tried in S6, from the pooled
         177 samples down to tongue only. Shows the sign is unstable and never
         significant -- the honest statement.
Panel B  The three genera that replicate: q values pooled vs oral vs tongue, showing
         the signal sharpens as the anatomy gets more homogeneous.
Panel C  Subset sizes, so the reader can see why tongue is the only usable subsite.

Run: python3 s6v5_plot.py
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

pooled = pd.read_csv(TAB / "s6v2_direction_check.csv")
v3 = pd.read_csv(TAB / "s6v3_stratified_tests.csv")
v5 = pd.read_csv(TAB / "s6v5_subsite_tests.csv")
v5s = json.loads((TAB / "s6v5_summary.json").read_text(encoding="utf-8"))
INCONS, CONS, GREY = "#d95f02", "#1f78b4", "0.72"

apply_figure_style(sizes=(9.6, 7.0, 6.4))
fig = plt.figure(figsize=(9.6, 3.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.15, 0.8], wspace=0.55)

# ---------------- panel a: Parvimonas across every grouping ----------------
ax = fig.add_subplot(gs[0, 0])
rows = []
p_pool = pooled.loc[pooled["genus"] == "Parvimonas"].iloc[0]
rows.append(("pooled (155/22)", float(p_pool["log2fc"]), float(p_pool["p_paired"])))
for sg in ["oral", "larynx"]:
    r = v3[(v3.genus == "Parvimonas") & (v3.site_group == sg)]
    if len(r):
        rows.append((f"{sg} (v3)", float(r["log2fc"].iloc[0]),
                     float(r["p_paired"].iloc[0]) if pd.notna(r["p_paired"].iloc[0]) else np.nan))
for sub in ["tongue_all", "tongue_NOS"]:
    r = v5[(v5.genus == "Parvimonas") & (v5.subset == sub)]
    if len(r):
        rows.append((f"{sub.replace('_', ' ')} (61/16)" if sub == "tongue_all"
                     else f"{sub.replace('_', ' ')} (54/15)",
                     float(r["log2fc"].iloc[0]), float(r["p_paired"].iloc[0])))
A = pd.DataFrame(rows, columns=["label", "lfc", "p"])
yy = np.arange(len(A))
cols = [CONS if v > 0 else INCONS for v in A["lfc"]]
ax.barh(yy, A["lfc"], color=cols, height=0.6, zorder=3)
ax.axvline(0, color="0.4", lw=0.9, zorder=4)
ax.set_yticks(yy)
ax.set_yticklabels(A["label"], fontsize=6.1)
for i, (v, p) in enumerate(zip(A["lfc"], A["p"])):
    txt = f"{v:+.2f}" + (f" p={p:.2f}" if pd.notna(p) else " n/a")
    ax.text(v + (0.06 if v >= 0 else -0.06), i, txt, fontsize=5.6, va="center",
            ha="left" if v >= 0 else "right", color="0.25")
ax.set_xlabel("Parvimonas log$_2$ (tumour / normal)", fontsize=6.8)
ax.set_title("Parvimonas: sign flips,\nnever significant",
             multialignment="center", loc="center", fontsize=7.4)
panel_letter(ax, "A", case="upper")
ax.margins(x=0.45, y=0.08)

# ---------------- panel b: the three replicating genera ----------------
ax = fig.add_subplot(gs[0, 1])
gen = ["Actinomyces", "Lancefieldella", "Veillonella"]
wide = []
for g in gen:
    q_pool = float(pooled.loc[pooled.genus == g, "q_paired"].iloc[0])
    q_oral = float(v3[(v3.genus == g) & (v3.site_group == "oral")]["q_paired"].iloc[0])
    q_tong = float(v5[(v5.genus == g) & (v5.subset == "tongue_all")]["q_paired"].iloc[0])
    wide.append({"genus": g, "pooled": q_pool, "oral": q_oral, "tongue": q_tong})
W = pd.DataFrame(wide)
x = np.arange(len(W)); w = 0.26
for k, (col, c, lab) in enumerate([("pooled", GREY, "pooled (177)"),
                                   ("oral", "#4a98c9", "oral (95)"),
                                   ("tongue", "#1f4e79", "tongue (77)")]):
    ax.bar(x + (k - 1) * w, -np.log10(W[col]), w, color=c, zorder=3, label=lab)
    for i, q in enumerate(W[col]):
        ax.text(i + (k - 1) * w, -np.log10(q) + 0.07, f"{q:.4f}" if q < 0.01 else f"{q:.3f}",
                ha="center", fontsize=5.0, color="0.25", rotation=90)
ax.axhline(-np.log10(0.05), color="#b2182b", lw=1.0, ls=(0, (4, 2)), zorder=4)
ax.text(len(W) - 0.5, -np.log10(0.05) + 0.1, "FDR 0.05", fontsize=5.9, color="#b2182b")
ax.set_xticks(x)
ax.set_xticklabels(W["genus"], fontsize=6.3, rotation=18, ha="right")
ax.set_ylabel("$-\\log_{10}$ q (paired Wilcoxon)", fontsize=6.8)
ax.set_title("Signal sharpens with\nhomogeneous anatomy",
             multialignment="center", loc="center", fontsize=7.4)
panel_letter(ax, "B", case="upper")
ax.legend(frameon=False, fontsize=5.7, loc="upper center", ncol=3, handlelength=1.0,
          labelspacing=0.2, columnspacing=1.0, borderaxespad=0.1)
ax.margins(y=0.28)

# ---------------- panel c: subset sizes ----------------
ax = fig.add_subplot(gs[0, 2])
S = pd.DataFrame([
    {"subset": "tongue\nall", **{k: v5s["subsets"]["tongue_all"][k]
                                 for k in ("n_tumour", "n_normal", "n_pairs")}},
    {"subset": "tongue\nNOS", **{k: v5s["subsets"]["tongue_NOS"][k]
                                 for k in ("n_tumour", "n_normal", "n_pairs")}},
    {"subset": "floor of\nmouth", **{k: v5s["subsets"]["floor_of_mouth"][k]
                                     for k in ("n_tumour", "n_normal", "n_pairs")}},
    {"subset": "oral\nothers", **{k: v5s["subsets"]["oral_others"][k]
                                  for k in ("n_tumour", "n_normal", "n_pairs")}},
])
xx = np.arange(len(S))
ax.bar(xx - 0.2, S["n_tumour"], 0.38, color="#7b3294", zorder=3, label="tumour")
ax.bar(xx + 0.2, S["n_normal"], 0.38, color="#a6dba0", zorder=3, label="normal")
for i, p in enumerate(S["n_pairs"]):
    ax.text(i, max(S["n_tumour"].iloc[i], S["n_normal"].iloc[i]) + 2.5, f"{int(p)} pr",
            ha="center", fontsize=5.6, color="0.3")
ax.set_xticks(xx)
ax.set_xticklabels([lbl.replace(" ", "\n", 1) for lbl in S["subset"]],
                   fontsize=5.4, rotation=0, ha="center")
ax.set_ylabel("samples", fontsize=6.8)
ax.set_title("Subsite sizes", multialignment="center", loc="center", fontsize=7.4)
panel_letter(ax, "C", case="upper")
ax.legend(frameon=False, fontsize=5.7, loc="upper right", handlelength=1.0,
          labelspacing=0.2)
ax.margins(y=0.20)

out = FIG / "s6-tongue-subsite-TCGA-HNSC.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print(f"file: {out.name} {out.stat().st_size:,} bytes")
try:
    panel_crops(fig)
    print("crops: .cache/figure-style/")
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
