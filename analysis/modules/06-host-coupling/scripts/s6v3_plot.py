#!/usr/bin/env python3
"""Figure for S6-v3: site-stratified analysis of TCGA-HNSC.

Three panels:
  a  Parvimonas across the pooled cohort and each site group. Shows the genus is NOT
     rescued by stratification: in the clean oral group the effect is flat.
  b  The three genera that do replicate, pooled vs oral-only, with FDR q values from
     both analyses side by side. Shows stratification STRENGTHENS the signal because
     the larynx samples point the other way.
  c  Testability per site group (tumour/normal counts), so the larynx/pharynx
     limitations are visible instead of hidden.

Run: python3 s6v3_plot.py
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

res = pd.read_csv(TAB / "s6v3_stratified_tests.csv")
pooled = pd.read_csv(TAB / "s6v2_direction_check.csv")
summary = json.loads((TAB / "s6v3_summary.json").read_text(encoding="utf-8"))

ORAL, ORALAMB, LAR = "#1f78b4", "#8da0cb", "#d95f02"
GREY = "#9e9e9e"

apply_figure_style(sizes=(9.8, 7.2, 6.4))
fig = plt.figure(figsize=(7.09, 5.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], hspace=0.55, wspace=0.38)

# ---------------- panel a: Parvimonas by site ----------------
ax = fig.add_subplot(gs[0, :])
pv = res[res["genus"] == "Parvimonas"].copy()
pv["label"] = pv["site_group"].str.replace("*", "\u2020", regex=False)
pv = pv.sort_values("log2fc")
yy = np.arange(len(pv))
cols = [ORAL if g == "oral" else ORALAMB if g.startswith("oral") else
        LAR if g == "larynx" else GREY for g in pv["site_group"]]
vals = pv["log2fc"].fillna(0)
testable = pv["log2fc"].notna().values
bars = ax.barh(yy, vals, color=[c if t else "white" for c, t in zip(cols, testable)],
               edgecolor=cols, linewidth=1.0, height=0.6, zorder=3)
ax.axvline(0, color="0.4", lw=0.9, zorder=4)
ax.set_yticks(yy)
ax.set_yticklabels([f"{l}\n(n={int(t)}T/{int(n)}N)"
                    for l, t, n in zip(pv["label"], pv["n_tumour"], pv["n_normal"])],
                   fontsize=6.2)
ax.set_xlabel("log$_2$ (tumour / normal)", fontsize=6.8)
ax.set_title("Parvimonas is not rescued\nby site stratification", multialignment="center", loc="center")
panel_letter(ax, "A", case="upper")
for i, (v, p, t) in enumerate(zip(pv["log2fc"], pv["p_paired"], testable)):
    if pd.notna(v):
        txt = f"{v:+.2f}" + (f"  p={p:.2f}" if pd.notna(p) else "")
        ax.text(v + (0.12 if v >= 0 else -0.12), i, txt, fontsize=5.6,
                va="center", ha="left" if v >= 0 else "right", color="0.25")
ax.margins(x=0.28, y=0.08)

# ---------------- panel b: pooled vs oral q for the replicating genera -------
ax = fig.add_subplot(gs[1, 0])
rep_gen = ["Actinomyces", "Lancefieldella", "Veillonella"]
rowsB = []
for g in rep_gen:
    qp = pooled.loc[pooled["genus"] == g, "q_paired"]
    qo = res.loc[(res["genus"] == g) & (res["site_group"] == "oral"), "q_paired"]
    lo = res.loc[(res["genus"] == g) & (res["site_group"] == "larynx"), "log2fc"]
    rowsB.append({"genus": g,
                  "pooled": float(qp.iloc[0]) if len(qp) and pd.notna(qp.iloc[0]) else np.nan,
                  "oral": float(qo.iloc[0]) if len(qo) and pd.notna(qo.iloc[0]) else np.nan,
                  "larynx_lfc": float(lo.iloc[0]) if len(lo) and pd.notna(lo.iloc[0]) else np.nan})
B = pd.DataFrame(rowsB)
x = np.arange(len(B))
w = 0.36
ax.bar(x - w / 2, -np.log10(B["pooled"]), w, color=GREY, zorder=3, label="pooled (177)")
ax.bar(x + w / 2, -np.log10(B["oral"]), w, color=ORAL, zorder=3, label="oral only (95)")
ax.axhline(-np.log10(0.05), color="#b2182b", lw=1.0, ls=(0, (4, 2)), zorder=4)
ax.text(len(B) - 0.55, -np.log10(0.05) + 0.12, "FDR 0.05", fontsize=5.9, color="#b2182b")
for i, (qp, qo) in enumerate(zip(B["pooled"], B["oral"])):
    ax.text(i - w / 2, -np.log10(qp) + 0.08, f"{qp:.3f}", ha="center", fontsize=5.5, color="0.25")
    ax.text(i + w / 2, -np.log10(qo) + 0.08, f"{qo:.4f}", ha="center", fontsize=5.5, color=ORAL)
ax.set_xticks(x)
ax.set_xticklabels(B["genus"], fontsize=6.4)
ax.set_ylabel("$-\\log_{10}$ q (paired Wilcoxon)", fontsize=6.8)
ax.set_title("Restricting to oral sites strengthens\nthe three replicating genera",
             multialignment="center", loc="center")
panel_letter(ax, "B", case="upper")
ax.set_ylim(0, 4.75)
ax.legend(frameon=False, fontsize=5.9, loc="upper right", handlelength=1.1,
          labelspacing=0.25, borderaxespad=0.2)
ax.margins(y=0.18)

# ---------------- panel c: testability ----------------
ax = fig.add_subplot(gs[1, 1])
counts = pd.crosstab(res.drop_duplicates("site_group")["site_group"], 1).iloc[:, 0]
tn = res.drop_duplicates("site_group").set_index("site_group")[["n_tumour", "n_normal", "n_pairs"]]
tn = tn.reindex(["oral", "oral*", "larynx", "pharynx"])
xx = np.arange(len(tn))
ax.bar(xx - 0.22, tn["n_tumour"], 0.42, color="#7b3294", zorder=3, label="tumour")
ax.bar(xx + 0.22, tn["n_normal"], 0.42, color="#a6dba0", zorder=3, label="normal")
for i, p in enumerate(tn["n_pairs"]):
    ax.text(i, max(tn["n_tumour"].iloc[i], tn["n_normal"].iloc[i]) + 3,
            f"{int(p)} pr", ha="center", fontsize=5.6, color="0.3")
ax.set_xticks(xx)
ax.set_xticklabels([s.replace("*", "\u2020") for s in tn.index], fontsize=6.2, rotation=25,
                   ha="right")
ax.set_ylabel("samples", fontsize=6.8)
ax.set_title("Testability by site", multialignment="center", loc="center")
panel_letter(ax, "C", case="upper")
ax.legend(frameon=False, fontsize=5.7, loc="upper right", handlelength=1.0,
          labelspacing=0.2)
ax.margins(y=0.22)

out = FIG / "s6-site-stratified-TCGA-HNSC.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print(f"file: {out.name} {out.stat().st_size:,} bytes")
try:
    panel_crops(fig)
    print("crops: .cache/figure-style/")
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
