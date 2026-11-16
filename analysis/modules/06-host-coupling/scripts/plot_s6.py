#!/usr/bin/env python3
"""Figure for S6 (host-side coupling, GSE227919).

Panel A  direction of the S5-replicating genera in RNA-seq-mined microbial data
         (log2 ratio Cancer/Control), coloured by whether it matches the 16S result.
         Effect SIZES are not comparable across the two measurements, so only the
         sign is scored -- that is stated on the axis.
Panel B  genus x host-program Spearman rho; no cell survives FDR, which is the
         point of the panel.
Panel C  composition of the RNA-seq microbial signal: how much of it is genuine
         oral community versus environmental/reagent background. This is what makes
         the null in (a) uninformative rather than a refutation.

Run: python3 s6_host_coupling.py && python3 s6_signal_composition.py
     python3 plot_s6.py
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
SKILL = _FIGSTYLE / "runtime.py"
sys.path.insert(0, str(SKILL.parent))
from runtime import apply_figure_style, panel_crops, panel_crops as _pc, panel_letter


def project_root() -> Path:
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results").is_dir() and (cand / "analysis").is_dir():
            return cand
    return here.parents[4]


ROOT = project_root()
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

UP = "#b2182b"      # agrees with the 16S direction
DOWN = "#2166ac"    # disagrees
GREY = "#9e9e9e"

dir_tab = pd.read_csv(TAB / "s6a_direction_check.csv")
grp = pd.read_csv(TAB / "s6a_group_abundance.csv", index_col=0)
corr = pd.read_csv(TAB / "s6b_genus_signature_corr.csv")
comp = json.loads((TAB / "s6b_signal_composition.json").read_text(encoding="utf-8"))
summary = json.loads((TAB / "s6_summary.json").read_text(encoding="utf-8"))

apply_figure_style(sizes=(10, 7.5, 6.5))
fig = plt.figure(figsize=(10.2, 4.0))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.25, 0.75], wspace=0.55,
                       top=0.80, bottom=0.20)

# ---------------- panel a: direction check -----------------------------
ax = fig.add_subplot(gs[0, 0])
# "not_detected" is a string, not NaN, so dropna alone would keep it
d = dir_tab[dir_tab["rna_direction"].isin(["up_in_tumour", "down_in_tumour"])].copy()
d = d[d["genus"].isin(grp.index)]
# log2 ratio Cancer/Control from mean relative abundance
eps = 1e-6
d["log2fc"] = np.log2((grp.loc[d["genus"], "Cancer"].values + eps) /
                      (grp.loc[d["genus"], "Control"].values + eps))
d = d.sort_values("log2fc")
y = np.arange(len(d))
colors = [UP if bool(a) else DOWN for a in d["direction_agrees"]]
ax.barh(y, d["log2fc"], color=colors, height=0.62, zorder=3)
ax.axvline(0, color="0.4", lw=0.9, zorder=4)
ax.set_yticks(y)
ax.set_yticklabels(d["genus"], fontsize=6.6)
ax.set_xlabel("log$_2$ (Cancer / Control)\nin RNA-seq-mined microbiome", fontsize=6.8)
n_agree = int(d["direction_agrees"].sum())
n_tot = len(d)
ax.set_title(f"Direction vs 16S: {n_agree}/{n_tot} agree\n"
             f"sign test p = {summary['s6a']['sign_test_p']:.2f}",
             multialignment="center", loc="center", fontsize=7.4)
panel_letter(ax, "A", case="upper")
ax.text(0.985, 0.02, "red = same as 16S\nblue = opposite",
        transform=ax.transAxes, fontsize=5.9, color="0.35", va="bottom", ha="right")
ax.margins(y=0.04)

# ---------------- panel b: coupling heatmap ----------------------------
ax = fig.add_subplot(gs[0, 1])
pivot = corr.pivot(index="genus", columns="signature", values="spearman_rho")
qmap = corr.pivot(index="genus", columns="signature", values="q_bh")
pivot = pivot.loc[d["genus"].values[::-1]]
order = ["cytotoxic", "ifn_gamma", "checkpoint", "pemt", "proliferation", "epithelial_diff"]
SHORT = {"cytotoxic": "cyto-\ntoxic", "ifn_gamma": "IFN-$\\gamma$",
         "checkpoint": "check-\npoint", "pemt": "p-EMT",
         "proliferation": "proli-\nfer.", "epithelial_diff": "epith.\ndiff."}
pivot = pivot[[c for c in order if c in pivot.columns]]
vmax = float(np.nanmax(np.abs(pivot.values)))
im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(np.arange(pivot.shape[1]))
ax.set_xticklabels([SHORT.get(c, c) for c in pivot.columns], fontsize=5.7)
ax.set_yticks(np.arange(pivot.shape[0]))
ax.set_yticklabels(pivot.index, fontsize=5.9)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cb.set_label("Spearman $\\rho$", fontsize=6.2)
cb.ax.tick_params(labelsize=5.6)
n_sig = int((qmap.values < 0.05).sum())
ax.set_title(f"Host-program coupling\n0/{corr.shape[0]} pass FDR, min q "
             f"{corr['q_bh'].min():.2f}",
             multialignment="center", loc="center", fontsize=7.4)
panel_letter(ax, "B", case="upper")

# ---------------- panel c: signal composition --------------------------
ax = fig.add_subplot(gs[0, 2])
shares = comp["share_pct"]
labels = [("oral_typical", "oral-typical", "#1f78b4"),
          ("environmental_background", "environmental /\nreagent background", "#d95f02"),
          ("other", "other", "0.72")]
left = 0.0
for key, lab, col in labels:
    val = shares.get(key, 0.0)
    ax.barh([0], [val], left=left, color=col, height=0.5, zorder=3)
    ax.text(left + val / 2, 0, f"{val:.0f}%", ha="center", va="center",
            fontsize=6.6, color="white", fontweight="bold")
    left += val
ax.set_yticks([])
ax.set_xlim(0, 100)
ax.set_xlabel("share of mean relative abundance (%)", fontsize=6.8)
ax.set_title("RNA-seq signal\ncomposition", multialignment="center",
             loc="center", fontsize=7.4)
panel_letter(ax, "C", case="upper")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in labels]
ax.legend(handles, [l.replace("\n", " ") for _, l, _ in labels], frameon=False,
          fontsize=5.8, loc="lower center", bbox_to_anchor=(0.5, -0.52), ncol=1,
          handlelength=1.0, labelspacing=0.25)
ax.set_ylim(-0.45, 0.45)

out = FIG / "s6-host-coupling-GSE227919.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print(f"file: {out.name} {out.stat().st_size:,} bytes")
try:
    panel_crops(fig)
    print("crops: .cache/figure-style/")
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
