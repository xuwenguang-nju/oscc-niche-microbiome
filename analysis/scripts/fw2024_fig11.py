#!/usr/bin/env python3
"""Draw Fig. 11: robustness across taxonomic ranks and whole-profile replication.

Cohort : FW2024-420, 20 paired OSCC tumour vs DISTANT normal mucosa (author-confirmed).
         The public comparator is the four-cohort random-effects meta-analysis.

A  number of taxa tested vs number passing FDR<0.05, at all six taxonomic ranks
B  phylum-level change, tumour minus normal (paired Wilcoxon, BH within rank)
C  every genus present in both the in-house cohort and the public meta-analysis:
   public pooled effect vs in-house log2 fold change

Input : results/tables/fw2024_420/rank_level_tests.csv
        results/tables/fw2024_420/whole_profile_vs_public.csv
Output: figures/main/F11-rank-and-profile-robustness.png      (176 mm, 330 dpi)
        figures/submission_tiff/F11-rank-and-profile-robustness.tif (600 dpi RGB LZW)

Run: python analysis/scripts/fw2024_fig11.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/tables/fw2024_420"
OUT_PNG = ROOT / "figures/main/F11-rank-and-profile-robustness.png"
OUT_TIF = ROOT / "figures/submission_tiff/F11-rank-and-profile-robustness.tif"
TARGET_MM = 176.0
RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def draw():
    rt = pd.read_csv(TAB / "rank_level_tests.csv", encoding="utf-8-sig")
    wp = pd.read_csv(TAB / "whole_profile_vs_public.csv", encoding="utf-8-sig")
    cnt = rt.groupby("rank").agg(tested=("taxon", "size"),
                                 fdr=("q", lambda s: (s < 0.05).sum())).reindex(RANKS)
    ph = rt[rt["rank"] == "Phylum"].copy()
    ph["delta_pct"] = (ph.mean_T - ph.mean_N) * 100
    PH = ph[(ph.mean_N > 0.005) | (ph.mean_T > 0.005)].sort_values("delta_pct")

    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.95))

    # A -- robustness across ranks
    ax = axes[0]
    y = np.arange(len(RANKS))
    ax.hlines(y, cnt.tested, cnt.fdr, color="#d5d5d5", lw=1.0, zorder=1)
    ax.scatter(cnt.tested, y, s=22, facecolor="white", edgecolor="#8a8a8a", lw=0.9,
               zorder=3, label="taxa tested")
    ax.scatter(cnt.fdr, y, s=26, color="C3", edgecolor="white", lw=0.4, zorder=4,
               label="FDR < 0.05")
    for i, (t_, f_) in enumerate(zip(cnt.tested, cnt.fdr)):
        ax.annotate(f"{int(f_)}", (f_, i), xytext=(-6, 0), textcoords="offset points",
                    ha="right", va="center", fontsize=6, color="C3")
    ax.set_yticks(y); ax.set_yticklabels(RANKS)
    ax.set_xscale("log"); ax.set_xlim(1.5, 900)
    ax.set_xlabel("number of taxa (log scale)")
    ax.set_title("Signal survives at every rank", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="lower right")

    # B -- phylum shift
    ax = axes[1]
    yy = np.arange(len(PH))
    cols = ["C3" if d > 0 else "C0" for d in PH.delta_pct]
    ax.barh(yy, PH.delta_pct, color=cols, height=0.62, zorder=2)
    ax.axvline(0, color="#555", lw=0.8, zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels([t.replace("_", " ") for t in PH.taxon], fontsize=6.0)
    for i, (d, q) in enumerate(zip(PH.delta_pct, PH.q)):
        star = "**" if q < 0.01 else ("*" if q < 0.05 else "")
        ax.annotate(f"{d:+.1f}{star}", (d, i), xytext=(5 if d > 0 else -5, 0),
                    textcoords="offset points", ha="left" if d > 0 else "right",
                    va="center", fontsize=5.8)
    ax.set_xlabel("$\\Delta$ relative abundance (points)")
    ax.set_title("Two phyla move in opposite directions", loc="left")
    ax.set_xlim(-19, 12)

    # C -- whole-profile replication
    ax = axes[2]
    hi = (wp.prev_N >= 0.5) & (wp.prev_T >= 0.5)
    ax.scatter(wp.pooled_effect[~hi], wp.log2FC[~hi], s=10, color="#b0b0b0",
               edgecolor="white", lw=0.3, zorder=2,
               label=f"prevalence < 50% (n={int((~hi).sum())})")
    ax.scatter(wp.pooled_effect[hi], wp.log2FC[hi], s=13, color="C0",
               edgecolor="white", lw=0.3, zorder=3,
               label=f"prevalent in both (n={int(hi.sum())})")
    lo = min(wp.pooled_effect.min(), wp.log2FC.min()) - 0.4
    hiL = max(wp.pooled_effect.max(), wp.log2FC.max()) + 0.4
    ax.plot([lo, hiL], [lo, hiL], ls=":", color="#444", lw=0.8, zorder=1)
    ax.axhline(0, color="#e2e2e2", lw=0.6, zorder=0)
    ax.axvline(0, color="#e2e2e2", lw=0.6, zorder=0)
    ax.set_xlim(lo, hiL); ax.set_ylim(lo, hiL)
    rho_a = stats.spearmanr(wp.pooled_effect, wp.log2FC).statistic
    rho_h = stats.spearmanr(wp.pooled_effect[hi], wp.log2FC[hi]).statistic
    ax.annotate(f"$\\rho$ = {rho_a:.2f} (n={len(wp)})", xy=(lo + 0.15, hiL * 0.92),
                fontsize=6, color="#666", va="top")
    ax.annotate(f"$\\rho$ = {rho_h:.2f} (n={int(hi.sum())})", xy=(lo + 0.15, hiL * 0.83),
                fontsize=6, color="C0", va="top")
    ax.set_xlabel("pooled effect, public meta-analysis")
    ax.set_ylabel("log$_2$FC here (tumour / normal)")
    ax.set_title("Direction agrees across 120 genera", loc="left")
    ax.legend(frameon=False, fontsize=5.8, loc="lower right", bbox_to_anchor=(1.02, -0.02))

    for a, L in zip(axes, "ABC"):
        a.text(-0.24, 1.18, L, transform=a.transAxes, fontsize=9, fontweight="bold",
               va="top", ha="left", clip_on=False)
    fig.tight_layout(pad=0.6)
    return fig


def main():
    fig = draw()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=330, bbox_inches="tight")
    with Image.open(OUT_PNG) as im:
        arr = im.convert("RGB").copy()
    d = round(arr.width / (TARGET_MM / 25.4))
    arr.save(OUT_PNG, format="PNG", dpi=(d, d))

    fig2 = draw()
    tmp = ROOT / ".cache/f11_600.png"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    fig2.savefig(tmp, dpi=600, bbox_inches="tight")
    with Image.open(tmp) as im:
        a2 = im.convert("RGB").copy()
    d2 = round(a2.width / (TARGET_MM / 25.4))
    OUT_TIF.parent.mkdir(parents=True, exist_ok=True)
    a2.save(OUT_TIF, format="TIFF", compression="tiff_lzw", dpi=(d2, d2))
    tmp.unlink(missing_ok=True)
    print(f"PNG {arr.size} {d} dpi = {arr.width/d*25.4:.1f} mm")
    print(f"TIF {a2.size} {d2} dpi = {a2.width/d2*25.4:.1f} mm {OUT_TIF.stat().st_size/1e6:.2f} MB")


if __name__ == "__main__":
    main()
