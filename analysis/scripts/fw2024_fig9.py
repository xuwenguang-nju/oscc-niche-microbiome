#!/usr/bin/env python3
"""Draw Fig. 9 from the down-sampling tables.

Kept separate from the analysis script so the figure can be re-styled without
re-running 210 subsampling replicates.

Input : results/tables/fw2024_420/downsampling_for_manuscript.csv
        results/derived/fw2024_420/genus_counts.csv
        results/tables/fw2024_420/downsampling_per_genus.csv
Cohort : FW2024-420, 20 paired OSCC tumour vs DISTANT normal mucosa (author-confirmed).
         No adjacent-mucosa sample exists, so this figure cannot address localisation.
Output: figures/main/F9-depth-sensitivity.png   (176 mm wide, 330 dpi)
        figures/submission_tiff/F9-depth-sensitivity.tif (600 dpi, RGB, LZW)

Run: python analysis/scripts/fw2024_fig9.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/tables/fw2024_420"
DER = ROOT / "results/derived/fw2024_420"
OUT_PNG = ROOT / "figures/main/F9-depth-sensitivity.png"
OUT_TIF = ROOT / "figures/submission_tiff/F9-depth-sensitivity.tif"

TARGET_MM = 180.0
FOC, MUT = "C0", "#b9b9b9"


def load():
    T = pd.read_csv(TAB / "downsampling_for_manuscript.csv", encoding="utf-8-sig").sort_values("depth")
    counts = pd.read_csv(DER / "genus_counts.csv", index_col=0, encoding="utf-8-sig")
    rel = counts.div(counts.sum(axis=0), axis=1)
    prev = (counts > 0).mean(axis=1)
    ab = rel.mean(axis=1) * 100
    PG = pd.read_csv(TAB / "downsampling_per_genus.csv", encoding="utf-8-sig")
    med = PG.groupby(["genus", "depth"])["p"].median().unstack()
    sig = med < 0.05
    keep = sig[20000]
    info = pd.DataFrame({"prev": prev * 100, "ab": ab})
    info["retained"] = keep
    info["retained_500"] = sig[500]
    return T, info


def draw():
    T, info = load()
    x = T.depth.values
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.0))

    ax = axes[0, 0]
    ax.plot(x, T.genera_detected, "-o", color=MUT, ms=3.2, lw=1.3, label="genera detected")
    ax.plot(x, T.nominal_p05, "-s", color=FOC, ms=3.2, lw=1.3, label="nominal p < 0.05")
    ax.plot(x, T.fdr05, "-^", color="C3", ms=3.2, lw=1.3, label="FDR < 0.05")
    ax.set_xscale("log")
    ax.set_xlabel("reads per sample (subsampled)")
    ax.set_ylabel("number of genera")
    ax.set_title("Discovery shrinks with depth", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="center left")

    ax = axes[0, 1]
    ax.plot(x, T.parvimonas_p, "-o", color=FOC, ms=3.2, lw=1.3, label="Parvimonas enrichment")
    ax.plot(x, T.anaerobe_score_p, "-s", color="C3", ms=3.2, lw=1.3, label="anaerobic score")
    ax.axhline(0.05, color="#444", lw=0.8, ls=":")
    ax.annotate("p = 0.05", xy=(x.max(), 0.05), xytext=(0, -9), textcoords="offset points",
                ha="right", va="top", fontsize=6, color="#444")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("reads per sample (subsampled)")
    ax.set_ylabel("paired Wilcoxon p")
    ax.set_title("Key signals are depth-robust", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="center left", bbox_to_anchor=(0.02, 0.55))

    ax = axes[1, 0]
    ax.plot(x, T.auc_prior_down, "-o", color="C3", ms=3.2, lw=1.3)
    ax.plot(x, T.auc_prior_up, "-s", color=FOC, ms=3.2, lw=1.3)
    ax.axhline(0.5, color="#444", lw=0.8, ls=":")
    ax.annotate("12 down genera, a priori (AUC 0.90)", xy=(x.max(), T.auc_prior_down.iloc[-1]),
                xytext=(-4, 5), textcoords="offset points", ha="right", va="bottom",
                fontsize=6, color="C3")
    ax.annotate("Parvimonas only, a priori (AUC 0.81)", xy=(x.max(), T.auc_prior_up.iloc[-1]),
                xytext=(-4, -9), textcoords="offset points", ha="right", va="top",
                fontsize=6, color=FOC)
    ax.annotate("chance", xy=(ax.get_xlim()[1], 0.5), xytext=(-4, 5),
                textcoords="offset points", ha="right", va="bottom", fontsize=6, color="#444")
    ax.set_xscale("log"); ax.set_ylim(0.45, 1.02)
    ax.set_xlabel("reads per sample (subsampled)")
    ax.set_ylabel("AUC (tumour vs normal)")
    ax.set_title("Consortium separates at any depth", loc="left")

    ax = axes[1, 1]
    lost = info[info.retained & ~info.retained_500]
    kept = info[info.retained & info.retained_500]
    ax.scatter(kept.prev, kept.ab, s=16, color=FOC, alpha=.85,
               label=f"retained at 500 (n={len(kept)})")
    ax.scatter(lost.prev, lost.ab, s=22, color="C3", marker="^",
               label=f"lost by 500 (n={len(lost)})")
    ax.set_yscale("log")
    ax.set_xlabel("prevalence (% of 40 samples)")
    ax.set_ylabel("mean relative abundance (%)")
    ax.set_title("Lost genera are the low-abundance ones", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="lower right")

    for a, L in zip(axes.ravel(), "ABCD"):
        a.text(-0.14, 1.10, L, transform=a.transAxes, fontsize=9, fontweight="bold",
               va="top", ha="left")
        a.margins(0.07)
    fig.tight_layout(pad=0.5)
    return fig


def main():
    fig = draw()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=330, bbox_inches="tight")
    with Image.open(OUT_PNG) as im:
        arr = im.convert("RGB").copy()
    dpi = round(arr.width / (TARGET_MM / 25.4))
    arr.save(OUT_PNG, format="PNG", dpi=(dpi, dpi))
    print(f"PNG {OUT_PNG.relative_to(ROOT)}  {arr.size}  {dpi} dpi "
          f"= {arr.width/dpi*25.4:.1f} mm")

    # submission TIFF: re-render the vectors at 600 dpi rather than upscaling
    fig2 = draw()
    tmp = ROOT / ".cache/f9_600.png"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    fig2.savefig(tmp, dpi=600, bbox_inches="tight")
    with Image.open(tmp) as im:
        a2 = im.convert("RGB").copy()
    dpi2 = round(a2.width / (TARGET_MM / 25.4))
    OUT_TIF.parent.mkdir(parents=True, exist_ok=True)
    a2.save(OUT_TIF, format="TIFF", compression="tiff_lzw", dpi=(dpi2, dpi2))
    tmp.unlink(missing_ok=True)
    print(f"TIF {OUT_TIF.relative_to(ROOT)}  {a2.size}  {dpi2} dpi "
          f"= {a2.width/dpi2*25.4:.1f} mm  {OUT_TIF.stat().st_size/1e6:.2f} MB")


if __name__ == "__main__":
    main()
