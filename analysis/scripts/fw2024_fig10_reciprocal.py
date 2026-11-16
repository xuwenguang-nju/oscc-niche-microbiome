#!/usr/bin/env python3
"""Draw Fig. 10: the reciprocal anaerobic:aerobic shift across all six cohorts.

This is the cross-cohort counterpart to Fig. 9C. Genus-level replication is tested only
in the four cohorts whose normal is peri-tumoural (Fig. 5); the two cohorts with distant
normal (PRJNA412445 and FW2024-420) are analysed here on the composite axis instead, so
that the comparison is like with like.

Input : results/tables/s5_reciprocal_shift_per_cohort.csv
        results/tables/s5_reciprocal_shift_pooled.csv
Output: figures/main/F10-reciprocal-shift-six-cohorts.png
        figures/submission_tiff/F10-reciprocal-shift-six-cohorts.tif
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/tables"
OUT_PNG = ROOT / "figures/main/F10-reciprocal-shift-six-cohorts.png"
OUT_TIF = ROOT / "figures/submission_tiff/F10-reciprocal-shift-six-cohorts.tif"
TARGET_MM = 176.0
FOC_T, FOC_N = "C3", "C0"


def draw():
    R = pd.read_csv(TAB / "s5_reciprocal_shift_per_cohort.csv", encoding="utf-8-sig")
    P = pd.read_csv(TAB / "s5_reciprocal_shift_pooled.csv", encoding="utf-8-sig").iloc[0]
    R2 = R.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(R2))
    hi = float((R2.lr_mean + 1.96 * R2.lr_se).max())
    lo = float((R2.lr_mean - 1.96 * R2.lr_se).min())

    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.85))

    ax = axes[0]
    ax.axvline(0, color="#555", lw=0.8, zorder=1)
    ax.errorbar(R2.lr_mean, y, xerr=1.96 * R2.lr_se, fmt="none", ecolor="#666",
                elinewidth=0.9, capsize=2.0, zorder=3)
    for i, nt in enumerate(R2.normal_type):
        ax.plot(R2.lr_mean.iloc[i], i, "o", ms=4.6,
                color=FOC_T if not str(nt).startswith("peri") else FOC_N, zorder=4)
    ax.set_yticks(y); ax.set_yticklabels(R2.study, fontsize=6.2)
    ax.set_xlim(lo - 0.5, hi + 0.45); ax.set_ylim(-0.7, len(R2) - 0.3)
    ax.set_xlabel("$\\Delta\\log_{10}$ (anaerobe : aerobe), tumour $-$ normal", fontsize=6.4)
    ax.set_title("Every cohort moves the same way", loc="left")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4.6, color=FOC_N,
                              label="peri-tumoural normal"),
                       Line2D([], [], marker="o", ls="", ms=4.6, color=FOC_T,
                              label="distant normal")],
              frameon=False, fontsize=5.2, loc="lower right", handletextpad=0.4)

    ax = axes[1]
    ax.barh([0], [P.pooled_log10_shift], height=0.40, color="#444", zorder=3)
    ax.errorbar(P.pooled_log10_shift, 0, xerr=1.96 * P.se, fmt="none", ecolor="#111",
                elinewidth=1.1, capsize=3, zorder=4)
    ax.axvline(0, color="#555", lw=0.8, ls=":", zorder=1)
    ax.set_yticks([0]); ax.set_yticklabels(["pooled (k=%d)" % int(P.k)], fontsize=6.2)
    ax.set_xlim(-0.9, hi + 0.45); ax.set_ylim(-0.7, 1.15)
    ax.set_xlabel("$\\Delta\\log_{10}$ (anaerobe : aerobe)", fontsize=6.4)
    ax.set_title("%d of %d cohorts positive" % (int(P.cohorts_positive), int(P.cohorts_total)),
                 loc="left")
    ax.annotate("+%.2f (95%% CI %.2f to %.2f)\n$p$ = %.1e,  $I^2$ = %.0f%%\n"
                "$\\approx$%.1f-fold higher in tumour"
                % (P.pooled_log10_shift, P.ci_low, P.ci_high, P.p, P.I2_pct,
                   10 ** P.pooled_log10_shift),
                xy=(0.03, 0.97), xycoords="axes fraction", fontsize=5.8, va="top")

    for a, L in zip(axes, "AB"):
        a.text(-0.20, 1.16, L, transform=a.transAxes, fontsize=9, fontweight="bold",
               va="top", ha="left", clip_on=False)
        a.tick_params(labelsize=6)
    fig.tight_layout(pad=0.55)
    return fig


def main():
    fig = draw()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=330, bbox_inches="tight")
    with Image.open(OUT_PNG) as im:
        arr = im.convert("RGB").copy()
    dpi = round(arr.width / (TARGET_MM / 25.4))
    arr.save(OUT_PNG, format="PNG", dpi=(dpi, dpi))
    print("PNG %s %s %d dpi = %.1f mm" % (OUT_PNG.relative_to(ROOT), arr.size, dpi,
                                          arr.width / dpi * 25.4))
    fig2 = draw()
    tmp = ROOT / ".cache/f10_600.png"; tmp.parent.mkdir(parents=True, exist_ok=True)
    fig2.savefig(tmp, dpi=600, bbox_inches="tight")
    with Image.open(tmp) as im:
        a2 = im.convert("RGB").copy()
    d2 = round(a2.width / (TARGET_MM / 25.4))
    OUT_TIF.parent.mkdir(parents=True, exist_ok=True)
    a2.save(OUT_TIF, format="TIFF", compression="tiff_lzw", dpi=(d2, d2))
    tmp.unlink(missing_ok=True)
    print("TIF %s %s %d dpi = %.1f mm %.2f MB" % (OUT_TIF.relative_to(ROOT), a2.size, d2,
          a2.width / d2 * 25.4, OUT_TIF.stat().st_size / 1e6))


if __name__ == "__main__":
    main()
