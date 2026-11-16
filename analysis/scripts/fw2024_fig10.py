#!/usr/bin/env python3
"""Draw Fig. 10 from the FW2024-420 cohort: three views of the same within-patient shift.

Cohort : 20 paired OSCC tumour vs DISTANT normal mucosa (author-confirmed).
         No adjacent-mucosa sample exists, so this figure tests tumour association only;
         it cannot address localisation.

A  paired PCoA (Bray-Curtis) with an arrow per patient, normal -> tumour
B  per-patient anaerobe : aerobe ratio (log10)
C  Fusobacterium nucleatum, species level, paired

Input : provider PCoA coordinates, results/tables/fw2024_420/consortium_ratio_per_patient.csv,
        results/tables/fw2024_420/fusobacterium_species_per_sample.csv
Output: figures/main/F10-primary-cohort-shift.png      (176 mm, 330 dpi)
        figures/submission_tiff/F10-primary-cohort-shift.tif (600 dpi RGB LZW)

Run: python analysis/scripts/fw2024_fig10.py
"""
from __future__ import annotations

import re
import zipfile
from io import StringIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/FW2024-420/results"
TAB = ROOT / "results/tables/fw2024_420"
OUT_PNG = ROOT / "figures/main/F10-primary-cohort-shift.png"
OUT_TIF = ROOT / "figures/submission_tiff/F10-primary-cohort-shift.tif"
TARGET_MM = 176.0
FOC_T, FOC_N = "C3", "C0"

pid = lambda s: re.sub(r"^[NTY]", "", s)
grp = lambda s: "T" if s[0] in "TY" else "N"


def rd_zip(z, member, **kw):
    with zipfile.ZipFile(RAW / z) as zf:
        b = zf.read(member)
    for enc in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            return pd.read_csv(StringIO(b.decode(enc)), sep="\t", **kw)
        except UnicodeDecodeError:
            continue
    raise RuntimeError(member)


def load():
    pcoa = rd_zip("05.beta_diversity.zip",
                  "05.beta_diversity/beta_PCoA/T-vs-N/T-vs-N.braycurtis.PCoA.xls")
    pcoa["pid"] = pcoa["sample"].map(pid)
    pcoa["grp"] = pcoa["sample"].map(grp)
    rat = pd.read_csv(TAB / "consortium_ratio_per_patient.csv",
                      encoding="utf-8-sig", index_col=0)
    fn = pd.read_csv(TAB / "fusobacterium_species_per_sample.csv", encoding="utf-8-sig")
    return pcoa, rat, fn


def draw():
    pcoa, rat, fn = load()
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.8))

    # A -- paired PCoA
    ax = axes[0]
    px = pcoa.pivot_table(index="pid", columns="grp", values="PCo1")
    py = pcoa.pivot_table(index="pid", columns="grp", values="PCo2")
    ax.scatter(px["N"], py["N"], s=15, color=FOC_N, edgecolor="white", lw=0.4,
               label="distant normal", zorder=3)
    for p_ in px.index:
        ax.annotate("", xy=(px.loc[p_, "T"], py.loc[p_, "T"]),
                    xytext=(px.loc[p_, "N"], py.loc[p_, "N"]),
                    arrowprops=dict(arrowstyle="-|>", color=FOC_T, lw=0.55,
                                    shrinkA=1.5, shrinkB=1.5, alpha=.9), zorder=2)
    ax.scatter(px["T"], py["T"], s=15, color=FOC_T, edgecolor="white", lw=0.4,
               label="tumour", zorder=4)
    ax.axvline(0, color="#dddddd", lw=0.7, zorder=0)
    ax.set_xlabel("PCoA1 (Bray–Curtis)")
    ax.set_ylabel("PCoA2")
    ax.set_title("Arrows run left to right, normal → tumour", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="lower right")

    # B -- ratio
    ax = axes[1]
    pr = rat.pivot_table(index="pid", columns="grp", values="logratio")
    for p_ in pr.index:
        ax.plot([0, 1], [pr.loc[p_, "N"], pr.loc[p_, "T"]],
                color="#c9c9c9", lw=0.7, zorder=1)
    ax.scatter(np.zeros(len(pr)), pr["N"], s=15, color=FOC_N, edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(np.ones(len(pr)), pr["T"], s=15, color=FOC_T, edgecolor="white", lw=0.4, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["distant\nnormal", "tumour"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("log$_{10}$ anaerobe : aerobe")
    ax.set_title("Ratio rises in %d/20 patients" % int((pr["T"] > pr["N"]).sum()), loc="left")

    # C -- F. nucleatum
    ax = axes[2]
    fv = fn.pivot_table(index="pid", columns="grp", values="fn")
    for p_ in fv.index:
        ax.plot([0, 1], [np.log10(fv.loc[p_, "N"] + 0.01), np.log10(fv.loc[p_, "T"] + 0.01)],
                color="#c9c9c9", lw=0.7, zorder=1)
    ax.scatter(np.zeros(len(fv)), np.log10(fv["N"] + 0.01), s=15, color=FOC_N,
               edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(np.ones(len(fv)), np.log10(fv["T"] + 0.01), s=15, color=FOC_T,
               edgecolor="white", lw=0.4, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["distant\nnormal", "tumour"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("log$_{10}$ $\\it{F.\\ nucleatum}$ (% of reads)")
    ax.set_title("$\\it{F.\\ nucleatum}$ up in %d/20, d = 0.84"
                 % int((fv["T"] > fv["N"]).sum()), loc="left")

    for a, L in zip(axes, "ABC"):
        # place the letter outside the axes, above the title's left edge
        a.text(-0.22, 1.22, L, transform=a.transAxes, fontsize=9,
               fontweight="bold", va="top", ha="left", clip_on=False)
        a.margins(0.08)
    fig.tight_layout(pad=0.5)
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
    tmp = ROOT / ".cache/f10_600.png"
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
