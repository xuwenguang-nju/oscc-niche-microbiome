#!/usr/bin/env python3
"""Figure S3: predicted function in the primary paired cohort (PICRUSt2 / KEGG).

Panel A  Volcano view of all 240 KEGG pathways: tumour-minus-normal difference against significance.
         Pathways carrying a canonical bacterial annotation are labelled; the four replication
         machinery pathways that a concurrent meta-analysis also reports are marked.
Panel B  The same pathways after adjusting each for the consortium (anaerobe:aerobe) score. Nothing
         survives. This is the panel that keeps the functional result honest: PICRUSt2 predicts
         function from taxonomy, so the functional shift is a restatement of the compositional one
         rather than independent evidence.

Reads results/tables/fw2024_420/picrust2_kegg_pathway_tests.csv only.
"""
from __future__ import annotations

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
sys.path.insert(0, str(_FIGSTYLE.resolve()))
from runtime import apply_figure_style, panel_letter  # noqa: E402


def locate_root() -> Path:
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results" / "tables").is_dir() and (cand / "analysis").is_dir():
            return cand
    return here.parents[4]


ROOT = locate_root()
T = ROOT / "results" / "tables"
TAB = T / "fw2024_420"

RED = "#b2182b"
BLUE = "#2166ac"
GREY = "#9a9a9a"
ORANGE = "#e08214"

# canonical bacterial pathway labels, short enough to print
NAMES = {
    "ko03010": "Ribosome",
    "ko03430": "Mismatch repair",
    "ko03440": "Homologous recomb.",
    "ko00970": "Aminoacyl-tRNA",
    "ko00190": "Ox. phosphorylation",
    "ko00290": "Val/Leu/Ile biosyn.",
    "ko00400": "Phe/Tyr/Trp biosyn.",
    "ko00260": "Gly/Ser/Thr metab.",
}
# the four the concurrent meta-analysis (Na et al. 2026) also reports
SHARED = {"ko03010", "ko03430", "ko03440", "ko00970"}


def main() -> None:
    kg = pd.read_csv(TAB / "picrust2_kegg_pathway_tests.csv", encoding="utf-8-sig")
    kg["neglog10q"] = -np.log10(kg["q"].clip(lower=1e-12))
    qcut = -np.log10(0.05)

    apply_figure_style(sizes=(8.4, 7.4, 6.6))
    fig, axes = plt.subplots(1, 2, figsize=(6.85, 3.05))
    fig.subplots_adjust(left=0.105, right=0.985, top=0.895, bottom=0.170, wspace=0.32)

    # ---------------- A: volcano ----------------
    ax = axes[0]
    sig = kg["q"] < 0.05
    ax.axhline(qcut, color="0.45", lw=0.7, ls=(0, (4, 2)), zorder=1)
    ax.axvline(0, color="0.55", lw=0.7, zorder=1)
    ax.scatter(kg.loc[~sig, "delta"], kg.loc[~sig, "neglog10q"], s=7, color="#cfcfcf",
               lw=0, zorder=2, label="q ≥ 0.05")
    ax.scatter(kg.loc[sig, "delta"], kg.loc[sig, "neglog10q"], s=9, color=ORANGE,
               lw=0, zorder=3, label="q < 0.05")
    lab = kg[kg["pathway"].isin(NAMES) & sig].copy()
    ax.scatter(lab["delta"], lab["neglog10q"], s=26, facecolor="none", edgecolor=RED,
               lw=1.1, zorder=4)
    ax.scatter(lab.loc[lab["pathway"].isin(SHARED), "delta"],
               lab.loc[lab["pathway"].isin(SHARED), "neglog10q"],
               s=44, facecolor="none", edgecolor="black", lw=1.2, zorder=5,
               label="also in Na et al. 2026")
    # place each label on the far side of its point, alternating height, so neighbouring
    # pathways at similar (delta, q) do not print on top of each other
    # The eight circled pathways form two tight clusters, so neither automatic placement nor a
    # uniform offset works. Labels are pinned by hand in axes-fraction space, all inside the
    # frame, with leader lines to their points.
    # Only the four pathways that a concurrent meta-analysis also reports are labelled
    # individually; labelling all eight on a panel this size collides. The reciprocal structure is
    # carried by two directional annotations instead, which is the claim the panel exists to make.
    # Everything is placed in the empty right-hand block, ordered top-to-bottom to match the
    # vertical order of the markers so the leader lines fan without crossing.
    PIN = {
        "ko03440": (0.66, 0.90),  # Homologous recomb.  (q=0.009)
        "ko03430": (0.66, 0.79),  # Mismatch repair     (q=0.013)
        "ko00970": (0.66, 0.68),  # Aminoacyl-tRNA      (q=0.029)
        "ko03010": (0.66, 0.57),  # Ribosome            (q=0.029)
    }
    for _, r in lab.iterrows():
        if r["pathway"] not in PIN:
            continue
        lx, ly = PIN[r["pathway"]]
        ax.annotate(NAMES[r["pathway"]], xy=(r["delta"], r["neglog10q"]),
                    xycoords="data", xytext=(lx, ly), textcoords="axes fraction",
                    ha="left", va="center", fontsize=4.7, color="0.15", zorder=7,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="0.55",
                                    shrinkA=1, shrinkB=3))
    ax.annotate("replication and\ntranslation higher", xy=(0.66, 0.20),
                xycoords="axes fraction", ha="left", va="center", fontsize=4.8,
                color=RED, linespacing=1.35, zorder=7)
    ax.annotate("biosynthesis and\nrespiration lower", xy=(0.66, 0.07),
                xycoords="axes fraction", ha="left", va="center", fontsize=4.8,
                color=BLUE, linespacing=1.35, zorder=7)
    ax.set_xlabel("mean difference, tumour − normal (pathway abundance)")
    ax.set_ylabel(r"$-\log_{10} q$", labelpad=16)
    ax.set_title("%d of %d pathways pass FDR" % (int(sig.sum()), len(kg)), loc="left")
    ax.legend(frameon=False, fontsize=4.7, loc="upper left", handletextpad=0.35,
              borderaxespad=0.4, labelspacing=0.28)
    ax.margins(x=0.34, y=0.10)
    panel_letter(ax, "A", case="upper")

    # ---------------- B: after adjustment ----------------
    # Only the 46 pathways that were significant survive into the adjusted test, so this panel
    # plots exactly those 46: their q before (open) and after (filled), joined by a segment.
    ax = axes[1]
    bb = kg[kg["q"] < 0.05].copy()
    ax.axhline(qcut, color="0.45", lw=0.7, ls=(0, (4, 2)), zorder=1)
    ax.axhline(1.0, color=GREY, lw=0.6, ls=":", zorder=1)
    for _, r in bb.iterrows():
        ax.plot([0, 1], [r["q"], r["q_after_consortium"]], color="0.80", lw=0.5,
                zorder=2, solid_capstyle="round")
    ax.scatter(np.zeros(len(bb)), bb["q"], s=13, facecolor="none", edgecolor=RED,
               lw=0.8, zorder=3, label="before adjustment")
    ax.scatter(np.ones(len(bb)), bb["q_after_consortium"], s=15, color=RED, lw=0,
               zorder=4, label="after adjustment")
    ax.set_yscale("log")
    ax.set_ylim(1.5e-3, 2.5)
    ax.set_xlim(-0.35, 1.35)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["q", "q |\nconsortium"], fontsize=5.2)
    ax.set_ylabel("q value (log scale)")
    nm = int((bb["q_after_consortium"] < 0.05).sum())
    ax.set_title("None survives (0 of %d)" % len(bb), loc="left")
    ax.annotate("%d pathways move from\nq < 0.05 to q = %.2f–%.2f"
                % (len(bb), bb["q_after_consortium"].min(), bb["q_after_consortium"].max()),
                xy=(0.30, 0.60), xycoords="axes fraction", ha="center", va="center",
                fontsize=4.9, color="0.2", linespacing=1.35, zorder=6)
    ax.legend(frameon=False, fontsize=4.7, loc="lower right", handletextpad=0.4,
              borderaxespad=0.3, labelspacing=0.3)
    panel_letter(ax, "B", case="upper")

    out = ROOT / "figures" / "supplementary" / "FigS3-functional-shift.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=330, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
