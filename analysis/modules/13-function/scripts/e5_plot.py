#!/usr/bin/env python3
"""Figure: functional potential of the six OSCC niches (PICRUSt2 / MetaCyc pathways).

Panel A  NSTI quality control. Two things on one axis: the per-ASV NSTI distribution and the
         share of community reads those ASVs carry. The point is that the ASVs PICRUSt2 is least
         sure about are also the rarest, so the prediction is not extrapolating from the bulk of
         the community.
Panel B  PCoA on pathway profiles, coloured by niche -- does function separate by niche at all?
Panel C  distance of each patient's tumour to each niche reference (leave-self-out), paired.
         This is the pathway-level counterpart of the S7 source-attribution result.

Reads e_*.csv / e2_qc_nsti.json only.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform

def locate_root() -> Path:
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results" / "tables").is_dir() and (cand / "analysis").is_dir():
            return cand
    return here.parents[4]


ROOT = locate_root()
T = ROOT / "results" / "tables"
FIG = ROOT / "figures"

import sys
# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
sys.path.insert(0, str(_FIGSTYLE.resolve()))
from runtime import apply_figure_style, panel_crops  # noqa: E402
from figlib import panel_letter_fixed as panel_letter  # noqa: E402

NICHE_ORDER = ["tumor", "cancer_surface", "paracancerous", "saliva",
               "tongue_coat", "contralateral_normal_mucosa"]
NICHE_LAB = {"tumor": "tumour", "cancer_surface": "cancer\nsurface",
             "paracancerous": "adjacent\nnormal", "saliva": "saliva",
             "tongue_coat": "tongue\ncoat",
             "contralateral_normal_mucosa": "contra-\nlateral"}
# the y-tick labels carry deliberate hyphen breaks for a narrow axis; the legend needs the plain
# names, because a space-joined hyphen ("contra- lateral") reads as two words
LEGEND_LAB = {"tumor": "tumour", "cancer_surface": "cancer surface",
              "paracancerous": "adjacent normal", "saliva": "saliva",
              "tongue_coat": "tongue coat",
              "contralateral_normal_mucosa": "contralateral mucosa"}
COL = {"tumor": "#b2182b", "cancer_surface": "#e08214", "paracancerous": "#2166ac",
       "saliva": "#4a8fc0", "tongue_coat": "#7b9ab5",
       "contralateral_normal_mucosa": "#8c8c8c"}


def pcoa(dist: np.ndarray):
    n = dist.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (dist ** 2) @ J
    w, v = np.linalg.eigh(B)
    o = np.argsort(w)[::-1]
    w, v = w[o], v[:, o]
    pos = w > 0
    return v[:, pos] * np.sqrt(w[pos]), w[pos] / w[pos].sum()


def main() -> None:
    qc = json.loads((T / "e2_qc_nsti.json").read_text())
    D = pd.read_csv(T / "e3_tumour_distances.csv", index_col=0)
    rel = pd.read_csv(T / "e_pathway_relabund.csv", index_col=0)
    man = pd.read_csv(T / "oscc_niche_sample_manifest.csv", encoding="utf-8-sig")
    man = man[man.study == "PRJNA866676"].drop_duplicates("run").set_index("run")
    common = [s for s in rel.columns if s in man.index]

    apply_figure_style(sizes=(8.4, 8.0, 7.0))
    # Three stacked panels made the figure 281 mm tall at 172 mm wide, above the ~240 mm A4 text
    # block and taller than every other figure in the set. Panels A and B are both short, so they
    # share one row and C keeps the full width it needs for five labelled rows. The canvas width is
    # unchanged (4.72 in) so the export scale, and therefore every font size, stays as it was.
    fig = plt.figure(figsize=(4.72, 4.48))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.94, 1.10], width_ratios=[1.0, 1.45],
                          hspace=0.30, wspace=0.42,
                          left=0.155, right=0.985, top=0.945, bottom=0.105)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]

    # ---- a NSTI QC ----
    # A 99.9 / 0.1 split is invisible on a linear bar chart, so plot on a log axis and put the
    # numbers in text. The message is "the ASVs PICRUSt2 is unsure about are also rare".
    ax = axes[0]
    med = qc["asv_nsti"]["median"]
    n50 = qc["asv_nsti"]["n_gt_0.5"]
    n_tot = qc["n_asv_placed"]
    vals = [100 * (n_tot - n50) / n_tot, 100 * n50 / n_tot]
    labs = ["NSTI\n≤ 0.5", "NSTI\n> 0.5"]
    cols = ["#2166ac", "#e08214"]
    ax.bar([0, 1], vals, 0.55, color=cols, edgecolor="0.25", lw=0.6)
    # A log axis was tried first and rejected: it makes 0.11% look comparable to 99.9%, which is
    # the opposite of the point. Linear keeps the second bar effectively invisible, and the value is
    # annotated so the bar is not mistaken for a rendering fault.
    ax.set_ylim(0, 136)
    for i, v in enumerate(vals):
        ax.text(i, v + 3, f"{v:.1f}%" if v >= 1 else f"{v:.2f}%", ha="center",
                fontsize=6.0, color="0.25")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labs, fontsize=6.0)
    ax.set_ylabel("% of placed ASVs", fontsize=6.6)
    ax.set_title("NSTI: uncertain ASVs are rare", fontsize=8.2, multialignment="center", loc="center")
    panel_letter(ax, "A", case="upper")
    ax.text(0.5, 0.985,
            f"median NSTI {med:.3f}\n(NSTI > 0.5 carries "
            f"{qc['read_mass_pct_nsti_gt_0.5']:.2f}% of reads)",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.3, color="0.25",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))

    # ---- b PCoA on pathways ----
    ax = axes[1]
    X = rel[common].T.values
    bc = squareform(pdist(X, metric="braycurtis"))
    co, evr = pcoa(bc)
    niche = man.loc[common, "niche"].values
    for n in NICHE_ORDER:
        sel = niche == n
        if sel.sum() == 0:
            continue
        # The label carries a deliberate hyphen break for the y tick ("contra-\nlateral"); in the
        # legend that hyphen plus the space read as "contra- lateral", so the legend gets its own
        # spelling rather than a mechanical newline-to-space replacement.
        ax.scatter(co[sel, 0], co[sel, 1], s=6.5, alpha=0.6, linewidths=0,
                   color=COL[n], label=LEGEND_LAB[n], zorder=3)
    ax.axhline(0, color="0.9", lw=0.6, zorder=1)
    ax.axvline(0, color="0.9", lw=0.6, zorder=1)
    ax.set_xlabel(f"PCoA1 ({evr[0]*100:.0f}%)", fontsize=6.6)
    ax.set_ylabel(f"PCoA2 ({evr[1]*100:.0f}%)", fontsize=6.6)
    ax.set_title("Pathway profiles separate by niche", fontsize=8.2, multialignment="center", loc="center")
    panel_letter(ax, "B", case="upper")
    pm = json.loads((T / "e2_permanova.json").read_text())["niche_permanova_pathway"]
    # Region counts over the 313 points: top-right box 0, top-left 2, bottom-left 28,
    # bottom-right 32. The legend therefore takes the empty top-right corner and this
    # note takes the near-empty top-left one.
    ax.text(0.02, 0.985,
            f"PERMANOVA F={pm['F']:.1f}, R²={pm['R2']:.2f}, p={pm['p']:.3f}",
            transform=ax.transAxes, fontsize=5.0, va="top", ha="left", color="0.25",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    # every in-axes corner holds points (bottom-right: 36 of 313), so the legend goes
    # legend sits inside the panel, in the one corner that holds no points at all.
    # The x-axis title keeps matplotlib's default offset: forcing it to -0.10 put it on
    # the same line as the "0.0" tick label and the two overlapped.
    ax.legend(frameon=False, fontsize=5.0, loc="upper right",
              handletextpad=0.2, columnspacing=0.6, labelspacing=0.12,
              borderaxespad=0.15, markerscale=1.1)

    # ---- c paired distances ----
    ax = axes[2]
    order = [c for c in D.columns if c != "tumor"]
    order = sorted(order, key=lambda c: D[c].mean())   # nearest first
    pos = np.arange(len(order))[::-1]
    for i, n in zip(pos, order):
        v = D[n].dropna().values
        ax.boxplot([v], positions=[i], widths=0.55, orientation='horizontal', patch_artist=True,
                   showfliers=False,
                   boxprops=dict(facecolor="white", edgecolor=COL[n], lw=0.8),
                   medianprops=dict(color=COL[n], lw=1.3),
                   whiskerprops=dict(color="0.45", lw=0.7),
                   capprops=dict(color="0.45", lw=0.7))
        ax.plot(v, np.random.default_rng(0).normal(i, 0.055, len(v)),
                "o", ms=1.9, color=COL[n], alpha=0.5, zorder=3)
    ax.set_yticks(pos)
    ax.set_yticklabels([NICHE_LAB[n] for n in order], fontsize=5.9)
    ax.set_xlabel("Bray–Curtis distance from the patient's tumour",
                  fontsize=6.6)
    ax.set_title("vs external sources, adjacent normal is nearest",
                 fontsize=8.2, multialignment="center", loc="center")
    panel_letter(ax, "C", case="upper")
    # The two notes used to share one line and overprinted each other. Split them and give the
    # block room below the lowest box.
    # The lowest row sits at y = 0 and its jitter reaches about -0.06, so the notes are placed in
    # the band below that rather than at a fraction of the axes, which put them on the saliva points.
    ax.set_ylim(-1.05, len(order) - 0.35)
    # The second line used to read "against saliva p = 0.023 (22 of 37)", which cannot be
    # reproduced from any table: no test in this cohort has n = 37 with that p, and the paired
    # saliva comparison that does exist is the backflow gradient over the 36 patients with both
    # samples (28 of 36, p = 2.8e-05; Results 3, Fig. 3B). The figure now states that instead.
    ax.text(0.015, -0.34, "excludes cancer surface: adjacent normal nearest in 20 of 37",
            transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.5, color="0.25")
    ax.text(0.015, -0.72, "saliva nearer the surface than the deep tumour in 28 of 36 (p < 0.001)",
            transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.5, color="0.25")


    p = FIG / "manuscript-v6" / "supplementary" / "FigS1-niche-functional.png"
    fig.savefig(p, dpi=330, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p.relative_to(ROOT), f"{p.stat().st_size/1024:.0f} KB")
    print(f"PCoA variance explained: {np.round(evr[:2]*100,1)}")


if __name__ == "__main__":
    main()
