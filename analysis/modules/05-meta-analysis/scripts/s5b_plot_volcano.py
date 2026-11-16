#!/usr/bin/env python3
"""S3 - volcano plot of the whole 615-genus meta-analysis.

Why this replaces a 615-row typeset table
-----------------------------------------
A 615 x 15 table is twenty pages nobody reads; F5-A only shows the 13 genera that
replicate, which hides how they were selected from 146 testable genera. One volcano
panel carries both: every genus at once, the 13 highlighted, and the cohort count
encoded as marker size so the reader can see how much independent evidence sits
behind each point.

The full table is still shipped, as a data file rather than typeset text.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "analysis" / "modules" / "14-figure-layout" / "scripts"))
# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
sys.path.insert(0, str(_FIGSTYLE))

exec(compile(open(_FIGSTYLE / "runtime.py", encoding="utf-8").read(),
              "runtime.py", "exec"))

TAB = ROOT / "results" / "tables"
OUT = ROOT / "figures" / "supplementary" / "S3-all-genera-volcano.png"

DOWN, UP, GREY = "#2166ac", "#b2182b", "0.72"


def main() -> None:
    apply_figure_style(sizes=(7.0, 6.2, 5.4))
    m = pd.read_csv(TAB / "s5_meta_analysis.csv")
    m = m[m.n_cohorts >= 2].copy()
    m["neglogq"] = -np.log10(m.pooled_q_bh.clip(lower=1e-6))

    rep = m[(m.n_cohorts >= 3) & (m.pooled_q_bh < 0.05) & m.direction_consistent]
    rest = m.drop(rep.index)

    fig, ax = plt.subplots(figsize=(7.09, 4.2))
    fig.subplots_adjust(left=0.10, right=0.98, top=0.90, bottom=0.13)

    # A label gutter on the left keeps the leader-line text off the data. The left
    # labels need ~0.30 of x-range of clear space, the single right label ~0.12.
    XMIN = m.pooled_effect.min() - 0.34
    XMAX = m.pooled_effect.max() + 0.10

    for d, frame in [("up", rep[rep.pooled_effect > 0]), ("down", rep[rep.pooled_effect < 0])]:
        ax.scatter(frame.pooled_effect, frame.neglogq,
                   s=8 + 11 * frame.n_cohorts, color=UP if d == "up" else DOWN,
                   edgecolor="white", lw=0.35, zorder=4,
                   label=f"replicates, higher in tumour ({len(frame)})" if d == "up"
                   else f"replicates, lower in tumour ({len(frame)})")

    ns = rest.n_cohorts
    for k, mk, lab in [(2, "o", "2 cohorts"), (3, "s", "3 cohorts"), (4, "^", "4 cohorts")]:
        f = rest[ns == k]
        if len(f):
            ax.scatter(f.pooled_effect, f.neglogq, s=9, marker=mk, color=GREY,
                       edgecolor="none", alpha=0.8, zorder=2, label=lab)

    thr = -np.log10(0.05)
    ax.axhline(thr, color="0.35", lw=0.7, ls=(0, (4, 2)), zorder=1)
    ax.axvline(0, color="0.55", lw=0.6, zorder=1)
    ax.text(ax.get_xlim()[1], thr + 0.05, " FDR 0.05", ha="right", va="bottom",
            fontsize=5.5, color="0.35")

    # Label the replicating genera. Four of them sit within a few hundredths of each
    # other on the x axis, so a stacked label column with leader lines is the only way
    # to keep both the text and the point identifiable.
    placed = {1: [], -1: []}
    for _, r in rep.sort_values("neglogq", ascending=False).iterrows():
        side = 1 if r.pooled_effect > 0 else -1
        x_txt = XMAX * 0.985 if side > 0 else XMIN * 0.985
        y = r.neglogq
        for py in placed[side]:
            if abs(y - py) < 0.15:
                y = py - 0.15
        placed[side].append(y)
        # left labels sit in a gutter; the one right-hand label goes just above its
        # point, because the gutter on that side would run off the axis
        if side > 0:
            ax.annotate(r.genus, xy=(r.pooled_effect, r.neglogq),
                        textcoords="offset points", xytext=(0, 7),
                        fontsize=5.4, color="0.12", zorder=6,
                        ha="center", va="bottom", fontweight="bold")
        else:
            ax.annotate(r.genus, xy=(r.pooled_effect, r.neglogq), xytext=(x_txt, y),
                        fontsize=5.0, color="0.12", zorder=6,
                        ha="left", va="center",
                        arrowprops=dict(arrowstyle="-", lw=0.35, color="0.55",
                                        shrinkA=0, shrinkB=1.5,
                                        connectionstyle="arc3,rad=0.0"))

    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(-0.15, m.neglogq.max() + 0.45)

    ax.set_xlabel("pooled log$_2$ fold change (tumour / normal)")
    ax.set_ylabel("$-\\log_{10}$ FDR q")
    ax.set_title(f"{len(m)} genera tested in at least two cohorts "
                 f"({len(rep)} replicate)", fontsize=7.0, loc="center", pad=5)
    ax.legend(frameon=False, fontsize=5.2, loc="lower right", ncol=2,
              handletextpad=0.3, columnspacing=1.0, labelspacing=0.25,
              scatterpoints=1, markerscale=0.8, borderaxespad=0.6)
    panel_letter(ax, "A", case="upper")

    fig.savefig(OUT, dpi=330, facecolor="white")
    print(f"wrote {OUT.relative_to(ROOT)}")
    print(f"  plotted {len(m)} genera, {len(rep)} replicating, "
          f"q<0.05 overall: {(m.pooled_q_bh < 0.05).sum()}")


if __name__ == "__main__":
    main()
