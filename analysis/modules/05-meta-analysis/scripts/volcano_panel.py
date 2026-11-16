#!/usr/bin/env python3
"""Volcano panel shared by the meta-analysis figure.

Why this is a function rather than a script
-------------------------------------------
The plot originally shipped as a standalone supplementary figure whose only panel was
this volcano, which left a supplement of one panel. It belongs inside the meta-analysis
figure: it shows the 146 testable genera and where the 13 that replicate sit among them,
which is the selection step the forest plot necessarily hides. Keeping it as a function
lets the main figure place it in a grid cell without duplicating the drawing code.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DOWN, UP, GREY = "#2166ac", "#b2182b", "0.72"
THR = -np.log10(0.05)


def draw_volcano(ax, letter: str | None = "A", panel_letter=None,
                 tab: Path | None = None, title: str | None = None) -> dict:
    """Draw the 615-genus meta-analysis volcano onto ``ax``.

    Marker size encodes how many cohorts tested each genus, so the eye separates a
    four-cohort result from a two-cohort one without a second legend.
    """
    if tab is None:
        tab = Path(__file__).resolve().parents[3] / "results" / "tables"
    m = pd.read_csv(tab / "s5_meta_analysis.csv")
    m = m[m.n_cohorts >= 2].copy()
    m["neglogq"] = -np.log10(m.pooled_q_bh.clip(lower=1e-6))

    # The pooled effect in this table is a natural-log ratio. Every other place it is quoted -
    # Table S1, Fig. 4L, Results - is in log2, so the axis is converted too; otherwise the panel
    # reads 0.76 for Parvimonas while the table says 1.09 and the axis claims log2.
    m["effect_l2"] = m.pooled_effect / np.log(2)

    rep = m[(m.n_cohorts >= 3) & (m.pooled_q_bh < 0.05) & m.direction_consistent]
    rest = m.drop(rep.index)

    # a label gutter on the left keeps leader-line text off the data
    xmin = m.effect_l2.min() - 0.95
    xmax = m.effect_l2.max() + 0.10

    for direction, frame in [("up", rep[rep.pooled_effect > 0]),
                             ("down", rep[rep.pooled_effect < 0])]:
        ax.scatter(frame.effect_l2, frame.neglogq,
                   s=8 + 11 * frame.n_cohorts,
                   color=UP if direction == "up" else DOWN,
                   edgecolor="white", lw=0.35, zorder=4,
                   label=(f"replicates, higher ({len(frame)})"
                          if direction == "up"
                          else f"replicates, lower ({len(frame)})"))

    for k, mk, lab in [(2, "o", "2 cohorts"), (3, "s", "3 cohorts"),
                       (4, "^", "4 cohorts")]:
        f = rest[rest.n_cohorts == k]
        if len(f):
            ax.scatter(f.effect_l2, f.neglogq, s=9, marker=mk, color=GREY,
                       edgecolor="none", alpha=0.8, zorder=2, label=lab)

    ax.axhline(THR, color="0.35", lw=0.7, ls=(0, (4, 2)), zorder=1)
    ax.axvline(0, color="0.55", lw=0.6, zorder=1)

    # Four genera sit within a few hundredths of each other on x, so the left-hand
    # labels are stacked into a gutter with leader lines; the one genus on the right
    # goes just above its point, since a right-hand gutter would run off the axis.
    # Twelve of the thirteen fall on the same side, and several sit within 0.1 of each other on q,
    # so a fixed 0.30 pitch cannot fit them: it either runs past the axis or overlaps. The labels
    # are therefore spread evenly over the span their own q values occupy, which keeps the reading
    # order (highest q at the top) while guaranteeing no two labels collide.
    neg = rep[rep.effect_l2 < 0].sort_values("neglogq", ascending=False)
    _n = len(neg)
    _pitch = 0.235
    # The stack has to satisfy three things at once: keep the top label at or below the highest
    # point, stay above the axis, and never have two labels closer than one text height. If the
    # points are too tightly bunched for all three, the pitch is reduced rather than the labels
    # being allowed off the axis.
    _lo_limit = 0.06
    _y0 = float(neg.neglogq.max())
    if _y0 - (_n - 1) * _pitch < _lo_limit:
        _y0 = _lo_limit + (_n - 1) * _pitch
        _pitch = (_y0 - _lo_limit) / max(_n - 1, 1)
    for k, (_, row) in enumerate(neg.iterrows()):
        y = _y0 - k * _pitch
        ax.annotate(row.genus, xy=(row.effect_l2, row.neglogq),
                    xytext=(xmin * 0.985, y), fontsize=5.0, color="0.12", zorder=6,
                    ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.35, color="0.55",
                                    shrinkA=0, shrinkB=1.5,
                                    connectionstyle="arc3,rad=0.0"))
    for _, row in rep[rep.effect_l2 > 0].iterrows():
        ax.annotate(row.genus, xy=(row.effect_l2, row.neglogq),
                    textcoords="offset points", xytext=(0, 7),
                    fontsize=5.4, color="0.12", zorder=6,
                    ha="center", va="bottom", fontweight="bold")

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.15, m.neglogq.max() + 0.45)
    ax.text(xmax, THR + 0.05, "FDR 0.05 ", ha="right", va="bottom",
            fontsize=5.5, color="0.35")

    ax.set_xlabel("pooled log$_2$ ratio (tumour / normal)")
    ax.set_ylabel("$-\\log_{10}$ FDR q")
    if title is None:
        title = f"{len(m)} genera tested in \u22652 cohorts ({len(rep)} replicate)"
    ax.set_title(title, fontsize=5.6, loc="center", pad=4)
    # The right half holds only Parvimonas (top) and a low grey cloud, so the
    # legend goes low on the right rather than beside the replicating points.
    # This panel is drawn narrow in the merged figure, where three legend groups plus
    # eleven leader-line labels do not fit inside the axes. The legend goes OUTSIDE,
    # under the panel, so it cannot cover the data or the labels.
    ax.legend(frameon=False, fontsize=3.9, loc="upper center",
              bbox_to_anchor=(0.5, -0.30), ncol=3,
              handletextpad=0.22, labelspacing=0.15, columnspacing=0.6,
              scatterpoints=1, markerscale=0.7, borderaxespad=0.0)

    if letter is not None:
        # Route through the fixed helper so this panel's letter lines up with the rest of
        # the left-hand column; the caller may still pass its own function, which we honour
        # only if it is the fixed one.
        from figlib import panel_letter_fixed
        panel_letter_fixed(ax, letter, case="upper")

    return {"plotted": len(m), "replicating": len(rep),
            "q_below_005": int((m.pooled_q_bh < 0.05).sum())}
