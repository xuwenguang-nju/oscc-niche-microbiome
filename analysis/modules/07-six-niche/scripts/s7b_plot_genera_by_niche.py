#!/usr/bin/env python3
"""S2 - the ten S5 genera across the six niches of PRJNA866676.

This figure previously had no script in the repository: it was rendered once from an
inline snippet, so it could not be re-rendered at a new width, could not follow the
panel-letter convention, and could not be audited. This file replaces that snippet and
reproduces the same content from `results/tables/s7c_s5_genera_by_niche.csv`.

Layout mirrors the other main/supplementary figures:
  * no figure-level title
  * standalone uppercase letter at the axes top-left
  * centred panel title (genus + S5 direction)
  * shared legend for the six niches, below the panels
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pathlib import Path

import sys

# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
RUNTIME = _FIGSTYLE / "runtime.py"
exec(compile(RUNTIME.read_text(encoding="utf-8"), str(RUNTIME), "exec"))


def project_root() -> Path:
    """Walk up until the project root (the directory holding figures/ and results/)."""
    env = None
    for cand in [Path(__file__).resolve(), *Path(__file__).resolve().parents]:
        if (cand / "results" / "tables").is_dir() and (cand / "figures").is_dir():
            return cand
    raise RuntimeError("project root not found")


ROOT = project_root()
TAB = ROOT / "results" / "tables"

NICHES = [
    ("tumor", "Tumour tissue", "#7b3294"),
    ("paracancerous", "Paracancerous", "#a6611a"),
    ("cancer_surface", "Cancer surface", "#d95f02"),
    ("contralateral_normal_mucosa", "Contralateral normal", "#1b7837"),
    ("saliva", "Saliva", "#2166ac"),
    ("tongue_coat", "Tongue coat", "#79c4e8"),
]

DIR = {"down_in_tumour": "\u2193", "up_in_tumour": "\u2191"}


def main() -> None:
    apply_figure_style(sizes=(7.0, 6.2, 5.4))
    df = pd.read_csv(TAB / "s7c_s5_genera_by_niche.csv")

    fig, axes = plt.subplots(2, 5, figsize=(7.09, 3.35))
    axes = axes.ravel()
    x = np.arange(len(NICHES))

    for k, (ax, (_, row)) in enumerate(zip(axes, df.iterrows())):
        vals = [float(row[c]) for c, _, _ in NICHES]
        ax.bar(x, vals, width=0.72,
               color=[col for _, _, col in NICHES],
               edgecolor="white", linewidth=0.4, zorder=2)
        ax.set_xticks([])
        ax.set_xlim(-0.6, len(NICHES) - 0.4)
        ax.set_ylim(0, max(vals) * 1.22 if max(vals) > 0 else 1)
        ax.tick_params(axis="x", length=0)
        ax.set_yticks(ax.get_yticks())
        if k % 5 == 0:
            ax.set_ylabel("rel. abund. (%)", fontsize=6.0)
        ax.set_title(f"{DIR.get(row['s5_direction'], '')} {row['genus']}",
                     fontsize=6.6, multialignment="center", pad=4,
                     color="#b2182b" if row["s5_direction"] == "up_in_tumour" else "#2166ac", loc="center")
        panel_letter(ax, chr(ord("A") + k), case="upper", dx=-0.14)

    handles = [plt.Rectangle((0, 0), 1, 1, fc=col, ec="none") for _, _, col in
               NICHES]
    fig.legend(handles, [lab for _, lab, _ in NICHES],
               loc="lower center", ncol=6, frameon=False, fontsize=6.0,
               handlelength=1.0, handletextpad=0.35, columnspacing=1.0,
               bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.055, right=0.995, top=0.93, bottom=0.13,
                        wspace=0.38, hspace=0.42)

    out = ROOT / "figures" / "supplementary" / "S2-s5-genera-by-niche.png"
    fig.savefig(out, dpi=330, facecolor="white")
    print(f"wrote {out}")
    try:
        panel_crops(fig)
    except Exception as e:  # noqa: BLE001
        print("panel_crops skipped:", e)


if __name__ == "__main__":
    main()
