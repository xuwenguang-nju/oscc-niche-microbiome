#!/usr/bin/env python3
"""Audit the supplementary figure scripts against the main-figure conventions.

Conventions (established on F1-F8):
  1. no figure-level title (suptitle) - titles live in the legend CSV
  2. panel letter standalone at the axis top-left, uppercase
  3. panel titles centred, with multialignment="center"
  4. every set_title / set_xlabel states its fontsize explicitly
  5. legends inside the axes (nothing hung below the axes)
  6. width <= 180 mm
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[4]
SUPP = ROOT / "figures" / "supplementary"

# figure -> source script, discovered from the savefig literals below
JOBS = {
    "S1-function-niches.png": "13-function/scripts/e5_plot.py",
    "S2-s5-genera-by-niche.png": "07-six-niche/scripts/s7_plot.py",
    "S3-single-cohort-attribution.png": "04-niche-analysis/scripts/plot_niche_attribution.py",
    "S4-host-coupling-panels.png": "06-host-coupling/scripts/s6v2_plot.py",
    "S5-site-stratified.png": "06-host-coupling/scripts/s6v3_plot.py",
    "S6-tongue-subsite.png": "06-host-coupling/scripts/s6v5_plot.py",
    "S7-rnaseq-signal-composition.png": "06-host-coupling/scripts/s6_signal_composition.py",
}


def top_margin(path: Path) -> tuple[int, int]:
    """Return (width_px, top_blank_px) for the rendered png."""
    a = np.array(Image.open(path).convert("L"))
    ink = (a < 250).any(axis=1)
    first = int(np.argmax(ink)) if ink.any() else a.shape[0]
    return a.shape[1], first


def audit(name: str, rel: str) -> dict:
    src_path = ROOT / "analysis" / "modules" / rel
    out = {"figure": name, "script": rel, "exists": src_path.is_file()}
    if not out["exists"]:
        return out

    s = src_path.read_text(encoding="utf-8")
    out["syntax_ok"] = True
    try:
        ast.parse(s)
    except SyntaxError as e:
        out["syntax_ok"] = False
        out["syntax_err"] = f"L{e.lineno}: {e.msg}"

    out["suptitle"] = s.count(".suptitle(")
    out["panel_letter"] = len(re.findall(r"panel_letter\(", s))
    out["title_left"] = len(re.findall(r'set_title\([^)]*loc="left"', s, re.S))
    out["title_multialign"] = s.count('multialignment="center"')
    out["title_total"] = len(re.findall(r"\.set_title\(", s)) + s.count("lbl(ax")

    # titles / xlabels without an explicit size fall back to the 8 pt rcParams base
    bare = []
    for m in re.finditer(r"\.set_(title|xlabel)\(((?:[^()]|\([^()]*\))*)\)", s, re.S):
        if "fontsize" not in m.group(2):
            bare.append(m.group(1))
    out["size_fallbacks"] = bare

    # legends hung below the axes push the canvas down
    out["legend_below"] = len(re.findall(r"bbox_to_anchor=\([^)]*-\d", s))

    png = SUPP / name
    if png.is_file():
        w, top = top_margin(png)
        out["width_mm"] = round(w / 330 * 25.4, 1)
        out["top_blank_px"] = top
        out["height_mm"] = round(Image.open(png).height / 330 * 25.4, 1)
    return out


def main() -> None:
    rows = [audit(n, r) for n, r in JOBS.items()]
    hdr = (f"{'supp figure':32} {'panels':>6} {'letter':>6} {'mlc':>4} {'sup':>4} "
           f"{'left':>4} {'bare':>4} {'leg↓':>4} {'mm':>6} {'top':>4}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        if not r.get("exists"):
            print(f"{r['figure']:32}  MISSING SCRIPT: {r['script']}")
            continue
        print(f"{r['figure']:32} {r['title_total']:>6} {r['panel_letter']:>6} "
              f"{r['title_multialign']:>4} {r['suptitle']:>4} {r['title_left']:>4} "
              f"{len(r['size_fallbacks']):>4} {r['legend_below']:>4} "
              f"{r.get('width_mm', 0):>6} {r.get('top_blank_px', 0):>4}")
    print("\nlegend: panels=titles found, letter=panel_letter calls, mlc=multialignment=center,")
    print("        sup=suptitle count, left=left-aligned titles, bare=titles/xlabels with no")
    print("        explicit fontsize, leg↓=legends anchored below the axes, top=top blank px")


if __name__ == "__main__":
    main()
