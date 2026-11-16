#!/usr/bin/env python3
"""Centre every panel title in the supplementary-figure scripts.

Root cause: `apply_figure_style()` installs left-flush titles by default (its docstring
says so: "left-flush regular-weight titles"). `multialignment="center"` only aligns the
*lines of a multi-line title* relative to each other -- it does not move the title block
itself. Without `loc="center"` the block stays flush left, so on a wide axes a short title
sits at the left edge while the panel letter sits at the top-left: they collide, and the
title looks left-aligned even though multialignment is centred.

Fix: give every set_title an explicit loc="center" (idempotent, and also overrides any
inherited "left").
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]

SCRIPTS = [
    "04-niche-analysis/scripts/plot_niche_attribution.py",
    "06-host-coupling/scripts/s6v2_plot.py",
    "06-host-coupling/scripts/s6v3_plot.py",
    "06-host-coupling/scripts/s6v5_plot.py",
    "06-host-coupling/scripts/plot_s6.py",
    "07-six-niche/scripts/s7b_plot_genera_by_niche.py",
    "13-function/scripts/e5_plot.py",
]


def add_center(src: str) -> tuple[str, int]:
    """Insert loc="center" into every set_title call that lacks a loc."""
    tree = ast.parse(src)
    offsets: list[int] = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "set_title"):
            if not any(k.arg == "loc" for k in node.keywords):
                offsets.append(node.end_col_offset + node.col_offset * 0 + 0
                               if False else node)
    # work on text: find each call's closing paren and insert before it
    edits: list[tuple[int, str]] = []
    lines = src.splitlines(keepends=True)
    # map absolute char offsets
    starts, pos = [], 0
    for ln in lines:
        starts.append(pos)
        pos += len(ln)

    def abspos(lineno: int, col: int) -> int:
        return starts[lineno - 1] + col

    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "set_title"):
            if any(k.arg == "loc" for k in node.keywords):
                continue
            end = abspos(node.end_lineno, node.end_col_offset)
            # end points just past ")"
            edits.append((end - 1, ', loc="center"'))

    for at, text in sorted(edits, reverse=True):
        src = src[:at] + text + src[at:]
    return src, len(edits)


def main() -> None:
    total = 0
    for rel in SCRIPTS:
        p = ROOT / "analysis" / "modules" / rel
        if not p.is_file():
            print(f"  MISSING {rel}")
            continue
        s = p.read_text(encoding="utf-8")
        out, n = add_center(s)
        if n:
            ast.parse(out)
            p.write_text(out, encoding="utf-8")
        total += n
        print(f"  {rel.split('/')[-1]:34} +{n} loc=center")
    print(f"\ntotal titles centred: {total}")

    # verify
    print("\nverification:")
    for rel in SCRIPTS:
        p = ROOT / "analysis" / "modules" / rel
        if not p.is_file():
            continue
        s = p.read_text(encoding="utf-8")
        t = ast.parse(s)
        bad = 0
        for node in ast.walk(t):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "set_title"):
                if not any(k.arg == "loc" for k in node.keywords):
                    bad += 1
        if bad:
            print(f"  !! {p.name}: {bad} titles still without loc")
    print("  (no lines above = every set_title has an explicit loc)")


if __name__ == "__main__":
    main()
