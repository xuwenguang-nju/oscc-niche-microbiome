#!/usr/bin/env python3
"""Bring e5_plot.py (supplementary S1) in line with the main-figure conventions.

Three defects, all visible in the rendered S1:
  * axis labels were 6.4 / 6.6 / 6.3 - same role, three sizes
  * panel B's legend hung below the axes (bbox_to_anchor y=-0.26), stretching the
    canvas downward and leaving dead space
  * panel C's note hung below the axes at y=-0.34 for the same reason
"""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "analysis" / "modules" / "13-function" / "scripts" / "e5_plot.py"


def main() -> None:
    s = SRC.read_text(encoding="utf-8")

    edits: list[tuple[str, str]] = [
        (
            'ax.set_ylabel("% of placed ASVs (log scale)", fontsize=6.4)',
            'ax.set_ylabel("% of placed ASVs (log scale)", fontsize=6.6)',
        ),
        (
            "fontsize=6.3)\n    ax.set_title(\"vs external sources, adjacent normal is nearest\",",
            "fontsize=6.6)\n    ax.set_title(\"vs external sources, adjacent normal is nearest\",",
        ),
        (
            'ax.legend(frameon=False, fontsize=6.0, loc="upper center",\n'
            "              bbox_to_anchor=(0.5, -0.26), ncol=3, handletextpad=0.25,",
            'ax.legend(frameon=False, fontsize=5.8, loc="lower right",\n'
            "              handletextpad=0.25,",
        ),
    ]
    for old, new in edits:
        if old not in s:
            raise SystemExit(f"NOT FOUND: {old[:70]!r}")
        if s.count(old) != 1:
            raise SystemExit(f"AMBIGUOUS ({s.count(old)}x): {old[:70]!r}")
        s = s.replace(old, new)

    # panel C note: move from below the axes to inside, bottom-left
    start = s.index("    ax.text(0.0, -0.34,")
    end = s.index("linespacing=1.3)", start) + len("linespacing=1.3)")
    note_new = (
        "    ax.text(0.02, 0.03,\n"
        '            "excludes cancer surface; adjacent normal nearest 20/37\\n"\n'
        '            "vs saliva p = 0.023 (22/37)",\n'
        '            transform=ax.transAxes, ha="left", va="bottom", fontsize=5.5,\n'
        '            color="0.25", linespacing=1.3)'
    )
    s = s[:start] + note_new + s[end:]

    SRC.write_text(s, encoding="utf-8")
    ast.parse(s)

    # verify
    ok = SRC.read_text(encoding="utf-8")
    import re
    xl = sorted(set(re.findall(r"set_xlabel\([^)]*?fontsize=([0-9.]+)", ok, re.S)))
    yl = sorted(set(re.findall(r"set_ylabel\([^)]*?fontsize=([0-9.]+)", ok, re.S)))
    tl = sorted(set(re.findall(r"set_title\([^)]*?fontsize=([0-9.]+)", ok, re.S)))
    below = re.findall(r"bbox_to_anchor=\([^)]*-\d[^)]*\)", ok)
    neg_y = re.findall(r"ax\.text\([^)]*-\d\.\d\d", ok)
    print(f"xlabel sizes {xl} | ylabel sizes {yl} | title sizes {tl}")
    print(f"legends below axes: {below or 'none'} | notes at negative y: {neg_y or 'none'}")


if __name__ == "__main__":
    main()
