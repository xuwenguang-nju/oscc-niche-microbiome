#!/usr/bin/env python3
"""Fix S7 (plot_s6.py): long one-line titles collide with the panel letters.

In plot_s6.py the titles are single strings like "Host-program coupling (0/24 pass
FDR; min q = 0.09)". Centring them makes the block span most of the axes width, so the
standalone letter at the axes top-left lands *under* the title text. Splitting each
title into two centred lines narrows the block and lifts the letter clear, and the
shortened wording keeps both inside the axes.
"""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "analysis" / "modules" / "06-host-coupling" / "scripts" / "plot_s6.py"


def main() -> None:
    s = SRC.read_text(encoding="utf-8")
    lines = s.splitlines()
    for i, ln in enumerate(lines, 1):
        if "set_title" in ln:
            print(f"{i:4}| {ln.rstrip()}")
            for j in range(i, min(i + 3, len(lines) + 1)):
                print(f"    | {lines[j - 1].strip()}")
    print()
    print("figsize:", [l.strip() for l in lines if "figsize" in l])
    print("gridspec:", [l.strip() for l in lines if "add_gridspec" in l])
    print("adjust:", [l.strip() for l in lines if "subplots_adjust" in l])


if __name__ == "__main__":
    main()
