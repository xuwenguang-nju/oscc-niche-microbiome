#!/usr/bin/env python3
"""Fix S6 (s6v5_plot.py): titles collide with letters; tick labels overlap.

Defects in the rendered S6:
  * panels A and B carry long two-line titles; centred they span the axes width and
    the standalone letter ends up underneath the title text
  * panel B's three genus names collide on the x axis
  * panel C's four subsite names collide on the x axis
"""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "analysis" / "modules" / "06-host-coupling" / "scripts" / "s6v5_plot.py"

EDITS = [
    ('ax.set_title("Parvimonas: sign flips with grouping,\\nnever significant", '
     'multialignment="center", loc="center")',
     'ax.set_title("Parvimonas: sign flips,\\nnever significant",\n'
     '             multialignment="center", loc="center", fontsize=7.4)'),
    ('ax.set_title("Signal sharpens as the anatomy\\ngets more homogeneous",\n'
     '             multialignment="center", loc="center")',
     'ax.set_title("Signal sharpens with\\nhomogeneous anatomy",\n'
     '             multialignment="center", loc="center", fontsize=7.4)'),
    ('ax.set_title("Subsite sizes", multialignment="center", loc="center")',
     'ax.set_title("Subsite sizes", multialignment="center", loc="center", fontsize=7.4)'),
]


def main() -> None:
    s = SRC.read_text(encoding="utf-8")
    for old, new in EDITS:
        if old not in s:
            print(f"  skip (not present as written): {old.splitlines()[0][:66]}")
            continue
        if s.count(old) != 1:
            raise SystemExit(f"AMBIGUOUS: {old[:60]!r}")
        s = s.replace(old, new)

    # rotate colliding x tick labels
    lines = s.splitlines()
    for i, ln in enumerate(lines):
        if "set_xticklabels" in ln and "rotation" not in ln:
            lines[i] = ln.rstrip().rstrip(")") + ', rotation=18, ha="right")'
    s = "\n".join(lines) + "\n"
    SRC.write_text(s, encoding="utf-8")
    ast.parse(s)

    for i, ln in enumerate(s.splitlines(), 1):
        if "set_xticklabels" in ln or "set_title" in ln:
            print(f"  {i:4}| {ln.strip()[:104]}")


if __name__ == "__main__":
    main()
