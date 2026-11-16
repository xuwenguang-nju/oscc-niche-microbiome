#!/usr/bin/env python3
"""Compress S4 (s6v2_plot.py) so the stacked 3-panel figure fits a journal page.

Two problems in the rendered S4:
  * total height 294 mm - taller than the ~240 mm usable page height
  * panel A's two-line title sits directly on top of the first bar, because the
    figure is tall and narrow and the title needs its own band

Fix: shorter canvas, tighter hspace, and give panel A explicit headroom via margins.
"""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "analysis" / "modules" / "06-host-coupling" / "scripts" / "s6v2_plot.py"


def main() -> None:
    s = SRC.read_text(encoding="utf-8")
    edits = [
        ("fig = plt.figure(figsize=(4.7, 8.4))",
         "fig = plt.figure(figsize=(4.7, 6.4))"),
        ("gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.62], hspace=0.62)",
         "gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.55], hspace=0.75)"),
        # panel A: reserve headroom so the 2-line title clears the top bar
        ("ax.margins(y=0.05)",
         "ax.margins(y=0.05)\nax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 0.9)"),
        # panel C legend was pinned to the far right corner; keep it inside the axes
        ('ax.legend(handles, ["Primary tumour", "Adjacent normal"], frameon=False, fontsize=5.9,',
         'ax.legend(handles, ["Primary tumour", "Adjacent normal"], frameon=False, fontsize=5.6,\n'
         '          loc="lower right",'),
    ]
    for old, new in edits:
        if old not in s:
            raise SystemExit(f"NOT FOUND: {old[:70]!r}")
        if s.count(old) != 1:
            raise SystemExit(f"AMBIGUOUS: {old[:70]!r}")
        s = s.replace(old, new)
    SRC.write_text(s, encoding="utf-8")
    ast.parse(s)
    print("s6v2_plot.py compressed: canvas 8.4 -> 6.4 in, hspace 0.62 -> 0.75, panel A headroom")


if __name__ == "__main__":
    main()
