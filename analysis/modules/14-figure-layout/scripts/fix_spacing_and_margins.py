#!/usr/bin/env python3
"""Four spacing/margin corrections requested after visual review.

1. F4  - row spacing (hspace) slightly too large -> tighten
2. F6-L - the provenance legend collided with the stacked bar -> give it real room
3. S2  - row spacing too large -> tighten
4. all - the white border around every delivered figure is too narrow -> widen

The border is controlled in one place (`figlib.pad_in`), so fixing it once applies to
every figure that goes through `render_at_width`.
"""

from __future__ import annotations

import ast
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[4]


def patch(rel: str, pairs: list[tuple[str, str]]) -> None:
    p = ROOT / "analysis" / "modules" / rel
    s = p.read_text(encoding="utf-8")
    for old, new in pairs:
        n = s.count(old)
        if n != 1:
            raise SystemExit(f"{rel}: expected 1 match, found {n}: {old[:70]!r}")
        s = s.replace(old, new)
    ast.parse(s)
    p.write_text(s, encoding="utf-8")
    print(f"  patched {pathlib.Path(rel).name}")


def main() -> None:
    # 1) F4: 2x5 grid was airy
    patch("08-replication/scripts/s8_plot.py",
          [("wspace=0.92, hspace=0.72", "wspace=0.92, hspace=0.56")])

    # 3) S2: 3x4 combined figure was airy
    patch("06-host-coupling/scripts/s_combined_host_site_platform.py",
          [("wspace=0.85, hspace=0.85", "wspace=0.85, hspace=0.60")])

    # 2) F6 panel L: legend sat on the bar; give the axes headroom and 3 columns
    patch("11-clinical-survival/scripts/s11_plot.py", [
        ('ax.set_ylim(-0.55, 0.55)', 'ax.set_ylim(-0.55, 1.45)'),
        ('''ax.legend([plt.Rectangle((0, 0), 1, 1, fc=COLC[b], ec="none") for b in shares.index],
          ["oral-typical", "environmental / reagent", "other"],
          frameon=False, fontsize=4.5, loc="upper center", ncol=1,
          bbox_to_anchor=(0.5, 0.86), handlelength=0.9, labelspacing=0.2)''',
         '''ax.legend([plt.Rectangle((0, 0), 1, 1, fc=COLC[b], ec="none") for b in shares.index],
          ["oral-typical", "background", "other"],
          frameon=False, fontsize=4.5, loc="upper left", ncol=3,
          bbox_to_anchor=(0.0, 1.02), handlelength=0.8, columnspacing=0.9)'''),
    ])

    # 4) wider white border on every delivered figure
    patch("14-figure-layout/scripts/figlib.py",
          [("pad_in: float = 0.045", "pad_in: float = 0.085")])

    print("\nall patches applied")


if __name__ == "__main__":
    main()
