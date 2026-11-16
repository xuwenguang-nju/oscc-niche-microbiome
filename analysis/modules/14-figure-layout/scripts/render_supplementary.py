#!/usr/bin/env python3
"""Re-render all seven supplementary figures at their target widths.

Widths follow the main-figure rule: 180 mm double column, except S1 which is a
single-column vertical stack at 120 mm.
"""

from __future__ import annotations

import contextlib
import importlib
import io
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "analysis" / "modules" / "14-figure-layout" / "scripts"))
import figlib  # noqa: E402

JOBS = [
    ("13-function/scripts/e5_plot.py", "figures/supplementary/S1-function-niches.png", 120),
    ("07-six-niche/scripts/s7b_plot_genera_by_niche.py", "figures/supplementary/S2-s5-genera-by-niche.png", 180),
    ("04-niche-analysis/scripts/plot_niche_attribution.py", "figures/supplementary/S3-single-cohort-attribution.png", 180),
    ("06-host-coupling/scripts/s6v2_plot.py", "figures/supplementary/S4-host-coupling-panels.png", 180),
    ("06-host-coupling/scripts/s6v3_plot.py", "figures/supplementary/S5-site-stratified.png", 180),
    ("06-host-coupling/scripts/s6v5_plot.py", "figures/supplementary/S6-tongue-subsite.png", 180),
    ("06-host-coupling/scripts/plot_s6.py", "figures/supplementary/S7-rnaseq-signal-composition.png", 180),
]


def main() -> int:
    importlib.reload(figlib)
    ok = fail = 0
    for rel, out, mm in JOBS:
        src = ROOT / "analysis" / "modules" / rel
        if not src.is_file():
            print(f"  MISSING {rel}")
            fail += 1
            continue
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                figlib.render_at_width(str(src), str(ROOT / out), mm=mm)
            im = Image.open(ROOT / out)
            w, h = im.width / 330 * 25.4, im.height / 330 * 25.4
            print(f"  {Path(out).name:38} {w:6.1f} x {h:6.1f} mm")
            ok += 1
        except Exception as e:  # noqa: BLE001
            print(f"  FAIL {Path(out).name:38} {type(e).__name__}: {e}")
            fail += 1
        finally:
            plt.close("all")
    print(f"\n{ok} re-rendered, {fail} failed")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
