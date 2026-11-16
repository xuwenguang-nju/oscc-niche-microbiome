#!/usr/bin/env python3
"""Re-render every main and supplementary figure after the S->main merge."""

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
    # ---- v6 set (six main + three supplementary). Files are named Fig<N>-<slug> so the name
    # states the figure number of THIS set; the previous F<N>-... set is archived under
    # figures/legacy/v5-figure-set/.
    ("09-figures/scripts/fig1_study_design.py",
     "figures/manuscript-v6/main/Fig1-design-and-cohort-overlap.png", 180),
    # Fig 2  identifiability           <- 03-simulation/scripts/plot_simulation_calibration.py
    # Fig 3  six-niche attribution     <- 07-six-niche/scripts/s7_plot.py
    # Fig 4  paired effect + depletion <- 08-replication/scripts/s8_plot.py + plot_forest.py (merge)
    # Fig 5  hospital paired cohort    <- scripts/fw2024_fig9_validation.py (moved forward)
    # Fig 6  platform + pan-cancer     <- s11_plot_platform.py + s12_plot.py (merge)
    # Fig S1 clinical (negative)       <- 11-clinical-survival/scripts/s11_plot.py
    # Fig S2 niche functional potential<- 13-function/scripts/e5_plot.py
    # Fig S3 TCGA host coupling        <- s_combined_host_site_platform.py
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
            flag = "  <-- over page height" if h > 240 else ""
            print(f"  {Path(out).name:40} {w:6.1f} x {h:6.1f} mm{flag}")
            ok += 1
        except Exception as e:  # noqa: BLE001
            print(f"  FAIL {Path(out).name:40} {type(e).__name__}: {e}")
            fail += 1
        finally:
            plt.close("all")
    print(f"\n{ok} rendered, {fail} failed")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
