#!/usr/bin/env python3
"""Print the panel letters and panel titles each figure actually renders.

Reads the plotting scripts (not the PNGs), so the letters and titles come from the same code
that produced the deliverable. Figure.savefig is disabled for the duration, otherwise the
scripts' own savefig calls overwrite the sized deliverables.

Use: compare this output with results/reports/figure-captions-jom.md and with any table or
checklist that points at a panel letter.
"""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis/modules/14-figure-layout/scripts"))
sys.path.insert(0, str(ROOT / "analysis/modules/05-meta-analysis/scripts"))

import matplotlib                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.figure              # noqa: E402
import matplotlib.pyplot as plt       # noqa: E402
from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: E402

import figlib                         # noqa: E402

SCRIPTS = {
    "Fig 1": "analysis/modules/09-figures/scripts/fig1_study_design.py",
    "Fig 2": "analysis/modules/03-simulation/scripts/plot_simulation_calibration.py",
    "Fig 3": "analysis/modules/07-six-niche/scripts/s7_plot.py",
    "Fig 4": "analysis/modules/08-replication/scripts/s8_plot.py",
    "Fig 5": None,   # re-reads its own saved PNG, so it cannot run with savefig disabled
    "Fig 6": "analysis/modules/11-clinical-survival/scripts/s11_plot_platform.py",
    "Fig 7": "analysis/modules/12-pancancer/scripts/s12_plot.py",
    "Fig S1": "analysis/modules/13-function/scripts/e5_plot.py",
    "Fig S2": "analysis/modules/05-meta-analysis/scripts/plot_forest.py",
    "Fig S3": "analysis/modules/11-clinical-survival/scripts/s11_plot.py",
    "Fig S4": "analysis/modules/06-host-coupling/scripts/s_combined_host_site_platform.py",
}


def panels(fig):
    """[(letter, title)] for the panels of one figure, ordered by letter."""
    if not isinstance(fig.canvas, FigureCanvasAgg):
        fig.set_canvas(FigureCanvasAgg(fig))
    fig.canvas.draw()
    out = []
    for ax in fig.axes:
        letter = None
        for t in ax.texts:
            s = t.get_text()
            if len(s) == 1 and s.isalpha() and t.get_fontweight() in ("bold", 700, "700"):
                letter = s.upper()
        # Several scripts set the panel title with loc="left" (or centre), so read all three.
        title = " ".join(x for x in (ax.get_title(), ax.get_title(loc="left"),
                                     ax.get_title(loc="right")) if x)
        if letter:
            out.append((letter, " ".join(title.split())))
    out.sort()
    return out


def main() -> None:
    real = matplotlib.figure.Figure.savefig
    matplotlib.figure.Figure.savefig = lambda *a, **k: None
    try:
        for label, rel in SCRIPTS.items():
            if rel is None:
                print(f"\n## {label}  skipped (script re-reads its own output; read the PNG)")
                continue
            plt.close("all")
            try:
                figs = figlib._run_capturing(ROOT / rel, 1.0)
            except Exception as exc:                              # noqa: BLE001
                print(f"\n## {label}  FAILED: {type(exc).__name__}: {exc}")
                continue
            found = []
            for f in figs:
                found.extend(panels(f))
            for f in figs:
                plt.close(f)
            letters = [p[0] for p in found]
            missing = [c for c in "ABCDEFGHIJKLMN" if c not in letters and
                       (not letters or c <= max(letters))]
            print(f"\n## {label}  ({len(found)} panels: {''.join(letters) or '-'})")
            if missing:
                print(f"   !! letters with no title/letter found: {''.join(missing)}")
            for letter, title in found:
                print(f"   {letter}  {title}")
    finally:
        matplotlib.figure.Figure.savefig = real


if __name__ == "__main__":
    main()
