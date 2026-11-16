#!/usr/bin/env python3
"""RETIRED — kept for reference only. Do not run.

This verified the v6 legend file (results/reports/figure-legends-v6.md) and the v6 Results draft,
both of which have since been superseded. It also called check_letter_alignment(), which re-renders
the figures at the renderer's default 330 dpi and overwrites the 600 dpi exports in
figures/manuscript-v6/ — the reason it is retired rather than left in the active folder.

What replaced it:
  analysis/scripts/verify_table_figure_numbering.py  numbering and first-citation order
  analysis/scripts/audit_citation_coverage.py        which panels each section actually cites
  analysis/scripts/verify_package_panels.py          captions against the package files
  analysis/scripts/export_figures.py                 the only supported way to write figure files
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]

FILES = {
    "Fig 1": "figures/manuscript-v6/main/Fig1-design-and-cohort-overlap.png",
    "Fig 2": "figures/manuscript-v6/main/Fig2-identifiability.png",
    "Fig 3": "figures/manuscript-v6/main/Fig3-six-niche-attribution.png",
    "Fig 4": "figures/manuscript-v6/main/Fig4-replication-and-meta-analysis.png",
    "Fig 5": "figures/manuscript-v6/main/Fig5-hospital-paired-cohort.png",
    "Fig 6": "figures/manuscript-v6/main/Fig6-platform-concordance.png",
    "Fig 7": "figures/manuscript-v6/main/Fig7-pan-cancer.png",
    "Fig S1": "figures/manuscript-v6/supplementary/FigS1-niche-functional.png",
    "Fig S2": "figures/manuscript-v6/supplementary/FigS2-cross-cohort-meta.png",
    "Fig S3": "figures/manuscript-v6/supplementary/FigS3-clinical-negative.png",
    "Fig S4": "figures/manuscript-v6/supplementary/FigS4-tcma-host-coupling.png",
}

PANEL_RE = re.compile(r"\(\*\*([A-U])\*\*\)")
CITE_RE = re.compile(r"Figs?\.\s*(S?\d+)")
SECTION_RE = re.compile(r"(?m)^### Figure .+$")


def check_letter_alignment() -> None:
    """Every figure's leftmost panel letter must sit at the same x, in inches."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    sys.path.insert(0, str(ROOT / "analysis/modules/14-figure-layout/scripts"))
    sys.path.insert(0, str(ROOT / "analysis/modules/05-meta-analysis/scripts"))
    import figlib

    scripts = {
        "Fig1": "analysis/modules/09-figures/scripts/fig1_study_design.py",
        "Fig2": "analysis/modules/03-simulation/scripts/plot_simulation_calibration.py",
        "Fig3": "analysis/modules/07-six-niche/scripts/s7_plot.py",
        "Fig4": "analysis/modules/08-replication/scripts/s8_plot.py",
        "Fig5": "analysis/scripts/fw2024_fig9_validation.py",
        "Fig6": "analysis/modules/11-clinical-survival/scripts/s11_plot_platform.py",
        "Fig7": "analysis/modules/12-pancancer/scripts/s12_plot.py",
        "S1": "analysis/modules/13-function/scripts/e5_plot.py",
        "S2": "analysis/modules/05-meta-analysis/scripts/plot_forest.py",
        "S3": "analysis/modules/11-clinical-survival/scripts/s11_plot.py",
        "S4": "analysis/modules/06-host-coupling/scripts/s_combined_host_site_platform.py",
    }
    xs = {}
    real_savefig = matplotlib.figure.Figure.savefig
    matplotlib.figure.Figure.savefig = lambda *a, **k: None    # never overwrite the deliverables
    try:
      for k, rel in scripts.items():
        plt.close("all")
        fig = figlib._run_capturing(ROOT / rel, 1.0)[0]
        if not isinstance(fig.canvas, FigureCanvasAgg):
            fig.set_canvas(FigureCanvasAgg(fig))
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        found = [(t.get_text(), t.get_window_extent(rend).x0 / fig.dpi)
                 for ax in fig.axes for t in ax.texts
                 if len(t.get_text()) == 1 and t.get_text().isupper()
                 and t.get_fontweight() in ("bold", 700, "700")]
        if found:
            xs[k] = min(found, key=lambda r: r[1])[1]
        plt.close(fig)
    finally:
        matplotlib.figure.Figure.savefig = real_savefig

    print("\nleftmost panel letter, inches from the figure edge:")
    for k, x in xs.items():
        print(f"  {k:6} {x:.3f}")
    if xs:
        spread = max(xs.values()) - min(xs.values())
        print(f"  spread {spread:.3f} in ->", "ALIGNED" if spread < 0.01 else "NOT ALIGNED")


def main() -> None:
    ok = True

    print("files:")
    for label, rel in FILES.items():
        present = (ROOT / rel).is_file()
        ok &= present
        print(f"  {'ok     ' if present else 'MISSING'} {label:7} {rel}")

    legend = (ROOT / "results/reports/figure-legends-v6.md").read_text(encoding="utf-8")
    heads = [(m.start(), m.group(0)) for m in SECTION_RE.finditer(legend)]
    print("\nlegend sections and the panels each describes:")
    for i, (start, head) in enumerate(heads):
        end = heads[i + 1][0] if i + 1 < len(heads) else len(legend)
        letters = "".join(dict.fromkeys(PANEL_RE.findall(legend[start:end])))
        print(f"  {head[:46]:48} {letters}")

    text = (ROOT / "results/reports/manuscript-results-v6.md").read_text(encoding="utf-8")
    body = text.split("## Results", 1)[1].split("## Notes for the authors")[0]
    order: list[str] = []
    for m in CITE_RE.finditer(body):
        if m.group(1) not in order:
            order.append(m.group(1))
    main = [x for x in order if not x.startswith("S")]
    supp = [x for x in order if x.startswith("S")]
    print("\nfirst-citation order:", " -> ".join(order))
    print("main figures  :", " -> ".join(main), " ascending:", main == ["1", "2", "3", "4", "5", "6", "7"])
    print("supplementary :", " -> ".join(supp), " ascending:", supp == ["S1", "S2", "S3", "S4"])

    print("\nALL FIGURE FILES PRESENT:", ok)
    check_letter_alignment()


if __name__ == "__main__":
    main()
