#!/usr/bin/env python3
"""Export every figure of the seven-figure set at a fixed physical width.

Why this exists: the per-figure scripts hard-code `savefig(..., dpi=330)` and their own
figsize, so the delivered file can drift in width. `figlib.render_at_width` fixes the width
but its dpi handling proved unreliable when called with a non-default dpi. This runner does
the one thing needed — iterate the figure scale until the ink width, plus margin, equals the
target, then save once at the requested dpi and crop deterministically.

Usage:
    python analysis/scripts/export_figures.py            # 330 dpi, review files
    python analysis/scripts/export_figures.py --dpi 600  # submission files
"""
from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "modules" / "14-figure-layout" / "scripts"))
SET = ROOT / "figures" / "manuscript-v6"
PAD_IN = 0.085
TOL_MM = 0.8

JOBS = {
    "main/Fig1-design-and-cohort-overlap": ("analysis/modules/09-figures/scripts/fig1_study_design.py", 172.3),
    "main/Fig2-identifiability": ("analysis/modules/03-simulation/scripts/plot_simulation_calibration.py", 171.6),
    "main/Fig3-six-niche-attribution": ("analysis/modules/07-six-niche/scripts/s7_plot.py", 172.6),
    "main/Fig4-replication-and-meta-analysis": ("analysis/modules/08-replication/scripts/s8_plot.py", 172.2),
    "main/Fig5-hospital-paired-cohort": ("analysis/scripts/fw2024_fig9_validation.py", 172.6),
    "main/Fig6-platform-concordance": ("analysis/modules/11-clinical-survival/scripts/s11_plot_platform.py", 171.7),
    "main/Fig7-pan-cancer": ("analysis/modules/12-pancancer/scripts/s12_plot.py", 172.5),
    "supplementary/FigS1-niche-functional": ("analysis/modules/13-function/scripts/e5_plot.py", 172.0),
    "supplementary/FigS2-cross-cohort-meta": ("analysis/modules/05-meta-analysis/scripts/plot_forest.py", 172.3),
    "supplementary/FigS3-clinical-negative": ("analysis/modules/11-clinical-survival/scripts/s11_plot.py", 171.8),
    "supplementary/FigS4-tcma-host-coupling": ("analysis/modules/06-host-coupling/scripts/s_combined_host_site_platform.py", 172.7),
}


def draw(script: Path, scale: float):
    """Run the script with a scaled figsize and return the figures it built."""
    held: list = []
    orig_fig, orig_save = plt.figure, matplotlib.figure.Figure.savefig

    def scaled_figure(*a, **k):
        fs = k.get("figsize")
        if fs:
            k["figsize"] = (fs[0] * scale, fs[1] * scale)
        return orig_fig(*a, **k)

    keep_writes = any(t in script.name for t in ("fw2024", "e5_plot", "s11_plot.py", "s_combined"))

    def swallow(self, *a, **k):
        held.append(self)                     # keep the figure
        if keep_writes:
            # these scripts read their own file back immediately after saving
            orig_save(self, *a, **k)

    plt.figure = scaled_figure
    matplotlib.figure.Figure.savefig = swallow
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            # run_name="__main__" is required: several scripts only draw inside
            # `if __name__ == "__main__":`
            import runpy
            runpy.run_path(str(script), run_name="__main__")
    finally:
        plt.figure, matplotlib.figure.Figure.savefig = orig_fig, orig_save
    if not held:
        # some scripts write through plt.savefig or read their own file back; in that case
        # keep whatever figure is still registered
        held = [plt.figure(n) for n in plt.get_fignums()]
    return held


def ink_crop(fig, dpi: int) -> Image.Image:
    tmp = ROOT / "analysis" / "out" / "_crop_tmp.png"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(tmp, dpi=dpi, facecolor="white")
    arr = np.asarray(Image.open(tmp).convert("RGB"))
    grey = arr.min(axis=2)
    rows, cols = np.where(grey < 250)
    pad = int(round(PAD_IN * dpi))
    y0, y1 = max(0, rows.min() - pad), min(arr.shape[0], rows.max() + 1 + pad)
    x0, x1 = max(0, cols.min() - pad), min(arr.shape[1], cols.max() + 1 + pad)
    tmp.unlink(missing_ok=True)
    return Image.fromarray(arr[y0:y1, x0:x1])


def export(slug: str, script: str, mm: float, dpi: int) -> str:
    path = ROOT / script
    if not path.is_file():
        return f"  SKIP {slug}: {script} not found"
    scale, best = 1.0, None
    for _ in range(8):
        figs = draw(path, scale)
        if not figs:
            return f"  FAIL {slug}: no figure object"
        fig = figs[-1]
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        if not isinstance(fig.canvas, FigureCanvasAgg):
            fig.set_canvas(FigureCanvasAgg(fig))
        fig.canvas.draw()
        w_in = (fig.get_tightbbox(fig.canvas.get_renderer()).width + 2 * PAD_IN)
        best = (fig, w_in)
        if abs(w_in * 25.4 - mm) <= TOL_MM:
            break
        plt.close(fig)
        scale *= (mm / 25.4) / w_in
    fig, w_in = best
    crop = ink_crop(fig, dpi)
    plt.close(fig)
    out_png = SET / f"{slug}.png"
    out_tif = SET / f"{slug}.tif"
    crop.save(out_png, dpi=(dpi, dpi))
    crop.save(out_tif, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
    got = crop.width / dpi * 25.4
    flag = "" if abs(got - mm) <= 1.5 else f"  <-- {got - mm:+.1f} mm"
    return (f"  {out_png.name:46}{crop.width:5}px {dpi}dpi {got:6.1f} mm"
            f"  {out_png.stat().st_size/1e6:5.2f} MB{flag}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=int, default=330)
    args = ap.parse_args()
    print(f"exporting {len(JOBS)} figures at {args.dpi} dpi")
    for slug, (script, mm) in JOBS.items():
        print(export(slug, script, mm, args.dpi))


if __name__ == "__main__":
    main()
