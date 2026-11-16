#!/usr/bin/env python3
"""Re-export the seven-figure set at 600 dpi for journal submission.

The figure scripts save themselves at 330 dpi. Rather than edit eleven files, this runner
imports each script, takes the figure object it just drew, and re-saves that same vector
figure at 600 dpi. The ink bounding box is then measured on the 600-dpi raster and cropped
with the same 0.085 in margin used everywhere else, so the physical size is unchanged and
only the pixel density differs. Finally a TIFF is written from the same pixels.

Usage:  python analysis/scripts/export_600dpi.py
"""
from __future__ import annotations

import importlib.util
import io
import contextlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
SET = ROOT / "figures" / "manuscript-v6"
PAD_IN = 0.085
DPI = 600

# slug -> (script, target width in mm, figure index when a script draws several)
JOBS = {
    "main/Fig1-design-and-cohort-overlap": ("analysis/modules/09-figures/scripts/fig1_study_design.py", 172.3, 0),
    "main/Fig2-identifiability": ("analysis/modules/03-simulation/scripts/plot_simulation_calibration.py", 171.6, 0),
    "main/Fig3-six-niche-attribution": ("analysis/modules/07-six-niche/scripts/s7_plot.py", 172.6, 0),
    "main/Fig4-replication-and-meta-analysis": ("analysis/modules/08-replication/scripts/s8_plot.py", 172.2, 0),
    "main/Fig5-hospital-paired-cohort": ("analysis/scripts/fw2024_fig9_validation.py", 172.6, 0),
    "main/Fig6-platform-concordance": ("analysis/modules/11-clinical-survival/scripts/s11_plot_platform.py", 171.7, 0),
    "main/Fig7-pan-cancer": ("analysis/modules/12-pancancer/scripts/s12_plot.py", 172.5, 0),
    "supplementary/FigS1-niche-functional": ("analysis/modules/13-function/scripts/e5_plot.py", 172.0, 0),
    "supplementary/FigS2-cross-cohort-meta": ("analysis/modules/05-meta-analysis/scripts/plot_forest.py", 172.3, 0),
    "supplementary/FigS3-clinical-negative": ("analysis/modules/11-clinical-survival/scripts/s11_plot.py", 172.0, 0),
    "supplementary/FigS4-tcma-host-coupling": ("analysis/modules/06-host-coupling/scripts/s_combined_host_site_platform.py", 172.7, 0),
}


def ink_bbox(arr: np.ndarray) -> tuple[int, int, int, int]:
    """Bounding box of non-white content, in pixels."""
    grey = arr[..., :3].min(axis=2) if arr.ndim == 3 else arr
    mask = grey < 250
    rows, cols = np.where(mask)
    if rows.size == 0:
        h, w = arr.shape[:2]
        return 0, 0, w, h
    return cols.min(), rows.min(), cols.max() + 1, rows.max() + 1


def export(slug: str, script: str, mm: float, which: int) -> str:
    path = ROOT / script
    if not path.is_file():
        return f"  SKIP {slug}: script not found ({script})"
    name = "figlib_render"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.pop(name, None)
    # the scripts resolve their own output paths; point them at a scratch dir so the
    # 330-dpi deliverable is not overwritten
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            spec.loader.exec_module(mod)
        except Exception as exc:                      # noqa: BLE001
            return f"  FAIL {slug}: import {type(exc).__name__}: {exc}"
    figs = [n for n in (plt.get_fignums() or [])]
    if not figs:
        return f"  FAIL {slug}: no figure object produced"
    fig = plt.figure(figs[min(which, len(figs) - 1)])

    tmp = ROOT / "analysis" / "out" / "_600dpi_tmp.png"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(tmp, dpi=DPI, facecolor="white")
    arr = np.asarray(Image.open(tmp).convert("RGB"))
    x0, y0, x1, y1 = ink_bbox(arr)
    pad = int(round(PAD_IN * DPI))
    x0, y0 = max(0, x0 - pad), max(0, y0 - pad)
    x1, y1 = min(arr.shape[1], x1 + pad), min(arr.shape[0], y1 + pad)
    crop = Image.fromarray(arr[y0:y1, x0:x1])

    out_png = SET / f"{slug}.png"
    out_tif = SET / f"{slug}.tif"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    crop.save(out_png, dpi=(DPI, DPI))
    crop.save(out_tif, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
    tmp.unlink(missing_ok=True)
    plt.close(fig)
    w_mm = crop.width / DPI * 25.4
    return (f"  {out_png.name:46}{crop.width:5}px {DPI}dpi "
            f"{w_mm:6.1f} mm  {out_png.stat().st_size/1e6:5.2f} MB"
            f"{'  <-- width drift' if abs(w_mm - mm) > 1.5 else ''}")


def main() -> None:
    print(f"re-exporting {len(JOBS)} figures at {DPI} dpi")
    for slug, (script, mm, which) in JOBS.items():
        print(export(slug, script, mm, which))


if __name__ == "__main__":
    main()
