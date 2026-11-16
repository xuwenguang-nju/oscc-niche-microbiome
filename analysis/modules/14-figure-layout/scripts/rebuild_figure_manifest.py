#!/usr/bin/env python3
"""Rebuild results/tables/figure_manifest.csv for the current (v6) figure set.

Why this exists: the manifest had been edited by hand through two renumberings and had drifted
from reality again — it still carried the previous five-figure file names (F1…F9, FigS1–S2) and
widths measured before the last re-layout, so any check that "panel letters match the manifest"
was checking the wrong set. It is derived data: regenerate it, do not patch it.

Sources of truth
----------------
geometry (px, dpi, mm) : measured from the delivered PNGs in figures/manuscript-v6
panel counts / letters : verified by rendering each figure from its own plotting script
                         (analysis/scripts/audit_caption_panels.py prints the letters and titles);
                         the counts below are those renders, not a guess from the code
titles                 : the title sentence of each caption in results/reports/figure-captions-jom.md

The full panel-by-panel legend text is *not* duplicated here. It lives in
results/reports/figure-captions-jom.md, which is what the journal receives.
"""
from __future__ import annotations

import csv
import pathlib

from PIL import Image

ROOT = pathlib.Path(__file__).resolve().parents[4]
OUT = ROOT / "results" / "tables" / "figure_manifest.csv"
MIRROR = ROOT / "submission" / "jom-2026" / "07-source-and-reproducibility" / "figure_manifest.csv"
FIGDIR = ROOT / "figures" / "manuscript-v6"

# (set, file, panels verified by render, title)
SPEC = [
    ("main", "Fig1-design-and-cohort-overlap.png", 6,
     "Seven OSCC datasets, what each can answer, and overlap with an existing pooled analysis"),
    ("main", "Fig2-identifiability.png", 9,
     "Simulation establishes which quantities can be estimated from these data, and which cannot"),
    ("main", "Fig3-six-niche-attribution.png", 9,
     "Adjacent tissue is the closest available source of the tumour community"),
    ("main", "Fig4-replication-and-meta-analysis.png", 13,
     "The source ranking replicates in a second cohort; the tumour-enriched signal resolves to "
     "one species"),
    ("main", "Fig5-hospital-paired-cohort.png", 14,
     "In an independent hospital sample the change is a reciprocal exchange, and it survives "
     "depth, rank and functional inference"),
    ("main", "Fig6-platform-concordance.png", 12,
     "A second measurement substrate supports the broad structure and not the leading taxa"),
    ("main", "Fig7-pan-cancer.png", 9,
     "Across cancer types the signal follows the anatomical site, not the tumour"),
    ("supplementary", "FigS1-niche-functional.png", 3,
     "Functional prediction and source distance in the six-niche cohort"),
    ("supplementary", "FigS2-cross-cohort-meta.png", 8,
     "Detail of the cross-cohort meta-analysis"),
    ("supplementary", "FigS3-clinical-negative.png", 11,
     "Clinical and survival analysis, all negative"),
    ("supplementary", "FigS4-tcma-host-coupling.png", 7,
     "Direction, site and platform detail behind Figure 6"),
]


def main() -> None:
    rows = []
    for setname, fname, n_panels, title in SPEC:
        png = FIGDIR / setname / fname
        if not png.is_file():
            raise SystemExit(f"missing deliverable: {png}")
        im = Image.open(png)
        dpi = float(im.info.get("dpi", (0, 0))[0] or 0)
        if not dpi:
            raise SystemExit(f"{fname} carries no dpi tag, so its physical size is unknown")
        letters = [chr(65 + i) for i in range(n_panels)]
        rows.append({
            "set": setname,
            "figure": fname,
            "panels": n_panels,
            "title": title,
            # full legend text is in results/reports/figure-captions-jom.md, not duplicated here
            "legend": "",
            "width_mm": f"{im.width / dpi * 25.4:.1f}",
            "height_mm": f"{im.height / dpi * 25.4:.1f}",
            "px": f"{im.width}x{im.height}",
            "panel_letters": ",".join(letters),
            "dpi": f"{dpi:.0f}",
            "tiff": "yes" if (FIGDIR / setname / fname).with_suffix(".tif").is_file() else "no",
        })

    for out in (OUT, MIRROR):
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    for r in rows:
        print(f"  {r['set'][:4]:5}{r['figure']:40}{r['panels']:>4}p  "
              f"{r['width_mm']:>6} x {r['height_mm']:>6} mm  {r['dpi']} dpi  {r['panel_letters']}")
    print(f"\nwrote {OUT.relative_to(ROOT)} and its package mirror ({len(rows)} figures)")


if __name__ == "__main__":
    main()
