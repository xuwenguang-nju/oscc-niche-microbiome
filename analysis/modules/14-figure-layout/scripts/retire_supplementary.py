#!/usr/bin/env python3
"""Retire the supplementary figures whose content moved into the main figures.

After the merge:
  * S2 (ten genera x six niches) -> main F4 panel J  (heatmap)
  * S3 (patient-level attribution) -> superseded by F4-A/D and F5-B, and it was drawn
    from the pre-host-removal tables, so its numbers no longer matched anything else
  * S5 (site stratified)    -> F6 panels J, K  (panel B folded into F7's scope)
  * S7 (signal composition)  -> F6 panel L
  * S6 (tongue subsite)      -> kept, but merged with S4 into one host-coupling figure

What remains:
  * S1  functional potential (PICRUST2)          - standalone, single column
  * S2  host coupling / site / tongue / platform  - built from s6v2+s6v3+s6v5+plot_s6
"""

from __future__ import annotations

import pathlib
import shutil

ROOT = pathlib.Path(__file__).resolve().parents[4]
SUPP = ROOT / "figures" / "supplementary"
LEGACY = ROOT / "figures" / "legacy"

RETIRE = [
    "S2-s5-genera-by-niche.png",
    "S3-single-cohort-attribution.png",
    "S5-site-stratified.png",
    "S6-tongue-subsite.png",
    "S7-rnaseq-signal-composition.png",
]
KEEP = ["S1-function-niches.png"]


def main() -> None:
    LEGACY.mkdir(parents=True, exist_ok=True)
    for name in RETIRE:
        src = SUPP / name
        if not src.is_file():
            print(f"  skip (absent): {name}")
            continue
        dst = LEGACY / name
        if dst.exists():
            dst.unlink()
        shutil.move(str(src), str(dst))
        print(f"  retired -> legacy/: {name}")
    for name in KEEP:
        print(f"  kept: {name}")
    print("\nsupplementary now:")
    for p in sorted(SUPP.glob("*.png")):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
