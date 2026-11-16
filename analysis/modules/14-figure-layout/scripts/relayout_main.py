#!/usr/bin/env python3
"""Re-lay out the main figures to an exact 180 mm double-column width.

Strategy (user-approved "A primarily, B as exception"):
  - main figures  : 180 mm double column, font sizes untouched
  - one exception : F9 functional panel -> 120 mm single column, vertical stack (see e5_plot.py)

Nothing here recomputes statistics; each source script is re-executed and only its figure
geometry changes. Output goes to figures/main/ so the deliverable set is unambiguous.

Run: python3 relayout_main.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import project_root, render_at_width  # noqa: E402

ROOT = project_root(Path(__file__).resolve().parent)

MAIN = ROOT / "figures" / "main"
MOD = ROOT / "analysis" / "modules"

# (script, figure index within that script, output name, width mm, what it shows)
JOBS = [
    (MOD / "09-figures/scripts/fig1_study_design.py", 0,
     "F1-study-design.png", 180,
     "cohort x niche coverage, pan-cancer check, analysis flow, evidence map"),
    (MOD / "03-simulation/scripts/plot_simulation_calibration.py", 0,
     "F2-identifiability.png", 180,
     "why source fractions are not reported; taxa detectable vs sequencing depth"),
    (MOD / "07-six-niche/scripts/s7_plot.py", 0,
     "F3-source-attribution.png", 180,
     "niche ranking, backflow gradient, PERMANOVA (PRJNA866676)"),
    (MOD / "08-replication/scripts/s8_plot.py", 0,
     "F4-replication-and-species.png", 180,
     "replication in a second cohort, P. micra, genus-unassigned blind spot"),
    (MOD / "05-meta-analysis/scripts/plot_forest.py", 0,
     "F5-cross-cohort-meta.png", 180,
     "13 replicating genera, absence of large effects, literature novelty"),
    (MOD / "11-clinical-survival/scripts/s11_plot.py", 0,
     "F6-clinical-prognosis.png", 180,
     "stage, KM, stage-adjusted hazard ratios (TCGA-HNSC)"),
    (MOD / "11-clinical-survival/scripts/s11_plot.py", 1,
     "F7-platform-concordance.png", 180,
     "TCMA vs amplicon 16S genus-level concordance"),
    (MOD / "12-pancancer/scripts/s12_plot.py", 0,
     "F8-pan-cancer.png", 180,
     "HNSC among five TCGA projects: PCoA, oral-gut axis, Parvimonas"),
]


def main() -> None:
    MAIN.mkdir(parents=True, exist_ok=True)
    rows = []
    for script, idx, name, mm, desc in JOBS:
        print(f"\n[{name}]  {script.name}  (figure #{idx}, target {mm} mm)")
        if not script.is_file():
            print("    MISSING SCRIPT, skipped")
            rows.append({"figure": name, "status": "missing script"})
            continue
        try:
            r = render_at_width(script, MAIN / name, mm=mm, which=idx)
            r.update({"figure": name, "status": "ok", "caption": desc,
                      "source": str(script.relative_to(ROOT))})
            rows.append(r)
        except Exception as e:                                   # noqa: BLE001
            print("    FAILED:", type(e).__name__, e)
            rows.append({"figure": name, "status": f"failed: {type(e).__name__}: {e}"})

    out = ROOT / "results" / "tables" / "figure_manifest.csv"
    import csv
    keys = ["figure", "status", "px", "mm", "kb", "caption", "source"]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nmanifest: {out.relative_to(ROOT)}")
    ok = [r for r in rows if r.get("status") == "ok"]
    print(f"laid out {len(ok)}/{len(JOBS)} main figures")
    for r in rows:
        name = r.get("figure", "?")
        print(f"  {name:34} {r.get('status'):10} {r.get('px',''):>12} "
              f"{r.get('mm','')} mm")


if __name__ == "__main__":
    main()
