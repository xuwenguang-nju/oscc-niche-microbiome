"""Summarise the Ion Torrent maxEE calibration into one table.

The calibration question: `maxEE = 2` is the Illumina default and is too strict
for Ion Torrent's higher per-base error. Raising it recovers reads, but does it
also manufacture diversity?

Answer recorded here: retention roughly doubles from maxEE 2 -> 4 while the ASV
count rises only ~9%. The recovered reads are therefore redundant copies of
sequences already detected, not new diversity. maxEE = 3 is adopted: it captures
most of the retention gain at a smaller risk of admitting error-derived ASVs.

Reads: analysis/out/dada2_calib/me*/PRJNA866676/run_status.json + qc.csv
Writes: results/tables/dada2_maxee_calibration.csv
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
CAL = ROOT / "analysis" / "out" / "dada2_calib"
TAB = ROOT / "results" / "tables"
# the maxEE = 2 run is the main pilot, which lives in results/dada2
BASELINE = ROOT / "results" / "dada2" / "PRJNA866676"

rows = []


def collect(label: str, study_dir: Path) -> None:
    st = json.loads((study_dir / "run_status.json").read_text(encoding="utf-8"))
    qc_path = study_dir / "qc.csv"
    ret, nonchim = [], []
    if qc_path.exists():
        with qc_path.open(encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                ret.append(int(r["reads_out_filter"]) / int(r["reads_in"]))
                nonchim.append(int(r["reads_nonchimera"]) / int(r["reads_in"]))
    rows.append({
        "maxEE": label,
        "n_samples": st["n_samples_kept"],
        "n_asv": st["n_asv_postchimera"],
        "median_filter_retention_pct": round(100 * sorted(ret)[len(ret) // 2], 1) if ret else "",
        "median_nonchimera_of_in_pct": round(100 * sorted(nonchim)[len(nonchim) // 2], 1) if nonchim else "",
        "elapsed_sec": st["elapsed_sec"],
    })


if BASELINE.exists():
    collect("2", BASELINE)
for me in ("3", "4"):
    d = CAL / f"me{me}" / "PRJNA866676"
    if d.exists():
        collect(me, d)

# derived: how much new diversity each step of maxEE buys
base_asv = next((r["n_asv"] for r in rows if r["maxEE"] == "2"), None)
for r in rows:
    r["asv_gain_vs_maxEE2_pct"] = (
        round(100 * (r["n_asv"] - base_asv) / base_asv, 1) if base_asv else "")

TAB.mkdir(parents=True, exist_ok=True)
out = TAB / "dada2_maxee_calibration.csv"
with out.open("w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print(f"{'maxEE':<6} {'samples':>7} {'ASV':>6} {'retention%':>11} "
      f"{'nonchim%':>9} {'ASV gain%':>10}")
for r in rows:
    print(f"{r['maxEE']:<6} {r['n_samples']:>7} {r['n_asv']:>6} "
          f"{r['median_filter_retention_pct']:>11} {r['median_nonchimera_of_in_pct']:>9} "
          f"{r['asv_gain_vs_maxEE2_pct']:>10}")
print(f"\nwrote {out.relative_to(ROOT)}")
