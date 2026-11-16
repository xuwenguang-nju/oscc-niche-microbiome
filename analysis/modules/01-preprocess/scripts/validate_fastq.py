"""Validate every downloaded FASTQ for integrity before it reaches DADA2.

Why this exists: a truncated or partially written .fastq.gz surfaces inside
DADA2 as the opaque `internal: buf != <newline>`, hours into a run, with no
indication of which file is at fault. It also happens silently after a
download retry. Checking first turns a 12-minute failure into a per-file report.

Two checks per file:
  1. gzip stream reads to the end without error (catches truncation)
  2. the decompressed line count is a multiple of 4 (catches a partial record)

Outputs
-------
results/tables/fastq_integrity.csv
analysis/out/fastq_integrity.txt

Usage
-----
python analysis/modules/01-preprocess/scripts/validate_fastq.py [--study NAME]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
RAW = ROOT / "data" / "raw" / "16s"
TAB = ROOT / "results" / "tables"
OUT = ROOT / "analysis" / "out"


def check(path: Path) -> dict:
    n_lines = 0
    n_records = 0
    err = ""
    try:
        with gzip.open(path, "rb") as fh:
            for _ in fh:
                n_lines += 1
        if n_lines % 4 != 0:
            err = f"line count {n_lines} not a multiple of 4"
        else:
            n_records = n_lines // 4
    except (EOFError, zlib.error, OSError) as e:
        err = f"{type(e).__name__}: {e}"
    return {
        "study": path.parent.name,
        "file": path.name,
        "bytes": path.stat().st_size,
        "reads": n_records,
        "status": "ok" if not err else "corrupt",
        "error": err,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default=None)
    ap.add_argument("--threads", type=int, default=6)
    args = ap.parse_args()

    studies = ([args.study] if args.study
               else sorted(p.name for p in RAW.iterdir() if p.is_dir()))
    files = [f for s in studies for f in sorted((RAW / s).glob("*.fastq.gz"))]
    print(f"checking {len(files)} files in {len(studies)} cohorts")

    rows = []
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futs = {pool.submit(check, f): f for f in files}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            rows.append(r)
            if r["status"] != "ok":
                print(f"  CORRUPT {r['study']}/{r['file']}: {r['error']}")
            if i % 100 == 0:
                print(f"  {i}/{len(files)}")

    rows.sort(key=lambda r: (r["study"], r["file"]))
    TAB.mkdir(parents=True, exist_ok=True)
    with (TAB / "fastq_integrity.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    bad = [r for r in rows if r["status"] != "ok"]
    lines = [f"checked {len(rows)} files; corrupt {len(bad)}"]
    for r in bad:
        lines.append(f"  {r['study']}/{r['file']}  bytes={r['bytes']}  {r['error']}")
    (OUT / "fastq_integrity.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nwrote {TAB / 'fastq_integrity.csv'}")


if __name__ == "__main__":
    main()
