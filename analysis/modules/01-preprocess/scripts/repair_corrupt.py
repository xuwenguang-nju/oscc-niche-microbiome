"""Delete files that fail gzip verification and re-fetch them.

Corrupt-but-right-sized files proved to be a recurring failure mode of the ENA
HTTPS mirror under load (7 seen so far across two cohorts). This wraps the
detect -> delete -> re-fetch cycle so it can be repeated until the cohort is
clean, instead of being handled ad hoc each time.

Usage
-----
python analysis/modules/01-preprocess/scripts/repair_corrupt.py --studies PRJNA1415268
"""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
MANIFEST = ROOT / "results" / "tables" / "oscc_niche_sample_manifest.csv"
RAW = ROOT / "data" / "raw" / "16s"
STATUS = ROOT / "data" / "processed" / "fetch_status.csv"
FETCH = ROOT / "analysis" / "modules" / "01-preprocess" / "scripts" / "fetch_reads.py"


def gzip_ok(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except Exception:
        return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", nargs="+", required=True)
    args = ap.parse_args()

    m = pd.read_csv(MANIFEST)
    m = m[m.study.isin(args.studies)]
    files = []
    for _, row in m.iterrows():
        for p in str(row["fastq_ftp"]).split(";"):
            p = p.strip()
            if p:
                files.append((row["study"], RAW / row["study"] / Path(p).name))

    print(f"checking {len(files)} files ...", flush=True)
    bad = [(s, f) for s, f in files if f.exists() and not gzip_ok(f)]
    if not bad:
        print("no corrupt files", flush=True)
        return 0

    print(f"corrupt: {len(bad)}", flush=True)
    for s, f in bad:
        print(f"  {s}/{f.name}  {f.stat().st_size / 1e6:.1f} MB", flush=True)
        f.unlink()

    # mark them failed so --retry-failed picks them up
    names = {f.name for _, f in bad}
    rows = []
    if STATUS.exists():
        for r in csv.DictReader(STATUS.open(encoding="utf-8-sig")):
            if r["file"] in names:
                r["status"] = "failed"
                r["message"] = "corrupt gzip (CRC), deleted for re-download"
            rows.append(r)
        with STATUS.open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=["study", "file", "status", "bytes", "message"])
            w.writeheader()
            w.writerows(rows)

    print("\nre-fetching ...", flush=True)
    cmd = [sys.executable, "-u", str(FETCH), "--retry-failed", "--threads", "4"]
    rc = subprocess.call(cmd)
    print(f"fetch exit={rc} (ignore: this env can return 0 on failure)", flush=True)

    still = [(s, f) for s, f in bad if not (f.exists() and gzip_ok(f))]
    if still:
        print(f"\nSTILL CORRUPT: {[f.name for _, f in still]}", flush=True)
        return 1
    print(f"\nrepaired all {len(bad)} files", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
