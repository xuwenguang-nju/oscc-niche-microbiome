"""Wait until every expected file for the given cohorts is on disk.

Why a dedicated waiter: the bulk ENA download runs at ~1.8 MB/s, so a cohort
takes tens of minutes. This environment cannot hold a long foreground poll, and
exit codes are not trustworthy (a killed native child can still leave 0), so
completion is decided by counting files that match the manifest and verifying
each one's gzip stream.

Usage
-----
python analysis/modules/01-preprocess/scripts/wait_for_fetch.py \
    --studies PRJNA1415268 PRJNA597251 --timeout 5400
"""

from __future__ import annotations

import argparse
import gzip
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
MANIFEST = ROOT / "results" / "tables" / "oscc_niche_sample_manifest.csv"
RAW = ROOT / "data" / "raw" / "16s"


def expected_files(studies: list[str]) -> dict[str, list[Path]]:
    m = pd.read_csv(MANIFEST)
    m = m[m.study.isin(studies)]
    out: dict[str, list[Path]] = {s: [] for s in studies}
    for _, row in m.iterrows():
        for p in str(row["fastq_ftp"]).split(";"):
            p = p.strip()
            if p:
                out[row["study"]].append(RAW / row["study"] / Path(p).name)
    return out


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
    ap.add_argument("--timeout", type=int, default=5400)
    ap.add_argument("--poll", type=int, default=60)
    ap.add_argument("--verify-gzip", action="store_true",
                    help="also decompress every file (slow but catches CRC damage)")
    args = ap.parse_args()

    exp = expected_files(args.studies)
    total = sum(len(v) for v in exp.values())
    print(f"expecting {total} files across {len(exp)} cohorts", flush=True)
    for s, v in exp.items():
        print(f"  {s}: {len(v)}", flush=True)

    t0 = time.time()
    last = -1
    while time.time() - t0 < args.timeout:
        present = {s: sum(1 for f in v if f.exists() and f.stat().st_size > 0)
                   for s, v in exp.items()}
        n = sum(present.values())
        if n != last:
            gb = sum(f.stat().st_size for v in exp.values() for f in v if f.exists()) / 1e9
            print(f"[{time.time() - t0:6.0f}s] {n}/{total} files, {gb:.2f} GB  "
                  + " ".join(f"{s}={present[s]}/{len(exp[s])}" for s in exp), flush=True)
            last = n
        if n == total:
            if args.verify_gzip:
                print("verifying gzip integrity ...", flush=True)
                bad = [f for v in exp.values() for f in v if not gzip_ok(f)]
                if bad:
                    print(f"CORRUPT: {[b.name for b in bad]}", flush=True)
                    return 1
                print("all files decompress cleanly", flush=True)
            print(f"\nCOMPLETE in {round(time.time() - t0)}s", flush=True)
            return 0
        time.sleep(args.poll)

    print(f"\nTIMEOUT after {args.timeout}s: {last}/{total} files present", flush=True)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
