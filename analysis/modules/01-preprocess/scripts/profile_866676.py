#!/usr/bin/env python3
"""Profile PRJNA866676 reads to choose DADA2 parameters empirically.

Ion Torrent single-end amplicons vary in length run to run, so truncLen and maxEE must
be measured rather than assumed. This samples the first N reads of a spread of runs and
reports:
  - length distribution (median, IQR, mode)
  - per-position mean quality for the first 420 positions
  - cumulative bases retained at candidate truncLen values

Run: python3 profile_866676.py [--per-run 2000] [--runs 12]
"""
from __future__ import annotations

import argparse
import gzip
import os
import statistics
from pathlib import Path

import numpy as np

def locate_root() -> Path:
    """Search upward for the directory holding the target FASTQ, so the script runs
    unchanged in the project tree and on the compute host."""
    env = os.environ.get("TCMA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "data" / "raw" / "16s" / "PRJNA866676").is_dir():
            return cand
    return here.parents[3]


ROOT = locate_root()
RAW = ROOT / "data" / "raw" / "16s" / "PRJNA866676"


def read_fastq(path: Path, limit: int):
    """Yield (seq, qual) pairs from a gzipped FASTQ, up to `limit` records."""
    n = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        while n < limit:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            fh.readline()
            qual = fh.readline().strip()
            if not seq:
                break
            yield seq, qual
            n += 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-run", type=int, default=2000)
    ap.add_argument("--runs", type=int, default=12)
    args = ap.parse_args()

    files = sorted(RAW.glob("*.fastq.gz"))
    print(f"fastq files: {len(files)}")
    if not files:
        raise SystemExit(f"no fastq under {RAW}")
    step = max(1, len(files) // args.runs)
    picks = files[::step][:args.runs]

    lens, qsum, qn = [], None, 0
    for p in picks:
        L = []
        for seq, qual in read_fastq(p, args.per_run):
            L.append(len(seq))
            q = np.frombuffer(qual.encode("ascii", "replace"), dtype=np.uint8).astype(np.int16) - 33
            q = q[:420]
            if qsum is None:
                qsum = np.zeros(420)
            qsum[:len(q)] += q
            qn += 1
        lens.extend(L)
        print(f"  {p.stem}: n={len(L)} median_len={int(statistics.median(L))} "
              f"min={min(L)} max={max(L)}", flush=True)

    arr = np.array(lens)
    mean_q = qsum / max(qn, 1)
    print("\n=== length distribution (all sampled reads) ===")
    print(f"  n={arr.size}  median={int(np.median(arr))}  "
          f"IQR=[{int(np.percentile(arr,25))}, {int(np.percentile(arr,75))}]  "
          f"max={int(arr.max())}")
    vals, counts = np.unique(arr, return_counts=True)
    top = np.argsort(counts)[::-1][:6]
    print("  modal lengths:", [(int(vals[i]), int(counts[i])) for i in top])

    print("\n=== mean quality by position (every 20 bp) ===")
    for pos in range(0, 420, 20):
        print(f"  {pos:>3}: {mean_q[pos]:5.2f}")

    print("\n=== bases retained at candidate truncLen ===")
    for cut in (200, 220, 240, 260, 280, 300, 350, 400):
        keep = float((arr >= cut).mean())
        bkeep = float(np.minimum(arr, cut).sum() / arr.sum())
        print(f"  truncLen={cut}: reads kept={100*keep:5.1f}%  bases kept={100*bkeep:5.1f}%")

    print("\n=== last position with mean Q >= 20 / 25 ===")
    for thr in (20, 25, 30):
        ok = np.where(mean_q >= thr)[0]
        print(f"  Q>={thr}: up to position {ok.max() if ok.size else 'n/a'}")


if __name__ == "__main__":
    main()
