"""Measure actual ENA throughput using the real URL from the project manifest.

The bulk download stalled with .part files stuck at 0 MB, which is consistent
with server-side throttling but must be measured. The earlier attempt to
diagnose this constructed the ENA path by hand and got a 404, so this version
reads the exact `fastq_ftp` path that the manifest already recorded.

Usage
-----
python analysis/modules/01-preprocess/scripts/diagnose_download.py --study PRJNA1415268
"""

from __future__ import annotations

import argparse
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
MANIFEST = ROOT / "results" / "tables" / "oscc_niche_sample_manifest.csv"
UA = {"User-Agent": "wisp-science/1.0 (download diagnostic)"}
MAX_BYTES = 20 << 20


def probe(url: str, timeout: int) -> dict:
    t0 = time.time()
    first, got, err = None, 0, ""
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=timeout) as r:
            while got < MAX_BYTES:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                if first is None:
                    first = time.time() - t0
                got += len(chunk)
    except Exception as e:
        err = f"{type(e).__name__}: {str(e)[:60]}"
    el = time.time() - t0
    return {"mb": round(got / 1e6, 2), "ttfb_s": None if first is None else round(first, 2),
            "elapsed_s": round(el, 1),
            "mb_per_s": round(got / 1e6 / el, 3) if el > 0 and got else 0.0,
            "error": err}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default="PRJNA1415268")
    ap.add_argument("--n", type=int, default=2)
    args = ap.parse_args()

    m = pd.read_csv(MANIFEST)
    sub = m[m.study == args.study].head(args.n)
    if sub.empty:
        raise SystemExit(f"no manifest rows for {args.study}")

    print(f"throughput probe: {args.study}, {len(sub)} file(s), "
          f"stop after {MAX_BYTES >> 20} MB each\n")
    hdr = f"{'run':<14}{'file':<24}{'timeout':>8}{'MB':>7}{'MB/s':>8}{'ttfb':>7}{'s':>7}  error"
    print(hdr)
    print("-" * len(hdr))
    for _, row in sub.iterrows():
        paths = [p.strip() for p in str(row["fastq_ftp"]).split(";") if p.strip()]
        if not paths:
            continue
        url = "https://" + paths[0]
        for timeout in (30, 60):
            r = probe(url, timeout)
            print(f"{row['run']:<14}{Path(paths[0]).name:<24}{timeout:>8}"
                  f"{r['mb']:>7.2f}{r['mb_per_s']:>8.3f}"
                  f"{str(r['ttfb_s']):>7}{r['elapsed_s']:>7.1f}  {r['error'][:40]}")

    print("\nreading the result:")
    print("  MB/s > 1.0                  -> link fine; re-run the bulk fetch")
    print("  MB/s < 0.1, ttfb large      -> server-side throttling; slow down threads")
    print("  'timed out' with 0 MB       -> connection accepted but starved; "
          "back off and retry later")


if __name__ == "__main__":
    main()
