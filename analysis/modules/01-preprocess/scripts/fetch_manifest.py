#!/usr/bin/env python3
"""Generic manifest-driven FASTQ fetcher (handles paired-end ';'-separated paths).

The per-cohort fetchers we already have each encode one study's quirks; this one takes a
study id from the shared manifest and handles:
  - fastq_ftp with or without an https:// scheme
  - single path, or two ';'-separated paths for paired-end
  - gzip integrity verification, with a tolerant recompress fallback for streams whose
    upstream gzip trailer is corrupt (happens on ENA for a small fraction of runs)

Run: TCMA_ROOT=<dir> python3 fetch_manifest.py --study PRJNA412445 [--threads 8]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
import time
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

UA = {"User-Agent": "wisp-science/1.0 (16S fetch)"}


def locate_root() -> Path:
    env = os.environ.get("TCMA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results" / "tables" / "oscc_niche_sample_manifest.csv").is_file():
            return cand
    return here.parents[3]


ROOT = locate_root()
MANIFEST = ROOT / "results" / "tables" / "oscc_niche_sample_manifest.csv"


def to_url(raw: str) -> str:
    raw = raw.strip()
    if raw.startswith(("http://", "https://", "ftp://")):
        return raw
    return "https://" + raw.lstrip("/")


def gzip_ok(p: Path) -> bool:
    try:
        with gzip.open(p, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except Exception:                                        # noqa: BLE001
        return False


def tolerant_recompress(src: Path, dst: Path) -> int:
    """Decompress as far as the stream allows, re-compress cleanly. Returns read count."""
    d = zlib.decompressobj(31)
    tmp = dst.with_suffix(".plain")
    with open(src, "rb") as fi, open(tmp, "wb") as fo:
        while True:
            chunk = fi.read(1 << 20)
            if not chunk:
                break
            try:
                fo.write(d.decompress(chunk))
            except zlib.error:
                break
        try:
            fo.write(d.flush())
        except Exception:                                    # noqa: BLE001
            pass
    n = 0
    with open(tmp, "rb") as fi, gzip.open(dst, "wb", compresslevel=6) as fo:
        while True:
            c = fi.read(1 << 20)
            if not c:
                break
            fo.write(c)
    with open(tmp, "rb") as fi:
        n = sum(1 for _ in fi) // 4
    tmp.unlink(missing_ok=True)
    return n


def fetch(url_raw: str, dest: Path, attempts: int = 4) -> tuple[str, str]:
    url = to_url(url_raw)
    if dest.exists() and dest.stat().st_size > 0 and gzip_ok(dest):
        return "cached", ""
    last = ""
    for i in range(attempts):
        part = dest.with_suffix(dest.suffix + ".part")
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=180) as r, open(part, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            if gzip_ok(part):
                part.replace(dest)
                return "ok", ""
            n = tolerant_recompress(part, dest)
            part.unlink(missing_ok=True)
            return "recovered", f"{n} reads after tolerant decompress"
        except Exception as e:                                # noqa: BLE001
            last = f"{type(e).__name__}: {str(e)[:70]}"
            part.unlink(missing_ok=True)
            time.sleep(2 * (i + 1))
    return "FAILED", last


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    with MANIFEST.open(encoding="utf-8-sig") as fh:
        rows = [r for r in csv.DictReader(fh) if r["study"] == args.study]
    if not rows:
        raise SystemExit(f"no manifest rows for {args.study}")

    dest_dir = ROOT / "data" / "raw" / "16s" / args.study
    dest_dir.mkdir(parents=True, exist_ok=True)
    jobs: list[tuple[str, str, Path]] = []
    for r in rows:
        paths = [p for p in r["fastq_ftp"].split(";") if p.strip()]
        for p in paths:
            name = p.rsplit("/", 1)[-1]
            jobs.append((r["run"], p, dest_dir / name))
    print(f"{args.study}: {len(rows)} runs -> {len(jobs)} files -> {dest_dir}", flush=True)

    t0 = time.time()
    done, failures, recovered = 0, [], []
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = {ex.submit(fetch, u, d): (r, u) for r, u, d in jobs}
        for f in as_completed(futs):
            run, url = futs[f]
            status, note = f.result()
            done += 1
            if status == "FAILED":
                failures.append((run, url.rsplit("/", 1)[-1], note))
            elif status == "recovered":
                recovered.append((run, url.rsplit("/", 1)[-1], note))
            if done % 25 == 0 or done == len(jobs):
                sz = sum(p.stat().st_size for p in dest_dir.glob("*.fastq.gz"))
                print(f"  {done}/{len(jobs)}  {sz/1e6:.0f} MB  {time.time()-t0:.0f}s  "
                      f"fail={len(failures)} recovered={len(recovered)}", flush=True)

    print(f"\ndone in {time.time()-t0:.0f}s; failures={len(failures)} "
          f"recovered={len(recovered)}", flush=True)
    for x in recovered:
        print("  recovered:", x, flush=True)
    for x in failures:
        print("  FAILED:", x, flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
