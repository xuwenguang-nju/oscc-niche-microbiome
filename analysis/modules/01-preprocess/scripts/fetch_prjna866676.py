#!/usr/bin/env python3
"""Download all PRJNA866676 runs (Ion Torrent, six ecological niches).

WHY THIS DATASET MATTERS
------------------------
PRJNA866676 is the only OSCC cohort in the project that samples SIX niches from the same
patients: tumour tissue, paracancerous tissue, contralateral normal mucosa, cancer
surface, saliva, and tongue coat. That is what the "which niche does the tumour
microbiome resemble" question needs, and it is also the only cohort that can test a
backflow gradient (cancer surface is in direct contact with the oral cavity; tumour
tissue is deep).

Runs are single-end Ion Torrent (~105k reads/run at ~428-441 bp).

Downloads are resumable-by-absence (skip existing non-empty files) and every file is
gzip-verified, because ENA produces occasional right-sized but corrupt streams.

Run: python3 fetch_prjna866676.py [--threads 8]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

def locate_root() -> Path:
    """Find the project root in either layout (project tree or compute-host workdir).

    A fixed parents[N] resolves to the home directory on the compute host, so search
    upward for the directory that actually holds the manifest instead.
    """
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
DEST = ROOT / "data" / "raw" / "16s" / "PRJNA866676"
STUDY = "PRJNA866676"
UA = {"User-Agent": "wisp-science/1.0 (16S fetch)"}


def gzip_ok(p: Path) -> bool:
    try:
        with gzip.open(p, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except Exception:                                         # noqa: BLE001
        return False


def to_url(raw: str) -> str:
    """ENA's fastq_ftp column has no scheme for some studies (this one), with it for
    others. Normalise to https, which the EBI mirror serves."""
    raw = raw.strip()
    if raw.startswith(("http://", "https://", "ftp://")):
        return raw
    return "https://" + raw.lstrip("/")


def fetch_one(run: str, url: str, attempts: int = 4) -> tuple[str, str]:
    url = to_url(url)
    dest = DEST / f"{run}.fastq.gz"
    if dest.exists() and dest.stat().st_size > 0 and gzip_ok(dest):
        return run, "cached"
    last = ""
    for i in range(attempts):
        part = dest.with_suffix(".part")
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=120) as r, open(part, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            if not gzip_ok(part):
                part.unlink(missing_ok=True)
                last = "gzip corrupted"
                time.sleep(2 * (i + 1))
                continue
            part.replace(dest)
            return run, "ok"
        except Exception as e:                                # noqa: BLE001
            last = f"{type(e).__name__}: {str(e)[:80]}"
            part.unlink(missing_ok=True)
            time.sleep(2 * (i + 1))
    return run, f"FAILED {last}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    DEST.mkdir(parents=True, exist_ok=True)
    # utf-8-sig: the manifest is written with a BOM, which otherwise makes the
    # first column name "\ufeffstudy" and silently yields zero rows.
    with MANIFEST.open(encoding="utf-8-sig") as fh:
        rows = [r for r in csv.DictReader(fh) if r["study"] == STUDY]
    print(f"runs to fetch: {len(rows)} -> {DEST}", flush=True)

    t0 = time.time()
    done = 0
    failures: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = {ex.submit(fetch_one, r["run"], r["fastq_ftp"]): r["run"] for r in rows}
        for f in as_completed(futs):
            run, status = f.result()
            done += 1
            if status.startswith("FAILED"):
                failures.append((run, status))
            if done % 25 == 0 or done == len(rows):
                sz = sum(p.stat().st_size for p in DEST.glob("*.fastq.gz"))
                print(f"  {done}/{len(rows)}  {sz/1e9:.2f} GB  "
                      f"{time.time()-t0:.0f}s  failures={len(failures)}", flush=True)

    print(f"\ndone in {time.time()-t0:.0f}s; failures: {len(failures)}", flush=True)
    for run, st in failures:
        print(f"  {run}  {st}", flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
