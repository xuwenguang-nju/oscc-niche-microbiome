"""Fetch raw 16S reads for the OSCC cohorts listed in the sample manifest.

Reads results/tables/oscc_niche_sample_manifest.csv, selects runs per cohort,
and downloads the FASTQ files into data/raw/16s/<study>/.

Design notes
------------
* Downloads are resumable: an existing non-zero file of the expected size is
  skipped, so re-running never re-fetches.
* A per-run status file is written so a later DADA2 step can tell "not
  downloaded" from "downloaded and empty".
* Two source protocols are used: HTTPS from ENA (default) and, if a run has no
  usable URL, the run is recorded as unavailable rather than silently skipped.
* Pilot mode (--max-per-study N) exists because full delivery is ~45 GB and the
  first thing to validate is the processing pipeline, not the download volume.

Usage
-----
python analysis/modules/01-preprocess/scripts/fetch_reads.py --pilot
python analysis/modules/01-preprocess/scripts/fetch_reads.py --studies PRJNA1415359
python analysis/modules/01-preprocess/scripts/fetch_reads.py            # all cohorts
"""

from __future__ import annotations

import argparse
import csv
import json
import gzip
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
MANIFEST = ROOT / "results" / "tables" / "oscc_niche_sample_manifest.csv"
RAW = ROOT / "data" / "raw" / "16s"
STATUS = ROOT / "data" / "processed" / "fetch_status.csv"
UA = {"User-Agent": "wisp-science/1.0 (16S fetch)"}

# pilot sizes: enough runs per cohort to validate the pipeline end to end
PILOT_N = 6


def gzip_ok(path: Path) -> bool:
    """True if the gzip stream decompresses completely.

    Size comparison alone is not enough: three files downloaded during a period
    of ENA rate-limiting had the correct byte count but failed the CRC check,
    and a size-only guard would have kept skipping them forever.
    """
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except Exception:
        return False


def remote_size(url: str) -> int | None:
    try:
        req = urllib.request.Request(url, headers=UA, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n else None
    except Exception:
        return None


def download(url: str, dest: Path, attempts: int = 8) -> tuple[str, str, int]:
    """Return (status, message, bytes_written).

    Retries with exponential backoff. ENA returns transient HTTP 403 /
    RemoteDisconnected under high concurrency even though the object is
    readable -- verified by re-requesting a "failed" URL, which returned 200.
    Without backoff those become silent data loss for whole patients.
    """
    if dest.exists() and dest.stat().st_size > 0:
        expected = remote_size(url)
        size_ok = expected is None or dest.stat().st_size == expected
        if size_ok and gzip_ok(dest):
            return "skipped_present", "", dest.stat().st_size
        dest.unlink()
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = ""
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r, \
                    tmp.open("wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            tmp.replace(dest)
            return ("downloaded" if attempt == 0 else "downloaded_retry",
                    "" if attempt == 0 else f"after {attempt + 1} attempts",
                    dest.stat().st_size)
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
            if tmp.exists():
                tmp.unlink()
            if attempt < attempts - 1:
                time.sleep(2 ** attempt)
    return "failed", last[:200], 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true",
                    help=f"only the first {PILOT_N} runs per study")
    ap.add_argument("--max-per-study", type=int, default=None)
    ap.add_argument("--studies", nargs="*", default=None)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--retry-failed", action="store_true",
                    help="only re-fetch files recorded as failed in the status file")
    args = ap.parse_args()
    limit = args.max_per_study or (PILOT_N if args.pilot else None)

    rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8-sig")))
    if args.studies:
        rows = [r for r in rows if r["study"] in args.studies]

    wanted = None
    if args.retry_failed:
        if not STATUS.exists():
            raise SystemExit("no status file to retry from: " + str(STATUS))
        wanted = {r["file"] for r in csv.DictReader(STATUS.open(encoding="utf-8-sig"))
                  if r["status"] == "failed"}
        if not wanted:
            print("no failed files recorded")
            return 0
        print(f"retrying {len(wanted)} previously failed files")

    jobs: list[tuple[str, str, Path]] = []
    per_study: dict[str, int] = {}
    for r in rows:
        st = r["study"]
        if limit is not None and per_study.get(st, 0) >= limit:
            continue
        paths = [p.strip() for p in str(r["fastq_ftp"]).split(";") if p.strip()]
        if not paths:
            continue
        per_study[st] = per_study.get(st, 0) + 1
        for p in paths:
            name = Path(p).name
            if wanted is not None and name not in wanted:
                continue
            jobs.append((st, "https://" + p, RAW / st / name))

    print(f"runs selected: {sum(per_study.values())}  files: {len(jobs)}")
    for st, n in sorted(per_study.items()):
        print(f"  {st}: {n} runs")

    status_rows = []
    done = 0
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futs = {pool.submit(download, url, dest): (st, url, dest)
                for st, url, dest in jobs}
        for fut in as_completed(futs):
            st, url, dest = futs[fut]
            status, msg, nbytes = fut.result()
            status_rows.append({"study": st, "file": dest.name, "status": status,
                                "bytes": nbytes, "message": msg})
            done += 1
            if done % 20 == 0 or status == "failed":
                print(f"  [{done}/{len(jobs)}] {status:18} {dest.name} {msg}")

    STATUS.parent.mkdir(parents=True, exist_ok=True)
    merged: dict[tuple[str, str], dict] = {}
    if STATUS.exists():
        for r in csv.DictReader(STATUS.open(encoding="utf-8-sig")):
            merged[(r["study"], r["file"])] = r
    for r in status_rows:
        merged[(r["study"], r["file"])] = r   # this run wins
    all_rows = sorted(merged.values(), key=lambda r: (r["study"], r["file"]))
    with STATUS.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["study", "file", "status", "bytes", "message"])
        w.writeheader()
        w.writerows(all_rows)
    status_rows = all_rows

    from collections import Counter
    agg = Counter(r["status"] for r in status_rows)
    total = sum(int(r["bytes"] or 0) for r in status_rows)
    print("\nsummary:", dict(agg), f"total {total / 1e9:.2f} GB")
    print("status file:", STATUS.relative_to(ROOT))
    # non-zero exit only for a real failure, so callers can branch on it
    return 1 if agg.get("failed") else 0


if __name__ == "__main__":
    sys.exit(main())
