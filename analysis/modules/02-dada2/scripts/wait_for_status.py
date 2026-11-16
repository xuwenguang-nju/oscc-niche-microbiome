"""Wait for a cohort's DADA2 run to produce run_status.json, then report.

Why: this environment cannot reliably hold a long foreground poll. The pipeline
writes run_status.json as its last action, so its appearance is the single
unambiguous completion signal (exit codes are not trustworthy here -- a crashed
native child can still leave exit 0).

Usage
-----
python analysis/modules/02-dada2/scripts/wait_for_status.py --study PRJNA1415359 --timeout 3600
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
STATUS = ROOT / "results" / "dada2"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True)
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--poll", type=int, default=30)
    args = ap.parse_args()

    target = STATUS / args.study / "run_status.json"
    t0 = time.time()
    last = -1
    while time.time() - t0 < args.timeout:
        # progress proxy: filtered file count
        fdir = STATUS / args.study / "filtered"
        n = len(list(fdir.glob("*.fastq.gz"))) if fdir.exists() else 0
        if n != last:
            print(f"[{time.time() - t0:6.0f}s] filtered files: {n}", flush=True)
            last = n
        if target.exists():
            time.sleep(3)   # let the writer finish
            try:
                st = json.loads(target.read_text(encoding="utf-8"))
            except Exception as e:
                print("status file not yet parseable:", e, flush=True)
                time.sleep(5)
                continue
            print("\n=== run_status ===", flush=True)
            print(json.dumps(st, indent=1, ensure_ascii=False), flush=True)
            print("\nelapsed", round(time.time() - t0), "s", flush=True)
            return 0 if st.get("status") == "ok" else 1
        time.sleep(args.poll)

    print(f"\nTIMEOUT after {args.timeout}s: {target} never appeared", flush=True)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
