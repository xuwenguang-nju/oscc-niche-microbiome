"""Watch a DADA2 run and report stage transitions without flooding the log.

Motivation: polling from the shell for a multi-hour step costs one tool call per
minute and tells us nothing new. This watcher samples the run every 30 s and
prints only when the observed stage changes or throughput stalls.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
STUDY = "PRJNA666746"
RD = ROOT / "results" / "dada2" / STUDY
LOG = ROOT / "analysis" / "out" / "dada2_PRJNA666746_ds2.log"


def counts() -> dict[str, int]:
    def n(p: Path) -> int:
        return len(list(p.glob("*"))) if p.exists() else 0

    return {
        "raw_ds": n(RD / "raw_ds"),
        "filtered": n(RD / "filtered"),
    }


def main() -> int:
    minutes = float(sys.argv[1]) if len(sys.argv) > 1 else 25.0
    deadline = time.time() + minutes * 60
    last = None
    last_change = time.time()
    while time.time() < deadline:
        st = counts()
        done = (RD / "run_status.json").exists()
        if st != last:
            print(f"[{time.strftime('%H:%M:%S')}] {st}", flush=True)
            if last is not None:
                moved = any(st[k] != last[k] for k in st)
                if moved:
                    last_change = time.time()
            last = st
        if done:
            print("RUN STATUS FILE PRESENT", flush=True)
            print((RD / "run_status.json").read_text(encoding="utf-8"), flush=True)
            return 0
        if time.time() - last_change > 900:
            print(f"[{time.strftime('%H:%M:%S')}] STALLED >15 min at {st}", flush=True)
            return 2
        time.sleep(30)
    print(f"window ended at {counts()}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
