#!/usr/bin/env python3
"""Profile a cohort's reads and check for 5' primer residue, before choosing DADA2 params.

MiSeq amplicon submissions differ on whether the sequencing primers were trimmed, and
getting that wrong silently truncates real sequence or leaves primer-derived errors. This
reports, per cohort:
  - read-length distribution (R1 and R2 separately)
  - mean quality by position
  - fraction of reads carrying an 8-mer from the V3-V4 primers (341F / 805R) within the
    first 30 bp, in both orientations

Run: TCMA_ROOT=<dir> python3 profile_cohort.py --study PRJNA412445 [--runs 8] [--per-run 2000]
"""
from __future__ import annotations

import argparse
import gzip
import os
import statistics
from pathlib import Path

import numpy as np

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
PRIMERS = {  # 8-bp 5' seeds of the standard V3-V4 primers
    "341F": "CCTACGGG",
    "805R": "GACTACHV",
}


def iupac_regex(seq: str) -> str:
    return "".join(IUPAC.get(c, c) for c in seq.upper())


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def locate_root() -> Path:
    env = os.environ.get("TCMA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "data" / "raw" / "16s").is_dir():
            return cand
    return here.parents[3]


ROOT = locate_root()


def read_fastq(path: Path, limit: int):
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
    import re
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True)
    ap.add_argument("--runs", type=int, default=8)
    ap.add_argument("--per-run", type=int, default=2000)
    args = ap.parse_args()

    raw = ROOT / "data" / "raw" / "16s" / args.study
    files = sorted(raw.glob("*.fastq.gz"))
    if not files:
        raise SystemExit(f"no fastq under {raw}")
    r1 = [f for f in files if f.name.endswith("_1.fastq.gz")]
    r2 = [f for f in files if f.name.endswith("_2.fastq.gz")]
    single = not r1
    print(f"{args.study}: {len(files)} files; paired={not single} "
          f"(R1={len(r1)}, R2={len(r2)})")

    regexes = {k: re.compile(iupac_regex(v)) for k, v in PRIMERS.items()}
    step = max(1, len(files) // args.runs)
    picks = files[::step][:args.runs]

    lens_all: dict[str, list[int]] = {}
    qsum: dict[str, np.ndarray | None] = {}
    qn: dict[str, int] = {}
    hits: dict[str, dict[str, int]] = {}

    for p in picks:
        mate = "R2" if p.name.endswith("_2.fastq.gz") else "R1"
        L = []
        for seq, qual in read_fastq(p, args.per_run):
            L.append(len(seq))
            q = np.frombuffer(qual.encode("ascii", "replace"), dtype=np.uint8).astype(np.int16) - 33
            q = q[:320]
            if qsum.get(mate) is None:
                qsum[mate] = np.zeros(320)
            qsum[mate][:len(q)] += q
            qn[mate] = qn.get(mate, 0) + 1
            head = seq[:30].upper()
            for name, rgx in regexes.items():
                d = hits.setdefault(mate, {})
                d[name] = d.get(name, 0) + (1 if rgx.search(head) else 0)
                rc = revcomp(head)
                d[name + "_rc"] = d.get(name + "_rc", 0) + (1 if rgx.search(rc) else 0)
        lens_all.setdefault(mate, []).extend(L)
        print(f"  {p.stem}: n={len(L)} median={int(statistics.median(L))} "
              f"min={min(L)} max={max(L)}", flush=True)

    print("\n=== length distribution ===")
    for mate, L in lens_all.items():
        a = np.array(L)
        print(f"  {mate}: n={a.size} median={int(np.median(a))} "
              f"IQR=[{int(np.percentile(a,25))},{int(np.percentile(a,75))}] max={a.max()}")

    print("\n=== mean quality by position (every 20 bp) ===")
    for mate in qsum:
        mq = qsum[mate] / max(qn[mate], 1)
        s = " ".join(f"{pos}:{mq[pos]:.0f}" for pos in range(0, 320, 20))
        print(f"  {mate}: {s}")
        for thr in (20, 25, 28):
            ok = np.where(mq >= thr)[0]
            print(f"     Q>={thr} up to position {ok.max() if ok.size else 'n/a'}")

    print("\n=== 5' primer seed hits (first 30 bp) ===")
    for mate, d in hits.items():
        n = qn.get(mate, 1)
        for k, v in sorted(d.items()):
            print(f"  {mate} {k:<10} {100*v/n:5.1f}%")

    print("\n=== bases retained at candidate truncLen ===")
    for mate, L in lens_all.items():
        a = np.array(L)
        for cut in (200, 220, 240, 250, 260, 280, 300):
            keep = float((a >= cut).mean())
            bkeep = float(np.minimum(a, cut).sum() / a.sum())
            print(f"  {mate} truncLen={cut}: reads kept={100*keep:5.1f}%  "
                  f"bases kept={100*bkeep:5.1f}%")


if __name__ == "__main__":
    main()
