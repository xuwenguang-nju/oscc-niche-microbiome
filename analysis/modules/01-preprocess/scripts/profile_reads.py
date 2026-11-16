"""Empirical read profiling for parameter selection.

Why this exists: truncLen / trimLeft / maxEE are the result-changing parameters
of the DADA2 step, and they must be chosen from the data, not from a template.
This script measures, per cohort:

  * read length distribution (which is what distinguishes an Ion Torrent
    single-end library from an Illumina paired library)
  * per-cycle mean Phred quality for R1 and R2
  * the cycle at which mean quality first drops below a threshold
  * the primer bases actually present at the 5' end (so trimLeft is justified)

Outputs
-------
results/tables/read_profile_by_cohort.csv     one row per cohort x read
results/tables/read_profile_summary.csv       chosen truncLen / trimLeft per cohort
analysis/out/read_profile_report.txt          human-readable reasoning

Usage
-----
python analysis/modules/01-preprocess/scripts/profile_reads.py
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
RAW = ROOT / "data" / "raw" / "16s"
TAB = ROOT / "results" / "tables"
OUT = ROOT / "analysis" / "out"
TAB.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

MAX_READS_PER_FILE = 20_000
Q_DROP = 30          # cycle at which mean quality falls below this
Q_HARD = 20          # never allow truncLen beyond a cycle below this


def parse(path: Path, max_reads: int = MAX_READS_PER_FILE):
    """Return (n, length_counts, per-cycle quality sums, per-cycle base counts)."""
    n = 0
    lengths: dict[int, int] = {}
    qsum: list[int] = []
    bases: list[dict[str, int]] = []
    with gzip.open(path, "rt", errors="replace") as fh:
        while n < max_reads:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            fh.readline()
            qual = fh.readline().strip()
            if not seq:
                continue
            n += 1
            lengths[len(seq)] = lengths.get(len(seq), 0) + 1
            for i, ch in enumerate(qual):
                q = ord(ch) - 33
                if i >= len(qsum):
                    qsum.append(0)
                    bases.append({})
                qsum[i] += q
                b = seq[i] if i < len(seq) else "N"
                bases[i][b] = bases[i].get(b, 0) + 1
    return n, lengths, qsum, bases


def main() -> None:
    rows, summary, report = [], [], []
    report.append("Empirical read profile (first 20,000 reads per file)\n" + "=" * 92)

    for study in sorted(p.name for p in RAW.iterdir() if p.is_dir()):
        # ENA names paired files <run>_1.fastq.gz / <run>_2.fastq.gz; single-end
        # libraries have a single <run>.fastq.gz
        pair = {"R1": ("*_1.fastq.gz", "_1.fastq.gz"),
                "R2": ("*_2.fastq.gz", "_2.fastq.gz")}
        single = [f for f in sorted(RAW.joinpath(study).glob("*.fastq.gz"))
                  if not f.name.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
        for read_tag in ("R1", "R2"):
            pat, suf = pair[read_tag]
            files = sorted(RAW.joinpath(study).glob(pat))
            if not files:
                if read_tag == "R2":
                    continue
                files = single
            if not files:
                continue
            use = files[:3]
            n_tot = 0
            lengths: dict[int, int] = {}
            qsum_acc: list[int] = []
            bases_acc: list[dict[str, int]] = []
            for f in use:
                n, lens, qs, bs = parse(f)
                n_tot += n
                for k, v in lens.items():
                    lengths[k] = lengths.get(k, 0) + v
                for i, v in enumerate(qs):
                    if i >= len(qsum_acc):
                        qsum_acc.append(0)
                        bases_acc.append({})
                    qsum_acc[i] += v
                for i, d in enumerate(bs):
                    if i >= len(bases_acc):
                        qsum_acc.append(0)
                        bases_acc.append({})
                    for b, c in d.items():
                        bases_acc[i][b] = bases_acc[i].get(b, 0) + c
            if not n_tot:
                continue
            mean_len = sum(k * v for k, v in lengths.items()) / sum(lengths.values())
            len_min, len_max = min(lengths), max(lengths)
            ncycles = len(qsum_acc)
            mean_q = [qsum_acc[i] / n_tot for i in range(ncycles)]

            # first cycle whose mean quality dips below Q_DROP
            drop = next((i for i, q in enumerate(mean_q) if q < Q_DROP), None)
            hard = next((i for i, q in enumerate(mean_q) if q < Q_HARD), None)
            last_ok = (hard - 1) if hard is not None else ncycles - 1
            # DADA2 truncLen must be > 0 and consistent with the overlap requirement
            trunc = max(0, min(last_ok, ncycles - 1))

            # dominant bases in the first 4 cycles (primer check)
            head = []
            for i in range(min(4, ncycles)):
                top = sorted(bases_acc[i].items(), key=lambda kv: -kv[1])[:2]
                head.append("".join(b for b, _ in top))

            rows.append({
                "study": study, "read": read_tag, "files_sampled": len(use),
                "reads_sampled": n_tot, "len_mean": round(mean_len, 1),
                "len_min": len_min, "len_max": len_max,
                "cycles": ncycles, "meanq_cycle0": round(mean_q[0], 1),
                "first_cycle_below_q30": drop, "first_cycle_below_q20": hard,
                "suggested_truncLen": trunc,
                "head_bases_c1_c4": "/".join(head),
            })
            summary.append({"study": study, "read": read_tag,
                            "suggested_truncLen": trunc,
                            "head_bases": "/".join(head),
                            "len_mean": round(mean_len, 1),
                            "len_min": len_min, "len_max": len_max})

            report.append(
                f"\n{study}  {read_tag}  ({len(use)} files, {n_tot:,} reads)")
            report.append(
                f"  length: mean {mean_len:.1f}  range {len_min}-{len_max}"
                f"  cycles {ncycles}")
            report.append(
                f"  mean Q cycle0 {mean_q[0]:.1f} | first cycle < Q30: {drop}"
                f" | first cycle < Q20: {hard} | suggested truncLen {trunc}")
            report.append(f"  5' dominant bases (cycles 1-4): {'/'.join(head)}")
            step = max(1, ncycles // 10)
            prof = ", ".join(f"{i}:{mean_q[i]:.0f}" for i in range(0, ncycles, step))
            report.append(f"  quality profile: {prof}")

    with (TAB / "read_profile_by_cohort.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with (TAB / "read_profile_summary.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    (OUT / "read_profile_report.txt").write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report))
    print(f"\nwrote {TAB / 'read_profile_by_cohort.csv'}")
    print(f"wrote {TAB / 'read_profile_summary.csv'}")


if __name__ == "__main__":
    main()
