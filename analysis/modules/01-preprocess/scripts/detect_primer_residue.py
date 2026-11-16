"""Measure per-cohort 5' primer residue to set DADA2 trimLeft empirically.

Why this matters: an un-trimmed primer produces mismatches at the 5' end, which
shows up as a low filterAndTrim retention rate and an inflated chimera rate. The
pilot run on PRJNA1415359 retained only 57% of reads at filtering, and 85% of R1
reads were found to begin with the 341F primer -- trimLeft = 0 was wrong for
that cohort and would have silently degraded every downstream comparison.

This script tests, for each cohort and each read direction, whether a known 16S
primer sits at the very start of the read, and reports the primer length that
must be trimmed. Cohorts whose reads are already trimmed get trimLeft = 0.

Outputs
-------
results/tables/primer_residue.csv      per cohort x read: fraction with primer at 5'
results/tables/dada2_trimlen.json      machine-readable trimLeft to use

Usage
-----
python analysis/modules/01-preprocess/scripts/detect_primer_residue.py
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
RAW = ROOT / "data" / "raw" / "16s"
TAB = ROOT / "results" / "tables"

N_READS = 3000
START_WINDOW = 26      # allow a few leading bases before the primer
THRESHOLD = 0.50       # >50% of reads carrying the primer => trim it

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}

# 16S primers that can appear at a read start, with their annealing lengths
PRIMERS = {
    "27F":  "AGAGTTTGATCMTGGCTCAG",
    "338F": "ACTCCTACGGGAGGCAGCAG",
    "341F": "CCTACGGGNGGCWGCAG",
    "515F": "GTGCCAGCMGCCGCGGTAA",
    "806R": "GACTACHVGGGTATCTAATCC",
    "926R": "CCGTCAATTCNTTTRAGT",
    "1492R": "TACGGYTACCTTGTTACGACTT",
}


def rc(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def to_regex(primer: str) -> re.Pattern:
    return re.compile("".join(IUPAC[c] for c in primer.upper()))


RX = {k: to_regex(v) for k, v in PRIMERS.items()}
RX_RC = {k: to_regex(rc(v)) for k, v in PRIMERS.items()}


def pick_files(study: str, tag: str) -> list[Path]:
    d = RAW / study
    if tag == "R1":
        f = sorted(d.glob("*_1.fastq.gz"))
        if f:
            return f
        return [x for x in sorted(d.glob("*.fastq.gz"))
                if not x.name.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
    return sorted(d.glob("*_2.fastq.gz"))


def scan(files: list[Path]) -> tuple[int, dict[str, int], dict[str, int]]:
    """Return (n_reads, primer_hits_forward, primer_hits_revcomp)."""
    fwd = {k: 0 for k in PRIMERS}
    rev = {k: 0 for k in PRIMERS}
    n = 0
    for f in files[:3]:
        with gzip.open(f, "rt", errors="replace") as fh:
            while n < N_READS:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline().strip()
                fh.readline()
                fh.readline()
                if not seq:
                    continue
                n += 1
                head = seq[:START_WINDOW]
                for k in PRIMERS:
                    if RX[k].match(head):
                        fwd[k] += 1
                    if RX_RC[k].match(head):
                        rev[k] += 1
    return n, fwd, rev


def main() -> None:
    rows, config = [], {}
    for study in sorted(p.name for p in RAW.iterdir() if p.is_dir()):
        cfg = {"trimLeft": {}, "evidence": {}}
        for tag in ("R1", "R2"):
            files = pick_files(study, tag)
            if not files:
                continue
            n, fwd, rev = scan(files)
            if not n:
                continue
            best = max(fwd.items(), key=lambda kv: kv[1])
            best_rc = max(rev.items(), key=lambda kv: kv[1])
            frac = best[1] / n
            frac_rc = best_rc[1] / n
            if frac >= THRESHOLD and frac >= frac_rc:
                primer, length, direction = best[0], len(PRIMERS[best[0]]), "forward"
            elif frac_rc >= THRESHOLD:
                primer, length, direction = best_rc[0], len(PRIMERS[best_rc[0]]), "revcomp"
            else:
                primer, length, direction = "", 0, "already_trimmed"
            cfg["trimLeft"][tag] = length
            cfg["evidence"][tag] = {
                "reads_tested": n, "primer": primer, "direction": direction,
                "fraction_with_primer_at_5prime": round(frac if direction != "revcomp" else frac_rc, 3),
                "all_hits_forward": {k: v for k, v in fwd.items() if v},
                "all_hits_revcomp": {k: v for k, v in rev.items() if v},
            }
            rows.append({"study": study, "read": tag, "reads_tested": n,
                         "primer": primer, "direction": direction,
                         "fraction_at_5prime": round(frac if direction != "revcomp" else frac_rc, 3),
                         "primer_length": length,
                         "trimLeft": length,
                         "top_forward": best[0] if best[1] else "",
                         "top_revcomp": best_rc[0] if best_rc[1] else ""})
        config[study] = cfg

    with (TAB / "primer_residue.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    (TAB / "dada2_trimlen.json").write_text(
        json.dumps(config, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"{'study':<14} {'read':<5} {'primer':<8} {'dir':<16} {'frac':<6} trimLeft")
    for r in rows:
        print(f"{r['study']:<14} {r['read']:<5} {r['primer'] or '-':<8} "
              f"{r['direction']:<16} {r['fraction_at_5prime']:<6} {r['trimLeft']}")
    print(f"\nwrote {TAB / 'primer_residue.csv'}")
    print(f"wrote {TAB / 'dada2_trimlen.json'}")


if __name__ == "__main__":
    main()
