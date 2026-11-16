#!/usr/bin/env python3
"""Download and inspect GSE227919 supplementary matrices.

GSE227919 is total RNA-seq of oral samples (healthy control, hyperkeratosis
non-reactive, leukoplakia with dysplasia, and OSCC). Because total RNA-seq retains
prokaryotic reads, the series deposits a host expression matrix AND a microbial
matrix for the SAME libraries -- which is what S6 needs to couple the two sides.

Downloads (small, supplementary only -- raw FASTQ not needed):
  GSE227919_2023_20_03_expr_mat_microbe.txt.gz   microbial abundance
  GSE227919_2023_20_03_expr_mat_raw.txt.gz       host raw counts
  GSE227919_2023_21_03_pml_annotations.csv.gz    sample annotations

Run: python3 fetch_gse227919.py
"""
from __future__ import annotations

import gzip
import io
import json
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DEST = ROOT / "data" / "raw" / "gse227919"
UA = {"User-Agent": "wisp-science/1.0 (S6 fetch)"}
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227919/suppl"
FILES = [
    "GSE227919_2023_20_03_expr_mat_microbe.txt.gz",
    "GSE227919_2023_20_03_expr_mat_raw.txt.gz",
    "GSE227919_2023_21_03_pml_annotations.csv.gz",
]


def fetch(name: str) -> Path:
    dest = DEST / name
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  cached: {name} ({dest.stat().st_size:,} B)")
        return dest
    url = f"{BASE}/{name}"
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    # integrity: the whole stream must decompress
    with gzip.open(dest, "rb") as fh:
        while fh.read(1 << 20):
            pass
    print(f"  {name}: {dest.stat().st_size:,} B  gzip OK")
    return dest


def peek(path: Path, n_lines: int = 3) -> dict:
    info = {"file": path.name, "bytes": path.stat().st_size}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        head = [next(fh) for _ in range(n_lines)]
    info["header"] = head[0].rstrip("\n")[:300]
    info["n_cols"] = head[0].count("\t") + 1 if "\t" in head[0] else head[0].count(",") + 1
    info["row_example"] = (head[1].rstrip("\n")[:200] if len(head) > 1 else "")
    return info


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    print("downloading GSE227919 supplementary ...")
    paths = [fetch(f) for f in FILES]

    report = {"series": "GSE227919", "files": [peek(p) for p in paths]}
    dst = ROOT / "results" / "tables" / "s6_gse227919_files.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")

    print("\n=== structure ===")
    for f in report["files"]:
        print(f"\n{f['file']}  ({f['bytes']:,} B, {f['n_cols']} cols)")
        print(f"  header: {f['header']}")
        print(f"  row1  : {f['row_example']}")
    print(f"\nwrote {dst}")


if __name__ == "__main__":
    main()
