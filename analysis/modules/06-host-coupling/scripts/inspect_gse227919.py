#!/usr/bin/env python3
"""Structural inspection of the GSE227919 matrices before analysis.

What has to be established before any statistics:
  - microbial matrix: are rows organism-level taxon strings, and what taxa/kingdoms?
    (row ids look like "ti|1969841|org|Proteus_phage_VB_PmiS-Isfahan|accession|NC_...")
  - host matrix: gene symbols, and whether values are raw counts or log-scaled
  - annotations: which column carries the diagnosis group, and what the levels are
  - sample overlap between the two matrices

Run: python3 inspect_gse227919.py
"""
from __future__ import annotations

import gzip
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "gse227919"
MIC = RAW / "GSE227919_2023_20_03_expr_mat_microbe.txt.gz"
HOST = RAW / "GSE227919_2023_20_03_expr_mat_raw.txt.gz"
ANN = RAW / "GSE227919_2023_21_03_pml_annotations.csv.gz"

TORG = re.compile(r"org\|([^|]+)\|")


def sniff_sep(path: Path) -> str:
    """GEOSeries matrices are inconsistently delimited; decide from the header."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
    n_tab = header.count("\t")
    n_sp = header.count(" ")
    sep = "\t" if n_tab > n_sp else r"\s+"
    return sep


def read_matrix(path: Path) -> pd.DataFrame:
    sep = sniff_sep(path)
    df = pd.read_csv(path, sep=sep, index_col=0)
    if df.shape[1] == 0:
        raise RuntimeError(f"{path.name}: parsed 0 columns (sep={sep!r})")
    return df


def genus_from_row(row_id: str) -> str | None:
    """Pull the organism name out of the Kraken-style row id and take its genus."""
    m = TORG.search(row_id)
    if not m:
        return None
    org = m.group(1).replace("_", " ")
    parts = org.split()
    if not parts:
        return None
    return parts[0]


def main() -> None:
    rep: dict = {}

    # ---- microbial matrix -------------------------------------------------
    # NOTE: this file is whitespace-separated (space), not tab-separated, whereas
    # the host matrix is tab-separated. Row ids contain "|" but no spaces, so a
    # whitespace split is unambiguous.
    rep["microbe_parse_sep"] = sniff_sep(MIC)
    mic = read_matrix(MIC)
    rep["microbe"] = {"shape": list(mic.shape),
                      "value_range": [float(mic.values.min()), float(mic.values.max())],
                      "n_negative": int((mic.values < 0).sum())}
    rows = list(mic.index)
    rep["microbe"]["row_id_examples"] = rows[:5]
    rep["microbe"]["n_with_org_field"] = sum(1 for r in rows if TORG.search(r))
    genera = [genus_from_row(r) for r in rows]
    rep["microbe"]["n_genus_parsed"] = sum(1 for g in genera if g)
    rep["microbe"]["top_genera_by_rows"] = Counter(
        g for g in genera if g).most_common(15)
    # how many rows are viruses / phages?
    rep["microbe"]["n_phage_rows"] = sum(1 for r in rows if "phage" in r.lower())

    # ---- host matrix ------------------------------------------------------
    rep["host_parse_sep"] = sniff_sep(HOST)
    host = read_matrix(HOST)
    rep["host"] = {"shape": list(host.shape),
                   "value_range": [float(host.values.min()), float(host.values.max())]}
    hv = host.values
    rep["host"]["looks_like_integer_counts"] = bool(
        abs(hv - hv.round()).max() < 1e-9)
    rep["host"]["looks_log_scaled"] = bool(hv.max() < 40 and hv.min() > -5
                                           and not rep["host"]["looks_like_integer_counts"])
    rep["host"]["gene_examples"] = list(host.index[:8])

    # ---- annotations ------------------------------------------------------
    ann = pd.read_csv(ANN, index_col=0)
    rep["annotations"] = {"shape": list(ann.shape), "columns": list(ann.columns)}
    # find the free-text diagnosis-like column
    for c in ann.columns:
        if ann[c].dtype == object:
            vals = ann[c].dropna().astype(str)
            joined = " ".join(vals).lower()
            if any(k in joined for k in ["leukoplakia", "hyperkeratosis", "carcinoma",
                                         "control", "normal", "dysplasia"]):
                rep["annotations"].setdefault("diagnosis_like_columns", {})[c] = \
                    vals.value_counts().head(10).to_dict()

    # ---- sample overlap ---------------------------------------------------
    mcols = set(mic.columns)
    hcols = set(host.columns)
    rep["overlap"] = {"n_microbe_samples": len(mcols), "n_host_samples": len(hcols),
                      "n_common": len(mcols & hcols),
                      "microbe_only_demo": sorted(mcols - hcols)[:8],
                      "host_only_demo": sorted(hcols - mcols)[:8]}

    dst = ROOT / "results" / "tables" / "s6_gse227919_structure.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(rep, indent=1, ensure_ascii=False, default=str),
                   encoding="utf-8")
    print(json.dumps(rep, indent=1, ensure_ascii=False, default=str)[:6000])
    print(f"\nwrote {dst}")


if __name__ == "__main__":
    main()
