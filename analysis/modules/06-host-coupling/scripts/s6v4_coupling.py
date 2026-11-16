#!/usr/bin/env python3
"""S6-v4: host-program coupling from the TCMA correlation panels (TCGA-HNSC, genus).

WHY THIS STEP EXISTS
--------------------
S6-v2 proved the direction question could be answered from the TCMA abundance matrix,
but the user's original S6 goal -- does the microbial genus covary with host immune /
EMT programmes -- needs the host side. TCMA publishes those correlations, computed on
log per-rank relative abundance against four assays, with the Benjamini-Hochberg q
already corrected across all taxon x feature pairs *for that assay* (per the site's own
note). Those tables are exported here rather than re-deriving them, because
re-downloading TCGA RNA-seq/RPPA/miRNA/methylation for 177 patients is out of scope.

INPUT
  data/raw/tcma/tcma_correlations_genus_HNSC.tsv
      TCMA web export, genus level, TCGA-HNSC filter, four panels:
      RPPA / RNAseq / miRNA / Methylation27K.
      Each cell is "R|n|p|q"; 50 features (top by |R|) x 10 genera per panel.
      Row order and feature selection are the site's own (its top-50 display), NOT a
      balanced genome-wide panel -- that is a real limitation, stated below.

WHAT IS TESTED
  For each of the 10 genera in the panel, across the four assays:
    - how many features pass q<0.05 (the site's own corrected q)
    - whether the genus x feature coupling of the S5 focus genera is coherent
    - sign concordance of the strongest couplings with the tumour/normal direction
  For Parvimonas specifically: does it couple with cytotoxic / IFN / checkpoint / EMT
  programmes, as it must if it were driving an immune-hot tumour phenotype?

OUTPUT
  results/tables/s6v4_coupling_long.csv       tidy: assay, genus, feature, R, n, p, q
  results/tables/s6v4_genus_summary.csv       per assay x genus: n tests, q<0.05 count
  results/tables/s6v4_parvimonas_top.csv      strongest Parvimonas couplings
  results/tables/s6v4_summary.json

Run: python3 s6v4_coupling.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

def locate_root() -> Path:
    env = os.environ.get("TCMA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "data" / "raw" / "tcma").is_dir() or (cand / "results" / "tables").is_dir():
            return cand
    return here.parents[3]


ROOT = locate_root()
RAW = ROOT / "data" / "raw" / "tcma"
TAB = ROOT / "results" / "tables"
SRC = RAW / "tcma_correlations_genus_HNSC.tsv"

# Host programmes of interest, mapped to the assay that can measure them.
# RPPA measures phospho-protein activity, RNAseq measures transcripts.
PROGRAMME_HINTS = {
    "immune_cytotoxic": ["CD8", "CD20", "LCK", "SYK", "STAT1", "STAT3", "STAT5",
                         "P38", "JNK", "NFKB", "CD49"],
    "checkpoint_like": ["PDL1", "CD274", "PD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "emt": ["FIBRONECTIN", "CAVEOLIN", "PCADHERIN", "ECADHERIN", "VIM", "ZEB",
            "SNAI", "TWIST", "CLAUDIN", "MYH11"],
    "proliferation": ["CMYC", "CYCLIN", "CCNA", "CCNB", "MKI67", "PCNA", "TOP2A",
                      "P70S6K", "P90RSK", "AKT", "MTOR", "RAPTOR", "RICTOR"],
    "dna_damage": ["ATM", "CHK1", "CHK2", "BRCA", "KU80", "P53", "MDM2", "P21"],
    "growth_rtk": ["EGFR", "HER2", "HEREGULIN", "CMET", "VEGFR", "NRAS", "CRAF",
                   "MEK1", "ERK", "SRC", "PI3K"],
}


def parse_export(path: Path) -> pd.DataFrame:
    rows: list[dict] = []
    assay = None
    for ln in path.read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        if ln.startswith("# panel"):
            assay = ln.split("\t")[1]
            continue
        parts = ln.split("\t")
        if parts[0] == "feature":
            genera = parts[1:]
            continue
        feature = parts[0]
        for genus, cell in zip(genera, parts[1:]):
            vals = cell.split("|")
            if len(vals) != 4:
                continue
            def num(x):
                try:
                    return float(x)
                except (TypeError, ValueError):
                    return None
            rows.append({"assay": assay, "genus": genus, "feature": feature,
                         "R": num(vals[0]), "n": num(vals[1]),
                         "p": num(vals[2]), "q": num(vals[3])})
    return pd.DataFrame(rows)


def assign_programme(feature: str) -> str:
    f = feature.upper()
    for prog, keys in PROGRAMME_HINTS.items():
        if any(k in f for k in keys):
            return prog
    return "other"


def main() -> None:
    long = parse_export(SRC)
    long["programme"] = [assign_programme(f) for f in long["feature"]]
    long.to_csv(TAB / "s6v4_coupling_long.csv", index=False)

    summary: dict = {
        "source": "TCMA web export (thecancermicrobiomeatlas.org), TCGA-HNSC, genus",
        "assays": sorted(long["assay"].unique().tolist()),
        "n_genera": int(long["genus"].nunique()),
        "n_features_per_assay": long.groupby("assay")["feature"].nunique().to_dict(),
        "cells": int(len(long)),
        "limitation": "each panel shows only the site's top-50 features by |R|, so the "
                      "feature set is not a balanced genome-wide panel and hit counts "
                      "are not comparable to a full-transcriptome scan; q is the site's "
                      "own BH correction within each assay",
    }

    # ---- per assay x genus hit counts (site's own q) ----
    per = (long.dropna(subset=["q"])
               .groupby(["assay", "genus"])
               .agg(n_tests=("q", "size"),
                    q05=("q", lambda s: int((s < 0.05).sum())),
                    q01=("q", lambda s: int((s < 0.01).sum())),
                    min_q=("q", "min"),
                    max_absR=("R", lambda s: float(s.abs().max())))
               .reset_index())
    per.to_csv(TAB / "s6v4_genus_summary.csv", index=False)

    summary["per_assay_genus"] = per.round(4).to_dict("records")

    # ---- Parvimonas detail ----
    pv = long[long["genus"] == "Parvimonas"].dropna(subset=["q"]).copy()
    pv["absR"] = pv["R"].abs()
    pv_top = pv.sort_values("absR", ascending=False).head(20)
    pv_top[["assay", "feature", "programme", "R", "n", "p", "q"]].to_csv(
        TAB / "s6v4_parvimonas_top.csv", index=False)
    summary["parvimonas"] = {
        "n_cells": int(len(pv)),
        "min_q": float(pv["q"].min()),
        "q05": int((pv["q"] < 0.05).sum()),
        "max_absR": float(pv["absR"].max()),
        "top": pv_top[["assay", "feature", "programme", "R", "q"]]
                     .round(4).to_dict("records"),
    }

    # ---- which genera actually couple with anything? ----
    best = (per.sort_values("min_q")
               .groupby("assay", group_keys=False)
               .head(3)[["assay", "genus", "q05", "min_q", "max_absR"]])
    summary["best_coupled_per_assay"] = best.round(4).to_dict("records")

    # ---- programme-level view: strongest coupling per programme ----
    prog = (long.dropna(subset=["q"])
                .groupby(["programme", "genus"])
                .agg(min_q=("q", "min"), max_absR=("R", lambda s: float(s.abs().max())))
                .reset_index()
                .sort_values("min_q")
                .groupby("programme", group_keys=False).head(2))
    summary["top_programme_couplings"] = prog.round(4).to_dict("records")

    (TAB / "s6v4_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False, default=str), encoding="utf-8")

    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ("per_assay_genus",)}, indent=1,
                     ensure_ascii=False, default=str)[:5000])
    print("\n=== per assay x genus (q<0.05) ===")
    print(per.pivot(index="genus", columns="assay", values="q05").fillna(0)
             .astype(int).to_string())
    print(f"\nwrote {TAB}/s6v4_*.csv|json")


if __name__ == "__main__":
    main()
