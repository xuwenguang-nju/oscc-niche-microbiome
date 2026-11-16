#!/usr/bin/env python3
"""S6-v2 step 1: build the TCMA(TCGA-HNSC) sample table and fetch metadata.

WHY THIS REPLACES THE GSE227919 APPROACH
----------------------------------------
The first S6 attempt used GSE227919 (total RNA-seq) because I wrongly concluded
that TCMA was unreachable -- I had probed the wrong hostname. TCMA
(thecancermicrobiomeatlas.org) is reachable, decontaminated, and pairs TCGA-HNSC
microbial counts with the SAME patients' RNA-seq / RPPA / miRNA / methylation.
That is strictly the better substrate for a host-coupling analysis, so S6 is redone
here and the GSE227919 result is kept only as a supplementary direction check.

INPUT
  data/raw/tcma/bacteria.sample.relabund.genus.txt
      TCMA web export, genus level, relative abundance.
      Rows = genera, columns = TCGA aliquot barcodes (plus a "name" column).

METADATA
  The web app's "Download metadata" button did not produce a file, so sample type
  and anatomic site are pulled from the public GDC API by barcode instead. TCGA
  barcode position 14-15 encodes the sample type: 01 = Primary solid Tumor,
  11 = Solid Tissue Normal.

OUTPUT
  results/tables/s6v2_tcma_genus.csv          samples x genera (relative abundance)
  results/tables/s6v2_tcma_samples.csv        per-sample metadata
  results/tables/s6v2_tcma_overview.json

Run: python3 s6v2_build_tcma.py
"""
from __future__ import annotations

import json
import os
import re
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

def locate_root() -> Path:
    """Find the project root in either layout (project tree or compute host work dir).

    A fixed parents[N] works in only one of the two layouts, so search upward for the
    directory that actually holds data/raw/tcma, then fall back to results/ + analysis/.
    TCMA_ROOT can override either.
    """
    env = os.environ.get("TCMA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "data" / "raw" / "tcma").is_dir() or (cand / "results").is_dir():
            return cand
    return here.parents[3]


ROOT = locate_root()
RAW = ROOT / "data" / "raw" / "tcma"
TAB = ROOT / "results" / "tables"
GENUS = RAW / "bacteria.sample.relabund.genus.txt"
UA = {"User-Agent": "wisp-science/1.0 (S6 TCMA)"}
GDC = "https://api.gdc.cancer.gov"

SAMPLE_TYPE = {"01": "Primary Tumor", "11": "Solid Tissue Normal",
               "02": "Recurrent Tumor", "06": "Metastatic"}


def gdc_samples(barcodes: list[str]) -> pd.DataFrame:
    """Look up sample metadata for TCGA aliquot barcodes via the GDC API."""
    # GDC accepts the 16-char sample barcode; derive it from the aliquot barcode.
    want = sorted({b[:16] for b in barcodes})
    rows: list[dict] = []
    chunk = 60
    for i in range(0, len(want), chunk):
        part = want[i:i + chunk]
        payload = {"filters": {"op": "in", "content": {
            "field": "samples.submitter_id", "value": part}},
            "fields": ("submitter_id,case.submitter_id,case.primary_site,"
                       "case.project.project_id,sample_type,tissue_type,"
                       "specimen_type,days_to_collection"),
            "format": "JSON", "size": len(part)}
        req = urllib.request.Request(
            f"{GDC}/samples", data=json.dumps(payload).encode(),
            headers={**UA, "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                hits = json.loads(r.read()).get("data", {}).get("hits", [])
        except Exception as e:                                # noqa: BLE001
            print(f"  GDC chunk {i} failed: {e}")
            hits = []
        for h in hits:
            case = h.get("case") or {}
            rows.append({
                "sample_barcode": h.get("submitter_id"),
                "case_id": case.get("submitter_id"),
                "primary_site": case.get("primary_site"),
                "project": (case.get("project") or {}).get("project_id"),
                "sample_type_gdc": h.get("sample_type"),
                "tissue_type": h.get("tissue_type"),
            })
        print(f"  GDC {min(i+chunk, len(want))}/{len(want)} sample barcodes queried")
    return pd.DataFrame(rows)


def main() -> None:
    if not GENUS.exists():
        raise SystemExit(f"missing {GENUS} -- export it from the TCMA web app "
                         f"(HNSC, genus level, Download data for selection)")
    # TCMA uses tab separators, genus labels may contain no tabs
    ab = pd.read_csv(GENUS, sep="\t", index_col=0)
    ab.index.name = "genus"
    print("TCMA matrix:", ab.shape, "(genera x samples)")
    samples = list(ab.columns)
    print("barcodes:", samples[:3], "...", samples[-1])

    # barcode -> sample type from the barcode itself (TCGA position 14-15)
    def stype(b: str) -> str:
        key = b.split("-")[3][:2] if len(b.split("-")) > 3 else ""
        return SAMPLE_TYPE.get(key, f"other({key})")

    meta = pd.DataFrame({"aliquot_barcode": samples})
    meta["sample_barcode"] = [b[:16] for b in samples]
    meta["case_id"] = [b[:12] for b in samples]
    meta["sample_type"] = [stype(b) for b in samples]
    meta["in_tcga_barcode"] = True

    gdc = gdc_samples(samples)
    if not gdc.empty:
        meta = meta.merge(gdc, on="sample_barcode", how="left",
                          suffixes=("", "_gdc"))
        agree = (meta["sample_type"] == meta["sample_type_gdc"]).mean()
        print(f"barcode-derived vs GDC sample_type agreement: {agree:.1%}")
    else:
        agree = None
        print("GDC metadata unavailable; relying on barcode-derived sample type")

    # samples x genera, relative abundance
    mat = ab.T
    mat.index.name = "aliquot_barcode"
    mat.to_csv(TAB / "s6v2_tcma_genus.csv")
    meta.to_csv(TAB / "s6v2_tcma_samples.csv", index=False)

    overview = {
        "source": "TCMA web export (thecancermicrobiomeatlas.org), v2 (2022)",
        "project": "TCGA-HNSC",
        "taxonomic_level": "genus",
        "n_samples": int(mat.shape[0]),
        "n_genera": int(mat.shape[1]),
        "sample_type_counts": meta["sample_type"].value_counts().to_dict(),
        "gdc_agreement": None if agree is None else round(float(agree), 4),
        "primary_site_counts": (meta["primary_site"].value_counts().head(12).to_dict()
                                if "primary_site" in meta else {}),
        "note": "relative abundance, rows sum to 1 per sample",
        "row_sum_check": {
            "min": round(float(mat.sum(axis=1).min()), 4),
            "max": round(float(mat.sum(axis=1).max()), 4)},
    }
    (TAB / "s6v2_tcma_overview.json").write_text(
        json.dumps(overview, indent=1, ensure_ascii=False, default=str),
        encoding="utf-8")
    print(json.dumps(overview, indent=1, ensure_ascii=False, default=str))
    print(f"\nwrote {TAB}/s6v2_tcma_genus.csv, s6v2_tcma_samples.csv, "
          f"s6v2_tcma_overview.json")


if __name__ == "__main__":
    main()
