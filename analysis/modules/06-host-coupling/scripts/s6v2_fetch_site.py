#!/usr/bin/env python3
"""Fetch per-sample anatomic site for TCGA-HNSC from the GDC API.

WHY
---
The user asked whether Parvimonas (the one genus that was up in tumour across all four
16S cohorts but came out DOWN in TCMA) is masked by pooling oral cavity with
pharynx/larynx. To stratify, every TCMA aliquot needs its anatomic site, which the TCMA
genus export does not carry.

The GDC /samples endpoint returns 404 on this deployment, but /cases works and carries
both `primary_site` and the nested `samples` list, so the mapping can be reconstructed
from there.

OUTPUT
  results/tables/s6v2_gdc_primary_site.csv   sample_barcode, case_id, primary_site,
                                             sample_type (GDC's own label)
Run: python3 s6v2_fetch_site.py
"""
from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

GDC = "https://api.gdc.cancer.gov/cases"
UA = {"User-Agent": "wisp-science/1.0 (S6 site lookup)"}


def locate_root() -> Path:
    env = os.environ.get("TCMA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results" / "tables").is_dir():
            return cand
    return here.parents[3]


ROOT = locate_root()
TAB = ROOT / "results" / "tables"


def post(url: str, payload: dict, timeout: int = 120) -> dict:
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={**UA, "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def main() -> None:
    rows: list[dict] = []
    page = 0
    size = 500
    while True:
        payload = {
            "filters": {"op": "in", "content": {
                "field": "project.project_id", "value": ["TCGA-HNSC"]}},
            "fields": "submitter_id,primary_site,samples.submitter_id,samples.sample_type",
            "format": "JSON", "size": size, "from": page * size,
            "sort": "submitter_id:asc"}
        d = post(GDC, payload)
        hits = d.get("data", {}).get("hits", [])
        if not hits:
            break
        for h in hits:
            case = h.get("submitter_id")
            site = h.get("primary_site")
            for s in (h.get("samples") or []):
                rows.append({"sample_barcode": s.get("submitter_id"),
                             "case_id": case, "primary_site": site,
                             "sample_type_gdc": s.get("sample_type")})
        total = d["data"]["pagination"]["total"]
        page += 1
        print(f"  cases fetched: {min(page * size, total)}/{total}")
        if page * size >= total:
            break

    df = pd.DataFrame(rows).drop_duplicates("sample_barcode")
    df.to_csv(TAB / "s6v2_gdc_primary_site.csv", index=False)
    print(f"\nsamples: {len(df)}  cases: {df['case_id'].nunique()}")
    print("\nprimary_site counts:")
    print(df["primary_site"].value_counts().to_string())
    print(f"\nwrote {TAB}/s6v2_gdc_primary_site.csv")


if __name__ == "__main__":
    main()
