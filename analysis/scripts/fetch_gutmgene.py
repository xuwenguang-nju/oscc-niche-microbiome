#!/usr/bin/env python3
"""下载 gutMGene v2.0 的六张关联表（人/鼠 × 菌-代谢物、菌-基因、代谢物-基因）。

来源：https://bio-computing.hrbmu.edu.cn/gutMGene2.0_api/search/searchdetaildata
（POST JSON，payload 由前端 search 视图构建，字段名取自其 chunk：dataset/species/gutMicrobe/gene/
 metaboliteNames/metaboliteChEBIOntology/subMetanameslist/datatype）

原始返回原样保存到 data/raw/gutmgene/（不改动），并转成 UTF-8 CSV 到 results/tables/。
用法：python analysis/scripts/fetch_gutmgene.py
"""
from __future__ import annotations

import csv
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

API = "https://bio-computing.hrbmu.edu.cn/gutMGene2.0_api/search/searchdetaildata"
RAW = Path("data/raw/gutmgene")
TAB = Path("results/tables")
DATASETS = ["microbe_metabolite", "microbe_gene", "metabolite_gene"]
SPECIES = ["human", "mouse"]


def fetch(dataset: str, species: str) -> dict:
    body = json.dumps(
        {
            "data": {
                "dataset": dataset,
                "species": species,
                "gutMicrobe": "",
                "gene": "",
                "metaboliteNames": "",
                "metaboliteChEBIOntology": "",
                "subMetanameslist": "",
                "datatype": "",
            }
        }
    ).encode()
    req = urllib.request.Request(
        API,
        data=body,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "Mozilla/5.0",
            "Origin": "https://bio-computing.hrbmu.edu.cn",
            "Referer": "https://bio-computing.hrbmu.edu.cn/gutmgene/",
        },
    )
    with urllib.request.urlopen(req, timeout=120) as fh:
        return json.loads(fh.read())


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    summary = []
    for species in SPECIES:
        for dataset in DATASETS:
            payload = fetch(dataset, species)
            rows = payload.get("resultSet", [])
            raw_path = RAW / f"gutmgene_v2_{species}_{dataset}.json"
            raw_path.write_text(
                json.dumps({"retrieved_utc": stamp, "api": API, "dataset": dataset,
                            "species": species, "payload": payload},
                           ensure_ascii=False, indent=1),
                encoding="utf-8",
            )
            if rows:
                cols = list(rows[0].keys())
                for r in rows:
                    for c in cols:
                        r.setdefault(c, None)
                csv_path = TAB / f"gutmgene_v2_{species}_{dataset}.csv"
                with csv_path.open("w", newline="", encoding="utf-8-sig") as fh:
                    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
                    w.writeheader()
                    w.writerows(rows)
            summary.append((species, dataset, len(rows), raw_path.name))
            print(f"{species:<6} {dataset:<20} rows={len(rows):<5} -> {raw_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
