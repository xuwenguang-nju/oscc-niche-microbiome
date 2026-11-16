"""Attach first-author surnames to the curated table (used for readable
inline citations) and emit data/processed/curated_authors.json."""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results" / "tables" / "intratumoral_bacteria_bioinformatics_curated.csv"
UA = {"User-Agent": "wisp-science/1.0 (literature sweep)"}

rows = list(csv.DictReader(TAB.open(encoding="utf-8-sig")))
ids = [r["pmid"] for r in rows]
out: dict[str, str] = {}
for i in range(0, len(ids), 150):
    q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids[i:i + 150]), "retmode": "json"})
    with urllib.request.urlopen(
            urllib.request.Request(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + q, headers=UA),
            timeout=120) as r:
        j = json.loads(r.read())["result"]
    for pid in j["uids"]:
        d = j[pid]
        name = d.get("sortfirstauthor") or ""
        out[pid] = name
    time.sleep(0.4)

(ROOT / "data" / "processed" / "curated_authors.json").write_text(
    json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print("authors:", len(out))
print(list(out.items())[:8])
