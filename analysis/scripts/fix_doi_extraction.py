"""Fix DOI extraction for the PubMed sweep records.

The first pass used .//ArticleId, which also matches DOIs inside the
reference list. The article DOI lives under PubmedData/ArticleIdList/ArticleId.
Re-fetch XML and repair data/processed/tumor_microbiome_pubmed_records.json
plus the derived CSV.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "processed"
JSON = OUT / "tumor_microbiome_pubmed_records.json"
CSV = OUT / "tumor_microbiome_pubmed_records.csv"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA = {"User-Agent": "wisp-science/1.0 (literature sweep)"}


def text(el) -> str:
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip() if el is not None else ""


recs = json.loads(JSON.read_text(encoding="utf-8"))
byid = {r["pmid"]: r for r in recs}
ids = list(byid)
fixed: dict[str, str] = {}

for i in range(0, len(ids), 150):
    batch = ids[i : i + 150]
    data = urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}
    ).encode()
    req = urllib.request.Request(EUTILS, data=data, headers=UA)
    with urllib.request.urlopen(req, timeout=180) as r:
        root = ET.fromstring(r.read())
    for art in root.findall(".//PubmedArticle"):
        pmid = text(art.find(".//MedlineCitation/PMID"))
        doi = ""
        aidl = art.find("PubmedData/ArticleIdList")
        if aidl is not None:
            for aid in aidl.findall("ArticleId"):
                if aid.get("IdType") == "doi":
                    doi = text(aid)
        fixed[pmid] = doi
    time.sleep(0.4)

changed = 0
for r in recs:
    new = fixed.get(r["pmid"])
    if new is not None and new != r["doi"]:
        changed += 1
        r["doi"] = new
print("records:", len(recs), "doi corrected:", changed, "missing DOI:", sum(1 for r in recs if not r["doi"]))

JSON.write_text(json.dumps(recs, ensure_ascii=False, indent=1), encoding="utf-8")
with CSV.open("w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(
        fh,
        fieldnames=["pmid", "year", "journal", "doi", "title", "pubtypes", "mesh", "country", "queries", "abstract"],
        extrasaction="ignore",
    )
    w.writeheader()
    for r in recs:
        rr = dict(r)
        rr["pubtypes"] = "; ".join(r["pubtypes"])
        rr["mesh"] = "; ".join(r["mesh"])
        rr["queries"] = "; ".join(r["queries"])
        w.writerow(rr)
