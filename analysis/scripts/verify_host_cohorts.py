"""Verify the host-side public cohort numbers quoted in the OSCC study design:
TCGA-HNSC case/sample counts (GDC API), CPTAC-HNSCC size, oral-microbiome GWAS
sample sizes, and OSCC GWAS resources.

Output: analysis/out/host_cohort_verification.txt
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "out"
UA = {"User-Agent": "wisp-science/1.0 (study design survey)"}
JSON_H = {"User-Agent": "wisp-science/1.0 (study design survey)", "Content-Type": "application/json"}
lines: list[str] = []


def post(url: str, payload: dict) -> dict:
    req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers=JSON_H, method="POST")
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def get(url: str):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
        return json.loads(r.read())


# ---- GDC: TCGA-HNSC counts ----
lines.append("=" * 100)
lines.append("GDC API :: TCGA-HNSC")
for prog in ("TCGA-HNSC",):
    try:
        cases = post("https://api.gdc.cancer.gov/cases",
                     {"filters": {"op": "in", "content": {"field": "project.project_id", "value": [prog]}},
                      "size": 0, "results_type": "count"})
        lines.append(f"  {prog}: cases = {cases['data']['pagination']['total']}")
    except Exception as e:
        lines.append(f"  {prog} cases FAILED {type(e).__name__} {e}")
    for dt in ("RNA-Seq", "WXS", "Methylation Array", "miRNA-Seq"):
        try:
            f = {"op": "and", "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": [prog]}},
                {"op": "in", "content": {"field": "data_type", "value": [dt]}},
                {"op": "in", "content": {"field": "access", "value": ["open"]}},
            ]}
            r = post("https://api.gdc.cancer.gov/files", {"filters": f, "size": 0, "results_type": "count"})
            lines.append(f"  {prog} open files {dt:<18} = {r['data']['pagination']['total']}")
        except Exception as e:
            lines.append(f"  {prog} {dt} FAILED {type(e).__name__} {e}")
        time.sleep(0.3)

# ---- PubMed: sizes of key resources ----
lines.append("\n" + "=" * 100)
lines.append("Key resource papers (title/abstract quoted from PubMed)")
import re
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PMIDS = ["33382980", "41650974", "42017663", "40901923", "42092351", "38299619",
         "39873508", "37409975", "41059290", "42063001", "33846399", "34873157",
         "34875929", "27749845", "41318233", "36769103"]
data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(PMIDS), "retmode": "xml"}).encode()
with urllib.request.urlopen(urllib.request.Request(EUTILS + "efetch.fcgi", data=data, headers=UA), timeout=180) as r:
    root = ET.fromstring(r.read())


def text(el):
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip() if el is not None else ""


for art in root.findall(".//PubmedArticle"):
    pmid = text(art.find(".//MedlineCitation/PMID"))
    title = text(art.find(".//ArticleTitle"))
    journal = text(art.find(".//Journal/ISOAbbreviation"))
    year = text(art.find(".//JournalIssue/PubDate/Year")) or text(art.find(".//JournalIssue/PubDate/MedlineDate"))[:4]
    doi = ""
    aidl = art.find("PubmedData/ArticleIdList")
    if aidl is not None:
        for aid in aidl.findall("ArticleId"):
            if aid.get("IdType") == "doi":
                doi = text(aid)
    ab = " ".join(text(a) for a in art.findall(".//Abstract/AbstractText"))
    lines.append(f"\n{pmid} | {year} | {journal} | {doi}")
    lines.append(f"  {title}")
    lines.append(f"  ABS: {ab[:1100]}")

(OUT / "host_cohort_verification.txt").write_text("\n".join(lines), encoding="utf-8")
print("ok")
