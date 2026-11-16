"""Confirm (a) public availability + sample counts for the candidate OSCC
tissue-microbiome projects, and (b) metadata for the key OSCC intratumoral
microbiome papers the user should position against.

Output: analysis/out/oscc_key_paper_and_release_check.txt
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "out"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "wisp-science/1.0 (study design survey)"}
lines: list[str] = []


def get(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
        return r.read()


def esearch(db, term, retmax=5):
    u = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmax": retmax, "retmode": "json"})
    j = json.loads(get(u))["esearchresult"]
    return j["idlist"], j["count"]


def text(el):
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip() if el is not None else ""


# ---- A. BioSample counts per candidate project (proxy for sample availability) ----
PROJ = ["PRJNA1415268", "PRJNA1415359", "PRJNA1430204", "PRJEB115488",
        "PRJDB35665", "PRJNA1471312", "PRJEB98225", "PRJNA1445929",
        "PRJNA1482781", "PRJNA1475993", "PRJNA1522555", "PRJNA1507331"]
lines.append("=" * 110)
lines.append("A. Sample/run counts per candidate project (NCBI SRA index)")
for acc in PROJ:
    try:
        sra_ids, sra_n = esearch("sra", f"{acc}[BioProject]", 1)
        bs_ids, bs_n = esearch("biosample", f"{acc}[BioProject]", 1)
        lines.append(f"  {acc:>12} | SRA runs={sra_n:>5} | BioSamples={bs_n:>5}")
    except Exception as e:
        lines.append(f"  {acc:>12} | FAILED {type(e).__name__} {e}")
    time.sleep(0.4)

# ---- B. newest submissions of OSCC tissue 16S (catch anything I missed) ----
lines.append("\n" + "=" * 110)
lines.append("B. All BioProjects matching 'oral squamous cell carcinoma' + 16S / tissue microbiome")
for term in ['(oral squamous cell carcinoma[All Fields] OR OSCC[All Fields]) AND (16S[All Fields] OR microbiome[All Fields] OR microbiota[All Fields])',
             '(head and neck squamous cell carcinoma[All Fields] OR HNSCC[All Fields]) AND (tissue microbiome[All Fields] OR intratumoral[All Fields] OR 16S[All Fields])']:
    try:
        ids, n = esearch("bioproject", term, 60)
        u = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode({"db": "bioproject", "id": ",".join(ids), "retmode": "json"})
        j = json.loads(get(u))["result"]
        lines.append(f"\n  term={term[:70]}... hits={n}")
        got = [(j[k].get("project_acc", ""), j[k].get("project_title", "")) for k in j.get("uids", [])]
        for a, t in got:
            if re.search(r"oral|OSCC|HNSCC|head and neck|tongue|buccal", t, re.I):
                lines.append(f"    {a:>13} | {t[:120]}")
    except Exception as e:
        lines.append(f"    FAILED {type(e).__name__} {e}")
    time.sleep(0.4)

# ---- C. key OSCC intratumoral microbiome papers ----
PMIDS = ["37409975", "41059290", "42063001", "39873508", "33846399",
         "28855542", "36445134", "36338051", "41629856", "27259999",
         "42013847", "41714823", "35862975", "33218162", "38299619",
         "42795648", "41234847"]
lines.append("\n" + "=" * 110)
lines.append("C. Key OSCC / HNSCC intratumoral- and oral-microbiome papers")
data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(PMIDS), "retmode": "xml"}).encode()
with urllib.request.urlopen(urllib.request.Request(EUTILS + "efetch.fcgi", data=data, headers=UA), timeout=180) as r:
    root = ET.fromstring(r.read())
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
    accs = sorted(set(re.findall(r"\b(GSE\d{4,7}|PRJ[A-Z]{1,2}\d{4,8}|HRA\d{6}|CRA\d{5,8})\b", title + " " + ab)))
    lines.append(f"\n{pmid} | {year} | {journal} | {doi}")
    lines.append(f"  {title}")
    lines.append(f"  ACC in abstract: {accs}")
    lines.append(f"  ABS: {ab[:700]}")

(OUT / "oscc_key_paper_and_release_check.txt").write_text("\n".join(lines), encoding="utf-8")
print("lines:", len(lines))
