"""Survey public resources that could support an OSCC intratumoral-microbiome
study: (1) GEO/SRA/BioProject datasets containing oral-cancer 16S or
metagenomic data, (2) PubMed papers on OSCC/oral-cancer microbiome that report
public data accessions, (3) bulk/single-cell OSCC host cohorts.

Output: analysis/out/oscc_public_resource_survey.json + .txt
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
OUT.mkdir(parents=True, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "wisp-science/1.0 (study design survey)"}


def get(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
        return r.read()


def esearch(db: str, term: str, retmax: int = 60) -> list[str]:
    u = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmax": retmax, "retmode": "json", "sort": "relevance"})
    return json.loads(get(u))["esearchresult"]["idlist"]


def esummary(db: str, ids: list[str]) -> list[dict]:
    out = []
    for i in range(0, len(ids), 100):
        u = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
            {"db": db, "id": ",".join(ids[i:i + 100]), "retmode": "json"})
        j = json.loads(get(u))["result"]
        for k in j.get("uids", []):
            out.append(j[k])
        time.sleep(0.35)
    return out


def text(el) -> str:
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip() if el is not None else ""


GEO_QUERIES = {
    "geo_oscc_microbiome": '(oral squamous cell carcinoma OR oral cancer OR head and neck squamous cell carcinoma) AND (microbiome OR microbiota OR 16S OR metagenom*) AND "Homo sapiens"[Organism]',
    "geo_oscc_16s_strict": '("oral squamous cell carcinoma"[All Fields] OR "oral cancer"[All Fields]) AND ("16S"[All Fields] OR "16S rRNA"[All Fields] OR "metagenome"[All Fields] OR "metagenomic"[All Fields])',
    "geo_hnsc_microbiome": '(head and neck[All Fields]) AND (microbiome[All Fields] OR microbiota[All Fields] OR bacteriome[All Fields]) AND (RNA-seq[All Fields] OR sequencing[All Fields])',
    "geo_oscc_scrna": '(oral squamous cell carcinoma[All Fields]) AND ("single cell"[All Fields] OR scRNA[All Fields])',
    "geo_oscc_spatial": '(oral squamous cell carcinoma[All Fields] OR oral cancer[All Fields]) AND (spatial transcriptomics[All Fields] OR "Visium"[All Fields] OR "10x Visium"[All Fields])',
    "bioproject_oscc_microbiome": '(oral squamous cell carcinoma[All Fields] OR oral cancer[All Fields]) AND (microbiome[All Fields] OR microbiota[All Fields] OR 16S[All Fields])',
    "sra_oscc_16s": '("oral squamous cell carcinoma"[All Fields] OR "oral cancer"[All Fields]) AND ("16S ribosomal RNA"[All Fields] OR "16S"[All Fields])',
    "bioproject_oral_16s": '(oral cavity[All Fields] OR oral[All Fields]) AND ("16S ribosomal RNA"[All Fields]) AND (carcinoma[All Fields] OR cancer[All Fields] OR tumor[All Fields])',
}

PUBMED_QUERIES = {
    "pm_oscc_16s": '("oral squamous cell carcinoma"[tiab] OR "oral cancer"[tiab] OR "oral cavity"[tiab]) AND ("16S"[tiab] OR "metagenom*"[tiab] OR "microbiome"[tiab] OR "microbiota"[tiab])',
    "pm_oscc_true_intratumoral": '("oral squamous cell carcinoma"[tiab] OR "oral cancer"[tiab] OR "head and neck squamous cell carcinoma"[tiab]) AND ("intratumoral microb*"[tiab] OR "intratumor microb*"[tiab] OR "tumor microbiome"[tiab] OR "tumor microbiota"[tiab] OR "tissue microb*"[tiab])',
    "pm_oscc_microbiome_prognosis": '("oral squamous cell carcinoma"[tiab] OR "oral cancer"[tiab] OR "head and neck"[tiab]) AND (microbiome[tiab] OR microbiota[tiab]) AND (prognos*[tiab] OR survival[tiab] OR signature[tiab] OR "risk model"[tiab] OR prognostic[tiab])',
}

ACC = re.compile(r"\b(GSE\d{4,7}|PRJNA\d{5,8}|PRJEB\d{4,7}|SRP\d{5,7}|ERP\d{5,7}|HRA\d{6}|OEP\d{6}|CRA\d{5,8}|EGAS\d{8}|phs\d{6}|SRR\d{5,9})\b")

report: dict = {"geo": {}, "pubmed": {}, "accessions": {}}
lines: list[str] = []

for name, term in GEO_QUERIES.items():
    try:
        ids = esearch("gds", term, 60)
        det = esummary("gds", ids)
        items = []
        for d in det:
            items.append({
                "id": d.get("accession") or d.get("gds", ""),
                "title": d.get("title", ""),
                "taxon": d.get("taxon", ""),
                "gdstype": d.get("gdstype", ""),
                "n_samples": d.get("n_samples", ""),
                "pdat": d.get("pdat", ""),
                "summary": (d.get("summary") or "")[:600],
            })
        report["geo"][name] = items
        lines.append(f"\n{'='*100}\nGEO/{name}  hits={len(ids)}\n")
        for it in items:
            lines.append(f'  {it["id"]:>12} | {it["pdat"]:<12} | n={it["n_samples"]:<4} | {it["gdstype"]:<28} | {it["title"][:120]}')
    except Exception as e:
        lines.append(f"\nGEO/{name} FAILED: {type(e).__name__} {e}")
    time.sleep(0.4)

for name, term in PUBMED_QUERIES.items():
    ids = esearch("pubmed", term, 80)
    recs = []
    for i in range(0, len(ids), 100):
        data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids[i:i + 100]), "retmode": "xml"}).encode()
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
            alltext = title + " " + ab
            accs = sorted(set(ACC.findall(alltext)))
            # only cell-line/animal/16S-of-tumor relevance filter
            recs.append({"pmid": pmid, "year": year, "journal": journal, "doi": doi,
                         "title": title, "accessions": accs, "abstract": ab})
            for a in accs:
                report["accessions"].setdefault(a, []).append(pmid)
        time.sleep(0.4)
    report["pubmed"][name] = recs
    lines.append(f"\n{'='*100}\nPubMed/{name}  hits={len(ids)}\n")
    for r in recs:
        acc = ("  ACC:" + ",".join(r["accessions"])) if r["accessions"] else ""
        lines.append(f'  {r["pmid"]} | {r["year"]} | {r["journal"]:<26} | {r["title"][:110]}{acc}')

(OUT / "oscc_public_resource_survey.json").write_text(json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")
(OUT / "oscc_public_resource_survey.txt").write_text("\n".join(lines), encoding="utf-8")
print("\n".join(lines))
