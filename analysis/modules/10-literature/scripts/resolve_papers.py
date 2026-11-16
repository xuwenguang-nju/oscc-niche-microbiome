#!/usr/bin/env python3
"""Resolve the papers behind a list of secondary (blog) summaries to PubMed records.

Why: the analysis plan was built from Chinese blog-style literature summaries. Those give a
title, journal and impact factor but not methods, statistics or supplementary data, so every
claim taken from them needs to be traced to the primary record before it is used to justify a
new analysis. This does that resolution and reports, per paper, whether an open-access full
text exists.

Run: python3 resolve_papers.py
"""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "wisp-science/1.0 (paper resolution)"}

# (label, query) -- queries are exact title fragments taken from the blog summaries
PAPERS = [
    ("01_HNSC_Fusobacterium_TCMA",
     "Multi-Omics Integration Reveals the Crucial Role of Fusobacterium"),
    ("02_ESCC_prognosis",
     "Intratumoral microbiome impacts immune infiltrates in tumor microenvironment and "
     "predicts prognosis in esophageal squamous cell carcinoma"),
    ("03_HCC_microbial_subtype",
     "Intratumoral microbial heterogeneity affected tumor immune microenvironment and "
     "determined clinical outcome of HBV-related HCC"),
    ("04_melanoma_CD8",
     "Intratumour microbiome associated with the infiltration of cytotoxic CD8 T cells and "
     "patient survival in cutaneous melanoma"),
    ("05_LUAD_TMEindex",
     "Tumor microenvironment signatures enhances lung adenocarcinoma prognosis prediction"),
    ("06_HCC_MRS",
     "Identification of a brand intratumor microbiome signature for predicting prognosis of "
     "hepatocellular carcinoma"),
    ("07_CRC_CMS",
     "Host-microbiota interactions contributing to the heterogeneous tumor microenvironment "
     "in colorectal cancer"),
    ("08_PDAC_Pseudomonas_metabolome",
     "intratumoral microbiota Pseudomonas amino acid metabolism pancreatic cancer"),
    ("09_PDAC_MRS_TCGA_ICGC",
     "Association of tumor microbiome with survival in resected early-stage PDAC"),
    ("10_mesothelioma",
     "intratumoral microbiome mesothelioma survival"),
    ("11_breast_bile_acid",
     "Elevated bile acid metabolism and microbiome are associated with suppressed cell "
     "proliferation and better survival in breast cancer"),
    ("12_TNBC_metabolites",
     "Intra-tumoral microbial community profiling and associated metabolites alterations of "
     "triple negative breast cancer"),
    ("13_ICC_GutMicrobes",
     "Intratumor microbiome features reveal antitumor potentials of intrahepatic "
     "cholangiocarcinoma"),
    ("14_LUAD_glycolysis_lactate",
     "intratumoral microbiome glycolysis lactate lung adenocarcinoma immune heterogeneity"),
    ("16_CRC_RNAseq_reconstruction",
     "Location and condition based reconstruction of colon cancer microbiome from human "
     "RNA sequencing data"),
]


def http_json(url: str) -> dict:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def esearch(term: str, retmax: int = 5) -> list[str]:
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax,
                                "retmode": "json", "sort": "relevance"})
    return http_json(EUTILS + "esearch.fcgi?" + q).get(
        "esearchresult", {}).get("idlist", []) or []


def esummary(pmids: list[str]) -> dict:
    if not pmids:
        return {}
    q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids), "retmode": "json"})
    return http_json(EUTILS + "esummary.fcgi?" + q).get("result", {})


def pmc_id(pmid: str) -> str:
    q = urllib.parse.urlencode({"dbfrom": "pubmed", "db": "pmc", "id": pmid,
                                "retmode": "json"})
    try:
        j = http_json(EUTILS + "elink.fcgi?" + q)
        for ls in j.get("linksets", []):
            for db in ls.get("linksetdbs", []):
                if db.get("linkname", "").startswith("pubmed_pmc"):
                    ids = db.get("links", [])
                    if ids:
                        return f"PMC{ids[0]}"
    except Exception:
        pass
    return ""


def main() -> None:
    rows = []
    for label, query in PAPERS:
        ids = esearch(query)
        pmid = journal = year = doi = ttl = pmc = ""
        if ids:
            pmid = ids[0]
            rec = esummary([pmid]).get(pmid, {})
            ttl = (rec.get("title") or "").rstrip(".")
            journal = rec.get("source", "")
            year = (rec.get("pubdate", "") or "")[:4]
            for aid in rec.get("articleids", []):
                if aid.get("idtype") == "doi":
                    doi = aid.get("value", "")
            pmc = pmc_id(pmid)
            time.sleep(0.4)
        flag = "FOUND" if pmid else "  --  "
        print(f"{flag} {label:32} pmid={pmid or '-':9} {year or '----'} "
              f"{journal[:26]:26} {pmc or '-'}")
        if ttl:
            print(f"       {ttl[:108]}")
        rows.append({"label": label, "pmid": pmid, "journal": journal, "year": year,
                     "doi": doi, "pmc": pmc, "title": ttl, "query": query})
        time.sleep(0.4)

    with open("resolved_papers.json", "w") as fh:
        json.dump(rows, fh, indent=1)

    found = [r for r in rows if r["pmid"]]
    free = [r for r in rows if r["pmc"]]
    print(f"\nresolved {len(found)}/{len(rows)}; open-access (PMC): {len(free)}")


if __name__ == "__main__":
    main()
