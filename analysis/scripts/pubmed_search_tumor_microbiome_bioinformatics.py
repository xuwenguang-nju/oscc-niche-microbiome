"""PubMed sweep for bioinformatics analyses of intratumoral bacteria / tumor
microbiome that reuse public data (TCGA, GEO, CPTAC, ICGC) and are primarily
dry-lab, optionally with minor wet-lab validation.

Outputs
-------
data/processed/tumor_microbiome_pubmed_records.json  full parsed records
data/processed/tumor_microbiome_pubmed_records.csv   flat table for screening
analysis/out/pubmed_query_counts.json                per-query hit counts
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
OUT.mkdir(parents=True, exist_ok=True)
(ROOT / "analysis" / "out").mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "wisp-science/1.0 (literature sweep)"}

QUERIES = {
    "A_intratumoral_bacteria": (
        '("intratumoral bacteria"[tiab] OR "intratumoural bacteria"[tiab] OR '
        '"intratumor bacteria"[tiab] OR "intratumoral microbiome"[tiab] OR '
        '"intratumoural microbiome"[tiab] OR "intratumoral microbiota"[tiab] OR '
        '"intratumor microbiota"[tiab] OR "tumor microbiome"[tiab] OR '
        '"tumour microbiome"[tiab] OR "tumor microbiota"[tiab] OR '
        '"cancer microbiome"[tiab] OR "tumor-resident bacteria"[tiab] OR '
        '"tissue-resident microbiota"[tiab] OR "intratumoral bacterial"[tiab])'
    ),
    "B_tcga_based": (
        '(TCGA[tiab] OR "The Cancer Genome Atlas"[tiab] OR CPTAC[tiab] OR '
        'ICGC[tiab] OR "public data"[tiab] OR "publicly available"[tiab]) AND '
        '(microbiome[tiab] OR microbiota[tiab] OR bacteri*[tiab] OR microb*[tiab]) '
        'AND (tumor[tiab] OR tumour[tiab] OR cancer[tiab])'
    ),
    "C_prognostic_signature": (
        '(microbiome[tiab] OR microbiota[tiab] OR bacteri*[tiab]) AND '
        '(TCGA[tiab] OR GEO[tiab] OR "Gene Expression Omnibus"[tiab]) AND '
        '(signature[tiab] OR "risk model"[tiab] OR "prognostic model"[tiab] OR '
        'prognosis[tiab] OR nomogram[tiab] OR "machine learning"[tiab] OR '
        '"deep learning"[tiab] OR "risk score"[tiab] OR biomarker[tiab]) AND '
        '(tumor[tiab] OR tumour[tiab] OR cancer[tiab])'
    ),
    "D_deconvolution": (
        '(bacterial[tiab] OR bacteria[tiab] OR microbial[tiab] OR microbiome[tiab]) AND '
        '(deconvolution[tiab] OR "read-based"[tiab] OR "in silico"[tiab] OR '
        'reanalysis[tiab] OR "re-analysis"[tiab] OR "metagenomic"[tiab] OR '
        '"microbial reads"[tiab] OR Kraken[tiab] OR "Centrifuge"[tiab] OR '
        '"PathSeq"[tiab] OR "GATK"[tiab] OR "Kraken2"[tiab]) AND '
        '(tumor[tiab] OR tumour[tiab] OR cancer[tiab] OR carcinoma[tiab])'
    ),
    "E_singlecell_spatial": (
        '(single-cell[tiab] OR "single cell"[tiab] OR spatial[tiab] OR '
        '"spatial transcriptomics"[tiab]) AND (bacteri*[tiab] OR microb*[tiab] OR '
        'microbiome[tiab] OR microbiota[tiab]) AND '
        '(tumor[tiab] OR tumour[tiab] OR cancer[tiab])'
    ),
    "F_pancancer_landscape": (
        '(pan-cancer[tiab] OR "pan cancer"[tiab] OR landscape[tiab]) AND '
        '(microbiome[tiab] OR microbiota[tiab] OR bacteri*[tiab]) AND '
        '(TCGA[tiab] OR "The Cancer Genome Atlas"[tiab] OR "public"[tiab])'
    ),
    "G_bacterial_associated_genes": (
        '(bacteri*[tiab] OR microbiome[tiab] OR microbiota[tiab]) AND '
        '("gene expression"[tiab] OR transcriptom*[tiab] OR "RNA sequencing"[tiab] OR '
        'RNA-seq[tiab] OR "WES"[tiab] OR "whole-exome"[tiab] OR "whole genome"[tiab]) AND '
        '(tumor[tiab] OR tumour[tiab] OR cancer[tiab]) AND '
        '(correlat*[tiab] OR associat*[tiab] OR interaction[tiab] OR crosstalk[tiab])'
    ),
}

MESH_MAJOR_KEEP = re.compile(
    r"neoplasms|microbiota|bacteria|bacterial infections|biomarkers|"
    r"computational biology|gene expression profiling|metagenomics",
    re.I,
)


def eutil(path: str, params: dict, html: bool = False) -> bytes:
    params = dict(params)
    url = EUTILS + path + "?" + urllib.parse.urlencode(params)
    data = urllib.parse.urlencode(params).encode() if html else None
    req = urllib.request.Request(url, data=data, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def esearch(term: str, retmax: int = 500) -> list[str]:
    out: list[str] = []
    retstart = 0
    while len(out) < retmax:
        raw = eutil(
            "esearch.fcgi",
            {
                "db": "pubmed",
                "term": term,
                "retmax": min(200, retmax - len(out)),
                "retstart": retstart,
                "retmode": "json",
                "sort": "relevance",
            },
        )
        js = json.loads(raw.decode())
        ids = js["esearchresult"]["idlist"]
        out.extend(ids)
        if len(ids) < 200:
            break
        retstart += 200
        time.sleep(0.4)
    return out


def text(el) -> str:
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip() if el is not None else ""


def parse_articles(xml_bytes: bytes) -> list[dict]:
    root = ET.fromstring(xml_bytes)
    recs = []
    for art in root.findall(".//PubmedArticle"):
        pmid = text(art.find(".//PMID"))
        title = text(art.find(".//ArticleTitle"))
        journal = text(art.find(".//Journal/Title"))
        iso = text(art.find(".//Journal/ISOAbbreviation"))
        year = (
            text(art.find(".//JournalIssue/PubDate/Year"))
            or (text(art.find(".//JournalIssue/PubDate/MedlineDate"))[:4])
        )
        doi = ""
        for aid in art.findall(".//ArticleId"):
            if aid.get("IdType") == "doi":
                doi = text(aid)
        abstract = " ".join(
            text(a) for a in art.findall(".//Abstract/AbstractText")
        )
        pts = [text(p) for p in art.findall(".//PublicationType")]
        mesh = [text(m) for m in art.findall(".//MeshHeading/DescriptorName")]
        affs = [text(a) for a in art.findall(".//AffiliationInfo/Affiliation")]
        recs.append(
            {
                "pmid": pmid,
                "year": year,
                "journal": iso or journal,
                "title": title,
                "abstract": abstract,
                "doi": doi,
                "pubtypes": pts,
                "mesh": mesh,
                "country": affs[0].rsplit(",", 1)[-1].strip() if affs else "",
            }
        )
    return recs


def fetch(ids: list[str]) -> list[dict]:
    recs = []
    for i in range(0, len(ids), 150):
        batch = ids[i : i + 150]
        data = urllib.parse.urlencode(
            {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}
        ).encode()
        req = urllib.request.Request(EUTILS + "efetch.fcgi", data=data, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as r:
            recs.extend(parse_articles(r.read()))
        time.sleep(0.4)
    return recs


def main() -> None:
    counts = {}
    idmap: dict[str, set[str]] = {}
    for name, term in QUERIES.items():
        ids = esearch(term, retmax=600)
        counts[name] = len(ids)
        print(f"{name}: {len(ids)}")
        for pid in ids:
            idmap.setdefault(pid, set()).add(name)
        time.sleep(0.4)

    allids = list(idmap)
    print("union:", len(allids))
    recs = fetch(allids)
    for r in recs:
        r["queries"] = sorted(idmap.get(r["pmid"], set()))
    recs.sort(key=lambda r: (r["year"] or "", r["title"]))

    (OUT / "tumor_microbiome_pubmed_records.json").write_text(
        json.dumps(recs, ensure_ascii=False, indent=1), encoding="utf-8"
    )
    with (OUT / "tumor_microbiome_pubmed_records.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "pmid", "year", "journal", "doi", "title", "pubtypes",
                "mesh", "country", "queries", "abstract",
            ],
            extrasaction="ignore",
        )
        w.writeheader()
        for r in recs:
            rr = dict(r)
            rr["pubtypes"] = "; ".join(r["pubtypes"])
            rr["mesh"] = "; ".join(r["mesh"])
            rr["queries"] = "; ".join(r["queries"])
            w.writerow(rr)

    (ROOT / "analysis" / "out" / "pubmed_query_counts.json").write_text(
        json.dumps({"counts": counts, "union": len(allids)}, indent=1),
        encoding="utf-8",
    )
    print("records:", len(recs))


if __name__ == "__main__":
    main()
