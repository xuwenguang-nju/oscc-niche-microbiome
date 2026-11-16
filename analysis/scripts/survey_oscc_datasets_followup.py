"""Targeted follow-up: find the concrete public datasets that an OSCC
intratumoral-microbiome study could reuse.

Covers: NCBI BioProject/SRA/BioSample, ENA portal, Europe PMC full-text
accession mentions in OSCC/oral-cancer microbiome papers, and the specific
key papers' data availability.

Output: analysis/out/oscc_dataset_followup.txt
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "out"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "wisp-science/1.0 (study design survey)"}
lines: list[str] = []


def get(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
        return r.read()


def esearch(db: str, term: str, retmax: int = 40) -> list[str]:
    u = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmax": retmax, "retmode": "json", "sort": "relevance"})
    return json.loads(get(u))["esearchresult"]["idlist"]


def esummary(db: str, ids: list[str], fields: list[str]) -> list[dict]:
    out = []
    for i in range(0, len(ids), 100):
        u = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
            {"db": db, "id": ",".join(ids[i:i + 100]), "retmode": "json"})
        j = json.loads(get(u))["result"]
        for k in j.get("uids", []):
            d = j[k]
            out.append({f: d.get(f, "") for f in fields})
        time.sleep(0.35)
    return out


# ---------- 1. NCBI BioProject / SRA / BioSample ----------
NCBI = {
    "bioproject_oral_cancer_16s": ("bioproject", '(oral cancer OR oral squamous cell carcinoma OR oral cavity) AND (16S OR microbiome OR microbiota OR metagenome)'),
    "bioproject_oscc_metagenome": ("bioproject", '(oral squamous cell carcinoma OR oral cancer) AND (metagenome OR metagenomic OR bacteriome OR "microbial community")'),
    "sra_oral_16s": ("sra", '(oral cancer[All Fields] OR oral squamous cell carcinoma[All Fields]) AND (16S[All Fields] OR metagenomic[All Fields])'),
    "biosample_oscc_tumor_microbiome": ("biosample", '(oral squamous cell carcinoma[All Fields] OR "oral cancer"[All Fields]) AND (16S[All Fields] OR microbiome[All Fields])'),
}
for name, (db, term) in NCBI.items():
    try:
        ids = esearch(db, term, 40)
        lines.append(f"\n{'#'*100}\nNCBI {db} :: {name} :: hits={len(ids)}")
        if db == "bioproject":
            for d in esummary(db, ids, ["project_acc", "project_title", "project_data_type", "organism_name", "submission_date"]):
                lines.append(f'  {d["project_acc"]:>12} | {d["submission_date"][:10]:<10} | {d["project_data_type"][:30]:<30} | {d["project_title"][:110]}')
        else:
            for d in esummary(db, ids, ["title", "organism", "createdate"]):
                lines.append(f'  {d["createdate"][:10]:<10} | {d["organism"][:34]:<34} | {d["title"][:110]}')
    except Exception as e:
        lines.append(f"\nNCBI {db} :: {name} FAILED: {type(e).__name__} {e}")
    time.sleep(0.4)

# ---------- 2. ENA portal ----------
ENA = {
    "ena_oral_cancer_16s": '(study_title="*oral cancer*" AND study_title="*16S*")',
    "ena_oral_microbiome": '(study_title="*oral*" AND study_title="*microbiome*" AND study_title="*cancer*")',
}
for name, q in ENA.items():
    try:
        u = ("https://www.ebi.ac.uk/ena/portal/api/search?result=study&fields=study_accession,study_title,"
             "first_public,tax_id&format=json&limit=100&query=" + urllib.parse.quote(q))
        rows = json.loads(get(u))
        lines.append(f"\n{'#'*100}\nENA {name} :: hits={len(rows)}")
        for r in rows[:60]:
            lines.append(f'  {r.get("study_accession",""):>12} | {r.get("first_public","")[:10]:<10} | {r.get("study_title","")[:120]}')
    except Exception as e:
        lines.append(f"\nENA {name} FAILED: {type(e).__name__} {e}")
    time.sleep(0.4)

# ---------- 3. Europe PMC full-text accession mentions ----------
EPMC = {
    "epmc_oscc_16s_accessions": '("oral squamous cell carcinoma" OR "oral cancer") AND ("16S" OR "microbiome") AND (ACCESSION_ID:"PRJ" OR ACCESSION_ID:"GSE" OR ACCESSION_ID:"SRP" OR ACCESSION_ID:"HRA")',
    "epmc_oscc_intratumoral": '"oral squamous cell carcinoma" AND ("intratumoral microbiota" OR "intratumoral microbiome" OR "tumor microbiota")',
}
for name, q in EPMC.items():
    try:
        u = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + urllib.parse.quote(q) +
             "&format=json&pageSize=60&resultType=core")
        j = json.loads(get(u))
        res = j.get("resultList", {}).get("result", [])
        lines.append(f"\n{'#'*100}\nEuropePMC {name} :: total={j.get('hitCount')} shown={len(res)}")
        for r in res:
            accs = []
            for a in (r.get("accessionNumbers") or []):
                accs.append(str(a))
            lines.append(f'  {r.get("pmid",""):>9} | {r.get("pubYear","")} | {r.get("journalInfo",{}).get("journal",{}).get("medlineAbbreviation",""):<28} | {(r.get("title") or "")[:100]}  ACC:{",".join(accs[:8])}')
    except Exception as e:
        lines.append(f"\nEuropePMC {name} FAILED: {type(e).__name__} {e}")
    time.sleep(0.4)

(OUT / "oscc_dataset_followup.txt").write_text("\n".join(lines), encoding="utf-8")
print("written", sum(len(x) for x in lines))
