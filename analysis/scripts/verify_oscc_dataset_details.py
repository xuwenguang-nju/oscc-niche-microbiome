"""Verify the specific public datasets an OSCC intratumoral-microbiome study
would reuse: sample counts / sample type / paired design for the candidate
oral-cancer 16S & metagenome projects; OSCC host-expression cohorts with
survival; oral-microbiome GWAS resources for Mendelian randomisation.

Output: analysis/out/oscc_dataset_details.txt
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


def jget(url: str):
    return json.loads(get(url))


def esearch(db, term, retmax=40):
    u = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmax": retmax, "retmode": "json", "sort": "relevance"})
    return jget(u)["esearchresult"]["idlist"]


def esummary(db, ids):
    out = []
    for i in range(0, len(ids), 100):
        u = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
            {"db": db, "id": ",".join(ids[i:i + 100]), "retmode": "json"})
        j = jget(u)["result"]
        for k in j.get("uids", []):
            out.append(j[k])
        time.sleep(0.35)
    return out


# ---------- A. candidate oral-cancer microbiome projects ----------
CAND = ["PRJNA1415268", "PRJNA1415359", "PRJNA1430204", "PRJNA1471312",
        "PRJNA1445929", "PRJNA1475993", "PRJNA1482781", "PRJNA1456678",
        "PRJNA1522555", "PRJNA1507331", "PRJNA1495532"]
lines.append("=" * 110)
lines.append("A. NCBI BioProject detail for candidate oral-cancer microbiome projects")
for acc in CAND:
    try:
        ids = esearch("bioproject", f"{acc}[Project Accession]")
        if not ids:
            lines.append(f"\n{acc}: NOT FOUND")
            continue
        d = esummary("bioproject", ids)[0]
        lines.append(f"\n{acc}  ({d.get('submission_date','')})")
        lines.append(f"  title    : {d.get('project_title','')}")
        lines.append(f"  data type: {d.get('project_data_type','')} | organism: {d.get('organism_name','')}")
        desc = re.sub(r"<[^>]+>", " ", d.get("project_description") or "")
        desc = re.sub(r"\s+", " ", desc).strip()
        if desc:
            lines.append(f"  desc     : {desc[:900]}")
        for rel in ("project_relevance", "relevance"):
            if d.get(rel):
                lines.append(f"  {rel}: {d[rel]}")
    except Exception as e:
        lines.append(f"\n{acc} FAILED: {type(e).__name__} {e}")
    time.sleep(0.4)

# ---------- B. ENA run counts + sample metadata for key studies ----------
lines.append("\n" + "=" * 110)
lines.append("B. ENA run counts / tissue source for key studies")
ENASTUDIES = ["PRJNA1415268", "PRJNA1415359", "PRJNA1430204", "PRJNA1471312",
              "PRJEB115488", "PRJEB98225", "PRJNA1445929", "PRJDB35665",
              "PRJDB29370", "PRJDB29373"]
for acc in ENASTUDIES:
    try:
        u = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + acc +
             "&result=read_run&fields=run_accession,library_strategy,library_source,scientific_name,"
             "description,instrument_platform,read_count,base_count,sample_title,sample_alias"
             "&format=json&limit=0")
        rows = json.loads(get(u))
        title = ""
        try:
            t = jget("https://www.ebi.ac.uk/ena/portal/api/search?result=study&fields=study_title,center_name,first_public"
                     "&format=json&limit=1&query=" + urllib.parse.quote(f'study_accession="{acc}"'))
            if t:
                title = f'{t[0].get("study_title","")} | {t[0].get("center_name","")} | {t[0].get("first_public","")}'
        except Exception:
            pass
        lines.append(f"\n{acc}  runs={len(rows)}")
        lines.append(f"  {title}")
        if rows:
            strat = {}
            src = {}
            descs = []
            for r in rows:
                strat[r.get("library_strategy", "")] = strat.get(r.get("library_strategy", ""), 0) + 1
                src[r.get("library_source", "")] = src.get(r.get("library_source", ""), 0) + 1
                if r.get("description") and len(descs) < 25:
                    descs.append(r["description"])
            lines.append(f"  strategy={strat} source={src}")
            if descs:
                lines.append("  sample descriptions (up to 25): " + " || ".join(descs))
    except Exception as e:
        lines.append(f"\n{acc} ENA FAILED: {type(e).__name__} {e}")
    time.sleep(0.4)

# ---------- C. OSCC host-expression cohorts in GEO with decent n ----------
lines.append("\n" + "=" * 110)
lines.append("C. GEO OSCC bulk expression cohorts (n>=40, human, tissue)")
try:
    ids = esearch("gds", '(oral squamous cell carcinoma OR oral cancer) AND "Homo sapiens"[Organism] AND '
                         '(expression profiling by array OR expression profiling by high throughput sequencing) '
                         'AND gse[Entry Type]', 80)
    det = esummary("gds", ids)
    rows = [{"acc": d.get("accession", ""), "n": d.get("n_samples", 0), "pdat": d.get("pdat", ""),
             "type": d.get("gdstype", ""), "title": d.get("title", "")} for d in det]
    rows = [r for r in rows if isinstance(r["n"], int) and r["n"] >= 40]
    rows.sort(key=lambda r: -r["n"])
    for r in rows:
        lines.append(f'  {r["acc"]:>12} | n={r["n"]:<5} | {r["pdat"]:<12} | {r["type"][:44]:<44} | {r["title"][:110]}')
except Exception as e:
    lines.append(f"  FAILED: {type(e).__name__} {e}")

# ---------- D. oral microbiome GWAS / MR resources ----------
lines.append("\n" + "=" * 110)
lines.append("D. Oral-microbiome GWAS and Mendelian-randomisation resources")
for name, term in {
    "gwas_oral_microbiome": '(oral microbiome[tiab] OR oral microbiota[tiab] OR salivary microbiome[tiab]) AND (genome-wide association[tiab] OR GWAS[tiab] OR heritab*[tiab])',
    "mr_oral_microbiome_cancer": '(Mendelian randomi*[tiab]) AND (oral microb*[tiab] OR salivary microb*[tiab] OR oral bacteri*[tiab] OR microbiome[tiab]) AND (oral cancer[tiab] OR head and neck[tiab] OR cancer[tiab])',
    "gwas_oscc": '(oral squamous cell carcinoma[tiab] OR oral cancer[tiab]) AND (genome-wide association[tiab] OR GWAS[tiab])',
}.items():
    try:
        ids = esearch("pubmed", term, 25)
        det = esummary("pubmed", ids)
        lines.append(f"\n  --- {name} hits={len(ids)}")
        for d in det:
            doi = ""
            for a in d.get("articleids", []):
                if a.get("idtype") == "doi":
                    doi = a.get("value", "")
            lines.append(f'  {d.get("uid",""):>9} | {d.get("pubdate","")[:4]} | {d.get("source",""):<26} | {doi:<34} | {(d.get("title") or "")[:100]}')
    except Exception as e:
        lines.append(f"  {name} FAILED: {type(e).__name__} {e}")
    time.sleep(0.4)

(OUT / "oscc_dataset_details.txt").write_text("\n".join(lines), encoding="utf-8")
print("lines:", len(lines))
