"""Round 3: sample counts + full descriptions for the tissue-oriented OSCC /
HNSCC microbiome projects discovered in round 2, so the study design can cite
concrete reused datasets.

Output: analysis/out/oscc_tissue_projects_detail.txt
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

PROJ = [
    "PRJNA1328350", "PRJNA822685", "PRJNA1203737", "PRJNA803155", "PRJNA866676",
    "PRJNA412445", "PRJNA415963", "PRJNA528843", "PRJEB87426", "PRJNA1321769",
    "PRJNA775681", "PRJNA775796", "PRJNA1366983", "PRJDB16847", "PRJNA1429063",
    "PRJNA947509", "PRJNA811533", "PRJNA870048", "PRJEB75386", "PRJNA1117634",
    "PRJNA666746", "PRJNA666891", "PRJNA597251", "PRJNA762320", "PRJNA735337",
    "PRJNA700849", "PRJNA647515", "PRJNA1246956", "PRJNA997108", "PRJEB75386",
]


def get(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
        return r.read()


def esearch(db, term, retmax=3):
    u = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmax": retmax, "retmode": "json"})
    j = json.loads(get(u))["esearchresult"]
    return j["idlist"], j["count"]


lines = ["=" * 115, "Tissue-oriented OSCC / HNSCC microbiome projects: size, type, description", "=" * 115]
for acc in dict.fromkeys(PROJ):
    try:
        sra_ids, sra_n = esearch("sra", f"{acc}[BioProject]", 1)
        ids, _ = esearch("bioproject", f"{acc}[Project Accession]", 1)
        title = desc = dtype = ""
        if ids:
            u = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
                {"db": "bioproject", "id": ",".join(ids), "retmode": "json"})
            j = json.loads(get(u))["result"]
            for k in j.get("uids", []):
                d = j[k]
                title = d.get("project_title", "")
                dtype = d.get("project_data_type", "")
                desk = d.get("project_description", "") or ""
                desk = re.sub(r"<[^>]+>", " ", desk.replace("&lt;", "<").replace("&gt;", ">"))
                desc = re.sub(r"\s+", " ", desk).strip()
        lines.append(f"\n{acc:>13} | runs={sra_n:>5} | {dtype[:26]:<26} | {title[:105]}")
        if desc:
            lines.append(f"    {desc[:1100]}")
    except Exception as e:
        lines.append(f"\n{acc:>13} | FAILED {type(e).__name__} {e}")
    time.sleep(0.4)

(OUT / "oscc_tissue_projects_detail.txt").write_text("\n".join(lines), encoding="utf-8")
print("ok", len(lines))
