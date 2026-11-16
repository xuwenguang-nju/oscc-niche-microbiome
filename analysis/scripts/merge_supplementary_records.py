"""Merge supplementary hits (found in targeted follow-up queries) into the
main PubMed record file so the curated table has a single source of truth."""

from __future__ import annotations

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
UA = {"User-Agent": "wisp-science/1.0 (literature sweep)"}

EXTRA = [
    "37811944", "40901921", "40901923", "37555750", "36179670", "38599211",
    "42092351", "40558028", "38195165", "42446235", "39982111", "41289425",
    "38556382", "36732743", "28410234", "39355268", "41823382", "40575164",
    "41480207", "41387030", "41552720", "42687557", "37292990", "38259095",
    "41647993", "41932249", "36937384", "28118849", "31639030", "39896561",
    "37074188", "35586863", "37465100", "33659025", "34422640", "34136388",
    "34741726", "36999009", "37188204", "41715092", "42482748", "42865545",
    "41608560", "40908816", "42623381", "40492114", "41650974", "39982111",
    "37074188", "36999042", "38421176", "38868430", "36884119", "36834564",
    "39173404", "39660865", "40606983", "40433663", "38515713", "39117066",
    "42436233", "39915534", "41142936", "40761254", "42491956", "40013793",
    "39565103", "40875911", "41373739", "40865517", "35782745", "35610202",
    "35972351", "38090530", "41438381", "39833403", "40475759", "39907412",
    "41395468", "41234866", "40029433", "40443528", "42729228", "40952000",
    "40677704", "40202312", "37260411", "35822725", "35122129", "36555172",
    "36362038", "37895136", "39497974", "38722731", "39072166", "36786568",
    "38073489", "38808774", "38442004", "38375207", "38721596", "39397906",
    "39531992", "41142936", "32214244", "33382980", "32257046", "26295344",
    "42017663", "41714823", "40314439", "38396294", "37131219", "37850841",
    "36189892", "42013847", "41885500", "42488410", "42560302", "42399737",
    "36702066", "40606983",
]


def text(el) -> str:
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip() if el is not None else ""


recs = json.loads(JSON.read_text(encoding="utf-8"))
have = {r["pmid"] for r in recs}
todo = [p for p in dict.fromkeys(EXTRA) if p not in have]
print("to fetch:", len(todo))

new = []
for i in range(0, len(todo), 100):
    batch = todo[i : i + 100]
    data = urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}
    ).encode()
    req = urllib.request.Request(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        data=data,
        headers=UA,
    )
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
        new.append(
            {
                "pmid": pmid,
                "year": text(art.find(".//JournalIssue/PubDate/Year"))
                or text(art.find(".//JournalIssue/PubDate/MedlineDate"))[:4],
                "journal": text(art.find(".//Journal/ISOAbbreviation"))
                or text(art.find(".//Journal/Title")),
                "title": text(art.find(".//ArticleTitle")),
                "abstract": " ".join(
                    text(a) for a in art.findall(".//Abstract/AbstractText")
                ),
                "doi": doi,
                "pubtypes": [text(p) for p in art.findall(".//PublicationType")],
                "mesh": [text(m) for m in art.findall(".//MeshHeading/DescriptorName")],
                "country": "",
                "queries": ["supplementary"],
            }
        )
    time.sleep(0.4)

print("new records:", len(new), "missing DOI:", sum(1 for r in new if not r["doi"]))
recs.extend(new)
JSON.write_text(json.dumps(recs, ensure_ascii=False, indent=1), encoding="utf-8")
print("total:", len(recs))
