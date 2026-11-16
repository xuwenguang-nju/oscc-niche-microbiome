"""Build the curated deliverable table for the request:
"PubMed 上关于瘤内菌的生信分析（尽量纯生信、复用 TCGA/GEO 等公共数据，
可含少量湿实验验证）的文章有哪些？"

Curated PMID set is assigned to themes by hand after reading title + abstract
evidence; the "data" and "study type" columns are derived mechanically from the
PubMed abstract text so they stay auditable.

Outputs
-------
results/tables/intratumoral_bacteria_bioinformatics_curated.csv
analysis/out/curation_evidence.txt   (abstract evidence strings per paper)
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RECS = {r["pmid"]: r for r in json.loads(
    (ROOT / "data" / "processed" / "tumor_microbiome_pubmed_records.json").read_text(encoding="utf-8"))}
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)
(ROOT / "analysis" / "out").mkdir(parents=True, exist_ok=True)

# theme -> ordered pmids
CURATED: dict[str, list[str]] = {
    "1 数据资源与分析工具": [
        "32214244", "33382980", "36179670", "38599211", "42092351", "32257046",
        "37850841", "38259095", "37131219", "38195165", "36189892", "42446235",
        "41650974", "39896561", "28118849", "31639030", "35610202", "35972351",
        "33659025", "40558028", "42482748",
    ],
    "2 方法学可靠性与争议": [
        "37811944", "37555750", "40901923", "40901921", "38396294", "40314439",
        "42017663", "41714823", "41289425",
    ],
    "3 泛癌/跨癌种生信分析": [
        "40952000", "39565103", "40875911", "41373739", "40865517", "35782745",
        "42013847", "41885500", "28410234", "41234847", "42560302", "41343929",
        "42865545",
    ],
    "4 肺癌/胸部肿瘤": [
        "37664112", "38375207", "40860241", "38721596", "41142936", "38090530",
        "39915534", "39845316", "39531992", "41562055", "41608560",
        "37074188", "41932249", "36999042",
    ],
    "5 消化道肿瘤": [
        "34778069", "35586863", "40677704", "40202312", "37260411", "36884119",
        "37465100", "35822725", "39961991", "38073489", "36786568", "39072166",
        "34422640", "39833403", "41438381", "42560302", "42488410", "42399737",
        "41823382", "40492114", "41472855", "32575865", "37114494", "37380815",
        "41143238", "36732743", "34136388", "32962112", "40013793", "42491956",
        "38556382", "41647993",
    ],
    "6 头颈/口腔/鼻咽/神经与泌尿生殖系统": [
        "35122129", "36555172", "36362038", "37895136", "39497974", "38722731",
        "39907412", "41395468", "40475759", "40932704", "39173404", "39660865",
        "40606983", "38808774", "39333148", "39751930", "40761254", "42436233",
        "39397906", "36975828", "40433663", "38515713", "39117066", "39982111",
        "40908816", "36999009", "38442004", "34741726", "38421176", "38868430",
        "41234866", "36834564", "40029433", "41361827", "39924593",
    ],
    "7 综述/系统评价/立场文章": [
        "40713785", "42409802", "41552720", "40575164", "41480207", "41387030",
        "42687557", "42270015", "42558585", "42835748", "38347100", "42132991",
        "34986765",
    ],
}

PUBLIC = {
    "TCGA": r"\bTCGA\b|cancer genome atlas",
    "TCMA": r"cancer microbiome atlas|\bTCMA\b",
    "GEO": r"gene expression omnibus|\bGEO\b",
    "CPTAC": r"\bCPTAC\b",
    "ICGC": r"\bICGC\b",
    "cBioPortal": r"cbioportal",
    "ORIEN": r"\bORIEN\b|Oncology Research Information Exchange",
    "100kGP": r"100,000 genomes|100,000 Genomes",
    "scRNA": r"single-?cell (RNA|transcriptom)|scRNA",
    "spatial": r"spatial transcriptom",
    "GDSC/药敏": r"GDSC|genomics of drug sensitivity",
    "免疫治疗队列": r"IMvigor|immunotherapy cohort|anti-PD-?1 cohort",
    "GTEx": r"\bGTEx\b",
}
REVIEW = re.compile(r"review|perspective|comment|editorial|systematic review|meta-analysis", re.I)
WETLAB = re.compile(
    r"\bwe (performed|conducted|cultured|isolated|collected|validated|stained|"
    r"confirmed|generated)\b|cell line|organoid|xenograft|western blot|"
    r"immunohistochemistr|in vivo|mice|mouse|qPCR|quantitative PCR|FISH|"
    r"flow cytometr|colony formation|transwell|16S rRNA (gene )?(amplicon )?sequencing"
    r"|knockdown|overexpress|siRNA|CCK-?8|EdU|wound healing",
    re.I,
)


def evidence(ab: str) -> str:
    sents = re.split(r"(?<=[.;])\s+", ab)
    keep = [s.strip() for s in sents if re.search(
        r"TCGA|cancer genome atlas|gene expression omnibus|\bGEO\b|CPTAC|\bICGC\b|"
        r"public|re-?analy|in silico|Expression Omnibus|database", s, re.I)]
    return " ".join(keep)


rows = []
evid_lines = []
seen = set()
for theme, pmids in CURATED.items():
    for p in pmids:
        if p in seen:
            continue
        seen.add(p)
        r = RECS.get(p)
        if r is None:
            print("MISSING RECORD", p)
            continue
        blob = r["title"] + " " + r["abstract"]
        data = [k for k, v in PUBLIC.items() if re.search(v, blob, re.I)]
        is_rev = bool(REVIEW.search(r["title"])) or any(
            re.search(r"review|meta-analysis", t, re.I) for t in r["pubtypes"])
        rows.append({
            "theme": theme,
            "pmid": p,
            "year": r["year"],
            "journal": r["journal"],
            "doi": r["doi"],
            "title": r["title"],
            "data_mentioned_in_abstract": "; ".join(data) or "未在摘要中明确",
            "article_type": "综述/评论" if is_rev else "研究论文",
            "wetlab_mention_in_abstract": "有提及" if WETLAB.search(r["abstract"]) else "未见",
            "pubtypes": "; ".join(r["pubtypes"]),
        })
        evid_lines.append(f"### {p} | {r['year']} | {r['journal']} | {r['doi']}\n{r['title']}\nEVID: {evidence(r['abstract'])[:700]}\n")

with (TAB / "intratumoral_bacteria_bioinformatics_curated.csv").open(
        "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

(ROOT / "analysis" / "out" / "curation_evidence.txt").write_text(
    "\n".join(evid_lines), encoding="utf-8")

missing_doi = [r["pmid"] for r in rows if not r["doi"]]
print("curated:", len(rows), "with DOI:", sum(1 for r in rows if r["doi"]),
      "missing DOI:", missing_doi)
from collections import Counter
print(Counter(r["theme"] for r in rows))
