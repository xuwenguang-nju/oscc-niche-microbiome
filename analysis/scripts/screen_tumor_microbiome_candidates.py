"""Screen the PubMed sweep for papers that are (mostly) dry-lab bioinformatics
analyses of intratumoral bacteria / tumor microbiome reusing public data
(TCGA, GEO, CPTAC, ICGC), optionally with minor wet-lab validation.

Writes:
  results/tables/intratumoral_microbiome_bioinformatics_candidates.csv
  analysis/out/screening_shortlist.md
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data" / "processed" / "tumor_microbiome_pubmed_records.json"
TAB = ROOT / "results" / "tables"
OUT = ROOT / "analysis" / "out"
TAB.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

recs = json.loads(SRC.read_text(encoding="utf-8"))

MICROBIOME = re.compile(
    r"intratumo(ur|r)al microb|intratumo(ur|r)al bacteri|intratumor bacteri|"
    r"tumo(ur)?r?[- ]?(associated |resident )?(microbiome|microbiota|bacteri)|"
    r"cancer microbiome|microbial (signature|community|composition|reads|abundance)|"
    r"tissue[- ]resident microb|oncobiome|bacteri(al|ome) (burden|load|signature)|"
    r"16S|metagenom|metataxonom|microbial deconvolution|bacterial reads",
    re.I,
)
PUBLIC_DATA = re.compile(
    r"TCGA|the cancer genome atlas|gene expression omnibus|\bGEO\b|CPTAC|"
    r"\bICGC\b|public(ly)? available|public (cohort|dataset|data|database|repositor)|"
    r"re-?analy|in silico|secondary analysis|existing dataset|"
    r"\bGTEx\b|\bTARGET\b|ArrayExpress|Sequence Read Archive|\bSRA\b|\bHRA\b|GSA\b",
    re.I,
)
DRYLAB = re.compile(
    r"bioinformatic|computational|machine learning|deep learning|random forest|"
    r"LASSO|WGCNA|deconvolution|transcriptom|RNA-?seq|single-?cell|"
    r"scRNA|risk (model|score|signature)|prognostic (model|signature|score)|"
    r"nomogram|immune infiltration|CIBERSORT|ssGSEA|GSVA|enrichment analysis|"
    r"differential(ly)? express|Mendelian randomi|GWAS|network analysis|"
    r"predictive model|nomogram|biomarker|co-?expression|NGS|sequencing data|"
    r"multi-?omics|integrat(ive|ed) analysis|microbial (qtl|QTL)|microbiome QTL",
    re.I,
)
WETLAB_HEAVY = re.compile(
    r"\bwe (isolated|performed|collected|cultured)|cell line|mouse|mice|in vivo|"
    r"xenograft|knockdown|western blot|immunohistochemis|in vitro|organoid|"
    r"flow cytometr|qPCR|FISH|animal model|patient cohort|prospectively|"
    r"16S rRNA (gene )?(amplicon )?sequencing of|we sequenced",
    re.I,
)
NONHUMAN_ONLY = re.compile(
    r"canine|dog|dogs|bovine|cattle|murine model only|chicken|pig |porcine|"
    r"zebrafish|drosophila",
    re.I,
)
REVIEW = re.compile(r"review|perspective|comment|editorial|minireview|hypothesis", re.I)

rows = []
for r in recs:
    blob = (r["title"] + " " + r["abstract"]).lower()
    t = r["title"].lower()
    if not MICROBIOME.search(blob):
        continue
    if not DRYLAB.search(blob):
        continue
    pub = bool(PUBLIC_DATA.search(blob))
    is_review = bool(
        any(REVIEW.search(p) for p in r["pubtypes"]) or REVIEW.search(t)
    )
    wet = len(set(WETLAB_HEAVY.findall(blob)))
    score = 0
    if "intratumo" in blob or "intratumor" in blob or "tumour microbiome" in blob or "tumor microbiome" in blob:
        score += 4
    if "microbiota" in blob or "microbiome" in blob:
        score += 2
    if pub:
        score += 5
    if re.search(r"\bTCGA\b|the cancer genome atlas", blob):
        score += 4
    if re.search(r"gene expression omnibus|\bGEO\b", blob):
        score += 2
    if DRYLAB.search(blob):
        score += 2
    if "single-cell" in blob or "single cell" in blob:
        score += 1
    if is_review:
        score -= 8
    if NONHUMAN_ONLY.search(t):
        score -= 6
    if not pub and wet >= 3:
        score -= 3
    rows.append(
        {
            "pmid": r["pmid"],
            "year": r["year"],
            "journal": r["journal"],
            "doi": r["doi"],
            "title": r["title"],
            "public_data": int(pub),
            "review": int(is_review),
            "total": score,
            "queries": "; ".join(r["queries"]),
            "abstract": r["abstract"],
        }
    )

rows.sort(key=lambda d: (-d["total"], d["year"], d["title"]))
with (TAB / "intratumoral_microbiome_bioinformatics_candidates.csv").open(
    "w", newline="", encoding="utf-8-sig"
) as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

with (OUT / "screening_shortlist.md").open("w", encoding="utf-8") as fh:
    for d in rows:
        if d["total"] < 6:
            continue
        fh.write(
            f"### {d['total']:>3} | {d['year']} | {d['journal']} | PMID {d['pmid']} | {d['doi']}\n"
            f"{d['title']}\n"
            f"[pubdata={d['public_data']} review={d['review']}] {d['abstract'][:800]}\n\n"
        )

print("candidates:", len(rows), "score>=6:", sum(1 for d in rows if d["total"] >= 6))
