#!/usr/bin/env python3
"""Extract what each paper actually did, from its own methods text.

The earlier analysis plan was built from Chinese blog summaries. This is the primary-source
version: for every retrieved full text it records, from the methods/results text itself,
  - cohort size and specimen type
  - sequencing platform / 16S region / data source
  - named software and R packages
and then builds a paper x analysis-category matrix.

Only what the text actually contains is recorded; absence means "not mentioned", which is
NOT the same as "not done". That distinction is reported explicitly because it bounds how the
matrix can be used.

Run: python3 extract_analyses.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

MD = Path("md")

CATEGORIES = {
    "alpha_diversity": ["shannon", "simpson", "chao1", "observed otus", "alpha divers",
                        "richness", "faith"],
    "beta_diversity": ["bray", "unifrac", "jaccard", "aitchison", "pcoa", "nmds",
                       "beta divers", "permanova", "anosim", "adonis"],
    "differential_taxa": ["lefse", "ancom", "aldex", "deseq", "edger", "limma",
                          "differential abundance", "differentially abundant",
                          "metastats", "maaslin"],
    "machine_clustering": ["nmf", "consensus clustering", "unsupervised clustering",
                           "k-means", "hierarchical clustering", "random forest",
                           "lasso", "elastic net", "xgboost", "svm"],
    "survival": ["cox", "kaplan", "log-rank", "hazard ratio", "nomogram", "c-index",
                 "concordance index", "survival curve"],
    "immune_infiltration": ["cibersort", "ssgsea", "gsva", "estimate", "mcp-counter",
                            "xcell", "timer", "quantiseq", "epic", "immucellai",
                            "immune infiltration", "deconvolution", "immunedeconv"],
    "pathway_enrichment": ["gsea", "kegg", "go enrichment", "reactome", "metascape",
                           "david", "gene ontology", "enrichment analysis"],
    "function_prediction": ["picrust", "tax4fun", "faprotax", "humann", "bugbase",
                            "functional prediction", "predicted metagenome"],
    "taxonomic_profiling_rnaseq": ["kraken", "metaphlan", "shogun", "bracken",
                                   "centrifuge", "pathoscope", "kallisto", "sigma",
                                   "rna-seq data", "microbial reads"],
    "batch_correction": ["conqur", "combat", "batch effect", "batch correction",
                         "quantile regression"],
    "decontamination": ["decontam", "contaminant", "negative control", "blank",
                        "reagent"],
    "host_microbe_network": ["mantel", "procrustes", "cytoscape", "correlation network",
                             "microbe-host", "host-microbe", "interaction network"],
    "multiomics": ["proteom", "metabolom", "methylat", "rppa", "transcriptom",
                   "single-cell", "scrna", "mirna", "lncrna", "cerna"],
    "wetlab": ["fish", "immunohistochem", "pcr", "cell line", "mouse", "xenograft",
               "in vitro", "western blot", "16s rrna gene sequencing of"],
    "clinical_association": ["stage", "grade", "lymph node", "hpv", "smoking",
                             "alcohol", "tnm", "clinicopathologic"],
    "novelty_literature_search": ["pubmed", "literature", "previously reported",
                                  "reported in", "been reported"],
}

PLATFORM = {
    "16S": ["16s"],
    "metagenomics": ["metagenom", "shotgun"],
    "RNA-seq(TCGA)": ["tcga", "rna-seq", "rna sequencing"],
    "metabolomics": ["metabolom", "lc-ms", "mass spectrom"],
    "scRNA-seq": ["single-cell", "scrna"],
}

SOFT = re.compile(
    r"\b(R package|R version|R software|QIIME2?|QIIME|mothur|vegan|phyloseq|"
    r"DADA2|USEARCH|VSEARCH|UCLUST|SILVA|Greengenes|GTDB|Kraken2?|MetaPhlAn[0-9]?|"
    r"SHOGUN|Bracken|HUMAnN[0-9]?|PICRUSt2?|Tax4Fun2?|FAPROTAX|"
    r"CIBERSORTx?|ssGSEA|GSVA|ESTIMATE|MCP-counter|xCell|TIMER[0-9.]*|quanTIseq|EPIC|"
    r"limma|edgeR|DESeq2?|LEfSe|ANCOM-?BC?|ALDEx2?|metagenomeSeq|MaAsLin2?|"
    r"survival|survminer|lifelines|glmnet|randomForest|WGCNA|NMF|ConsensusClusterPlus|"
    r"Cytoscape|STRING|Metascape|DAVID|clusterProfiler|GSEA|GSVA|fgsea|"
    r"ConQuR|ComBat|decontam|SourceTracker|FEAST|microbiomeHD|"
    r"STAR|HISAT2|featureCounts|salmon|kallisto|Bowtie2?|BWA|"
    r"CellChat|Seurat|Scanpy|inferCNV|GEPIA|cBioPortal|TIMER)\b"
)


def load(label_prefix: str) -> str | None:
    hits = sorted(MD.glob(f"{label_prefix}_*.md"))
    return hits[0].read_text(encoding="utf-8", errors="replace") if hits else None


def cohort_size(txt: str) -> str:
    pats = [
        r"(\d{2,5})\s*(?:paired\s+)?(?:tumou?r|tumor)\s*(?:and|/)\s*"
        r"(?:adjacent|normal|paracancerous|non-?tumou?r)[^.]{0,40}?(\d{2,5})",
        r"a\s+total\s+of\s+(\d{1,5})\s+(?:patients|samples|specimens|cases)",
        r"(\d{1,5})\s+(?:paired\s+)?(?:tumou?r|tumor)\s+(?:and|/)\s+"
        r"(?:adjacent|normal|paracancerous)",
        r"included\s+(\d{1,5})\s+(?:patients|samples|cases)",
        r"(\d{1,5})\s+(?:OSCC|HNSCC|ESCC|HCC|PDAC|CRC|LUAD|TNBC|ICC)\s+"
        r"(?:patients|samples|specimens)",
    ]
    for p in pats:
        m = re.search(p, txt, re.I)
        if m:
            return " / ".join(g for g in m.groups() if g) or m.group(0)[:40]
    return ""


def main() -> None:
    labels = sorted({p.name.split("_PMC")[0] for p in MD.glob("*.md")})
    rows = []
    for lab in labels:
        txt = load(lab.split("_m")[0] if False else lab) or ""
        pass
    # simpler: iterate actual files
    rows = []
    for p in sorted(MD.glob("*.md")):
        lab = re.sub(r"_\d+_PMC\d+$", "", p.stem)
        txt = p.read_text(encoding="utf-8", errors="replace")
        low = txt.lower()
        body = "**NO BODY" not in txt
        rec = {"label": lab, "has_body": body, "chars": len(txt)}
        rec["cohort"] = cohort_size(txt)
        rec["platforms"] = [k for k, v in PLATFORM.items() if any(s in low for s in v)]
        rec["software"] = sorted({m.group(0) for m in SOFT.finditer(txt)})
        rec["categories"] = [k for k, v in CATEGORIES.items()
                             if any(s in low for s in v)]
        rows.append(rec)

    json.dump(rows, open("paper_analysis_matrix.json", "w"), indent=1)

    cats = list(CATEGORIES)
    print(f"{'paper':32} " + " ".join(f"{c[:6]:>6}" for c in cats))
    for r in rows:
        if not r["has_body"]:
            print(f"{r['label']:32} (abstract only)")
            continue
        marks = " ".join(f"{'  X   ' if c in r['categories'] else '   .  '}" for c in cats)
        print(f"{r['label'][:32]:32} {marks}")

    print("\n=== per-category: how many of the "
          f"{sum(1 for r in rows if r['has_body'])} full-text papers do it ===")
    nb = [r for r in rows if r["has_body"]]
    for c in cats:
        n = sum(1 for r in nb if c in r["categories"])
        who = ", ".join(r["label"].split("_")[0] for r in nb if c in r["categories"])
        print(f"  {c:28} {n:2}/{len(nb)}   {who}")

    print("\n=== cohort descriptions found ===")
    for r in nb:
        print(f"  {r['label'][:30]:30} {r['cohort'] or '(not parsed)'}")

    print("\n=== platforms ===")
    for r in nb:
        print(f"  {r['label'][:30]:30} {','.join(r['platforms'])}")


if __name__ == "__main__":
    main()
