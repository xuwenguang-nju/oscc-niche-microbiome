#!/usr/bin/env python3
"""Convert PMC JATS XML to compact, readable text for close reading.

Goal is not a pretty PDF but a file a reader (human or model) can scan for what each study
actually did: design, n, platform, the named analysis steps, and the reported statistics.
Figures/tables are collapsed to their captions, because the analysis modules are what matter
for extending our own study.

Per paper it writes <label>.md with abstract, section headings and their paragraphs, and a
one-screen "digest" at the top listing which analysis method names were detected.

Run: python3 to_text.py
"""
from __future__ import annotations

import re
from pathlib import Path
import xml.etree.ElementTree as ET

# method vocabulary to scan for -- these are the modules a reader would need to reproduce
VOCAB = [
    "PERMANOVA", "PCoA", "NMDS", "Shannon", "Simpson", "Chao1", "LEfSe", "ANCOM",
    "ALDEx2", "DESeq2", "edgeR", "limma", "Wilcoxon", "Kruskal", "Spearman", "Pearson",
    "Mantel", "Procrustes", "CLR", "Bray", "UniFrac", "Jaccard",
    "CIBERSORT", "CIBERSORTx", "ssGSEA", "GSVA", "ESTIMATE", "MCP-counter", "xCell",
    "TIMER", "quanTIseq", "EPIC", "ImmuCellAI", "deconvolution",
    "Cox", "Kaplan", "log-rank", "LASSO", "Ridge", "Elastic net", "random forest",
    "nomogram", "ROC", "AUC", "C-index", "concordance", "multivariable", "hazard ratio",
    "WGCNA", "NMF", "consensus clustering", "k-means", "hierarchical clustering",
    "GSEA", "KEGG", "GO enrichment", "Reactome", "GSVA", "Metascape", "DAVID",
    "PICRUSt", "PICRUSt2", "Tax4Fun", "FAPROTAX", "HUMAnN", "Kraken", "MetaPhlAn",
    "SHOGUN", "Kraken2", "Bracken", "centrifuge", "PathoScope", "ConQuR", "ComBat",
    "decontam", "microbiomeHD", "SourceTracker", "FEAST",
    "Cytoscape", "MCODE", "STRING", "ceRNA", "miRcode", "starBase", "GDCRNATools",
    "scRNA", "single-cell", "FISH", "immunohistochemistry", "TMA", "multiplex",
    "16S", "metagenom", "metatranscriptom", "metabolom", "LC-MS", "RPPA", "methylation",
    "HPV", "stage", "grade", "smoking", "alcohol",
]

SRC = Path("texts")
OUT = Path("md")
OUT.mkdir(exist_ok=True)


def strip_ns(tag: str) -> str:
    return tag.split("}")[-1]


def text_of(el) -> str:
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip()


def walk(el, depth=0, buf=None):
    if buf is None:
        buf = []
    for child in el:
        tag = strip_ns(child.tag)
        if tag == "sec":
            title = child.find("./{*}title")
            if title is None:
                title = child.find("./title")
            if title is not None:
                buf.append(f"\n{'#' * min(depth + 2, 6)} {text_of(title)}\n")
            walk(child, depth + 1, buf)
        elif tag == "p":
            t = text_of(child)
            if t:
                buf.append(t + "\n")
        elif tag in ("fig", "table-wrap"):
            cap = child.find("./{*}caption")
            if cap is None:
                cap = child.find("./caption")
            if cap is not None:
                buf.append(f"[{tag} caption] {text_of(cap)[:400]}\n")
        elif tag in ("title", "abstract"):
            if tag == "abstract":
                t = text_of(child)
                if t:
                    buf.append(t + "\n")
        else:
            walk(child, depth, buf)
    return buf


def main() -> None:
    rows = []
    for x in sorted(SRC.glob("*.xml")):
        label = x.name.split("_")[0] + "_" + x.name.split("_")[1]
        try:
            root = ET.parse(x).getroot()
        except Exception as e:
            print(f"  ERR {x.name}: {e}")
            continue

        title_el = root.find(".//{*}article-title")
        if title_el is None:
            title_el = root.find(".//article-title")
        title = text_of(title_el) if title_el is not None else x.stem

        body = root.find(".//{*}body")
        if body is None:
            body = root.find(".//body")

        buf = [f"# {title}\n", f"source: {x.name}\n"]
        if body is None:
            buf.append("\n**NO BODY: this PMC record is abstract-only.**\n")
            walk(root, 0, buf)
        else:
            walk(body, 0, buf)

        full = "\n".join(buf)
        low = full.lower()
        hits = sorted({v for v in VOCAB if v.lower() in low},
                      key=lambda v: low.index(v.lower()))
        digest = ("## method vocabulary detected\n"
                  + ", ".join(hits) + "\n\n")
        dst = OUT / (x.stem + ".md")
        dst.write_text(full[:200] + "\n\n" + digest + "\n---\n\n" + full,
                       encoding="utf-8")
        rows.append((dst.name, len(full), len(hits), body is not None))
        print(f"  {dst.name:58} {len(full)//1000:4}K chars  {len(hits):2} methods  "
              f"body={'yes' if body is not None else 'NO'}")

    print(f"\nconverted {len(rows)}; with body: {sum(1 for r in rows if r[3])}")


if __name__ == "__main__":
    main()
