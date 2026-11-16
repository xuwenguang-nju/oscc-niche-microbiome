#!/usr/bin/env python3
"""评估 gutMGene v2.0 能否与 OSCC 瘤内菌研究（方向三）联合。

输入：data/raw/gutmgene/gutmgene_v2_*.json（由 analysis/scripts/fetch_gutmgene.py 下载）
输出：
  results/tables/gutmgene_oscc_taxon_coverage.csv   每个 OSCC 相关类群在 gutMGene 中的证据覆盖
  results/tables/gutmgene_oscc_candidate_network.csv 类群 -> 代谢物 -> 宿主基因 三元链（可做预注册候选集）
  results/tables/gutmgene_human_gene_hgnc_check.csv  238 个人源基因的 HGNC 符号校验

用法：python analysis/scripts/analyse_gutmgene_oscc_linkage.py
"""
from __future__ import annotations

import csv
import json
import urllib.parse
import urllib.request
from pathlib import Path

RAW = Path("data/raw/gutmgene")
TAB = Path("results/tables")

# OSCC / 口腔相关类群（属级，名称包含匹配；gutMGene 的分类级不统一，故用包含匹配）
OSCC_GENERA = [
    "Fusobacterium", "Porphyromonas", "Prevotella", "Streptococcus", "Veillonella",
    "Parvimonas", "Peptostreptococcus", "Capnocytophaga", "Treponema", "Actinomyces",
    "Rothia", "Neisseria", "Haemophilus", "Gemella", "Eikenella", "Aggregatibacter",
    "Tannerella", "Alloprevotella", "Selenomonas", "Leptotrichia", "Solobacterium",
    "Filifactor", "Olsenella", "Dialister", "Megasphaera", "Campylobacter",
    "Stomatobaculum", "Lactobacillus", "Bifidobacterium",
]

DATASETS = {"microbe_metabolite": "菌-代谢物", "microbe_gene": "菌-基因", "metabolite_gene": "代谢物-基因"}


def load(species: str, dataset: str) -> list[dict]:
    path = RAW / f"gutmgene_v2_{species}_{dataset}.json"
    return json.loads(path.read_text(encoding="utf-8"))["payload"].get("resultSet", [])


def is_causal(mode: str | None) -> bool:
    return bool(mode) and "causally" in str(mode)


def in_panel(name: str | None) -> bool:
    if not name:
        return False
    low = str(name).lower()
    return any(g.lower() in low for g in OSCC_GENERA)


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def hgnc_approved(symbol: str) -> bool | None:
    q = urllib.parse.quote(f'symbol:"{symbol}" AND status:Approved')
    req = urllib.request.Request(
        f"https://rest.genenames.org/search/{q}",
        headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=20) as fh:
            return json.loads(fh.read()).get("response", {}).get("numFound", 0) > 0
    except Exception:
        return None


def main() -> int:
    mm = load("human", "microbe_metabolite")
    mg = load("human", "microbe_gene")
    gg = load("human", "metabolite_gene")
    m_mg = load("mouse", "microbe_gene")

    mm_panel = [r for r in mm if in_panel(r.get("gut_microbiota"))]
    mg_panel = [r for r in mg if in_panel(r.get("gut_microbiota"))]

    # 代谢物 -> 基因 索引（人）
    met2genes: dict[str, list[dict]] = {}
    for r in gg:
        met = r.get("metabolite")
        if met and r.get("gene"):
            met2genes.setdefault(met, []).append(r)

    # 三元链：类群 -> 代谢物 -> 基因
    chains: list[dict] = []
    for r in mm_panel:
        met = r.get("metabolite")
        for g in met2genes.get(met, []):
            chains.append(
                {
                    "taxon": r.get("gut_microbiota"),
                    "taxon_rank": r.get("rank"),
                    "metabolite": met,
                    "metabolite_chebi": r.get("metabolite_chebi"),
                    "gene": g.get("gene"),
                    "taxon_metabolite_evidence": r.get("associative_modes"),
                    "metabolite_gene_evidence": g.get("associative_modes"),
                    "taxon_metabolite_pmid": r.get("pmid"),
                    "metabolite_gene_pmid": g.get("pmid"),
                }
            )
    write_csv(
        TAB / "gutmgene_oscc_candidate_network.csv",
        chains,
        ["taxon", "taxon_rank", "metabolite", "metabolite_chebi", "gene",
         "taxon_metabolite_evidence", "metabolite_gene_evidence",
         "taxon_metabolite_pmid", "metabolite_gene_pmid"],
    )

    # 每个类群的可达基因（直接 + 经代谢物）
    coverage: list[dict] = []
    taxa = sorted({r.get("gut_microbiota") for r in mm_panel + mg_panel if r.get("gut_microbiota")})
    for t in taxa:
        direct = [r for r in mg_panel if r.get("gut_microbiota") == t and r.get("gene")]
        mets = {r.get("metabolite") for r in mm_panel if r.get("gut_microbiota") == t and r.get("metabolite")}
        via = {g.get("gene") for m in mets for g in met2genes.get(m, []) if g.get("gene")}
        coverage.append(
            {
                "taxon": t,
                "n_microbe_gene_edges": len(direct),
                "n_causal_microbe_gene_edges": sum(1 for r in direct if is_causal(r.get("associative_modes"))),
                "n_metabolites": len(mets),
                "n_genes_via_metabolites": len(via),
                "n_genes_total": len({r.get("gene") for r in direct} | via),
                "direct_genes": ";".join(sorted({r.get("gene") for r in direct})),
                "metabolites": ";".join(sorted(m for m in mets if m))[:400],
            }
        )
    coverage.sort(key=lambda x: (-x["n_genes_total"], -x["n_metabolites"]))
    write_csv(
        TAB / "gutmgene_oscc_taxon_coverage.csv",
        coverage,
        ["taxon", "n_microbe_gene_edges", "n_causal_microbe_gene_edges", "n_metabolites",
         "n_genes_via_metabolites", "n_genes_total", "direct_genes", "metabolites"],
    )

    # HGNC 校验（人源基因并集）
    genes = sorted({r.get("gene") for r in mg if r.get("gene")}
                   | {r.get("gene") for r in gg if r.get("gene")})
    hgnc_rows = []
    for s in genes:
        ok = hgnc_approved(s)
        hgnc_rows.append({"gene": s, "hgnc_approved": {True: "yes", False: "no", None: "unavailable"}[ok]})
    write_csv(TAB / "gutmgene_human_gene_hgnc_check.csv", hgnc_rows, ["gene", "hgnc_approved"])

    uniq_genes_mg = len({r.get("gene") for r in mg if r.get("gene")})
    uniq_genes_gg = len({r.get("gene") for r in gg if r.get("gene")})
    uniq_mg = len({r.get("gut_microbiota") for r in mg if r.get("gut_microbiota")})
    reach = {c["gene"] for c in chains} | {r.get("gene") for r in mg_panel if r.get("gene")}
    bad = [r["gene"] for r in hgnc_rows if r["hgnc_approved"] == "no"]

    print(f"human 菌-基因：{len(mg)} 边 / {uniq_mg} 类群 / {uniq_genes_mg} 基因")
    print(f"human 代谢物-基因：{len(gg)} 边 / {uniq_genes_gg} 基因")
    print(f"OSCC 面板：菌-代谢物 {len(mm_panel)} 边、菌-基因 {len(mg_panel)} 边")
    print(f"三元链 {len(chains)} 条 -> {TAB/'gutmgene_oscc_candidate_network.csv'}")
    print(f"面板可达基因 {len(reach)} 个；类群覆盖表 {len(coverage)} 行")
    print(f"HGNC 校验：{len(genes)} 个符号，未通过 {len(bad)} 个 -> {bad}")
    print(f"小鼠菌-基因表（规模对照）：{len(m_mg)} 边 / "
          f"{len({r.get('gut_microbiota') for r in m_mg})} 类群 / {len({r.get('gene') for r in m_mg})} 基因")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
