#!/usr/bin/env python3
"""S6 discovery: what host-side data actually exists for OSCC/HNSCC coupling.

The four 16S cohorts we analysed carry no host transcriptome, so S6 cannot
correlate our own genus abundances with expression. It has to stand on public data
where microbial signal and host expression come from the SAME samples. This script
only probes availability -- it downloads nothing large.

Sources probed
  1. TCMA (The Cancer Microbiome Atlas): decontaminated microbial abundances mined
     from TCGA RNA-seq. If reachable, TCGA-HNSC gives both sides per patient.
  2. GDC API: TCGA-HNSC case/sample counts and available data types.
  3. GEO: OSCC/HNSCC studies that might contain BOTH 16S and RNA-seq.

Run: python3 discover_s6_sources.py
Out: prints a report; writes results/tables/s6_data_availability.json
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "tables"
UA = {"User-Agent": "wisp-science/1.0 (S6 data discovery)"}


def get(url: str, timeout: int = 30, headers: dict | None = None) -> tuple[int, bytes]:
    req = urllib.request.Request(url, headers={**UA, **(headers or {})})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, b""
    except Exception as e:                                    # noqa: BLE001
        return -1, str(e).encode()


def probe_tcma(report: dict) -> None:
    """TCMA hosts decontaminated TCGA microbiome matrices."""
    cands = {
        "galaxy_tcma": "https://tcma.galaxyproject.org/",
        "cell_chom_doi": "https://doi.org/10.1016/j.chom.2020.12.001",
    }
    res: dict = {}
    for k, u in cands.items():
        code, body = get(u, timeout=25)
        res[k] = {"code": code, "note": body[:100].decode("utf-8", "replace") if code != 200 else "OK"}

    # figshare search API (TCMA figure data / supplementary matrices live there)
    try:
        req = urllib.request.Request(
            "https://api.figshare.com/v2/articles/search",
            data=json.dumps({"search_for": "Cancer Microbiome Atlas",
                             "page_size": 10}).encode(),
            headers={**UA, "Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as r:
            hits = json.loads(r.read())
        res["figshare"] = {"n": len(hits),
                           "titles": [h.get("title", "")[:100] for h in hits[:8]]}
    except Exception as e:                                    # noqa: BLE001
        res["figshare"] = {"error": str(e)[:120]}

    # Zenodo
    code, body = get("https://zenodo.org/api/records?q=%22Cancer%20Microbiome%20Atlas%22&size=10")
    if code == 200:
        try:
            h = json.loads(body).get("hits", {}).get("hits", [])
            res["zenodo"] = {"n": len(h),
                             "titles": [x.get("metadata", {}).get("title", "")[:100] for x in h]}
        except Exception as e:                                # noqa: BLE001
            res["zenodo"] = {"error": str(e)[:120]}
    else:
        res["zenodo"] = {"code": code}
    report["tcma"] = res


def probe_gdc(report: dict) -> None:
    """TCGA-HNSC totals via the public GDC API."""
    base = "https://api.gdc.cancer.gov"
    code, body = get(f"{base}/projects/HNSC?expand=summary,summary.experimental_strategies"
                     f"&fields=project_id,name,summary.case_count")
    out: dict = {"project_query_code": code}
    if code == 200:
        try:
            d = json.loads(body)["data"]
            out["project"] = d.get("project_id")
            out["case_count"] = d.get("summary", {}).get("case_count")
            out["experimental_strategies"] = [
                {"strategy": s.get("experimental_strategy"), "n": s.get("case_count")}
                for s in (d.get("summary", {}).get("experimental_strategies") or [])]
        except Exception as e:                                # noqa: BLE001
            out["parse_error"] = str(e)[:120]

    flt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-HNSC"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}}]}
    c2, b2 = get(f"{base}/files?filters={urllib.parse.quote(json.dumps(flt))}&size=0")
    if c2 == 200:
        try:
            out["rnaseq_star_counts_files"] = json.loads(b2)["data"]["pagination"]["total"]
        except Exception:                                     # noqa: BLE001
            pass
    report["tcga_hnsc"] = out


def probe_geo(report: dict) -> None:
    """GEO studies that may pair 16S with host RNA-seq in OSCC/HNSCC."""
    queries = {
        "oscc_16s": "(oral squamous cell carcinoma) AND (16S OR microbiome)",
        "oscc_16s_and_rnaseq": "(oral squamous cell carcinoma) AND 16S AND (RNA-seq OR transcriptome)",
        "hnsc_microbiome_rnaseq": "(head and neck) AND microbiome AND RNA-seq",
    }
    res: dict = {}
    for name, term in queries.items():
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds"
               f"&retmax=25&retmode=json&term={urllib.parse.quote(term)}")
        code, body = get(url)
        if code != 200:
            res[name] = {"code": code}
            continue
        try:
            er = json.loads(body)["esearchresult"]
            ids, total = er["idlist"], er["count"]
        except Exception:                                     # noqa: BLE001
            res[name] = {"parse_error": True}
            continue
        titles = []
        if ids:
            c2, b2 = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds"
                         f"&retmode=json&id={','.join(ids)}")
            if c2 == 200:
                j = json.loads(b2)["result"]
                for u in j["uids"]:
                    r = j[u]
                    titles.append({"gse": r.get("accession"), "n_samples": r.get("n_samples"),
                                   "title": (r.get("title") or "")[:110]})
        res[name] = {"count": total, "top": titles[:12]}
    report["geo"] = res


def main() -> None:
    report: dict = {}
    probe_tcma(report)
    probe_gdc(report)
    probe_geo(report)

    OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT / "s6_data_availability.json"
    dst.write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")

    print("=" * 78)
    print("TCMA")
    for k, v in report["tcma"].items():
        print(f"  {k}: {v}")
    print("\nTCGA-HNSC")
    for k, v in report["tcga_hnsc"].items():
        print(f"  {k}: {v}")
    print("\nGEO")
    for k, v in report["geo"].items():
        print(f"  {k}: count={v.get('count')}")
        for t in v.get("top", [])[:6]:
            print(f"     {t['gse']}  n={t['n_samples']}  {t['title']}")
    print(f"\nwrote {dst}")


if __name__ == "__main__":
    main()
