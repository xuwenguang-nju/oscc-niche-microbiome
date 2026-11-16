#!/usr/bin/env python3
"""Fix two bad PubMed resolutions, then fetch open-access full texts from PMC.

Two of the blog-derived queries resolved badly:
  - 01 matched a glioma multi-omics paper instead of the HNSC Fusobacterium paper
  - 07 failed to resolve at all
Both get tighter, title-exact queries here. Every hit is then verified by comparing the
returned title against the expected title fragment, so a wrong match cannot slip through
silently (which is exactly what happened on the first pass).

Full texts are pulled as JATS XML from PMC (open-access subset only). Nothing is fetched for
papers without a PMC id: those would need a subscription and are left alone.

Run: python3 fetch_fulltext.py
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "wisp-science/1.0 (paper fulltext)"}


def http(url: str, timeout: int = 90) -> bytes:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def http_json(url: str, timeout: int = 90) -> dict:
    return json.loads(http(url, timeout))


def esearch(term: str, retmax: int = 5) -> list[str]:
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax,
                                "retmode": "json", "sort": "relevance"})
    return http_json(EUTILS + "esearch.fcgi?" + q).get(
        "esearchresult", {}).get("idlist", []) or []


def esummary(pmid: str) -> dict:
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    return http_json(EUTILS + "esummary.fcgi?" + q).get("result", {}).get(pmid, {})


def pmc_id(pmid: str) -> str:
    q = urllib.parse.urlencode({"dbfrom": "pubmed", "db": "pmc", "id": pmid,
                                "retmode": "json"})
    try:
        j = http_json(EUTILS + "elink.fcgi?" + q)
        for ls in j.get("linksets", []):
            for db in ls.get("linksetdbs", []):
                if db.get("linkname", "").startswith("pubmed_pmc"):
                    ids = db.get("links", [])
                    if ids:
                        return f"PMC{ids[0]}"
    except Exception:
        pass
    return ""


# label -> (title-fragment that MUST appear in the returned title, query)
RETRY = {
    "01_HNSC_Fusobacterium_TCMA": (
        "Fusobacterium",
        'Multi-omics integration reveals the crucial role of Fusobacterium '
        '[Title] AND head and neck[Title]',
    ),
    "07_CRC_CMS": (
        "Host-microbiota interactions",
        '"Host-microbiota interactions contributing to the heterogeneous tumor '
        'microenvironment in colorectal cancer"[Title]',
    ),
}

# expected title fragment for post-hoc verification of EVERY paper
EXPECT = {
    "01_HNSC_Fusobacterium_TCMA": "Fusobacterium",
    "02_ESCC_prognosis": "esophageal squamous",
    "03_HCC_microbial_subtype": "HBV",
    "04_melanoma_CD8": "melanoma",
    "05_LUAD_TMEindex": "lung adenocarcinoma",
    "06_HCC_MRS": "hepatocellular carcinoma",
    "08_PDAC_Pseudomonas_metabolome": "pancreatic",
    "09_PDAC_MRS_TCGA_ICGC": "PDAC",
    "10_mesothelioma": "mesothelioma",
    "11_breast_bile_acid": "bile acid",
    "12_TNBC_metabolites": "TNBC",
    "13_ICC_GutMicrobes": "cholangiocarcinoma",
    "14_LUAD_glycolysis_lactate": "glycolysis",
    "16_CRC_RNAseq_reconstruction": "colon cancer",
}

SECTIONS = ("abstract", "introduction", "results", "discussion", "methods")


def main() -> None:
    rows = json.load(open("resolved_papers.json"))

    # ---- 1. retry the two bad ones -------------------------------------
    by_label = {r["label"]: r for r in rows}
    for label, (must, query) in RETRY.items():
        print(f"\n[retry] {label}")
        ids = esearch(query)
        chosen = None
        for pmid in ids:
            rec = esummary(pmid)
            ttl = (rec.get("title") or "")
            if must.lower() in ttl.lower():
                chosen = (pmid, rec)
                break
            time.sleep(0.3)
        if not chosen:
            print("   no verified match")
            continue
        pmid, rec = chosen
        doi = ""
        for aid in rec.get("articleids", []):
            if aid.get("idtype") == "doi":
                doi = aid.get("value", "")
        r = by_label[label]
        r.update({"pmid": pmid, "journal": rec.get("source", ""),
                  "year": (rec.get("pubdate", "") or "")[:4], "doi": doi,
                  "title": (rec.get("title") or "").rstrip("."),
                  "pmc": pmc_id(pmid)})
        print(f"   -> {pmid} {r['journal']} {r['year']} {r['pmc']}")
        print(f"      {r['title'][:100]}")
        time.sleep(0.4)

    # ---- 2. verify every row against the expected title fragment -------
    print("\n=== verification ===")
    dropped = []
    for r in rows:
        exp = EXPECT.get(r["label"])
        if not exp:
            continue
        ok = exp.lower() in (r.get("title") or "").lower()
        if not ok:
            print(f"  MISMATCH {r['label']}: expected '{exp}' not in title")
            dropped.append(r["label"])
    if not dropped:
        print("  all resolved titles match their expected topic")

    json.dump(rows, open("resolved_papers.json", "w"), indent=1)

    # ---- 3. fetch open-access full text --------------------------------
    print("\n=== fetching PMC full text (JATS XML) ===")
    import os
    from pathlib import Path
    out = Path("texts")
    out.mkdir(exist_ok=True)

    got = 0
    for r in rows:
        if not r.get("pmc"):
            print(f"  skip  {r['label']:32} (no PMC id)")
            continue
        pmcid = r["pmc"]
        dst = out / f"{r['label']}_{r['pmid']}_{pmcid}.xml"
        if dst.exists() and dst.stat().st_size > 5000:
            got += 1
            continue
        q = urllib.parse.urlencode({"db": "pmc", "id": pmcid, "retmode": "xml"})
        try:
            data = http(EUTILS + "efetch.fcgi?" + q, timeout=120)
            if len(data) < 2000:
                print(f"  thin  {r['label']:32} {pmcid} got {len(data)} bytes")
                continue
            dst.write_bytes(data)
            # crude section presence check so a stub record is visible
            txt = data.decode("utf-8", "replace")
            secs = [s for s in SECTIONS if f"<sec" in txt and s in txt.lower()]
            print(f"  ok    {r['label']:32} {pmcid} {len(data)//1024:5} KB  "
                  f"sections~{len(secs)}")
            got += 1
            time.sleep(0.5)
        except Exception as e:
            print(f"  ERR   {r['label']:32} {pmcid} {type(e).__name__}: {str(e)[:50]}")

    print(f"\nfull texts retrieved: {got}")


if __name__ == "__main__":
    main()
