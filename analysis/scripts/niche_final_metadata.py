"""Resolve the last open metadata questions for the OSCC source-attribution
design, via NCBI SRA (ENA exposes no attributes for PRJNA1415359):

1. PRJNA1415359: per-patient niche completeness (isolation_source is populated
   in SRA: saliva / Preoperative saliva / Postoperative saliva / tumor tissue /
   normal tissue).
2. Re-confirm PRJNA1415268 pairing (47 T+N pairs, Seoul).
3. 16S hypervariable region per core dataset, where the submitter states it.

Output: analysis/out/niche_final_metadata.txt
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "out"
UA = {"User-Agent": "wisp-science/1.0 (feasibility check)"}


def get(url: str) -> bytes:
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=180) as r:
                return r.read()
        except Exception:
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("failed: " + url)


def sra_attrs(acc: str) -> list[dict]:
    q = urllib.parse.urlencode({"db": "sra", "term": f"{acc}[BioProject]",
                                "retmax": 400, "retmode": "json"})
    ids = json.loads(get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + q))[
        "esearchresult"]["idlist"]
    out = []
    for i in range(0, len(ids), 60):
        d = urllib.parse.urlencode({"db": "sra", "id": ",".join(ids[i:i + 60]),
                                    "retmode": "xml"}).encode()
        req = urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                                     data=d, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as r:
            root = ET.fromstring(r.read())
        for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
            s = pkg.find(".//SAMPLE")
            a = {}
            if s is not None:
                for x in s.findall(".//SAMPLE_ATTRIBUTE"):
                    a[x.findtext("TAG") or ""] = x.findtext("VALUE") or ""
                a["_lib"] = (pkg.findtext(".//LIBRARY_STRATEGY") or "")
            out.append(a)
        time.sleep(0.6)
    return out


lines: list[str] = []

# ---------- PRJNA1415359 ----------
at = sra_attrs("PRJNA1415359")
lines.append("=" * 110)
lines.append(f"PRJNA1415359  packages={len(at)}")
lines.append(f"  isolation_source: {Counter(a.get('isolation_source', '') for a in at).most_common()}")
suf = Counter()
for a in at:
    sid = a.get("source_material_id", "")
    m = re.match(r"^(AP\d{1,3})-(.+)$", sid)
    suf[(m.group(2) if m else sid[:12], a.get("isolation_source", ""))] += 1
lines.append("  source_material_id suffix x isolation_source:")
for k, v in suf.most_common(25):
    lines.append(f"    {str(k):<52} n={v}")

NICHE = {
    "Preoperative saliva": "S_pre_唾液",
    "Postoperative saliva": "S_post_唾液",
    "saliva": "S_唾液",
    "tumor tissue": "T_肿瘤",
    "normal tissue": "N_正常组织",
}
by = defaultdict(set)
for a in at:
    sid = a.get("source_material_id", "")
    m = re.match(r"^(AP\d{1,3})", sid)
    if m:
        by[m.group(1)].add(NICHE.get(a.get("isolation_source", ""), a.get("isolation_source", "?")))
combo = Counter(tuple(sorted(v)) for v in by.values())
lines.append(f"  patients parsed={len(by)}")
lines.append("  niche-combination distribution:")
for k, v in combo.most_common(12):
    lines.append(f"    n={v:<3} {k}")
n_tissue = sum(1 for v in by.values() if "T_肿瘤" in v and "N_正常组织" in v)
n_tissue_sal = sum(1 for v in by.values()
                   if "T_肿瘤" in v and "N_正常组织" in v and any("唾液" in x for x in v))
lines.append(f"  含肿瘤+正常组织         : {n_tissue}")
lines.append(f"  含肿瘤+正常组织+任一唾液 : {n_tissue_sal}")

# ---------- PRJNA1415268 ----------
at2 = sra_attrs("PRJNA1415268")
by2 = defaultdict(set)
for a in at2:
    m = re.match(r"^([TN])(\d{1,3})$", a.get("source_material_id", "") or "")
    if m:
        by2[m.group(2)].add(m.group(1))
pairs = [p for p, s in by2.items() if s == {"T", "N"}]
lines.append("\n" + "=" * 110)
lines.append(f"PRJNA1415268  packages={len(at2)}  patients={len(by2)}  "
             f"完整 tumor+normal 配对={len(pairs)}")
lines.append(f"  isolation_source: {Counter(a.get('isolation_source','') for a in at2).most_common()}")

# ---------- 16S region statements ----------
lines.append("\n" + "=" * 110)
lines.append("16S 区段（仅记录提交方在元数据中明示的内容，不推断）")
REG = re.compile(r"\b(V[1-9](?:-V?[1-9])?)\b|16S[^,;.]{0,40}", re.I)
for acc in ["PRJNA866676", "PRJNA412445", "PRJNA666746", "PRJNA597251",
            "PRJNA1415359", "PRJNA1415268", "PRJNA1203737"]:
    try:
        rr = json.loads(get("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + acc +
                            "&result=read_run&fields=sample_title,sample_alias,isolate,"
                            "isolation_source,library_construction_protocol,description"
                            "&format=json&limit=0"))
        blob = " ".join(
            str(r.get(k) or "") for r in rr
            for k in ("library_construction_protocol", "description", "isolate",
                      "isolation_source", "sample_title"))
        hits = Counter(m.group(0).strip() for m in REG.finditer(blob) if m.group(0))
        lines.append(f"  {acc:<13} region-ish tokens: {dict(list(hits.most_common(6)))}")
    except Exception as e:
        lines.append(f"  {acc} FAILED {type(e).__name__} {e}")
    time.sleep(0.4)

(OUT / "niche_final_metadata.txt").write_text("\n".join(lines), encoding="utf-8")
print("\n".join(lines))
