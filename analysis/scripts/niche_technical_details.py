"""Round 3 of the feasibility check: the details that decide whether source
attribution is technically valid.

1. PRJNA1415359 tissue-sample alias format (the 161 unparsed runs).
2. PRJNA1415268 (94 runs, paired tumor/adjacent) sample metadata via NCBI SRA,
   because ENA exposes no BioSample attributes for it.
3. Sequencing region + read depth per core dataset: source tracking needs
   comparable depth and a comparable 16S region WITHIN a dataset (sources and
   sink must share the protocol).
4. Whether the core dataset has negative controls / blanks.

Output: analysis/out/niche_technical_details.txt
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
lines: list[str] = []


def get(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=180) as r:
        return r.read()


# ---------- 1. PRJNA1415359 alias formats ----------
F = ("run_accession,sample_alias,sample_title,isolate,isolation_source,"
     "library_construction_protocol,library_strategy,library_source,"
     "read_count,fastq_bytes,instrument_model")
rows = json.loads(get("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1415359"
                      "&result=read_run&fields=" + F + "&format=json&limit=0"))
lines.append("=" * 112)
lines.append("1. PRJNA1415359 — alias 形态与样本类型")
lines.append("=" * 112)
pat = Counter()
for r in rows:
    a = (r.get("sample_alias") or "").strip()
    a2 = re.sub(r"\d+", "#", a)
    pat[(a2, (r.get("isolation_source") or "").strip())] += 1
for (a2, lab), n in pat.most_common(20):
    lines.append(f"  {n:>4}  pattern={a2:<22} isolation_source={lab!r}")
lines.append(f"  library_strategy: {Counter(r.get('library_strategy') for r in rows)}")
prot = list(dict.fromkeys((r.get("library_construction_protocol") or "") for r in rows))
lines.append(f"  library_construction_protocol (distinct): {prot[:4]}")

# ---------- 2. PRJNA1415268 via NCBI SRA ----------
lines.append("\n" + "=" * 112)
lines.append("2. PRJNA1415268 — 通过 NCBI SRA 取样本元数据")
lines.append("=" * 112)
try:
    q = urllib.parse.urlencode({"db": "sra", "term": "PRJNA1415268[BioProject]",
                                "retmax": 200, "retmode": "json"})
    ids = json.loads(get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + q))[
        "esearchresult"]["idlist"]
    lines.append(f"  SRA UIDs: {len(ids)}")
    got = []
    for i in range(0, len(ids), 100):
        d = urllib.parse.urlencode({"db": "sra", "id": ",".join(ids[i:i + 100]),
                                    "retmode": "xml"}).encode()
        req = urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                                     data=d, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as r:
            root = ET.fromstring(r.read())
        for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
            s = pkg.find(".//SAMPLE")
            title = (s.findtext("TITLE") or "") if s is not None else ""
            attrs = {}
            if s is not None:
                for a in s.findall(".//SAMPLE_ATTRIBUTE"):
                    attrs[a.findtext("TAG") or ""] = a.findtext("VALUE") or ""
            got.append({"title": title, "attrs": attrs})
        time.sleep(0.4)
    lines.append(f"  packages parsed: {len(got)}")
    tag = Counter()
    for g in got:
        for k in g["attrs"]:
            tag[k] += 1
    lines.append(f"  attribute tags: {dict(tag)}")
    for g in got[:8]:
        lines.append(f"    title={g['title'][:40]!r} attrs={ {k: v[:40] for k, v in g['attrs'].items()} }")
except Exception as e:
    lines.append(f"  FAILED {type(e).__name__} {e}")

# ---------- 3. depth + region per dataset ----------
lines.append("\n" + "=" * 112)
lines.append("3. 核心数据集的读段深度（reads/run）与测序策略")
lines.append("=" * 112)
CORE = ["PRJNA866676", "PRJNA412445", "PRJNA666746", "PRJNA597251",
        "PRJNA666891", "PRJNA1203737"]
for acc in CORE:
    try:
        rr = json.loads(get("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + acc +
                            "&result=read_run&fields=run_accession,read_count,library_strategy,"
                            "library_source,instrument_model,fastq_bytes&format=json&limit=0"))
        rc = sorted(int(x.get("read_count") or 0) for x in rr)
        if not rc:
            continue
        q = lambda p: rc[min(len(rc) - 1, int(p * len(rc)))]
        lines.append(f"  {acc:<13} n={len(rc):>4}  reads min={rc[0]:>9,}  p25={q(.25):>9,}  "
                     f"median={q(.5):>9,}  p75={q(.75):>9,}  max={rc[-1]:>10,}")
        lines.append(f"                 strategy={Counter(x.get('library_strategy') for x in rr)}")
        lines.append(f"                 platform={Counter(x.get('instrument_model') for x in rr)}")
    except Exception as e:
        lines.append(f"  {acc} FAILED {type(e).__name__} {e}")
    time.sleep(0.35)

# ---------- 4. negative controls ----------
lines.append("\n" + "=" * 112)
lines.append("4. 是否存在阴性对照 / 空白（关键字扫描 sample_title / alias / isolate）")
lines.append("=" * 112)
KW = re.compile(r"blank|negative|control\s*blank|extraction|nuclease|"
                r"pbs|water|mock|kit", re.I)
for acc in CORE:
    try:
        rr = json.loads(get("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + acc +
                            "&result=read_run&fields=sample_alias,sample_title,isolate,"
                            "isolation_source&format=json&limit=0"))
        hits = [x for x in rr if KW.search(" ".join(str(x.get(k) or "") for k in
                                                   ("sample_alias", "sample_title", "isolate",
                                                    "isolation_source")))]
        lines.append(f"  {acc:<13} 命中阴性对照关键词的 run 数 = {len(hits)}")
        for h in hits[:5]:
            lines.append(f"      {h}")
    except Exception as e:
        lines.append(f"  {acc} FAILED {type(e).__name__} {e}")
    time.sleep(0.35)

(OUT / "niche_technical_details.txt").write_text("\n".join(lines), encoding="utf-8")
print("\n".join(lines))
