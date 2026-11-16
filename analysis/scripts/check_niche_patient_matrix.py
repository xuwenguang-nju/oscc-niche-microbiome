"""Decisive feasibility check: can tumor / adjacent-normal / contralateral
normal mucosa / saliva / (subgingival) plaque be MATCHED WITHIN PATIENT in the
public OSCC 16S projects?

Pulls per-run niche fields (isolate, isolation_source, host_body_site,
tissue_type, sample_title, sample_alias) for each candidate study, normalises
them into niche classes, extracts a patient identifier with a set of explicit
patterns, and reports the niche x patient completeness matrix.

No value is invented: every niche string printed is the submitter's own text.

Output: analysis/out/niche_patient_matrix.txt + .json
"""

from __future__ import annotations

import json
import re
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "out"
UA = {"User-Agent": "wisp-science/1.0 (feasibility check)"}
FIELDS = ("run_accession,sample_accession,sample_title,sample_alias,sample_description,"
          "isolate,isolation_source,host_body_site,tissue_type,experiment_title,read_count")

STUDIES = {
    "PRJNA866676": "OSCC 6 生态位（济南，山东大学口腔医院）",
    "PRJNA412445": "OSCC 5 生态位含龈下菌斑（四川大学），25 例",
    "PRJNA1203737": "OSCC 组织+唾液+口腔拭子（上海九院）",
    "PRJNA666891": "OSCC 舌 肿瘤/非肿瘤+唾液（华沙）",
    "PRJNA597251": "OSCC 24 例癌灶+配对癌旁（青岛）",
    "PRJNA666746": "OSCC 50 例恶性灶+配对癌旁（印度）",
    "PRJNA1415359": "OSCC 肿瘤+癌旁+术前术后唾液",
    "PRJNA1415268": "OSCC 配对肿瘤+癌旁",
}

# niche classification: ordered rules, first match wins
NICHE_RULES = [
    ("T_肿瘤", r"tumor|tumour|cancer(?!\s*surface)|carcinom|malign|lesion|carcinoma"),
    ("CS_癌表", r"cancer\s*surface|surface|ulcer"),
    ("P_癌旁", r"paracancer|peri-?tumou?r|adjacent|tumor.?adjacent|para-?carcinom"),
    ("N_对侧正常黏膜", r"contralateral|control\s*side|opposite|normal\s*mucosa|"
                    r"matched\s*control|healthy\s*mucosa"),
    ("S_唾液", r"saliva|spit|oral\s*rinse|rinse"),
    ("SW_拭子", r"swab"),
    ("TC_舌苔", r"tongue\s*coat|coated\s*tongue|tongue\s*coating"),
    ("DP_龈下菌斑", r"plaque|subgingiv|gingival"),
]
COMP = [(k, re.compile(v, re.I)) for k, v in NICHE_RULES]

# patient-id patterns, tried in order on the concatenated sample labels
PID_PATTERNS = [
    re.compile(r"(?:patient|pt|subject|case)[\s_\-:]*(\d{1,3})", re.I),
    re.compile(r"\bP[\s_\-]?(\d{1,3})(?=[\s_\-]|$|[A-Z])"),
    re.compile(r"\b(?:OSCC|HNSCC|TSCC|OPC|OC)[\s_\-]?(\d{1,3})\b", re.I),
    re.compile(r"\b(?:Tissue|Saliva|Swab|Sample)[A-Za-z]{0,3}[\s_\-]?(\d{1,3})\b", re.I),
    re.compile(r"\b(?:Normal|Tumor|Tumour|Cancer)[\s_\-]?(\d{1,3})\b", re.I),
    re.compile(r"\bN(\d{1,3})[_\s\-]?(?:Matched|Control)?", re.I),
    re.compile(r"(\d{1,3})"),
]


def get(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=180) as r:
        return r.read()


def classify(text: str) -> str:
    for name, pat in COMP:
        if pat.search(text):
            return name
    return "OTHER"


def patient_id(text: str) -> str:
    for pat in PID_PATTERNS:
        m = pat.search(text)
        if m:
            return m.group(1)
    return ""


report: dict = {}
lines: list[str] = []

for acc, label in STUDIES.items():
    u = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + acc +
         "&result=read_run&fields=" + FIELDS + "&format=json&limit=0")
    try:
        rows = json.loads(get(u))
    except Exception as e:
        lines.append(f"\n{acc}: ENA FAILED {type(e).__name__} {e}")
        continue
    recs = []
    for r in rows:
        # the richest niche-bearing string available for this sample
        niche_text = " | ".join(
            str(r.get(k) or "") for k in
            ("isolate", "isolation_source", "host_body_site", "tissue_type",
             "sample_title", "sample_alias", "sample_description"))
        id_text = " ".join(str(r.get(k) or "") for k in
                           ("sample_title", "sample_alias", "isolate", "isolation_source"))
        recs.append({
            "run": r.get("run_accession", ""),
            "sample": r.get("sample_accession", ""),
            "niche": classify(niche_text),
            "niche_text": niche_text[:160],
            "pid": patient_id(id_text),
            "reads": r.get("read_count", ""),
        })
    niche_counts = Counter(x["niche"] for x in recs)
    by_pid = defaultdict(set)
    for x in recs:
        if x["pid"]:
            by_pid[x["pid"]].add(x["niche"])
    multi = {p: sorted(v) for p, v in by_pid.items() if len(v) >= 2}
    with_tumor = [p for p, v in by_pid.items() if "T_肿瘤" in v]
    tumor_plus_saliva = [p for p, v in by_pid.items()
                         if "T_肿瘤" in v and "S_唾液" in v]
    tumor_plus_normal = [p for p, v in by_pid.items()
                         if "T_肿瘤" in v and ("N_对侧正常黏膜" in v or "P_癌旁" in v)]
    tumor_saliva_normal = [p for p, v in by_pid.items()
                           if {"T_肿瘤", "S_唾液"} <= v and
                           ({"N_对侧正常黏膜", "P_癌旁"} & v)]
    report[acc] = {
        "label": label, "runs": len(recs),
        "niche_counts": dict(niche_counts),
        "n_patients_with_id": len(by_pid),
        "n_patients_multiple_niches": len(multi),
        "n_patients_with_tumor": len(with_tumor),
        "n_tumor_plus_saliva": len(tumor_plus_saliva),
        "n_tumor_plus_normal": len(tumor_plus_normal),
        "n_tumor_plus_saliva_plus_normal": len(tumor_saliva_normal),
        "examples_multiple_niche_patients": dict(list(multi.items())[:6]),
    }
    lines.append(f"\n{'='*112}\n{acc} :: {label}   runs={len(recs)}")
    lines.append(f"  niche distribution: {dict(niche_counts)}")
    lines.append(f"  patients with a parsed ID           : {len(by_pid)}")
    lines.append(f"  ... with >=2 distinct niches        : {len(multi)}")
    lines.append(f"  ... with tumor                     : {len(with_tumor)}")
    lines.append(f"  ... with tumor + saliva            : {len(tumor_plus_saliva)}")
    lines.append(f"  ... with tumor + (normal|adjacent) : {len(tumor_plus_normal)}")
    lines.append(f"  ... with tumor + saliva + normal   : {len(tumor_saliva_normal)}")
    if multi:
        lines.append("  example matched patients (pid -> niches):")
        for p, v in list(multi.items())[:8]:
            lines.append(f"    {p} -> {sorted(v)}")
    lines.append("  raw niche strings (up to 12 distinct):")
    for s in list(dict.fromkeys(x["niche_text"] for x in recs))[:12]:
        lines.append(f"    [{classify(s)}] {s}")
    time.sleep(0.4)

(OUT / "niche_patient_matrix.txt").write_text("\n".join(lines), encoding="utf-8")
(OUT / "niche_patient_matrix.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")
print("\n".join(lines))
