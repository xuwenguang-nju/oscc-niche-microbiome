"""Feasibility check for a niche-resolved source-attribution study of the
intratumoral microbiome in OSCC, using public data only.

Question answered here: do public OSCC oral-microbiome projects actually
release SAMPLE-LEVEL metadata that (a) labels the anatomical niche and
(b) carries a patient identifier, so that tumor / saliva / contralateral
normal mucosa / subgingival plaque can be matched within the same patient?

Method: pull every read_run record in ENA for each candidate study, plus the
linked BioSample attributes, and count how many distinct sample attributes are
populated. Everything printed is a live index response.

Output: analysis/out/niche_metadata_census.txt  (+ .json)
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "out"
OUT.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "wisp-science/1.0 (feasibility check)"}

# Projects that could contain multiple oral niches for OSCC / HNSCC
STUDIES = {
    "PRJNA866676": "OSCC 六个口腔生态位（肿瘤/癌旁/癌表/对侧正常黏膜/唾液/舌苔）",
    "PRJNA412445": "OSCC 五生态位（龈下菌斑/病灶侧/对侧/唾液/瘤内），25 例",
    "PRJNA1203737": "OSCC 组织+唾液+口腔拭子",
    "PRJNA666891": "OSCC 肿瘤/非肿瘤组织+唾液",
    "PRJNA1415359": "OSCC 肿瘤+癌旁+术前术后唾液",
    "PRJNA1415268": "OSCC 配对肿瘤+癌旁",
    "PRJNA666746": "OSCC 50 例恶性灶+配对癌旁",
    "PRJNA597251": "OSCC 24 例癌灶+配对癌旁",
    "PRJNA822685": "OSCC 口腔菌群+宿主转录组",
    "PRJNA1475993": "唾液：健康/牙结石/牙周炎/OSCC",
    "PRJNA1482781": "唾液 OSCC+健康大队列",
    "PRJEB87426": "HPV+/- HNSCC 肿瘤",
}

TISSUE_HINTS = re.compile(
    r"tumor|tumour|cancer|carcinoma|malign|adjacent|paracancer|normal|control|"
    r"saliva|swab|plaque|gingiv|mucosa|tongue|coat|lesion|biopsy|tissue",
    re.I,
)


def get(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=180) as r:
        return r.read()


report: dict = {}
lines: list[str] = []

for acc, label in STUDIES.items():
    entry = {"label": label, "runs": 0, "attributes": {}, "niche_like_values": {},
             "error": None}
    try:
        u = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + acc +
             "&result=read_run&fields=run_accession,sample_accession,sample_title,"
             "sample_alias,description,experiment_title,library_strategy,library_source,"
             "instrument_model,read_count,fastq_ftp&format=json&limit=0")
        rows = json.loads(get(u))
        entry["runs"] = len(rows)
        attr_cover = Counter()
        values_by_attr = defaultdict(Counter)
        for r in rows:
            for k, v in r.items():
                if v not in ("", None):
                    attr_cover[k] += 1
            for field in ("sample_title", "sample_alias", "description", "experiment_title"):
                v = (r.get(field) or "").strip()
                if v:
                    values_by_attr[field][v] += 1
        entry["attributes"] = dict(attr_cover)
        for field, ctr in values_by_attr.items():
            sample_vals = [k for k in ctr if TISSUE_HINTS.search(k)]
            entry["niche_like_values"][field] = {
                "n_distinct_values": len(ctr),
                "n_value_strings_matching_niche_words": len(sample_vals),
                "examples": sample_vals[:25],
            }
    except Exception as e:
        entry["error"] = f"{type(e).__name__}: {e}"
    report[acc] = entry
    time.sleep(0.4)

lines.append("=" * 118)
lines.append("A. run-level (ENA read_run) metadata coverage")
lines.append("=" * 118)
for acc, e in report.items():
    lines.append(f"\n{acc}  runs={e['runs']}  :: {e['label']}")
    if e["error"]:
        lines.append(f"  ERROR {e['error']}")
        continue
    for field, info in e["niche_like_values"].items():
        if info["n_value_strings_matching_niche_words"] == 0:
            continue
        lines.append(f"  [{field}] {info['n_distinct_values']} distinct values, "
                     f"{info['n_value_strings_matching_niche_words']} contain niche words")
        for ex in info["examples"]:
            lines.append(f"      - {ex[:110]}")

# ---- B. BioSample attribute fields (the real determinant of pairing) ----
lines.append("\n" + "=" * 118)
lines.append("B. BioSample attribute schema per study (does it carry niche + host/subject id?)")
lines.append("=" * 118)
BS = ["PRJNA866676", "PRJNA412445", "PRJNA1203737", "PRJNA666891",
      "PRJNA1415359", "PRJNA1415268", "PRJNA666746", "PRJNA597251"]
bs_report = {}
for acc in BS:
    try:
        u = ("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run"
             "&fields=sample_accession&format=json&limit=0&query="
             + urllib.parse.quote(f'study_accession="{acc}"'))
        rows = json.loads(get(u))
        sacc = sorted({r["sample_accession"] for r in rows if r.get("sample_accession")})
        fields_seen = Counter()
        example_attrs = None
        # inspect up to 8 samples per study
        for s in sacc[:8]:
            try:
                su = ("https://www.ebi.ac.uk/biosamples/samples/" + s + ".json")
                j = json.loads(get(su))
                ch = j.get("characteristics", {})
                for k in ch:
                    fields_seen[k] += 1
                if example_attrs is None:
                    example_attrs = {k: [x.get("text") for x in v] for k, v in ch.items()}
            except Exception:
                pass
            time.sleep(0.15)
        bs_report[acc] = {"n_samples": len(sacc), "fields": dict(fields_seen),
                          "example": example_attrs}
        lines.append(f"\n{acc}  samples={len(sacc)}")
        lines.append(f"  attribute fields seen: {sorted(fields_seen)}")
        if example_attrs:
            for k, v in list(example_attrs.items())[:18]:
                lines.append(f"      {k} = {v}")
    except Exception as e:
        lines.append(f"\n{acc}  BioSample lookup FAILED {type(e).__name__} {e}")
    time.sleep(0.4)

(OUT / "niche_metadata_census.txt").write_text("\n".join(lines), encoding="utf-8")
(OUT / "niche_metadata_census.json").write_text(
    json.dumps({"run_level": report, "biosample_level": bs_report}, ensure_ascii=False, indent=1),
    encoding="utf-8")
print("\n".join(lines))
