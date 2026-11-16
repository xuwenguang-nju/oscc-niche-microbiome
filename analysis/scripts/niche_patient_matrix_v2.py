"""Definitive niche x patient feasibility matrix for the OSCC source-tracking
design, using each study's OWN naming convention (decoded from the observed
aliases) instead of a generic regex.

Decoded conventions (validated by cross-tabulating alias prefix against the
submitter's own isolate / isolation_source label, printed in section A):

  PRJNA866676  (C|CR|CS|N|S|T)<patient>          -> 6 niches, one numbering
  PRJNA412445  <patient>(a|b|c|d|e)              -> 5 niches, one numbering
  PRJNA666746  title "OSCC <n>" / "Normal <n>"   -> tumor + matched control
  PRJNA597251  CT<n> / PT<n>                     -> tumor + paracancerous
  PRJNA1415359 AP<n>-<...>                       -> tissue / pre- & post-op saliva
  PRJNA666891  P<patient>_TT|_SAL|_HT            -> tumor / saliva / healthy tissue

Output: analysis/out/niche_patient_matrix_v2.txt + .json
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
FIELDS = ("run_accession,sample_accession,sample_title,sample_alias,isolate,"
          "isolation_source,read_count,fastq_bytes")


def get(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=180) as r:
        return r.read()


def runs(acc: str) -> list[dict]:
    u = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + acc +
         "&result=read_run&fields=" + FIELDS + "&format=json&limit=0")
    return json.loads(get(u))


def show_prefix_table(acc, rows, alias_re, key_label="prefix"):
    """Cross-tabulate alias prefix against the submitter's own niche label."""
    tab = defaultdict(Counter)
    for r in rows:
        m = alias_re.search(r.get("sample_alias", "") or "")
        if not m:
            continue
        pref = m.group(1)
        label = (r.get("isolate") or r.get("isolation_source") or
                 r.get("sample_title") or "").strip()
        tab[pref][label] += 1
    return {p: dict(c) for p, c in sorted(tab.items())}


lines: list[str] = []
report: dict = {}

# ============================ 1. PRJNA866676 ============================
rows = runs("PRJNA866676")
ALIAS_866676 = re.compile(r"^([A-Z]{1,2})(\d{1,3})$")
tab = show_prefix_table("PRJNA866676", rows, ALIAS_866676)
by_patient = defaultdict(dict)
for r in rows:
    m = ALIAS_866676.match((r.get("sample_alias") or "").strip())
    if not m:
        continue
    by_patient[m.group(2)][m.group(1)] = r.get("read_count", "")
niche_names = {p: list(v)[0] for p, v in tab.items()}
full6 = [p for p, d in by_patient.items() if len(d) >= 6]
has_tumor = [p for p, d in by_patient.items() if "C" in d]
tumor_saliva = [p for p, d in by_patient.items() if "C" in d and "S" in d]
tumor_saliva_normal = [p for p, d in by_patient.items()
                       if {"C", "S"} <= d.keys() and ({"N", "CR"} & d.keys())]
tumor_saliva_normal_plaque = [p for p, d in by_patient.items()
                              if {"C", "S"} <= d.keys() and ({"N", "CR"} & d.keys())
                              and "T" in d]
report["PRJNA866676"] = {
    "runs": len(rows),
    "alias_prefix_to_submitter_label": niche_names,
    "prefix_counts": {p: sum(c.values()) for p, c in tab.items()},
    "n_patients": len(by_patient),
    "n_patients_all6": len(full6),
    "n_patients_with_tumor": len(has_tumor),
    "n_patients_tumor_plus_saliva": len(tumor_saliva),
    "n_patients_tumor_saliva_normal": len(tumor_saliva_normal),
    "n_patients_tumor_saliva_normal_tonguecoat": len(tumor_saliva_normal_plaque),
}
lines.append("=" * 116)
lines.append("1. PRJNA866676 (山东大学口腔医院, 济南) — 别名前缀 = 生态位, 数字 = 患者")
lines.append("=" * 116)
lines.append("  prefix -> submitter's own label (count):")
for p, c in tab.items():
    lines.append(f"    {p:<3} -> {c}")
lines.append(f"  runs={len(rows)}  distinct patient numbers={len(by_patient)}")
lines.append(f"  患者拥有全部 6 个生态位 : {len(full6)}")
lines.append(f"  含肿瘤(C)              : {len(has_tumor)}")
lines.append(f"  含肿瘤 + 唾液(C+S)      : {len(tumor_saliva)}")
lines.append(f"  含肿瘤 + 唾液 + 正常黏膜(C+S+N/CR): {len(tumor_saliva_normal)}")
lines.append(f"  ... 再加舌苔(T)         : {len(tumor_saliva_normal_plaque)}")
ex = [(p, sorted(d)) for p, d in by_patient.items() if len(d) >= 5][:5]
lines.append(f"  示例（患者号 -> 已有生态位）: {ex}")

# ============================ 2. PRJNA412445 ============================
rows = runs("PRJNA412445")
ALIAS_412445 = re.compile(r"^(\d{1,3})([a-e])$")
tab2 = defaultdict(Counter)
by_patient2 = defaultdict(dict)
for r in rows:
    m = ALIAS_412445.match((r.get("sample_alias") or "").strip())
    if not m:
        continue
    pid, niche = m.group(1), m.group(2)
    by_patient2[pid][niche] = r.get("read_count", "")
    tab2[niche][(r.get("isolation_source") or "").strip()] += 1
report["PRJNA412445"] = {
    "runs": len(rows),
    "letter_to_submitter_label": {k: dict(v) for k, v in sorted(tab2.items())},
    "n_patients": len(by_patient2),
    "n_patients_4plus": sum(1 for d in by_patient2.values() if len(d) >= 4),
    "n_patients_has_e_tumor": sum(1 for d in by_patient2.values() if "e" in d),
    "n_patients_e_plus_d_saliva": sum(1 for d in by_patient2.values()
                                      if {"e", "d"} <= d.keys()),
    "n_patients_e_d_a_plaque": sum(1 for d in by_patient2.values()
                                   if {"e", "d", "a"} <= d.keys()),
    "n_patients_e_d_a_b_c": sum(1 for d in by_patient2.values() if len(d) >= 5),
}
lines.append("\n" + "=" * 116)
lines.append("2. PRJNA412445 (四川大学) — 别名 = 患者号 + 字母(生态位)")
lines.append("=" * 116)
lines.append("  letter -> submitter's own isolation_source (count):")
for k, v in sorted(tab2.items()):
    lines.append(f"    {k} -> {dict(v)}")
lines.append(f"  runs={len(rows)}  distinct patients={len(by_patient2)}")
lines.append(f"  患者拥有 >=4 个生态位            : {report['PRJNA412445']['n_patients_4plus']}")
lines.append(f"  含 e(tumor)                     : {report['PRJNA412445']['n_patients_has_e_tumor']}")
lines.append(f"  含 e(tumor) + d(saliva)         : {report['PRJNA412445']['n_patients_e_plus_d_saliva']}")
lines.append(f"  含 e + d + a(plaque)            : {report['PRJNA412445']['n_patients_e_d_a_plaque']}")
lines.append(f"  含全部 5 个 (a,b,c,d,e)         : {report['PRJNA412445']['n_patients_e_d_a_b_c']}")
lines.append(f"  示例: {[(p, sorted(d)) for p, d in list(by_patient2.items())[:6]]}")

# ============================ 3. PRJNA666746 ============================
rows = runs("PRJNA666746")
pat = defaultdict(set)
for r in rows:
    m = re.search(r"(OSCC|Normal)\s*(\d{1,3})", r.get("sample_title") or "")
    if m:
        pat[m.group(2)].add(m.group(1))
pairs = [p for p, s in pat.items() if s == {"OSCC", "Normal"}]
report["PRJNA666746"] = {"runs": len(rows), "n_patients": len(pat),
                         "n_complete_pairs": len(pairs)}
lines.append("\n" + "=" * 116)
lines.append("3. PRJNA666746 (印度) — title 'OSCC <n>' / 'Normal <n>'")
lines.append("=" * 116)
lines.append(f"  runs={len(rows)}  patients={len(pat)}  完整配对(OSCC+Normal)={len(pairs)}")

# ============================ 4. PRJNA597251 ============================
rows = runs("PRJNA597251")
pat = defaultdict(set)
for r in rows:
    m = re.match(r"^(CT|PT)(\d{1,3})$", (r.get("sample_alias") or "").strip())
    if m:
        pat[m.group(2)].add(m.group(1))
pairs = [p for p, s in pat.items() if s == {"CT", "PT"}]
report["PRJNA597251"] = {"runs": len(rows), "n_patients": len(pat),
                         "n_complete_pairs": len(pairs)}
lines.append("\n" + "=" * 116)
lines.append("4. PRJNA597251 (青岛) — CT<n>=癌灶 / PT<n>=配对癌旁")
lines.append("=" * 116)
lines.append(f"  runs={len(rows)}  patients={len(pat)}  完整配对={len(pairs)}")

# ============================ 5. PRJNA1415359 ============================
rows = runs("PRJNA1415359")
suffix = Counter()
for r in rows:
    a = (r.get("sample_alias") or "").strip()
    label = (r.get("isolation_source") or "").strip()
    m = re.match(r"^(AP\d{1,3})-(.+)$", a)
    suffix[(m.group(2) if m else a[:14], label)] += 1
by_pat = defaultdict(set)
for r in rows:
    a = (r.get("sample_alias") or "").strip()
    m = re.match(r"^(AP\d{1,3})-(.+)$", a)
    if m:
        by_pat[m.group(1)].add(m.group(2))
n_tissue_pat = sum(1 for v in by_pat.values() if any("t" in x for x in v))
report["PRJNA1415359"] = {
    "runs": len(rows),
    "suffix_to_label": {f"{k[0]} | {k[1]}": v for k, v in sorted(suffix.items())},
    "n_patients": len(by_pat),
    "n_patients_with_tissue_and_saliva": sum(
        1 for v in by_pat.values()
        if any("s" in x for x in v) and any("t" in x for x in v)),
    "example": {k: sorted(v) for k, v in list(by_pat.items())[:5]},
}
lines.append("\n" + "=" * 116)
lines.append("5. PRJNA1415359 — 别名 = AP<患者号>-<样本类型后缀>")
lines.append("=" * 116)
lines.append("  suffix | isolation_source (count):")
for (s, lab), n in sorted(suffix.items())[:20]:
    lines.append(f"    {s:<16} | {lab:<18} n={n}")
lines.append(f"  runs={len(rows)}  patients={len(by_pat)}")
lines.append(f"  同时有组织与唾液的患者: {report['PRJNA1415359']['n_patients_with_tissue_and_saliva']}")
lines.append(f"  示例: {report['PRJNA1415359']['example']}")

# ============================ 6. PRJNA666891 ============================
rows = runs("PRJNA666891")
tab6 = defaultdict(Counter)
by_pat6 = defaultdict(set)
for r in rows:
    a = (r.get("sample_alias") or "").strip()
    m = re.match(r"^(P\d{4})_(\w+)$", a)
    if m:
        tab6[m.group(2)][""] += 1
        by_pat6[m.group(1)].add(m.group(2))
reads = [int(r.get("read_count") or 0) for r in rows]
report["PRJNA666891"] = {
    "runs": len(rows),
    "suffixes": {k: sum(v.values()) for k, v in tab6.items()},
    "n_patients": len(by_pat6),
    "n_patients_with_3_types": sum(1 for v in by_pat6.values() if len(v) >= 3),
    "median_reads": sorted(reads)[len(reads) // 2] if reads else 0,
    "max_reads": max(reads) if reads else 0,
}
lines.append("\n" + "=" * 116)
lines.append("6. PRJNA666891 (华沙) — P<患者号>_TT(肿瘤)/_SAL(唾液)/_HT(健康组织)")
lines.append("=" * 116)
lines.append(f"  runs={len(rows)}  suffix counts={dict(tab6)}")
lines.append(f"  patients={len(by_pat6)}  拥有 3 种样本的患者={report['PRJNA666891']['n_patients_with_3_types']}")
lines.append(f"  read_count 中位数={report['PRJNA666891']['median_reads']} 最大={report['PRJNA666891']['max_reads']}"
             "   <-- 注意: 读段数极低, 说明这不是常规深度 16S, 需谨慎")

(OUT / "niche_patient_matrix_v2.txt").write_text("\n".join(lines), encoding="utf-8")
(OUT / "niche_patient_matrix_v2.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")
print("\n".join(lines))
