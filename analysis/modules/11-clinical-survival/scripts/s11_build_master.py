#!/usr/bin/env python3
"""Build the master clinical + microbiome table for the survival / clinical analyses.

Sources
  GDC          TCGA-HNSC clinical fields. Two gotchas found by probing:
                 - `primary_site` is a TOP-LEVEL case field, not `diagnoses.primary_site`
                 - `demographic.gender` is empty in this API; use `demographic.sex_at_birth`
               GDC omits null fields entirely, so every field is read defensively.
  cBioPortal   HPV status, which GDC does not carry for HNSC: HPV_STATUS_ISH / HPV_STATUS_P16
               from study `hnsc_tcga`.
  TCMA         genus-level relative abundance for 177 TCGA-HNSC aliquots (S6 deliverable).

Derived per sample: Shannon, genus richness, Parvimonas relative abundance.
Survival: OS time = days_to_death if dead, else days_to_last_follow_up; event = dead.

Writes master.csv with one row per TCMA sample that has clinical data.

Run: python3 s11_build_master.py
"""
from __future__ import annotations

import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

UA = {"User-Agent": "wisp-science/1.0", "Content-Type": "application/json",
      "Accept": "application/json"}
GDC = "https://api.gdc.cancer.gov"
CBIO = "https://www.cbioportal.org/api"

TCMA_CANDIDATES = [
    Path("s6v2_tcma_genus.csv"),                      # uploaded next to this script
    Path("/home/v3724/oscc-mb/work/clin/s6v2_tcma_genus.csv"),
    Path("/home/v3724/results/tables/s6v2_tcma_genus.csv"),
]

GDC_FIELDS = [
    "submitter_id", "primary_site",
    "demographic.vital_status", "demographic.days_to_death",
    "demographic.age_at_index", "demographic.sex_at_birth", "demographic.race",
    "diagnoses.ajcc_pathologic_stage", "diagnoses.ajcc_pathologic_t",
    "diagnoses.ajcc_pathologic_n", "diagnoses.tumor_grade",
    "diagnoses.primary_diagnosis", "diagnoses.days_to_last_follow_up",
    "exposures.tobacco_smoking_status", "exposures.pack_years_smoked",
    "exposures.alcohol_history",
]


def gdc_post(path: str, payload: dict, timeout: int = 120) -> dict:
    req = urllib.request.Request(GDC + path, data=json.dumps(payload).encode(),
                                 headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def dig(rec, dotted):
    cur = rec
    for part in dotted.split("."):
        if isinstance(cur, list):
            cur = cur[0] if cur else None
        if not isinstance(cur, dict):
            return None
        cur = cur.get(part)
    return cur


def fetch_gdc() -> pd.DataFrame:
    flt = {"op": "in", "content": {"field": "project.project_id",
                                   "value": ["TCGA-HNSC"]}}
    rows, frm = [], 0
    while True:
        j = gdc_post("/cases", {"filters": flt, "fields": ",".join(GDC_FIELDS),
                                "size": 200, "from": frm, "format": "JSON"})
        hits = j["data"]["hits"]
        if not hits:
            break
        for h in hits:
            r = {"case_id": h.get("submitter_id")}
            for f in GDC_FIELDS:
                key = f.split(".")[-1]
                if key == "case_id":
                    continue
                r[key] = dig(h, f)
            rows.append(r)
        frm += len(hits)
        if frm >= j["data"]["pagination"]["total"]:
            break
        time.sleep(0.3)
    df = pd.DataFrame(rows)
    print(f"GDC cases: {len(df)}")
    return df


def fetch_cbio_hpv() -> pd.DataFrame:
    """Patient-level HPV status from cBioPortal hnsc_tcga."""
    url = f"{CBIO}/studies/hnsc_tcga/clinical-data?clinicalDataType=PATIENT"
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as r:
            recs = json.loads(r.read())
    except Exception as e:
        print(f"  cBioPortal HPV fetch failed: {type(e).__name__} {str(e)[:80]}")
        return pd.DataFrame(columns=["case_id", "hpv_ish", "hpv_p16"])

    keep = {"HPV_STATUS_ISH": "hpv_ish", "HPV_STATUS_P16": "hpv_p16"}
    rows = []
    for r in recs:
        aid = r.get("clinicalAttributeId")
        if aid in keep:
            rows.append({"case_id": r.get("patientId"), keep[aid]: r.get("value")})
    if not rows:
        print("  cBioPortal: no HPV attributes in response")
        return pd.DataFrame(columns=["case_id", "hpv_ish", "hpv_p16"])
    df = (pd.DataFrame(rows).groupby("case_id").first().reset_index())
    print(f"cBioPortal HPV rows: {len(df)} | "
          f"ish filled {df['hpv_ish'].notna().sum()} | p16 filled {df['hpv_p16'].notna().sum()}")
    return df


def shannon(row) -> float:
    v = row[row > 0]
    tot = v.sum()
    if tot <= 0:
        return float("nan")
    p = v / tot
    return float(-(p * p.map(math.log)).sum())


def main() -> None:
    gdc = fetch_gdc()
    hpv = fetch_cbio_hpv()

    tcma = None
    for c in TCMA_CANDIDATES:
        if c.is_file():
            tcma = pd.read_csv(c, index_col=0)
            print(f"TCMA table: {c}  {tcma.shape}")
            break
    if tcma is None:
        raise SystemExit("TCMA genus table not found; expected results/tables/s6v2_tcma_genus.csv")

    # per-sample microbiome metrics
    metrics = pd.DataFrame(index=tcma.index)
    metrics["shannon"] = tcma.apply(shannon, axis=1)
    metrics["richness"] = (tcma > 0).sum(axis=1)
    rel = tcma.div(tcma.sum(axis=1).replace(0, float("nan")), axis=0)
    metrics["parvimonas_rel"] = rel["Parvimonas"] if "Parvimonas" in rel.columns else float("nan")
    metrics["case_id"] = [str(s)[:12] for s in metrics.index]
    metrics["sample_barcode"] = metrics.index

    m = metrics.merge(gdc, on="case_id", how="left")
    m = m.merge(hpv, on="case_id", how="left")

    # ---- survival derivation -------------------------------------------
    vital = m["vital_status"].astype(str).str.lower()
    dtd = pd.to_numeric(m["days_to_death"], errors="coerce")
    dtf = pd.to_numeric(m["days_to_last_follow_up"], errors="coerce")
    m["os_event"] = (vital == "dead").astype(int)
    m["os_time"] = dtd.where(m["os_event"] == 1, dtf)
    m.loc[m["os_time"].isna() & (m["os_event"] == 0), "os_time"] = dtf

    # ---- tidy categoricals ---------------------------------------------
    # NOTE: no trailing \b here. GDC writes subdivisions ("Stage IVA"), and \b after the
    # roman numeral fails because the next char is a letter -- that silently dropped all
    # 74 stage-IV cases on the first pass. Alternation order keeps IV before III/II/I.
    stage_raw = m["ajcc_pathologic_stage"].astype(str)
    m["stage_roman"] = stage_raw.str.extract(r"Stage\s+(IV|III|II|I)")[0]
    m["stage"] = m["stage_roman"].map({"I": 1, "II": 2, "III": 3, "IV": 4})
    m["stage_letter"] = stage_raw.str.extract(r"Stage\s+IV|Stage\s+I{1,3}([ABC])")[0]
    m["grade_num"] = pd.to_numeric(
        m["tumor_grade"].astype(str).str.extract(r"G(\d)")[0], errors="coerce")
    m["age"] = pd.to_numeric(m["age_at_index"], errors="coerce")
    m["sex"] = m["sex_at_birth"]
    m["pack_years"] = pd.to_numeric(m["pack_years_smoked"], errors="coerce")
    smoke = m["tobacco_smoking_status"].astype(str)
    m["smoker"] = smoke.map(lambda s: (
        "never" if "non-smoker" in s.lower()
        else ("current" if s.lower().startswith("current smoker")
              else ("former" if "reformed" in s.lower() else None))))
    # HPV: prefer p16 (used clinically), fall back to ISH
    def hpv_clean(v):
        if not isinstance(v, str) or not v.strip():
            return None
        s = v.strip().lower()
        if s.startswith(("positive", "pos", "+")):
            return "positive"
        if s.startswith(("negative", "neg", "-")):
            return "negative"
        return None
    m["hpv_p16c"] = m.get("hpv_p16", pd.Series(index=m.index)).map(hpv_clean)
    m["hpv_ishc"] = m.get("hpv_ish", pd.Series(index=m.index)).map(hpv_clean)
    m["hpv"] = m["hpv_p16c"].fillna(m["hpv_ishc"])

    site = m["primary_site"].astype(str)
    m["site_group"] = site.map(lambda s: (
        "oral" if any(k in s.lower() for k in (
            "tongue", "floor of mouth", "buccal", "palate", "lip",
            "oral cavity", "gingiv", "alveolar"))
        else ("pharynx" if any(k in s.lower() for k in ("pharynx", "tonsil", "oropharynx"))
              else ("larynx" if "larynx" in s.lower() else "other"))))

    out = Path("master.csv")
    m.to_csv(out, index=False)
    print(f"\nwrote {out}: {m.shape[0]} samples x {m.shape[1]} columns")

    print("\n=== coverage ===")
    for c in ["shannon", "richness", "parvimonas_rel", "os_time", "os_event", "stage",
              "grade_num", "age", "sex", "smoker", "pack_years", "hpv", "site_group"]:
        if c in m:
            nn = m[c].notna().sum()
            print(f"  {c:16} {nn:4}/{len(m)}  {100*nn/len(m):5.1f}%")

    print("\n=== site_group ===")
    print(m["site_group"].value_counts().to_string())
    print("\n=== stage ===")
    print(m["stage"].value_counts(dropna=False).to_string())
    print("\n=== hpv ===")
    print(m["hpv"].value_counts(dropna=False).to_string())
    print("\n=== smoker ===")
    print(m["smoker"].value_counts(dropna=False).to_string())
    print(f"\nOS: {m['os_time'].notna().sum()} with time, {int(m['os_event'].sum())} deaths")


if __name__ == "__main__":
    main()
