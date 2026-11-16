#!/usr/bin/env python3
"""Discover which clinical fields GDC actually exposes for TCGA-HNSC.

Why: the analysis plan needs stage, grade, HPV status, smoking and alcohol. GDC's field set
is not uniform across projects, and asking for a field that does not exist silently returns
nulls, which would look like missing data rather than a wrong request. So probe first and
report the fill rate of every candidate field.

Run: python3 probe_clinical.py
"""
from __future__ import annotations

import json
import urllib.request

UA = {"User-Agent": "wisp-science/1.0", "Content-Type": "application/json"}
GDC = "https://api.gdc.cancer.gov"


def post(path: str, payload: dict, timeout: int = 90) -> dict:
    req = urllib.request.Request(GDC + path, data=json.dumps(payload).encode(), headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def main() -> None:
    # ---- 1. dump the FULL nested record for one case to see every key ----
    j = post("/cases", {
        "filters": {"op": "in", "content": {"field": "project.project_id",
                                            "value": ["TCGA-HNSC"]}},
        "size": 1, "format": "JSON",
    })
    rec = j["data"]["hits"][0]
    print("=== full record structure (1 case) ===")
    print(json.dumps(rec, indent=1)[:3000])

    # ---- 2. fill-rate for every candidate field across all cases ----
    candidates = [
        "submitter_id",
        "demographic.vital_status", "demographic.days_to_death",
        "demographic.age_at_index", "demographic.gender",
        "demographic.race", "demographic.ethnicity",
        "diagnoses.ajcc_pathologic_stage",
        "diagnoses.ajcc_pathologic_t", "diagnoses.ajcc_pathologic_n",
        "diagnoses.ajcc_pathologic_m",
        "diagnoses.tumor_grade", "diagnoses.primary_diagnosis",
        "diagnoses.morphology", "diagnoses.primary_site",
        "diagnoses.days_to_last_follow_up",
        "diagnoses.days_to_last_known_disease_status",
        "diagnoses.prior_malignancy", "diagnoses.Year_of_diagnosis",
        "exposures.cigarettes_per_day", "exposures.pack_years_smoked",
        "exposures.alcohol_history", "exposures.alcohol_intensity",
        "exposures.years_smoked", "exposures.tobacco_smoking_status",
    ]
    j = post("/cases", {
        "filters": {"op": "in", "content": {"field": "project.project_id",
                                            "value": ["TCGA-HNSC"]}},
        "fields": ",".join(candidates), "size": 300, "format": "JSON",
    })
    hits = j["data"]["hits"]
    total = len(hits)
    print(f"\n=== fill rate over {total} cases ===")

    def walk(d, dotted):
        cur = d
        for part in dotted.split("."):
            if isinstance(cur, list):
                if not cur:
                    return None
                cur = cur[0]
            if not isinstance(cur, dict):
                return None
            cur = cur.get(part)
            if cur is None:
                return None
        return cur

    for f in candidates:
        n = sum(1 for h in hits if walk(h, f) not in (None, "", "not reported"))
        bar = "#" * int(40 * n / total)
        print(f"  {f:44} {n:3}/{total} {100*n/total:5.1f}%  {bar}")


if __name__ == "__main__":
    main()
