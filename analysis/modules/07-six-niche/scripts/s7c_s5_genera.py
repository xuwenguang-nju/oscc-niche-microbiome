#!/usr/bin/env python3
"""S7c: do the S5-replicating genera behave consistently in the six-niche cohort?

WHY
---
PRJNA866676 is the only cohort in the project with a tumour/normal contrast AND a full
set of oral source niches. S5 established that ten genera replicate in direction across
four independent cohorts. This script asks two things in the new cohort, within patient:

  (a) tumour vs paracancerous  -- the closest thing to S5's tumour/adjacent contrast
  (b) tumour vs each oral source niche -- where does each genus actually live?

The second question is the one no other cohort can answer, and it explains any genus
whose S5 direction does not hold: a genus that is simply more abundant in saliva will
look "up in tumour" in a cohort whose normal tissue is saliva-contaminated, and vice
versa. That is a mechanism, not an excuse.

Genus names are matched through the same synonym map used for TCMA, because SILVA and
older references differ on split genera.

INPUT  genus table + manifest
OUTPUT results/tables/s7c_s5_genera_by_niche.csv, s7c_s5_paired_tests.csv, s7c_summary.json
Run: TCMA_ROOT=<dir> python3 s7c_s5_genera.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SYNONYM = {"Cutibacterium": "Propionibacterium",
           "Schaalia": "Actinobaculum",
           "Lancefieldella": "Atopobium"}
SOURCES = ["paracancerous", "cancer_surface", "contralateral_normal_mucosa",
           "saliva", "tongue_coat"]
CONTAMINANTS = {
    "Stenotrophomonas", "Pseudomonas", "Bradyrhizobium", "Acinetobacter",
    "Sphingomonas", "Microbacterium", "Corynebacterium", "Flavobacterium",
    "Delftia", "Comamonas", "Methylobacterium", "Ralstonia", "Burkholderia",
    "Phyllobacterium", "Bosea", "Sphingobacterium", "Chryseobacterium",
    "Variovorax", "Cupriavidus", "Aquabacterium", "Methylophilus", "Massilia",
    "Hymenobacter", "Sphingopyxis", "Herbaspirillum", "Brevundimonas",
    "Novosphingobium", "Devosia", "Rhizobium", "Mesorhizobium", "Aeromonas",
}


def locate_root() -> Path:
    env = os.environ.get("TCMA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results" / "tables").is_dir():
            return cand
    return here.parents[3]


ROOT = locate_root()
TAB = ROOT / "results" / "tables"


def bh(p):
    p = np.asarray(p, float)
    n = p.size
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    gtab = pd.read_csv(ROOT / "results" / "dada2" / "PRJNA866676"
                       / "genus_table_silva.csv", index_col=0)
    man = pd.read_csv(TAB / "oscc_niche_sample_manifest.csv", encoding="utf-8-sig")
    man = man[man["study"] == "PRJNA866676"][["run", "niche", "patient_id"]].set_index("run")
    meta = man.loc[man.index.intersection(gtab.index)]
    rel = gtab.loc[meta.index]
    rel = rel.div(rel.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    s5 = pd.read_csv(TAB / "s5_meta_analysis.csv")
    k = int(s5["n_cohorts"].max())
    rep = s5[(s5["n_cohorts"] == k) & (s5["pooled_q_bh"] < 0.05) & s5["direction_consistent"]]
    rep = rep[["genus", "pooled_effect"]].copy()
    rep["s5_direction"] = np.where(rep["pooled_effect"] > 0, "up_in_tumour", "down_in_tumour")

    genera = set(rel.columns)
    resolved = {}
    for g in rep["genus"]:
        if g in genera:
            resolved[g] = g
        elif SYNONYM.get(g) in genera:
            resolved[g] = SYNONYM[g]
        else:
            resolved[g] = None

    # --- (a) mean abundance by niche, decontaminated composition ------------
    drop = [c for c in rel.columns if c in CONTAMINANTS]
    rel_d = rel.drop(columns=drop)
    rel_d = rel_d.div(rel_d.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    tmp = rel_d.copy()
    tmp["__niche"] = meta["niche"].reindex(tmp.index)
    niche_mean = tmp.groupby("__niche").mean()

    rows = []
    for g in rep["genus"]:
        col = resolved[g]
        rec = {"genus": g, "resolved_column": col,
               "s5_direction": rep.loc[rep["genus"] == g, "s5_direction"].iloc[0]}
        for n in ["tumor"] + SOURCES:
            rec[n] = round(100 * float(niche_mean.loc[n, col]), 3) if col and n in niche_mean.index else np.nan
        rows.append(rec)
    byniche = pd.DataFrame(rows)
    byniche.to_csv(TAB / "s7c_s5_genera_by_niche.csv", index=False)

    # --- (b) within-patient paired tests -----------------------------------
    tests = []
    for g in rep["genus"]:
        col = resolved[g]
        if col is None:
            continue
        for src in SOURCES:
            pairs = []
            for p, sub in meta.groupby("patient_id"):
                t = sub.index[sub["niche"] == "tumor"]
                s = sub.index[sub["niche"] == src]
                if len(t) == 0 or len(s) == 0:
                    continue
                pairs.append(float(rel_d.loc[t, col].mean() - rel_d.loc[s, col].mean()))
            if len(pairs) < 5 or np.all(np.array(pairs) == 0):
                continue
            w = stats.wilcoxon(pairs)
            tests.append({"genus": g, "contrast": f"tumor_vs_{src}", "n_pairs": len(pairs),
                          "median_delta_rel_abundance": float(np.median(pairs)),
                          "mean_delta": float(np.mean(pairs)),
                          "p": float(w.pvalue),
                          "direction": "up_in_tumour" if np.median(pairs) > 0 else "down_in_tumour"})
    tt = pd.DataFrame(tests)
    if not tt.empty:
        tt["q_bh"] = bh(tt["p"].values).round(4)
        tt = tt.sort_values(["contrast", "p"])
    tt.to_csv(TAB / "s7c_s5_paired_tests.csv", index=False)

    # --- (c) does the S5 direction reproduce vs paracancerous? -------------
    verdict = []
    for g in rep["genus"]:
        s5dir = rep.loc[rep["genus"] == g, "s5_direction"].iloc[0]
        r = tt[(tt["genus"] == g) & (tt["contrast"] == "tumor_vs_paracancerous")]
        if r.empty:
            verdict.append({"genus": g, "s5_direction": s5dir, "observed": "not_testable"})
            continue
        r = r.iloc[0]
        verdict.append({"genus": g, "s5_direction": s5dir,
                        "observed": r["direction"], "agrees": r["direction"] == s5dir,
                        "median_delta": round(float(r["median_delta_rel_abundance"]), 5),
                        "p": round(float(r["p"]), 4)})
    vt = pd.DataFrame(verdict)
    vt.to_csv(TAB / "s7c_s5_direction_check.csv", index=False)
    testable = vt[vt["observed"] != "not_testable"]
    n_agree = int(testable["agrees"].sum()) if "agrees" in testable else 0

    summary = {
        "n_s5_genera": int(len(rep)),
        "resolved_in_82": int(sum(1 for v in resolved.values() if v)),
        "direction_check_vs_paracancerous": {
            "testable": int(len(testable)), "agrees": n_agree,
            "detail": vt.to_dict("records")},
        "parvimonas_by_niche_pct": {
            k2: round(100 * float(niche_mean.loc[k2, resolved["Parvimonas"]]), 3)
            for k2 in ["tumor"] + SOURCES if resolved.get("Parvimonas")},
        "note": "abundances computed after dropping reagent-contaminant genera and "
                "renormalising; paired tests are within patient",
    }
    (TAB / "s7c_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False, default=str), encoding="utf-8")

    print("=== S5 genera mean rel abundance (%) by niche, decontaminated ===")
    print(byniche.to_string(index=False))
    print("\n=== direction check vs paracancerous ===")
    print(vt.to_string(index=False))
    print(f"\nagreement: {n_agree}/{len(testable)}")
    print("\n=== tumour vs each source niche (|median delta| >= 0.001, p<0.05) ===")
    if not tt.empty:
        sig = tt[(tt["q_bh"] < 0.05)].copy()
        print(sig[["genus", "contrast", "n_pairs", "median_delta_rel_abundance",
                   "direction", "q_bh"]].to_string(index=False))
    print(f"\nwrote {TAB}/s7c_*.csv|json")


if __name__ == "__main__":
    main()
