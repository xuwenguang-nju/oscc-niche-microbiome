#!/usr/bin/env python3
"""S7b: decontamination sensitivity analysis for the six-niche conclusions.

WHY THIS IS NECESSARY
---------------------
The first pass over PRJNA866676 produced a clean-looking story (tumour closest to
paracancerous, far from saliva; strong depth gradient). But the per-niche contaminant
burden is wildly uneven:

    paracancerous   47.5%  (Stenotrophomonas alone 42.6%)
    tumor           13.2%  (Stenotrophomonas 8.9%)
    contralateral    7.9%  (Stenotrophomonas 6.8%)
    cancer_surface   8.8%
    saliva           0.5%
    tongue_coat      0.8%

Saliva and tongue coat are clean; the tissue niches are not. Two tissue niches sharing
a reagent contaminant will look "similar" for reasons that have nothing to do with
biology, so the tumour-to-paracancerous proximity and the tumour-to-contralateral
distance could both be artefacts of contamination.

No negative or blank controls are deposited for this study, so decontam-style
frequency-based identification is not possible. The defensible alternative is to drop
genera that are established reagent/environment contaminant taxa, renormalise, and ask
whether each conclusion survives. If a conclusion flips, it was contamination.

INPUTS  results/tables/s7_* (from s7_niche_analysis.py) + the genus table
OUTPUTS results/tables/s7b_*  and a comparison table
Run: TCMA_ROOT=<dir> python3 s7b_decontam_sensitivity.py
"""
from __future__ import annotations

import json
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

NICHE_ORDER = ["tumor", "paracancerous", "cancer_surface",
               "contralateral_normal_mucosa", "saliva", "tongue_coat"]
SOURCES = ["paracancerous", "cancer_surface", "contralateral_normal_mucosa",
           "saliva", "tongue_coat"]
# Established reagent / environment background genera in low-biomass 16S work.
# Deliberately conservative: only taxa repeatedly reported as kit or water contaminants,
# so that removing them cannot manufacture the result we want.
CONTAMINANTS = {
    "Stenotrophomonas", "Pseudomonas", "Bradyrhizobium", "Acinetobacter",
    "Sphingomonas", "Microbacterium", "Corynebacterium", "Flavobacterium",
    "Delftia", "Comamonas", "Methylobacterium", "Ralstonia", "Burkholderia",
    "Phyllobacterium", "Bosea", "Sphingobacterium", "Chryseobacterium",
    "Variovorax", "Cupriavidus", "Aquabacterium", "Methylophilus", "Massilia",
    "Hymenobacter", "Sphingopyxis", "Herbaspirillum", "Brevundimonas",
    "Novosphingobium", "Devosia", "Rhizobium", "Mesorhizobium", "Aeromonas",
}
MIN_PREVALENCE = 0.20
N_PERM = 2000


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


def bray_curtis(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    num = np.abs(A[:, None, :] - B[None, :, :]).sum(axis=2)
    den = (A[:, None, :] + B[None, :, :]).sum(axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = num / den
    return np.where(den > 0, d, np.nan)


def permanova(D, groups, strata, n_perm=N_PERM, seed=42):
    rng = np.random.default_rng(seed)
    n = D.shape[0]
    uniq = np.unique(groups)
    ss_total = (D[np.triu_indices(n, 1)] ** 2).sum() / n

    def pseudo_f(g):
        ss_within = 0.0
        for lev in uniq:
            idx = np.where(g == lev)[0]
            if len(idx) < 2:
                continue
            sub = D[np.ix_(idx, idx)]
            ss_within += (sub[np.triu_indices(len(idx), 1)] ** 2).sum() / len(idx)
        a = len([l for l in uniq if (g == l).sum() > 0])
        if a < 2 or (n - a) <= 0 or ss_within <= 0:
            return np.nan
        return ((ss_total - ss_within) / (a - 1)) / (ss_within / (n - a))

    F = pseudo_f(groups)
    null = []
    for _ in range(n_perm):
        g = groups.copy()
        for s in np.unique(strata):
            idx = np.where(strata == s)[0]
            g[idx] = rng.permutation(g[idx])
        null.append(pseudo_f(g))
    null = np.array(null)
    return {"pseudo_F": float(F), "p": float((np.sum(null >= F) + 1) / (n_perm + 1))}


def run_analysis(rel: pd.DataFrame, meta: pd.DataFrame, tag: str) -> dict:
    D = bray_curtis(rel.values, rel.values)
    groups = meta["niche"].values
    strata = meta["patient_id"].values
    out: dict = {"tag": tag, "n_genera": int(rel.shape[1])}
    out["permanova"] = permanova(D, groups, strata)

    # source ranking within patient
    rows = []
    for p, sub in meta.groupby("patient_id"):
        if "tumor" not in set(sub["niche"]):
            continue
        t_idx = sub.index[sub["niche"] == "tumor"]
        for t in t_idx:
            ti = meta.index.get_loc(t)
            for src in SOURCES:
                for s in sub.index[sub["niche"] == src]:
                    rows.append({"patient_id": p, "source_niche": src,
                                 "bray_curtis": float(D[ti, meta.index.get_loc(s)])})
    src = pd.DataFrame(rows)
    per_pat = src.groupby(["patient_id", "source_niche"])["bray_curtis"].mean().reset_index()
    rankings = {}
    for n, g in per_pat.groupby("source_niche"):
        ranks = []
        for p, sub in per_pat.groupby("patient_id"):
            if n not in set(sub["source_niche"]):
                continue
            sub = sub.sort_values("bray_curtis")
            ranks.append(int(np.where(sub["source_niche"].values == n)[0][0]) + 1)
        rankings[n] = {"mean_rank": float(np.mean(ranks)), "n": len(ranks)}
    out["ranking"] = rankings

    # per-niche paired test vs all other sources
    tests = []
    for n in per_pat["source_niche"].unique():
        d_this = per_pat[per_pat["source_niche"] == n].set_index("patient_id")["bray_curtis"]
        d_other = (per_pat[per_pat["source_niche"] != n]
                   .groupby("patient_id")["bray_curtis"].mean())
        common = d_this.index.intersection(d_other.index)
        if len(common) < 5:
            continue
        w = stats.wilcoxon(d_this[common], d_other[common])
        tests.append({"source_niche": n, "n_patients": int(len(common)),
                      "delta": float(d_this[common].mean() - d_other[common].mean()),
                      "p": float(w.pvalue)})
    out["tests"] = tests

    # backflow
    grad = []
    for p, sub in meta.groupby("patient_id"):
        have = set(sub["niche"])
        if not {"tumor", "cancer_surface", "saliva"}.issubset(have):
            continue

        def md(n1, n2):
            a = sub.index[sub["niche"] == n1]; b = sub.index[sub["niche"] == n2]
            v = [D[meta.index.get_loc(x), meta.index.get_loc(y)] for x in a for y in b]
            return float(np.mean(v)) if v else np.nan

        grad.append({"patient_id": p, "dt_s": md("tumor", "saliva"),
                     "ds_s": md("cancer_surface", "saliva")})
    gdf = pd.DataFrame(grad).dropna()
    if len(gdf) >= 5:
        w = stats.wilcoxon(gdf["dt_s"], gdf["ds_s"])
        out["backflow"] = {
            "n": int(len(gdf)),
            "median_d_tumour_saliva": round(float(gdf["dt_s"].median()), 4),
            "median_d_surface_saliva": round(float(gdf["ds_s"].median()), 4),
            "frac_surface_closer": round(float((gdf["ds_s"] < gdf["dt_s"]).mean()), 3),
            "p": float(w.pvalue)}
    return out


def main() -> None:
    gtab = pd.read_csv(ROOT / "results" / "dada2" / "PRJNA866676"
                       / "genus_table_silva.csv", index_col=0)
    man = pd.read_csv(TAB / "oscc_niche_sample_manifest.csv", encoding="utf-8-sig")
    man = man[man["study"] == "PRJNA866676"][["run", "niche", "patient_id"]].set_index("run")
    meta = man.loc[man.index.intersection(gtab.index)]

    prev = (gtab > 0).mean(axis=0)
    keep = prev[prev >= MIN_PREVALENCE].index
    rel_all = gtab.loc[meta.index, keep]
    rel_all = rel_all.div(rel_all.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    dropped = [c for c in rel_all.columns if c in CONTAMINANTS]
    rel_decon = rel_all.drop(columns=dropped)
    rel_decon = rel_decon.div(rel_decon.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    print(f"genera used: {rel_all.shape[1]}; dropped as contaminants: {len(dropped)}")
    print("dropped:", ", ".join(sorted(dropped)))

    before = run_analysis(rel_all, meta, "before_decontam")
    after = run_analysis(rel_decon, meta, "after_decontam")

    # contaminant burden per niche, for the record
    burden = []
    contam_cols = [c for c in rel_all.columns if c in CONTAMINANTS]
    tmp = rel_all.copy()
    tmp["__niche"] = meta["niche"].reindex(tmp.index)
    for n, sub in tmp.groupby("__niche"):
        burden.append({"niche": n,
                       "contaminant_pct": round(100 * sub[contam_cols].sum(axis=1).mean(), 1) if contam_cols else 0.0,
                       "n": int(len(sub))})
    pd.DataFrame(burden).to_csv(TAB / "s7b_contaminant_burden.csv", index=False)

    cmp_rows = []
    for tag, r in [("before", before), ("after", after)]:
        row = {"analysis": tag, "n_genera": r["n_genera"],
               "permanova_F": round(r["permanova"]["pseudo_F"], 3),
               "permanova_p": round(r["permanova"]["p"], 4)}
        for n, v in r["ranking"].items():
            row[f"rank_{n}"] = round(v["mean_rank"], 2)
        for t in r["tests"]:
            row[f"delta_{t['source_niche']}"] = round(t["delta"], 4)
            row[f"p_{t['source_niche']}"] = round(t["p"], 4)
        b = r.get("backflow", {})
        row["backflow_median_d_tumour_saliva"] = b.get("median_d_tumour_saliva")
        row["backflow_median_d_surface_saliva"] = b.get("median_d_surface_saliva")
        row["backflow_frac_surface_closer"] = b.get("frac_surface_closer")
        row["backflow_p"] = b.get("p")
        cmp_rows.append(row)
    comp = pd.DataFrame(cmp_rows)
    comp.to_csv(TAB / "s7b_before_after.csv", index=False)

    # do the conclusions survive?
    verdict = {}
    b_tests = {t["source_niche"]: t for t in before["tests"]}
    a_tests = {t["source_niche"]: t for t in after["tests"]}
    for n in SOURCES:
        if n in b_tests and n in a_tests:
            verdict[n] = {
                "before_delta": round(b_tests[n]["delta"], 4),
                "before_p": round(b_tests[n]["p"], 4),
                "after_delta": round(a_tests[n]["delta"], 4),
                "after_p": round(a_tests[n]["p"], 4),
                "sign_flips": bool(np.sign(b_tests[n]["delta"]) != np.sign(a_tests[n]["delta"])),
                "significant_after": bool(a_tests[n]["p"] < 0.05)}
    summary = {
        "contaminant_genera_dropped": sorted(dropped),
        "n_genera_before": int(rel_all.shape[1]), "n_genera_after": int(rel_decon.shape[1]),
        "contaminant_burden_pct": {b["niche"]: b["contaminant_pct"] for b in burden},
        "before": before, "after": after, "per_source_verdict": verdict,
        "backflow_survives": bool(
            before.get("backflow", {}).get("p", 1) < 0.05
            and after.get("backflow", {}).get("p", 1) < 0.05),
    }
    (TAB / "s7b_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False, default=str), encoding="utf-8")

    print("\n=== contaminant burden by niche (%) ===")
    print(pd.DataFrame(burden).to_string(index=False))
    print("\n=== before vs after ===")
    print(comp.T.to_string())
    print("\n=== per-source verdict ===")
    for n, v in verdict.items():
        flip = "FLIPS" if v["sign_flips"] else "stable"
        print(f"  {n:<30} before delta={v['before_delta']:+.4f} p={v['before_p']:.4f} | "
              f"after delta={v['after_delta']:+.4f} p={v['after_p']:.4f}  {flip}")
    print(f"\nwrote {TAB}/s7b_*.csv|json")


if __name__ == "__main__":
    main()
