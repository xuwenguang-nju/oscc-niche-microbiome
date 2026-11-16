#!/usr/bin/env python3
"""S7: six-niche analysis of PRJNA866676 — which ecological niche does the tumour
microbiome actually resemble?

WHY THIS COHORT IS THE PAYLOAD OF THE WHOLE PROJECT
--------------------------------------------------
Every other OSCC cohort we analysed has only tumour and one adjacent tissue. PRJNA866676
samples up to SIX niches in the SAME patients (26 patients have all six):

    tumor                     deep tumour tissue
    paracancerous             adjacent tissue
    contralateral_normal_mucosa   the patient's own healthy mucosa, other side
    cancer_surface            the tumour surface, in direct contact with the oral cavity
    saliva                    whole oral cavity
    tongue_coat               heavy oral biofilm reservoir

That makes questions answerable that nothing else in the project can answer:
  Q1  Is the tumour microbiome distinguishable from the surrounding oral niches at all?
  Q2  Which niche is it closest to -- i.e. what is the likely SOURCE?
  Q3  Is there a BACKFLOW GRADIENT: cancer surface (exposed) more saliva-like than deep
      tumour (shielded)? If yes, the tumour community is plausibly acquired from the
      oral cavity rather than being tumour-intrinsic.

METHOD NOTES (carried over from S3)
  - Bray-Curtis on genus relative abundance. The S3 simulation established that a
    QUANTITATIVE mixing fraction is not identifiable from these data, so this script
    reports distances and RANKINGS only, never "X% of the tumour community is saliva".
  - Within-patient contrasts are paired, so all comparisons of interest use
    patient-blocked permutation, not pooled group tests.
  - Only genera present in >=20% of samples are used, so a niche is not "separated" by
    single-read noise taxa.

INPUTS
  results/dada2/PRJNA866676/genus_table_silva.csv
  results/tables/oscc_niche_sample_manifest.csv

OUTPUTS (results/tables/)
  s7_niche_metadata.csv  s7_niche_composition.csv  s7_discriminative_genera.csv
  s7_permanova.csv  s7_pairwise_niche_distances.csv
  s7_paired_distances.csv  s7_source_ranking.csv  s7_source_ranking_tests.csv
  s7_backflow_gradient.csv  s7_summary.json
Run: TCMA_ROOT=<dir> python3 s7_niche_analysis.py
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
MIN_PREVALENCE = 0.20
N_PERM = 2000
SEED = 42


def bray_curtis(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Row-wise Bray-Curtis between two composition matrices with equal columns."""
    num = np.abs(A[:, None, :] - B[None, :, :]).sum(axis=2)
    den = (A[:, None, :] + B[None, :, :]).sum(axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = num / den
    return np.where(den > 0, d, np.nan)


def permanova(D: np.ndarray, groups: np.ndarray, strata: np.ndarray | None,
              n_perm: int = N_PERM, seed: int = SEED) -> dict:
    """PERMANOVA on a distance matrix, with optional within-stratum permutation.

    Blocking by patient matters: the six niches come from the same people, so an
    unrestricted permutation tests the wrong null and inflates significance.
    """
    rng = np.random.default_rng(seed)
    n = D.shape[0]
    uniq = np.unique(groups)
    iu = np.triu_indices(n, 1)
    ss_total = (D[iu] ** 2).sum() / n

    def pseudo_f(g: np.ndarray) -> float:
        ss_within = 0.0
        for lev in uniq:
            idx = np.where(g == lev)[0]
            if len(idx) < 2:
                continue
            sub = D[np.ix_(idx, idx)]
            ss_within += (sub[np.triu_indices(len(idx), 1)] ** 2).sum() / len(idx)
        a = len([lev for lev in uniq if (g == lev).sum() > 0])
        ss_among = ss_total - ss_within
        if a < 2 or (n - a) <= 0 or ss_within <= 0:
            return np.nan
        return (ss_among / (a - 1)) / (ss_within / (n - a))

    F_obs = pseudo_f(groups)
    perms = []
    for _ in range(n_perm):
        if strata is None:
            perms.append(rng.permutation(groups))
        else:
            g = groups.copy()
            for s in np.unique(strata):
                idx = np.where(strata == s)[0]
                g[idx] = rng.permutation(g[idx])
            perms.append(g)
    null = np.array([pseudo_f(g) for g in perms])
    p = (np.sum(null >= F_obs) + 1) / (n_perm + 1)
    a = len(uniq)
    r2 = np.nan if np.isnan(F_obs) else (F_obs * (a - 1)) / (F_obs * (a - 1) + (n - a))
    return {"pseudo_F": float(F_obs), "p": float(p), "n_perm": n_perm,
            "n": int(n), "r2": float(r2)}


def main() -> None:
    gtab_path = ROOT / "results" / "dada2" / "PRJNA866676" / "genus_table_silva.csv"
    if not gtab_path.is_file():
        raise SystemExit(f"missing {gtab_path}")
    gtab = pd.read_csv(gtab_path, index_col=0)

    man = pd.read_csv(TAB / "oscc_niche_sample_manifest.csv", encoding="utf-8-sig")
    man = man[man["study"] == "PRJNA866676"][["run", "niche", "patient_id"]]
    meta = man.set_index("run")
    meta = meta.loc[meta.index.intersection(gtab.index)]
    rel = gtab.loc[meta.index]
    print(f"samples with both genus table and metadata: {len(meta)}")

    prev = (rel > 0).mean(axis=0)
    keep = prev[prev >= MIN_PREVALENCE].index
    rel = rel[keep]
    rel = rel.div(rel.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    print(f"genera kept (prevalence >= {MIN_PREVALENCE}): {len(keep)} of {gtab.shape[1]}")

    meta.to_csv(TAB / "s7_niche_metadata.csv")
    summary: dict = {
        "cohort": "PRJNA866676", "platform": "Ion Torrent, single-end",
        "n_samples": int(len(meta)), "n_patients": int(meta["patient_id"].nunique()),
        "n_genera_used": int(len(keep)), "min_prevalence": MIN_PREVALENCE,
        "niche_counts": meta["niche"].value_counts().to_dict(),
        "scope": "distances and rankings only; quantitative mixing fractions are NOT "
                 "identifiable from these data (see S3)",
    }

    # Group SAMPLES by niche. Doing this on rel.T with a bare Series does NOT work:
    # groupby aligns the Series on the frame's index (genera), so every group comes out
    # empty and you get a silently empty table rather than an error. Attach the niche as
    # a column of the sample-indexed frame instead.
    tmp = rel.copy()
    tmp["__niche"] = meta["niche"].reindex(tmp.index)
    if tmp["__niche"].isna().any():
        raise SystemExit("niche labels missing for some samples")
    comp = tmp.groupby("__niche").mean()                  # niches x genera
    comp = comp.loc[[n for n in NICHE_ORDER if n in comp.index]]
    if len(comp.index) < 2:
        raise SystemExit(f"only {len(comp.index)} niches present: {list(comp.index)}")
    comp.index.name = "niche"
    comp.to_csv(TAB / "s7_niche_composition.csv")

    disc = []
    for n in comp.index:
        others = [c for c in comp.index if c != n]
        d = comp.loc[n] - comp.loc[others].mean(axis=0)
        for g, v in pd.concat([d.nlargest(8), d.nsmallest(4)]).items():
            disc.append({"niche": n, "genus": g, "delta_rel_abundance": float(v),
                         "mean_in_niche": float(comp.loc[n, g]),
                         "mean_elsewhere": float(comp.loc[others, g].mean())})
    pd.DataFrame(disc).to_csv(TAB / "s7_discriminative_genera.csv", index=False)

    # ---- Q1 -----------------------------------------------------------------
    D = bray_curtis(rel.values, rel.values)
    groups = meta["niche"].values
    strata = meta["patient_id"].values
    perm = permanova(D, groups, strata)
    perm_free = permanova(D, groups, None, n_perm=500)
    pd.DataFrame([{"test": "niche (patient-blocked)", **perm},
                  {"test": "niche (unblocked)", **perm_free}]
                 ).to_csv(TAB / "s7_permanova.csv", index=False)
    summary["permanova"] = {"blocked": perm, "unblocked": perm_free}

    pair_rows = []
    for a, b in combinations([n for n in NICHE_ORDER if n in comp.index], 2):
        ia = np.where(groups == a)[0]
        ib = np.where(groups == b)[0]
        sub = D[np.ix_(ia, ib)]
        pair_rows.append({"niche_a": a, "niche_b": b,
                          "mean_bray_curtis": float(np.nanmean(sub)),
                          "median": float(np.nanmedian(sub))})
    pairs = pd.DataFrame(pair_rows).sort_values("mean_bray_curtis")
    pairs.to_csv(TAB / "s7_pairwise_niche_distances.csv", index=False)
    summary["closest_niche_pairs"] = pairs.head(4).to_dict("records")
    summary["farthest_niche_pairs"] = pairs.tail(4).to_dict("records")

    # ---- Q2 -----------------------------------------------------------------
    pats = meta.groupby("patient_id")["niche"].apply(set)
    with_tumour = [p for p, s in pats.items() if "tumor" in s]
    rows = []
    for p in with_tumour:
        t_idx = meta.index[(meta["patient_id"] == p) & (meta["niche"] == "tumor")]
        for t in t_idx:
            ti = meta.index.get_loc(t)
            for src in SOURCES:
                s_idx = meta.index[(meta["patient_id"] == p) & (meta["niche"] == src)]
                for s in s_idx:
                    rows.append({"patient_id": p, "tumor_run": t, "source_niche": src,
                                 "source_run": s,
                                 "bray_curtis": float(D[ti, meta.index.get_loc(s)])})
    src = pd.DataFrame(rows)
    if src.empty:
        raise SystemExit("no within-patient tumour-source pairs")
    src.to_csv(TAB / "s7_paired_distances.csv", index=False)

    per_pat = (src.groupby(["patient_id", "source_niche"])["bray_curtis"]
                  .mean().reset_index())
    rank_rows = []
    for p, sub in per_pat.groupby("patient_id"):
        sub = sub.sort_values("bray_curtis")
        for r, n in enumerate(sub["source_niche"], start=1):
            rank_rows.append({"patient_id": p, "source_niche": n, "rank": r,
                              "bray_curtis": float(
                                  sub.loc[sub["source_niche"] == n, "bray_curtis"].iloc[0])})
    ranks = pd.DataFrame(rank_rows)
    ranks.to_csv(TAB / "s7_source_ranking.csv", index=False)

    rank_summary = (ranks.groupby("source_niche")["rank"]
                         .agg(["size", "mean", "median"]).reset_index()
                         .rename(columns={"size": "n_patients"}).sort_values("mean"))

    test_rows = []
    for n in rank_summary["source_niche"]:
        d_this = per_pat.loc[per_pat["source_niche"] == n].set_index("patient_id")["bray_curtis"]
        d_other = (per_pat.loc[per_pat["source_niche"] != n]
                        .groupby("patient_id")["bray_curtis"].mean())
        common = d_this.index.intersection(d_other.index)
        if len(common) < 5:
            continue
        w = stats.wilcoxon(d_this[common], d_other[common], alternative="two-sided")
        test_rows.append({"source_niche": n, "n_patients": int(len(common)),
                          "mean_dist_to_tumour": float(d_this[common].mean()),
                          "mean_dist_to_others": float(d_other[common].mean()),
                          "delta": float(d_this[common].mean() - d_other[common].mean()),
                          "p_vs_others": float(w.pvalue)})
    tests = pd.DataFrame(test_rows)
    tests.to_csv(TAB / "s7_source_ranking_tests.csv", index=False)
    summary["source_attribution"] = {
        "n_patients_with_tumour": int(len(with_tumour)),
        "ranking_by_mean_rank": rank_summary.round(3).to_dict("records"),
        "tests": tests.round(4).to_dict("records")}

    # ---- Q3 backflow --------------------------------------------------------
    grad_rows = []
    for p in with_tumour:
        idx = meta.index[meta["patient_id"] == p]
        have = meta.loc[idx, "niche"]
        if not {"tumor", "cancer_surface", "saliva"}.issubset(set(have)):
            continue

        def mean_d(n1, n2):
            a = idx[have == n1]; b = idx[have == n2]
            vals = [D[meta.index.get_loc(x), meta.index.get_loc(y)]
                    for x in a for y in b]
            return float(np.mean(vals)) if vals else np.nan

        grad_rows.append({
            "patient_id": p,
            "d_tumour_saliva": mean_d("tumor", "saliva"),
            "d_surface_saliva": mean_d("cancer_surface", "saliva"),
            "d_tumour_surface": mean_d("tumor", "cancer_surface"),
            "d_tumour_paracancerous": mean_d("tumor", "paracancerous"),
            "d_surface_paracancerous": mean_d("cancer_surface", "paracancerous"),
            "d_tumour_contralateral": mean_d("tumor", "contralateral_normal_mucosa"),
            "d_surface_contralateral": mean_d("cancer_surface",
                                              "contralateral_normal_mucosa")})
    grad = pd.DataFrame(grad_rows).dropna(subset=["d_tumour_saliva", "d_surface_saliva"])
    if len(grad) >= 5:
        grad["saliva_closer_to_surface"] = (grad["d_surface_saliva"]
                                            < grad["d_tumour_saliva"])
        w = stats.wilcoxon(grad["d_tumour_saliva"], grad["d_surface_saliva"],
                           alternative="two-sided")
        grad.to_csv(TAB / "s7_backflow_gradient.csv", index=False)
        summary["backflow"] = {
            "n_patients": int(len(grad)),
            "median_d_tumour_saliva": round(float(grad["d_tumour_saliva"].median()), 4),
            "median_d_surface_saliva": round(float(grad["d_surface_saliva"].median()), 4),
            "frac_patients_surface_closer_to_saliva":
                round(float(grad["saliva_closer_to_surface"].mean()), 3),
            "wilcoxon_p": float(w.pvalue),
            "interpretation": "surface closer to saliva than tumour = consistent with an "
                              "oral-cavity source whose penetration decreases with depth"}
    else:
        summary["backflow"] = {"n_patients": int(len(grad)),
                               "note": "too few patients for a paired test"}

    (TAB / "s7_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False, default=str), encoding="utf-8")

    print("\n=== niche counts ===")
    for k, v in summary["niche_counts"].items():
        print(f"  {k:<30} {v}")
    print("\n=== Q1 PERMANOVA (patient-blocked) ===")
    print(f"  pseudo-F={perm['pseudo_F']:.3f}  p={perm['p']:.4f}  R2={perm['r2']:.3f}"
          f"   (unblocked p={perm_free['p']:.4f})")
    print("\n=== closest niche pairs ===")
    print(pairs.head(6).round(4).to_string(index=False))
    print("\n=== Q2 source ranking (1 = closest to tumour) ===")
    print(rank_summary.round(3).to_string(index=False))
    print("\n=== Q2 per-niche tests ===")
    print(tests.round(4).to_string(index=False))
    b = summary.get("backflow", {})
    if "wilcoxon_p" in b:
        print(f"\n=== Q3 backflow ===\n  n={b['n_patients']}  "
              f"d(tumour,saliva)={b['median_d_tumour_saliva']}  "
              f"d(surface,saliva)={b['median_d_surface_saliva']}  "
              f"frac surface closer={b['frac_patients_surface_closer_to_saliva']}  "
              f"p={b['wilcoxon_p']:.4f}")
    print(f"\nwrote {TAB}/s7_*.csv|json")


if __name__ == "__main__":
    main()
