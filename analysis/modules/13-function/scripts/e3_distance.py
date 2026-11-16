#!/usr/bin/env python3
"""E3: distance-level comparison of the tumour functional profile to each niche.

Why a second, distance-based pass
---------------------------------
E2's ranking counted how often each niche was "nearest", but that puts `cancer_surface` at the
top (20/37). That is not a contradiction of S7 -- cancer_surface IS the tumour surface, i.e. the
same lesion, so of course it resembles the tumour. A count-based ranking conflates "same
structure" with "external source".

The informative question is comparative: for the same patient, is the tumour functionally
CLOSER to adjacent normal mucosa than to saliva / tongue coat / contralateral mucosa?
That is a paired test on distances and it mirrors S7's genus-level claim directly, so the two
analyses can be compared term by term.

Method
------
For each patient i and niche n, the reference profile is the mean pathway profile of niche n
across all OTHER patients (leave-self-out), so patient i cannot match themselves.
Distance = Bray-Curtis between the patient's tumour profile and that reference.
Pairs compared with Wilcoxon signed-rank; effect size reported as median difference and as
the fraction of patients in which the comparison goes the hypothesised way.

Run: python3 e3_distance.py
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import braycurtis

OUT = Path("/home/v3724/oscc-mb/work/picrust/PRJNA866676/out")
RES = Path("/home/v3724/results/tables")
MANIFEST = Path("/home/v3724/oscc-mb/work/s8/results/tables/"
                "oscc_niche_sample_manifest.csv")
STUDY = "PRJNA866676"

# niches that are anatomically external to the tumour mass
EXTERNAL = ["paracancerous", "saliva", "tongue_coat", "contralateral_normal_mucosa"]


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for r, i in enumerate(o[::-1]):
        k = n - r
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return np.clip(q, 0, 1)


def main() -> None:
    P = pd.read_csv(OUT / "pathways_out/path_abun_unstrat.tsv.gz", sep="\t", index_col=0)
    rel = P.div(P.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)

    man = pd.read_csv(MANIFEST, encoding="utf-8-sig")
    man = man[man.study == STUDY].drop_duplicates("run").set_index("run")
    common = [s for s in rel.columns if s in man.index]
    rel, meta = rel[common], man.loc[common]

    mat = rel.T.values
    niche = meta.niche.values
    pid = meta.patient_id.values

    # per-patient per-niche mean profile
    prof = {}
    for p_ in pd.unique(pid):
        for n_ in np.unique(niche):
            sel = (pid == p_) & (niche == n_)
            if sel.sum():
                prof[(p_, n_)] = mat[sel].mean(axis=0)

    patients = [p_ for p_ in pd.unique(pid) if (p_, "tumor") in prof]
    print(f"patients with tumour: {len(patients)}")

    # distances with leave-self-out references
    recs = []
    for p_ in patients:
        row = {"patient_id": p_}
        for n_ in np.unique(niche):
            others = [prof[(q, n_)] for q in patients if q != p_ and (q, n_) in prof]
            if len(others) >= 3:
                row[n_] = float(braycurtis(prof[(p_, "tumor")],
                                           np.mean(np.vstack(others), axis=0)))
        recs.append(row)
    D = pd.DataFrame(recs).set_index("patient_id")
    D.to_csv(RES / "e3_tumour_distances.csv")
    print("\n=== mean Bray-Curtis distance from tumour, by niche (leave-self-out) ===")
    cols = [c for c in D.columns if c != "tumor"]
    for c in sorted(cols, key=lambda x: D[x].mean()):
        print(f"  {c:32} n={D[c].notna().sum():3}  mean={D[c].mean():.4f}")

    # paired comparisons among external niches + cancer_surface
    order = ["cancer_surface"] + EXTERNAL
    order = [o for o in order if o in D.columns]
    rows = []
    for a, b in combinations(order, 2):
        sub = D[[a, b]].dropna()
        if len(sub) < 8:
            continue
        d = sub[a] - sub[b]
        try:
            _, pv = stats.wilcoxon(sub[a], sub[b])
        except ValueError:
            pv = 1.0
        rows.append({"closer_to": a if d.median() < 0 else b,
                     "pair": f"{a} vs {b}", "n": len(sub),
                     "median_d_a_minus_b": round(float(d.median()), 4),
                     "frac_a_closer": round(float((d < 0).mean()), 3),
                     "p": float(pv)})
    T = pd.DataFrame(rows)
    if len(T):
        T["q"] = bh(T.p.values)
        T = T.sort_values("p")
        T.to_csv(RES / "e3_distance_paired_tests.csv", index=False)
        print("\n=== paired distance comparisons (negative median => first is closer) ===")
        print(T.round(4).to_string(index=False))

    # the S7-parallel claim, stated explicitly
    if "paracancerous" in D.columns and "saliva" in D.columns:
        sub = D[["paracancerous", "saliva"]].dropna()
        d = sub["paracancerous"] - sub["saliva"]
        _, pv = stats.wilcoxon(sub["paracancerous"], sub["saliva"])
        print(f"\n=== S7-parallel test ===")
        print(f"  tumour closer to adjacent normal than to saliva in "
              f"{(d < 0).sum()}/{len(sub)} patients ({100*(d<0).mean():.1f}%)")
        print(f"  median difference {d.median():+.4f}, Wilcoxon p={pv:.4g}")

    # external-source-only ranking (cancer_surface excluded: it is the same lesion)
    ext = [c for c in EXTERNAL if c in D.columns]
    sub = D[ext].dropna()
    if len(sub):
        nearest = sub.idxmin(axis=1)
        print(f"\n=== nearest EXTERNAL source (cancer_surface excluded), n={len(sub)} ===")
        print("  " + nearest.value_counts().to_string().replace("\n", "\n  "))
        k = int((nearest == "paracancerous").sum())
        pb = stats.binomtest(k, len(sub), 1 / len(ext)).pvalue
        print(f"  adjacent normal nearest in {k}/{len(sub)} (chance 1/{len(ext)}), "
              f"binom p={pb:.4g}")

    (RES / "e3_summary.json").write_text(json.dumps({
        "n_patients_tumour": len(patients),
        "mean_distance_by_niche": {c: round(float(D[c].mean()), 4) for c in cols},
        "external_niches": ext,
    }, indent=1), encoding="utf-8")
    print("\ndone")


if __name__ == "__main__":
    main()
