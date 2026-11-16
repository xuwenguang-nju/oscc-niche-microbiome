#!/usr/bin/env python3
"""S8a: replicate the S7 niche-attribution result in two independent cohorts.

S7 (PRJNA866676) found that the tumour microbiome is closest to the paracancerous tissue,
and that the surface-exposed sample is closer to saliva than the deep tumour is. Two
questions follow, and only one is fully answerable with existing data:

  (1) Does "the tissue next to the tumour is the closest source" replicate?  -> testable
      in PRJNA412445 (tumour + lesion-side mucosa + contralateral mucosa + saliva + plaque)
      and in PRJNA1415359 (tumour + adjacent normal + saliva).

  (2) Does the surface->saliva gradient replicate?  -> NOT directly testable. A direct
      replication needs a cohort that sampled BOTH a surface and a deep tumour site in the
      same patient; PRJNA866676 is the only such cohort here. What we can test is a weaker
      proxy: is the surface-exposed mucosa (lesion side) more saliva-like than the
      contralateral mucosa? Reported separately and labelled as a proxy, not as a
      replication.

Statistics: per patient, rank each candidate source by distance to that patient's tumour
(rank 1 = closest). Under the null of no preference the mean rank is (k+1)/2, so the mean
rank is tested against that with a Wilcoxon signed-rank test, alongside the per-source
distances. Bray-Curtis and Aitchison are both reported because they weight rare taxa
differently.

Run: python3 s8a_replicate_ranking.py --genus-dir data/genus --out results/tables
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import braycurtis

COHORTS = {
    "PRJNA412445": {
        "tumour": "tumor",
        "sources": ["lesion_side_mucosa", "contralateral_normal_mucosa",
                    "saliva", "subgingival_plaque"],
        "surface_proxy": ("lesion_side_mucosa", "contralateral_normal_mucosa"),
    },
    "PRJNA1415359": {
        "tumour": "tumor",
        "sources": ["adjacent_normal_tissue", "preop_saliva"],
        "surface_proxy": None,
    },
}


def clr(mat: np.ndarray, pseudo: float = 1e-6) -> np.ndarray:
    lg = np.log(mat + pseudo)
    return lg - lg.mean(axis=1, keepdims=True)


def distances_to_ref(genus: pd.DataFrame, meta: pd.DataFrame, study: str,
                     tumour_niche: str) -> pd.DataFrame:
    sub = genus.loc[genus.index.intersection(meta.index)]
    rel = sub.div(sub.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    rel = rel.loc[:, (rel > 0).mean(axis=0) >= 0.10]
    cl = pd.DataFrame(clr(rel.values), index=rel.index, columns=rel.columns)

    m = meta.loc[rel.index]
    rows = []
    for pid, grp in m.groupby("patient_id"):
        if tumour_niche not in set(grp["niche"]):
            continue
        for tr in grp.index[grp["niche"] == tumour_niche]:
            for src in sorted(set(grp["niche"])):
                if src == tumour_niche:
                    continue
                for sr in grp.index[grp["niche"] == src]:
                    rows.append({
                        "study": study, "patient_id": pid, "source": src,
                        "bray": float(braycurtis(rel.loc[tr].values, rel.loc[sr].values)),
                        "aitchison": float(np.linalg.norm(
                            cl.loc[tr].values - cl.loc[sr].values)),
                    })
    return pd.DataFrame(rows)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def rank_test(long: pd.DataFrame, metric: str, sources: list[str]) -> pd.DataFrame:
    d = long[long.source.isin(sources)].groupby(
        ["patient_id", "source"])[metric].mean().unstack()
    d = d[[s for s in sources if s in d.columns]].dropna()
    if d.empty:
        return pd.DataFrame()
    ranks = d.rank(axis=1, method="average")
    k = d.shape[1]
    null = (k + 1) / 2
    rows = []
    for s in d.columns:
        r = ranks[s]
        try:
            p_rank = stats.wilcoxon(r - null).pvalue
        except ValueError:
            p_rank = np.nan
        rows.append({
            "study": long["study"].iloc[0], "metric": metric, "source": s,
            "n_patients": int(len(r)),
            "mean_rank": round(float(r.mean()), 3),
            "frac_rank1": round(float((r == 1).mean()), 3),
            "p_rank_vs_null": p_rank,
            "median_distance": round(float(d[s].median()), 4),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus-dir", default="data/genus")
    ap.add_argument("--manifest", default="results/tables/oscc_niche_sample_manifest.csv")
    ap.add_argument("--out", default="results/tables")
    args = ap.parse_args()

    genus_dir, out = Path(args.genus_dir), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(args.manifest, encoding="utf-8-sig")

    all_long, all_rank, all_proxy = [], [], []
    for study, cfg in COHORTS.items():
        f = genus_dir / study / "genus_table.csv"
        if not f.is_file():
            print(f"[skip] {study}: no genus table at {f}")
            continue
        genus = pd.read_csv(f, index_col=0)
        # Run ids appear as SRRxxxx_1_F_filt / SRRxxxx_F_filt / SRRxxxx; keep only the SRR id.
        norm = [re.match(r"(SRR\d+|ERR\d+|DRR\d+)", str(i)) for i in genus.index]
        genus.index = [m.group(1) if m else str(i) for i, m in zip(genus.index, norm)]
        meta = man[man.study == study].drop_duplicates("run").set_index("run")
        common = genus.index.intersection(meta.index)
        print(f"{study}: {len(common)} samples with genus + metadata "
              f"(tumour={int((meta.loc[common, 'niche'] == cfg['tumour']).sum())})")
        if not len(common):
            continue

        long = distances_to_ref(genus, meta, study, cfg["tumour"])
        if long.empty:
            print(f"[skip] {study}: no paired tumour/source structure")
            continue
        all_long.append(long)

        for metric in ("bray", "aitchison"):
            r = rank_test(long, metric, cfg["sources"])
            if not r.empty:
                all_rank.append(r)

        sp = cfg["surface_proxy"]
        if sp:
            d = long[long.source.isin(sp)].groupby(
                ["patient_id", "source"])["aitchison"].mean().unstack()
            d = d[[s for s in sp if s in d.columns]].dropna()
            if len(d) >= 5 and d.shape[1] == 2:
                diff = d[sp[0]] - d[sp[1]]
                try:
                    p = stats.wilcoxon(diff).pvalue
                except ValueError:
                    p = np.nan
                all_proxy.append({
                    "study": study,
                    "test": f"d(tumour->{sp[0]}) - d(tumour->{sp[1]})",
                    "metric": "aitchison", "n_pairs": int(len(d)),
                    "median_diff": round(float(diff.median()), 4),
                    "frac_negative": round(float((diff < 0).mean()), 3),
                    "p": p,
                })

    if all_long:
        L = pd.concat(all_long, ignore_index=True)
        L.to_csv(out / "s8a_distances_long.csv", index=False)
        summ = (L.groupby(["study", "source"])
                 .agg(n=("patient_id", "nunique"),
                      bray_median=("bray", "median"),
                      aitchison_median=("aitchison", "median"))
                 .round(4).reset_index())
        summ.to_csv(out / "s8a_source_distance_summary.csv", index=False)
        print("\n=== median Aitchison distance to tumour, by cohort ===")
        print(summ.pivot(index="source", columns="study",
                         values="aitchison_median").to_string())

    if all_rank:
        R = pd.concat(all_rank, ignore_index=True)
        R["q_within_study_metric"] = np.nan
        for (st, me), idx in R.groupby(["study", "metric"]).groups.items():
            R.loc[idx, "q_within_study_metric"] = bh(
                R.loc[idx, "p_rank_vs_null"].values)
        R = R.round({"p_rank_vs_null": 6, "q_within_study_metric": 4,
                     "mean_rank": 3, "median_distance": 4})
        R.to_csv(out / "s8a_rank_tests.csv", index=False)
        print("\n=== rank tests (rank 1 = closest to tumour) ===")
        print(R[["study", "metric", "source", "n_patients", "mean_rank",
                 "frac_rank1", "p_rank_vs_null", "q_within_study_metric"]]
              .to_string(index=False))

    if all_proxy:
        P = pd.DataFrame(all_proxy)
        P.to_csv(out / "s8a_surface_proxy.csv", index=False)
        print("\n=== surface-exposure proxy (not a replication) ===")
        print(P.to_string(index=False))


if __name__ == "__main__":
    main()
