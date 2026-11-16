#!/usr/bin/env python3
"""A2: clinical associations of the intratumoral microbiome (TCGA-HNSC / TCMA).

Question set, mirroring what 13/13 and 12/13 of the retrieved papers do:
  does microbial diversity or Parvimonas abundance vary with sex, anatomical site, stage,
  grade, smoking, or HPV status?

Statistics are chosen per variable type rather than applied uniformly:
  two groups            -> Mann-Whitney U (Wilcoxon rank-sum)
  three or more groups  -> Kruskal-Wallis
  ordinal (stage 1-4, grade 1-3) -> also Spearman rho as a trend test, which is more
                          appropriate than an unordered omnibus test
Everything is BH-corrected within each microbiome metric.

Caveats recorded in the output: stage is 75% complete, HPV only 16%, so HPV is reported as
exploratory and not used to gate any conclusion.

Run: python3 s11_clinical_assoc.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

METRICS = ["shannon", "richness", "parvimonas_rel"]


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    r = pv[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    q[ok] = out
    return q


def describe(d: pd.DataFrame, metric: str, var: str) -> str:
    g = d.groupby(var)[metric]
    return "; ".join(f"{k} (n={int(v)}) med={m:.3f}"
                     for k, v, m in zip(g.size().index, g.size().values,
                                        g.median().reindex(g.size().index).values))


def main() -> None:
    d = pd.read_csv("master.csv")
    # one row per patient: average metrics across multiple aliquots of the same case
    agg = {m: "mean" for m in METRICS}
    for c in ["sex", "site_group", "stage", "grade_num", "smoker", "hpv", "age"]:
        agg[c] = "first"
    p = d.groupby("case_id").agg(agg).reset_index()
    print(f"samples {len(d)} -> patients {len(p)}")

    rows = []
    ordinal = {"stage", "grade_num"}

    for metric in METRICS:
        for var in ["sex", "site_group", "stage", "grade_num", "smoker", "hpv"]:
            sub = p[[metric, var]].dropna()
            if sub[var].nunique() < 2:
                continue
            groups = [g[metric].values for _, g in sub.groupby(var) if len(g) >= 3]
            if len(groups) < 2:
                continue
            if len(groups) == 2:
                stat, pval = stats.mannwhitneyu(groups[0], groups[1],
                                                alternative="two-sided")
                test = "Mann-Whitney U"
            else:
                stat, pval = stats.kruskal(*groups)
                test = "Kruskal-Wallis"
            rec = {
                "metric": metric, "variable": var, "test": test,
                "n": len(sub), "n_groups": len(groups),
                "statistic": round(float(stat), 4),
                "p": float(pval),
                "group_summary": describe(sub, metric, var),
            }
            if var in ordinal:
                rho, prho = stats.spearmanr(sub[var], sub[metric])
                rec["spearman_rho"] = round(float(rho), 4)
                rec["spearman_p"] = float(prho)
            rows.append(rec)

    R = pd.DataFrame(rows)
    R["q_bh"] = np.nan
    for metric, idx in R.groupby("metric").groups.items():
        R.loc[idx, "q_bh"] = bh(R.loc[idx, "p"].values)
    R = R.sort_values(["metric", "p"])
    R.to_csv("s11_clinical_assoc.csv", index=False)

    print("\n=== clinical associations (BH within metric) ===")
    show = R[["metric", "variable", "test", "n", "p", "q_bh"]].copy()
    show["p"] = show["p"].map(lambda x: f"{x:.4f}")
    show["q_bh"] = show["q_bh"].map(lambda x: f"{x:.3f}")
    print(show.to_string(index=False))

    print("\n=== ordinal trend tests ===")
    tr = R.dropna(subset=["spearman_rho"])[
        ["metric", "variable", "spearman_rho", "spearman_p"]]
    print(tr.to_string(index=False))

    print("\n=== group summaries ===")
    for _, r in R.iterrows():
        print(f"  {r.metric} ~ {r.variable}: {r.group_summary}")

    print("\n=== coverage (patients) ===")
    for v in ["sex", "site_group", "stage", "grade_num", "smoker", "hpv"]:
        print(f"  {v:12} {p[v].notna().sum():3}/{len(p)}")


if __name__ == "__main__":
    main()
