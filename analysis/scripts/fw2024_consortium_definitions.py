#!/usr/bin/env python3
"""Compute the anaerobic:aerobic composite shift under three genus-set definitions.

Motivation: the original genus sets used for this axis were selected from the differential
result of the same cohort, which makes any AUC computed from them circular, and three genera
were mis-assigned by oxygen tolerance. The axis is therefore recomputed under (i) the original
sets, for transparency only, (ii) all genera with an unambiguous oxygen-tolerance assignment,
and (iii) the same with contested genera removed. All three are reported.

Output: results/tables/s5_consortium_three_definitions_{per_cohort,pooled,inhouse}.csv
"""
from __future__ import annotations
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

import importlib.util
ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "down", ROOT / "analysis/scripts/fw2024_downsampling.py")
down = importlib.util.module_from_spec(spec); spec.loader.exec_module(down)

spec2 = importlib.util.spec_from_file_location(
    "ccm", ROOT / "analysis/modules/05-meta-analysis/scripts/cross_cohort_meta.py")
ccm = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(ccm)

TAB = ROOT / "results/tables"
STUDIES = ["PRJNA1415359", "PRJNA1415268", "PRJNA597251", "PRJNA666746", "PRJNA412445"]
pid = lambda s: re.sub(r"^[NTY]", "", s)
grp = lambda s: "T" if s[0] in "TY" else "N"


def inhouse(A, E):
    c = pd.read_csv(ROOT / "results/derived/fw2024_420/genus_counts.csv",
                    index_col=0, encoding="utf-8-sig")
    rel = c.div(c.sum(axis=0), axis=1)
    a = [g for g in A if g in rel.index]; e = [g for g in E if g in rel.index]
    d = pd.DataFrame({"ana": rel.reindex(a).sum(), "aer": rel.reindex(e).sum()})
    d["logratio"] = np.log10((d.ana + 1e-9) / (d.aer + 1e-9))
    d["grp"] = [grp(s) for s in d.index]; d["pid"] = [pid(s) for s in d.index]
    p = d.pivot_table(index="pid", columns="grp",
                      values=["ana", "aer", "logratio"]).dropna()
    return p, len(a), len(e)


def cohort(study, A, E):
    g = ccm.load_cohort(study)
    if g is None:
        return None
    c, m = g
    e = ccm.paired_effect(c, m)
    if e.empty:
        return None
    nn = e.normal_niche_used.iloc[0]
    R = c.div(c.sum(axis=1), axis=0)
    a = [x for x in A if x in R.columns]; ee = [x for x in E if x in R.columns]
    d = pd.DataFrame({"ana": R[a].sum(axis=1), "aer": R[ee].sum(axis=1)})
    d["logratio"] = np.log10((d.ana + 1e-9) / (d.aer + 1e-9))
    d["niche"] = m.reindex(d.index)["niche"].values
    d["pid"] = m.reindex(d.index)["patient_id"].values
    p = d[d.niche.isin(["tumor", nn])].pivot_table(
        index="pid", columns="niche", values="logratio").dropna()
    if len(p) < 5:
        return None
    dif = (p["tumor"] - p[nn]).values
    return dict(n=len(p), effect=float(dif.mean()),
                se=float(dif.std(ddof=1) / np.sqrt(len(dif))),
                p=stats.wilcoxon(dif).pvalue)


def main():
    per, pooled, inh = [], [], []
    for tag, (A, E) in down.CONSORTIUM_DEFINITIONS.items():
        p, na, ne = inhouse(A, E)
        lr = p[("logratio", "T")] - p[("logratio", "N")]
        y = np.r_[np.zeros(len(p)), np.ones(len(p))]
        s = np.r_[p[("logratio", "N")].values, p[("logratio", "T")].values]
        inh.append(dict(variant=tag, n_ana=na, n_aer=ne, n=len(p),
                        ana_N=p[("ana", "N")].median() * 100,
                        ana_T=p[("ana", "T")].median() * 100,
                        aer_N=p[("aer", "N")].median() * 100,
                        aer_T=p[("aer", "T")].median() * 100,
                        ana_up=int((p[("ana", "T")] > p[("ana", "N")]).sum()),
                        aer_dn=int((p[("aer", "T")] < p[("aer", "N")]).sum()),
                        p_ana=stats.wilcoxon(p[("ana", "T")], p[("ana", "N")]).pvalue,
                        p_aer=stats.wilcoxon(p[("aer", "T")], p[("aer", "N")]).pvalue,
                        lr_up=int((lr > 0).sum()),
                        p_lr=stats.wilcoxon(lr).pvalue,
                        auc=roc_auc_score(y, s)))
        for study in STUDIES:
            r = cohort(study, A, E)
            if r:
                per.append(dict(variant=tag, study=study, **r))
        per.append(dict(variant=tag, study="hospital cohort", n=len(lr),
                        effect=float(lr.mean()),
                        se=float(lr.std(ddof=1) / np.sqrt(len(lr))),
                        p=stats.wilcoxon(lr).pvalue))
        T = pd.DataFrame([r for r in per if r["variant"] == tag])
        pl = ccm.dersimonian_laird(T.effect.values, T.se.values ** 2)
        pooled.append(dict(variant=tag, n_ana=na, n_aer=ne, k=len(T),
                           n_pos=int((T.effect > 0).sum()),
                           pooled=pl["pooled"], lo=pl["pooled"] - 1.96 * pl["se"],
                           hi=pl["pooled"] + 1.96 * pl["se"],
                           I2=pl["I2"], p=pl["p_pooled"]))
        print(f"{tag:16} k={len(T)} positive={int((T.effect>0).sum())} "
              f"pooled={pl['pooled']:+.3f} I2={pl['I2']:.0f}% p={pl['p_pooled']:.1e}")

    pd.DataFrame(per).to_csv(TAB / "s5_consortium_three_definitions_per_cohort.csv",
                             index=False, encoding="utf-8-sig")
    pd.DataFrame(pooled).to_csv(TAB / "s5_consortium_three_definitions_pooled.csv",
                                index=False, encoding="utf-8-sig")
    pd.DataFrame(inh).to_csv(TAB / "s5_consortium_three_definitions_inhouse.csv",
                             index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()
