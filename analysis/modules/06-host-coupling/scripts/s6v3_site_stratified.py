#!/usr/bin/env python3
"""S6-v3: stratify TCGA-HNSC by anatomic site and re-test the S5 genera.

WHY
---
Pooled TCMA analysis (S6-v2) had Parvimonas going DOWN in tumour (-0.33, q=1.0), the
opposite of the four 16S cohorts where it was up in every one. Two hypotheses were
offered: (a) TCMA is insensitive to low-abundance obligate anaerobes, (b) pooling oral
cavity with pharynx/larynx masks a site-specific signal. This script tests (b) by
splitting the cohort by primary_site group.

SITE GROUPS
  oral       tongue (base + NOS), floor of mouth, gum, palate, lip, other mouth,
             other ill-defined lip/oral cavity/pharynx  -> oral cavity
  pharynx    tonsil, oropharynx, hypopharynx
  larynx     larynx
The GDC "Other and ill-defined sites in lip, oral cavity and pharynx" bucket is not
resolvable to a subsite, so it is kept in `oral` but its ambiguity is reported.

TESTS (per group, per genus)
  unpaired  Mann-Whitney tumour vs adjacent normal within the group
  paired    Wilcoxon on within-patient differences, where a patient contributes both
Only groups with >= 5 normals can be tested unpaired; paired needs >= 5 pairs.

OUTPUT
  results/tables/s6v3_site_group_abundance.csv
  results/tables/s6v3_stratified_tests.csv
  results/tables/s6v3_summary.json
Run: TCMA_ROOT=<dir> python3 s6v3_site_stratified.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SITE_GROUP = {
    "Base of tongue": "oral",
    "Other and unspecified parts of tongue": "oral",
    "Floor of mouth": "oral",
    "Gum": "oral",
    "Palate": "oral",
    "Lip": "oral",
    "Other and unspecified parts of mouth": "oral",
    "Other and ill-defined sites in lip, oral cavity and pharynx": "oral*",
    "Tonsil": "pharynx",
    "Oropharynx": "pharynx",
    "Hypopharynx": "pharynx",
    "Larynx": "larynx",
}
SYNONYM = {"Cutibacterium": "Propionibacterium",
           "Schaalia": "Actinobaculum",
           "Lancefieldella": "Atopobium"}


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


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def s5_expectations() -> pd.DataFrame:
    m = pd.read_csv(TAB / "s5_meta_analysis.csv")
    k = int(m["n_cohorts"].max())
    rep = m[(m["n_cohorts"] == k) & (m["pooled_q_bh"] < 0.05) & m["direction_consistent"]]
    rep = rep[["genus", "pooled_effect", "pooled_q_bh"]].copy()
    rep["s5_direction"] = np.where(rep["pooled_effect"] > 0, "up_in_tumour",
                                   "down_in_tumour")
    return rep


def main() -> None:
    ab = pd.read_csv(TAB / "s6v2_tcma_genus.csv", index_col=0)
    meta = pd.read_csv(TAB / "s6v2_tcma_samples.csv")
    site = pd.read_csv(TAB / "s6v2_gdc_primary_site.csv", dtype=str)
    meta = meta.merge(site[["sample_barcode", "primary_site"]],
                      on="sample_barcode", how="left")
    meta["site_group"] = meta["primary_site"].map(SITE_GROUP)
    meta = meta.set_index("aliquot_barcode")
    ab = ab.loc[meta.index]

    summary: dict = {
        "dataset": "TCMA v2 / TCGA-HNSC, stratified by GDC primary_site",
        "site_group_counts": pd.crosstab(meta["site_group"], meta["sample_type"])
                               .to_dict(),
        "site_detail": pd.crosstab(meta["primary_site"], meta["sample_type"])
                         .to_dict(),
        "note": "site group marked * contains GDC's ambiguous 'other and ill-defined "
                "sites in lip, oral cavity and pharynx' bucket",
    }

    expct = s5_expectations()
    genera = list(ab.columns)
    resolve: dict[str, str | None] = {}
    for g in expct["genus"]:
        if g in genera:
            resolve[g] = g
        elif SYNONYM.get(g) in genera:
            resolve[g] = SYNONYM[g]
        else:
            resolve[g] = None

    # group-level abundance table
    grp_rows = []
    for sg in sorted(meta["site_group"].dropna().unique()):
        idx = meta.index[meta["site_group"] == sg]
        mt = meta.loc[idx]
        for st in ["Primary Tumor", "Solid Tissue Normal"]:
            sel = mt.index[mt["sample_type"] == st]
            if len(sel) == 0:
                continue
            means = ab.loc[sel].mean()
            for g in expct["genus"]:
                col = resolve[g]
                grp_rows.append({"site_group": sg, "sample_type": st,
                                 "n": len(sel), "genus": g,
                                 "mean_rel_abundance":
                                     float(means[col]) if col else np.nan})
    grp = pd.DataFrame(grp_rows)
    grp.to_csv(TAB / "s6v3_site_group_abundance.csv", index=False)

    # stratified tests
    rows = []
    for sg in sorted(meta["site_group"].dropna().unique()):
        idx = meta.index[meta["site_group"] == sg]
        mt = meta.loc[idx]
        tum = mt.index[mt["sample_type"] == "Primary Tumor"]
        nor = mt.index[mt["sample_type"] == "Solid Tissue Normal"]
        both = sorted(set(mt.loc[tum, "case_id"]) & set(mt.loc[nor, "case_id"]))
        for g in expct["genus"]:
            col = resolve[g]
            rec = {"site_group": sg, "genus": g, "tcma_column": col,
                   "n_tumour": len(tum), "n_normal": len(nor), "n_pairs": len(both),
                   "s5_direction": expct.loc[expct["genus"] == g, "s5_direction"].iloc[0]}
            if col is None or len(nor) == 0:
                rec.update({"mean_tumour": np.nan, "mean_normal": np.nan,
                            "log2fc": np.nan, "direction": "not_testable",
                            "p_unpaired": np.nan, "p_paired": np.nan})
                rows.append(rec)
                continue
            a = ab.loc[tum, col].astype(float)
            b = ab.loc[nor, col].astype(float)
            eps = 1e-9
            rec["mean_tumour"] = float(a.mean())
            rec["mean_normal"] = float(b.mean())
            rec["log2fc"] = float(np.log2((a.mean() + eps) / (b.mean() + eps)))
            rec["direction"] = "up_in_tumour" if rec["log2fc"] > 0 else "down_in_tumour"
            p_un = np.nan
            if len(tum) >= 3 and len(nor) >= 3:
                try:
                    p_un = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
                except ValueError:
                    p_un = np.nan
            rec["p_unpaired"] = p_un
            p_pair = np.nan
            if len(both) >= 5:
                av = a.groupby(meta.loc[a.index, "case_id"]).mean()
                bv = b.groupby(meta.loc[b.index, "case_id"]).mean()
                d = np.array([av[p] - bv[p] for p in both if p in av.index and p in bv.index])
                if len(d) >= 5 and np.any(d != 0):
                    p_pair = float(stats.wilcoxon(d, alternative="two-sided").pvalue)
            rec["p_paired"] = p_pair
            rec["direction_agrees"] = rec["direction"] == rec["s5_direction"]
            rows.append(rec)

    res = pd.DataFrame(rows)
    # FDR within each site group (separate families per group)
    res["q_unpaired"] = np.nan
    res["q_paired"] = np.nan
    for sg, sub in res.groupby("site_group"):
        m1 = sub["p_unpaired"].notna()
        if m1.any():
            res.loc[sub.index[m1], "q_unpaired"] = bh(sub.loc[m1, "p_unpaired"].values).round(4)
        m2 = sub["p_paired"].notna()
        if m2.any():
            res.loc[sub.index[m2], "q_paired"] = bh(sub.loc[m2, "p_paired"].values).round(4)
    res.to_csv(TAB / "s6v3_stratified_tests.csv", index=False)

    # focus on Parvimonas
    pv = res[res["genus"] == "Parvimonas"][
        ["site_group", "n_tumour", "n_normal", "n_pairs", "mean_tumour", "mean_normal",
         "log2fc", "direction", "p_unpaired", "p_paired"]]
    summary["parvimonas_by_site"] = pv.round(4).to_dict("records")

    # which site groups still reproduce the S5 direction?
    agree = {}
    for sg, sub in res.groupby("site_group"):
        t = sub[sub["direction"].isin(["up_in_tumour", "down_in_tumour"])]
        agree[sg] = {"testable_genera": int(len(t)),
                     "direction_agrees": int(t["direction_agrees"].sum()),
                     "q05_paired": int((sub["q_paired"] < 0.05).sum())}
    summary["direction_by_site"] = agree

    (TAB / "s6v3_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False, default=str), encoding="utf-8")

    print("=== site group x sample type ===")
    print(pd.crosstab(meta["site_group"], meta["sample_type"]).to_string())
    print("\n=== Parvimonas by site ===")
    print(pv.round(4).to_string(index=False))
    print("\n=== direction agreement by site ===")
    for k, v in agree.items():
        print(f"  {k:<9} testable={v['testable_genera']:>2}  agree={v['direction_agrees']:>2}"
              f"  q<0.05(paired)={v['q05_paired']}")
    print(f"\nwrote {TAB}/s6v3_*.csv|json")


if __name__ == "__main__":
    main()
