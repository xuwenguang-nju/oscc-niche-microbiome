#!/usr/bin/env python3
"""S6-v5: restrict to tongue only, and test finer oral subsites.

WHY
---
S6-v3 grouped the oral cavity by area. The user asked for a finer cut. Tongue is the
largest single subsite in TCGA-HNSC (61 tumour / 16 normal / 16 paired patients), and it
is the most anatomically homogeneous oral site available, so it is the cleanest test of
whether Parvimonas' tumour enrichment exists anywhere in the oral cavity of TCGA-HNSC.

Two tongue readings are reported deliberately:
  tongue_all   GDC "Base of tongue" + "Other and unspecified parts of tongue"
  tongue_NOS   only "Other and unspecified parts of tongue"
Base of tongue is technically the posterior tongue/oropharynx border and has only 1
normal, so it cannot be tested on its own; including it with NOS is also reported so the
choice is visible rather than hidden.

Also computed: floor of mouth (the next-largest subsite with any normals).

OUTPUT
  results/tables/s6v5_subsite_abundance.csv
  results/tables/s6v5_subsite_tests.csv
  results/tables/s6v5_summary.json
Run: TCMA_ROOT=<dir> python3 s6v5_tongue.py
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

SUBSETS = {
    "tongue_all": ["Base of tongue", "Other and unspecified parts of tongue"],
    "tongue_NOS": ["Other and unspecified parts of tongue"],
    "floor_of_mouth": ["Floor of mouth"],
    "oral_others": ["Gum", "Palate", "Lip", "Other and unspecified parts of mouth"],
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


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def s5_genera() -> pd.DataFrame:
    m = pd.read_csv(TAB / "s5_meta_analysis.csv")
    k = int(m["n_cohorts"].max())
    rep = m[(m["n_cohorts"] == k) & (m["pooled_q_bh"] < 0.05) & m["direction_consistent"]]
    rep = rep[["genus", "pooled_effect"]].copy()
    rep["s5_direction"] = np.where(rep["pooled_effect"] > 0, "up_in_tumour",
                                   "down_in_tumour")
    return rep


def main() -> None:
    ab = pd.read_csv(TAB / "s6v2_tcma_genus.csv", index_col=0)
    meta = pd.read_csv(TAB / "s6v2_tcma_samples.csv")
    site = pd.read_csv(TAB / "s6v2_gdc_primary_site.csv", dtype=str)
    meta = meta.merge(site[["sample_barcode", "primary_site"]],
                      on="sample_barcode", how="left").set_index("aliquot_barcode")
    ab = ab.loc[meta.index]

    expct = s5_genera()
    genera = set(ab.columns)
    resolve = {}
    for g in expct["genus"]:
        resolve[g] = g if g in genera else SYNONYM.get(g)

    summary: dict = {"subsets": {}}
    ab_rows, test_rows = [], []

    for name, sites in SUBSETS.items():
        idx = meta.index[meta["primary_site"].isin(sites)]
        mt = meta.loc[idx]
        tum = mt.index[mt["sample_type"] == "Primary Tumor"]
        nor = mt.index[mt["sample_type"] == "Solid Tissue Normal"]
        both = sorted(set(mt.loc[tum, "case_id"]) & set(mt.loc[nor, "case_id"]))
        summary["subsets"][name] = {"sites": sites, "n_tumour": len(tum),
                                    "n_normal": len(nor), "n_pairs": len(both)}
        if len(tum) == 0:
            continue
        for g in expct["genus"]:
            col = resolve[g]
            rec = {"subset": name, "genus": g, "tcma_column": col,
                   "n_tumour": len(tum), "n_normal": len(nor), "n_pairs": len(both),
                   "s5_direction": expct.loc[expct["genus"] == g, "s5_direction"].iloc[0]}
            if col is None or len(nor) == 0:
                rec.update({"mean_tumour": np.nan, "mean_normal": np.nan,
                            "log2fc": np.nan, "direction": "not_testable",
                            "p_unpaired": np.nan, "p_paired": np.nan})
                test_rows.append(rec)
                continue
            a = ab.loc[tum, col].astype(float)
            b = ab.loc[nor, col].astype(float)
            eps = 1e-9
            rec["mean_tumour"] = float(a.mean())
            rec["mean_normal"] = float(b.mean())
            rec["pct_detected_tumour"] = round(100 * float((a > 0).mean()), 1)
            rec["pct_detected_normal"] = round(100 * float((b > 0).mean()), 1)
            rec["log2fc"] = float(np.log2((a.mean() + eps) / (b.mean() + eps)))
            rec["direction"] = "up_in_tumour" if rec["log2fc"] > 0 else "down_in_tumour"
            ab_rows.append({"subset": name, "genus": g, "sample_type": "tumour",
                            "n": len(tum), "mean": rec["mean_tumour"]})
            ab_rows.append({"subset": name, "genus": g, "sample_type": "normal",
                            "n": len(nor), "mean": rec["mean_normal"]})
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
            test_rows.append(rec)

    res = pd.DataFrame(test_rows)
    res["q_unpaired"] = np.nan
    res["q_paired"] = np.nan
    for sub, s in res.groupby("subset"):
        m1 = s["p_unpaired"].notna()
        if m1.any():
            res.loc[s.index[m1], "q_unpaired"] = bh(s.loc[m1, "p_unpaired"].values).round(4)
        m2 = s["p_paired"].notna()
        if m2.any():
            res.loc[s.index[m2], "q_paired"] = bh(s.loc[m2, "p_paired"].values).round(4)
    res.to_csv(TAB / "s6v5_subsite_tests.csv", index=False)
    pd.DataFrame(ab_rows).to_csv(TAB / "s6v5_subsite_abundance.csv", index=False)

    pv = res[res["genus"] == "Parvimonas"][
        ["subset", "n_tumour", "n_normal", "n_pairs", "mean_tumour", "mean_normal",
         "pct_detected_tumour", "pct_detected_normal", "log2fc", "direction",
         "p_unpaired", "q_unpaired", "p_paired", "q_paired"]]
    summary["parvimonas"] = pv.round(4).to_dict("records")

    for sub in res["subset"].unique():
        s = res[res["subset"] == sub]
        t = s[s["direction"].isin(["up_in_tumour", "down_in_tumour"])]
        summary["subsets"][sub].update({
            "testable_genera": int(len(t)),
            "direction_agrees": int(t["direction_agrees"].sum()) if len(t) else 0,
            "q05_paired": int((s["q_paired"] < 0.05).sum())})

    (TAB / "s6v5_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False, default=str), encoding="utf-8")

    print("=== subset sizes ===")
    for k, v in summary["subsets"].items():
        print(f"  {k:<16} T={v['n_tumour']:>3} N={v['n_normal']:>3} pairs={v['n_pairs']:>2}"
              f"  testable={v.get('testable_genera','-')} agree={v.get('direction_agrees','-')}"
              f"  q<0.05={v.get('q05_paired','-')}")
    print("\n=== Parvimonas by subsite ===")
    print(pv.round(4).to_string(index=False))
    print("\n=== tongue_all full table ===")
    tt = res[res["subset"] == "tongue_all"].sort_values("p_paired")
    print(tt[["genus", "s5_direction", "direction", "log2fc", "n_pairs",
              "p_paired", "q_paired", "p_unpaired", "q_unpaired"]]
          .round(4).to_string(index=False))
    print(f"\nwrote {TAB}/s6v5_*.csv|json")


if __name__ == "__main__":
    main()
