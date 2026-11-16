#!/usr/bin/env python3
"""Regenerate Table S1 with the effect in log2, as the figures and the text now use.

Why this exists
---------------
The pooled genus effects are natural-log ratios (`np.log` in cross_cohort_meta.paired_effect), but
Table S1 labelled the column "Pooled log2FC" and Fig. 4L labelled its axis "mean log-ratio". The two
therefore read as a factor of ln 2 apart, and Fig. 5J drew a ln quantity against a log2 one with an
identity line between them.

The written convention is now log2 everywhere. This converts the pooled effect, its standard error
and its interval to log2 and writes the column name that matches, so the table, the panel and the
Methods cannot disagree. The underlying analysis is untouched: conversion is multiplication by
1/ln 2, applied to the estimate and its interval so the interval stays symmetric.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
TAB = ROOT / "results/tables"
L2 = np.log(2)

# Prior OSCC reporting is NOT transcribed here. It comes from s11_novelty_assessment.csv, which
# s11_novelty.py derives from the ORMCKB curation. A hand-copied dict previously sat in this file
# and disagreed with that table for two genera (Schaalia 4 vs 0, Cutibacterium 6 vs 0), which put
# "five genera with no prior report" into a caption while the source table said three.


def main() -> None:
    pooled = pd.read_csv(TAB / "s5_pooled_dz_and_logratio.csv", encoding="utf-8-sig")
    rep = pd.read_csv(TAB / "s5_replicating_genera.csv", encoding="utf-8-sig") \
        if (TAB / "s5_replicating_genera.csv").is_file() else None

    # the thirteen replicating genera, ordered by pooled effect
    if rep is not None and "genus" in rep.columns:
        names = rep.sort_values("pooled_effect")["genus"].tolist()
    else:
        names = ["Actinomyces", "Bifidobacterium", "Veillonella", "Staphylococcus",
                 "Alloprevotella", "Cardiobacterium", "Lancefieldella", "Paracoccus",
                 "Actinobacillus", "Schaalia", "Scardovia", "Cutibacterium", "Parvimonas"]
    d = pooled[pooled.genus.isin(names)].copy()
    d["order"] = d.genus.map({g: i for i, g in enumerate(names)})
    d = d.sort_values("order")

    out = pd.DataFrame({
        "Genus": d.genus.values,
        "Cohorts tested (n)": d.k.values,
        "Pooled log2 ratio": (d.pooled_logratio / L2).round(3).values,
        "SE": (d.pooled_logratio_se / L2).round(3).values,
        "95% CI": [f"{lo / L2:+.2f} to {hi / L2:+.2f}" for lo, hi in
                   zip(d.pooled_logratio - 1.96 * d.pooled_logratio_se,
                       d.pooled_logratio + 1.96 * d.pooled_logratio_se)],
        "FDR q": d.pooled_logratio.map(lambda v: "").values,
        "I2 (%)": d.I2_logratio.round(1).values,
        "Largest |d_z|": d.max_abs_dz.round(3).values,
    })
    # FDR q lives in the meta-analysis table, keyed per genus
    meta = pd.read_csv(TAB / "s5_meta_analysis.csv", encoding="utf-8-sig")
    q = meta.set_index("genus")["pooled_q_bh"]        # BH across pooled genera, the tabulated q
    out["FDR q"] = out.Genus.map(q).map(lambda v: f"{v:.4g}" if pd.notna(v) else "")
    nov = pd.read_csv(TAB / "s11_novelty_assessment.csv", encoding="utf-8-sig").set_index("genus")
    papers = out.Genus.map(nov["oscc_papers"])
    direction = out.Genus.map(nov["lit_majority"])
    out["Prior OSCC reports (n)"] = papers.fillna(0).astype(int).values
    out["Literature direction"] = [
        "none reported" if (pd.isna(d) or p == 0) else str(d)
        for p, d in zip(papers, direction)]
    out["Agrees with literature"] = [
        "n/a" if p == 0 or pd.isna(d) else
        ("conflicting" if str(d) == "conflict" else ("yes" if a else "no"))
        for p, d, a in zip(papers, direction, out.Genus.map(nov["direction_agrees"]))]

    out.to_csv(TAB / "tableS1_replicating_genera.csv", index=False, encoding="utf-8-sig")
    print(out.to_string(index=False))
    print(f"\nwrote {TAB / 'tableS1_replicating_genera.csv'} ({len(out)} genera, log2)")


if __name__ == "__main__":
    main()
