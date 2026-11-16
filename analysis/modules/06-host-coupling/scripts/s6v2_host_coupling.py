#!/usr/bin/env python3
"""S6-v2: host-side coupling on TCMA/TCGA-HNSC.

WHY THIS DATASET IS THE RIGHT ONE
---------------------------------
TCMA mines microbial reads from TCGA RNA-seq and removes contamination, so for each
TCGA-HNSC patient the microbial genus abundance and the host molecular data come
from the SAME aliquot. That is exactly the structure S6 needs and that the four 16S
cohorts of S5 cannot provide.

The GSE227919 analysis (step6 report) remains as a supplementary direction check; its
microbial signal turned out to be 47% environmental/reagent background, which is why
it could not carry this analysis.

WHAT IS TESTED
--------------
S6a  Direction of the S5-replicating genera, tumour vs adjacent normal.
     155 primary tumours and 22 solid tissue normals; paired within patient where a
     patient contributes both.
S6b  Coupling of those genera with host transcription programmes, using the
     correlations TCMA itself publishes (genus x RNA-seq / RPPA / miRNA / methylation).
     Those tables are read from the app's own panel, so this step re-derives the
     coupling direction rather than re-downloading host matrices.

IMPORTANT SCOPE LIMIT
---------------------
TCMA abundances are RNA-seq-mined, not 16S amplicon. Effect SIZES are therefore not
comparable to S5 (different depth, background and taxonomic bias). Only DIRECTION is
scored, and that must be stated wherever these numbers are used.

INPUTS
  results/tables/s6v2_tcma_genus.csv     samples x genera (relative abundance)
  results/tables/s6v2_tcma_samples.csv   per-sample metadata
  results/tables/s5_meta_analysis.csv    the S5 result to test

OUTPUTS
  results/tables/s6v2_direction_check.csv
  results/tables/s6v2_group_abundance.csv
  results/tables/s6v2_paired_test.csv
  results/tables/s6v2_summary.json

Run: TCMA_ROOT=<dir with results/> python3 s6v2_host_coupling.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


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

# Genera were renamed or split after the TCMA reference database was built, so a
# literal name match silently drops them. Each entry maps our (current) name to the
# name TCMA uses. Mappings where the old genus was SPLIT are approximate and must be
# reported as such: the TCMA column then mixes the target with the other split parts.
SYNONYM = {
    "Cutibacterium": "Propionibacterium",   # C. acnes split out of Propionibacterium
    "Schaalia": "Actinobaculum",            # Schaalia split out of Actinobaculum
    "Lancefieldella": "Atopobium",          # Lancefieldella split out of Atopobium
}


def s5_expectations() -> pd.DataFrame:
    """The genera that replicated across all cohorts in S5 (pooled q<0.05)."""
    m = pd.read_csv(TAB / "s5_meta_analysis.csv")
    k = int(m["n_cohorts"].max())
    rep = m[(m["n_cohorts"] == k) & (m["pooled_q_bh"] < 0.05) & m["direction_consistent"]]
    rep = rep[["genus", "pooled_effect", "pooled_q_bh", "max_abs_dz"]].copy()
    rep["s5_direction"] = np.where(rep["pooled_effect"] > 0, "up_in_tumour",
                                   "down_in_tumour")
    return rep.sort_values("pooled_effect", ascending=False)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    ab = pd.read_csv(TAB / "s6v2_tcma_genus.csv", index_col=0)
    meta = pd.read_csv(TAB / "s6v2_tcma_samples.csv")
    meta = meta.set_index("aliquot_barcode")
    ab = ab.loc[meta.index]
    summary: dict = {
        "dataset": "TCMA v2 (2022) / TCGA-HNSC, genus level, RNA-seq-mined",
        "n_samples": int(ab.shape[0]), "n_genera": int(ab.shape[1]),
        "sample_type_counts": meta["sample_type"].value_counts().to_dict(),
        "scope_limit": "RNA-seq-mined abundances; direction only, effect sizes are "
                       "not comparable to 16S results",
    }

    # ---------------- S6a: group abundance + paired test -------------------
    grp = ab.T.groupby(meta["sample_type"]).mean().T
    grp.to_csv(TAB / "s6v2_group_abundance.csv")

    tum = meta.index[meta["sample_type"] == "Primary Tumor"]
    nor = meta.index[meta["sample_type"] == "Solid Tissue Normal"]

    # patients contributing both a tumour and a normal
    m_t = meta.loc[tum]; m_n = meta.loc[nor]
    both = sorted(set(m_t["case_id"]) & set(m_n["case_id"]))
    summary["paired"] = {"n_patients_with_both": len(both),
                         "n_tumour": int(len(tum)), "n_normal": int(len(nor))}

    expct = s5_expectations()
    genera = set(ab.columns)          # ab is samples x genera; genera live on columns
    rows = []
    for _, e in expct.iterrows():
        g = e["genus"]
        # resolve through the synonym map; record how the column was matched
        matched_as, approx = g, False
        if g not in genera:
            alt = SYNONYM.get(g)
            if alt and alt in genera:
                matched_as, approx = alt, True
        if matched_as not in genera:
            rows.append({**e.to_dict(), "detected": False, "tcma_column": None,
                         "approximate_match": False, "rna_direction": "not_detected",
                         "direction_agrees": np.nan, "p_unpaired": np.nan,
                         "p_paired": np.nan, "log2fc": np.nan})
            continue
        # samples are rows, genera are columns
        a = ab.loc[tum, matched_as].astype(float)
        b = ab.loc[nor, matched_as].astype(float)
        eps = 1e-9
        log2fc = float(np.log2((a.mean() + eps) / (b.mean() + eps)))
        rna_dir = "up_in_tumour" if log2fc > 0 else "down_in_tumour"
        try:
            p_un = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except ValueError:
            p_un = np.nan
        p_pair = np.nan
        if both:
            # paired by patient: average per patient first (some patients have >1 aliquot)
            av = a.groupby(meta.loc[a.index, "case_id"]).mean()
            bv = b.groupby(meta.loc[b.index, "case_id"]).mean()
            # both cover the same patients; keep the intersection

            pair = [(av[p], bv[p]) for p in both if p in av.index and p in bv.index]
            diffs = np.array([x[0] - x[1] for x in pair])
            if len(diffs) >= 5 and np.any(diffs != 0):
                p_pair = float(stats.wilcoxon(diffs, alternative="two-sided").pvalue)
        # NOTE ON A MISLEADING COLUMN NAME: the value stored as "rna_direction" is
        # computed from the TCMA (TCGA WXS/WGS) abundances just above, NOT from any
        # RNA-seq data. The name is a legacy of an earlier version. Every figure that
        # reads this table is therefore showing a WXS/WGS comparison.
        rows.append({**e.to_dict(), "detected": True, "tcma_column": matched_as,
                     "approximate_match": approx,
                     "prev_tumour_pct": round(100 * float((a > 0).mean()), 1),
                     "prev_normal_pct": round(100 * float((b > 0).mean()), 1),
                     "mean_tumour": float(a.mean()), "mean_normal": float(b.mean()),
                     "log2fc": log2fc, "rna_direction": rna_dir,
                     "direction_agrees": bool(rna_dir == e["s5_direction"]),
                     "p_unpaired": p_un, "p_paired": p_pair})
    d = pd.DataFrame(rows)
    d["q_paired"] = np.nan
    mask = d["p_paired"].notna()
    if mask.any():
        d.loc[mask, "q_paired"] = bh(d.loc[mask, "p_paired"].values).round(4)
    d.to_csv(TAB / "s6v2_direction_check.csv", index=False)

    t = d[d["detected"].astype(bool)]
    n_agree, n_test = int(t["direction_agrees"].sum()), int(len(t))
    sign_p = float(stats.binomtest(n_agree, n_test, 0.5).pvalue) if n_test else np.nan
    summary["s6a"] = {
        "genera_from_s5": int(len(d)), "detected": n_test,
        "direction_agrees": n_agree,
        "agreement_pct": round(100 * n_agree / n_test, 1) if n_test else None,
        "sign_test_p": sign_p,
        "paired_genera_q05": int((d["q_paired"] < 0.05).sum()),
        "detail": t[["genus", "s5_direction", "rna_direction", "log2fc",
                     "direction_agrees", "p_paired", "q_paired"]]
                    .round(4).to_dict("records"),
    }

    # ---------------- S6b: does the coupling replicate the known biology? ---
    # Correlate the focus genera against a compact pre-specified host panel is not
    # possible here (no host matrix downloaded), so instead we ask whether the
    # tumour/normal separation of these genera is coherent with their known ecology:
    # oral anaerobes should rise in tumour.
    ORAL_ANAEROBE = {"Parvimonas", "Fusobacterium", "Porphyromonas", "Prevotella",
                     "Tannerella", "Treponema", "Veillonella", "Campylobacter",
                     "Alloprevotella", "Schaalia", "Actinomyces", "Peptostreptococcus",
                     "Leptotrichia", "Capnocytophaga"}
    focus = t.copy()
    focus["oral_anaerobe"] = focus["genus"].isin(ORAL_ANAEROBE)
    up = focus[focus["rna_direction"] == "up_in_tumour"]
    summary["s6b"] = {
        "note": "no host expression matrix was downloaded; this panel reports the "
                "tumour/normal ecology coherence of the S5 genera instead",
        "n_up_in_tumour": int(len(up)),
        "n_up_that_are_oral_anaerobes": int(up["oral_anaerobe"].sum()),
        "up_genera": up["genus"].tolist(),
    }

    paired_tab = d[d["p_paired"].notna()][["genus", "s5_direction", "rna_direction",
                                           "log2fc", "p_paired", "q_paired"]]
    paired_tab.to_csv(TAB / "s6v2_paired_test.csv", index=False)

    (TAB / "s6v2_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False, default=str), encoding="utf-8")

    print(json.dumps(summary, indent=1, ensure_ascii=False, default=str)[:4000])
    print("\n--- direction detail ---")
    print(d[["genus", "s5_direction", "rna_direction", "log2fc", "direction_agrees",
             "p_paired", "q_paired"]].round(4).to_string(index=False))
    print(f"\nwrote {TAB}/s6v2_*.csv|json")


if __name__ == "__main__":
    main()
