#!/usr/bin/env python3
"""S6 -- host-side coupling for the OSCC intratumoral-microbiome project.

WHY THIS DATASET
----------------
The four 16S cohorts used in S5 carry no host transcriptome, so S6 cannot couple
microbial abundance with expression in those patients. GSE227919 is total RNA-seq of
oral samples -- healthy control, hyperkeratosis non-reactive (HkNR), leukoplakia with
dysplasia, and OSCC -- and total RNA-seq retains prokaryotic reads, so ONE library
yields both the microbial profile and the host transcriptome for the SAME sample.

WHAT S6 CAN AND CANNOT CLAIM
----------------------------
CAN: check the DIRECTION of the S5-replicating genera on an independent, completely
     different measurement (RNA-seq-mined microbial reads instead of 16S amplicons),
     and test whether those genera covary with host immune / EMT programs.
CANNOT: replicate S5 effect sizes. RNA-seq microbial reads differ from 16S in depth
     (orders of magnitude lower), poly-A/host-DNA background, and taxonomic bias.
     n=9 OSCC is small. Every number here is a direction/coupling check, not an
     effect-size estimate.

INPUTS
  data/raw/gse227919/GSE227919_2023_20_03_expr_mat_microbe.txt.gz   (space-delimited)
  data/raw/gse227919/GSE227919_2023_20_03_expr_mat_raw.txt.gz       (space-delimited)
  data/raw/gse227919/GSE227919_2023_21_03_pml_annotations.csv.gz
  results/tables/s5_meta_analysis.csv        (uploaded from the project)

OUTPUTS (results/tables/)
  s6a_direction_check.csv        per-genus direction vs the S5 expectation
  s6a_group_abundance.csv        genus relative abundance by Class
  s6b_genus_signature_corr.csv   Spearman rho, genus x host signature
  s6_host_signature_scores.csv   per-sample signature scores
  s6_summary.json

Run: python3 s6_host_coupling.py
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "gse227919"
TAB = ROOT / "results" / "tables"
MIC = RAW / "GSE227919_2023_20_03_expr_mat_microbe.txt.gz"
HOST = RAW / "GSE227919_2023_20_03_expr_mat_raw.txt.gz"
ANN = RAW / "GSE227919_2023_21_03_pml_annotations.csv.gz"
S5 = TAB / "s5_meta_analysis.csv"

TORG = re.compile(r"org\|([^|]+)\|")
# Non-bacterial rows that total-RNA-seq microbial profiling picks up. They are
# excluded so that relative abundance is computed within prokaryotes, matching the
# host-removal rule used for the 16S genus tables.
#
# The patterns must be precise. An earlier version used a bare "homo", which
# silently deleted Stenotrophomonas (it contains "homo" inside "homon") and threw
# away 40% of the reads. Match host/fungal/viral markers explicitly instead.
EXCLUDE = re.compile(
    r"homo_sapiens|human|phage|virus|viral|eukaryot|"
    r"candida|malassezia|saccharomyces|aspergillus|pneumocystis|"
    r"toxoplasma|trichomonas|plasmodium|trypanosoma",
    re.I)

# Host programs to couple against. Kept small and pre-specified (not selected after
# looking at the data) to avoid fishing.
SIGNATURES = {
    "cytotoxic": ["CD8A", "GZMA", "GZMB", "PRF1", "NKG7", "GNLY"],
    "ifn_gamma": ["IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1"],
    "checkpoint": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "pemt": ["CDH2", "VIM", "ZEB1", "SNAI2", "TWIST1", "FN1"],
    "proliferation": ["MKI67", "TOP2A", "PCNA", "CCNB1"],
    "epithelial_diff": ["CDH1", "KRT5", "KRT13", "KRT14", "IVL"],
}


def sniff_sep(path: Path) -> str:
    """GEO matrices are inconsistently delimited; decide from the header line."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
    return "\t" if header.count("\t") > header.count(" ") else r"\s+"


def read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sniff_sep(path), index_col=0)
    if df.shape[1] == 0:
        raise RuntimeError(f"{path.name}: parsed 0 columns")
    return df


def norm_sample(x: str) -> str:
    """'MK-10' (microbe) and 'MK_10' (host/annotations) must map to one id."""
    return re.sub(r"[\s\-]+", "_", str(x).strip())


def genus_of(row_id: str) -> str | None:
    m = TORG.search(row_id)
    if not m:
        return None
    parts = m.group(1).replace("_", " ").split()
    return parts[0] if parts else None


def bacterial_genus_matrix(mic: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Sum organism-level counts to genus, drop non-prokaryotic rows, normalise."""
    rows = pd.DataFrame({"row": mic.index})
    rows["genus"] = [genus_of(r) for r in rows["row"]]
    rows["excluded"] = rows["row"].str.contains(EXCLUDE)

    kept = rows.loc[rows["genus"].notna() & ~rows["excluded"]]
    dropped = rows.loc[rows["genus"].isna() | rows["excluded"]]
    counts = mic.loc[kept["row"].values]
    counts.index = kept["genus"].values
    genus = counts.groupby(level=0).sum()

    rel = genus.div(genus.sum(axis=0), axis=1)
    meta = {
        "taxa_rows_total": int(mic.shape[0]),
        "taxa_rows_kept": int(kept.shape[0]),
        "taxa_rows_dropped_nonbacterial_or_unparsed": int(dropped.shape[0]),
        "n_genera": int(genus.shape[0]),
        "reads_total": float(mic.values.sum()),
        "reads_kept": float(counts.values.sum()),
        "reads_kept_pct": round(100 * counts.values.sum() / mic.values.sum(), 2),
        "dropped_examples": dropped["row"].head(8).tolist(),
    }
    return rel, meta


def load_annotations() -> pd.DataFrame:
    ann = pd.read_csv(ANN, index_col=0)
    ann.index = [norm_sample(i) for i in ann.index]
    if "Class" not in ann.columns:
        raise RuntimeError("expected a 'Class' column in the annotations")
    return ann


def s5_expectations() -> pd.DataFrame:
    """The genera that replicated in S5 (all four cohorts, pooled q<0.05)."""
    m = pd.read_csv(S5)
    k = int(m["n_cohorts"].max())
    rep = m[(m["n_cohorts"] == k) & (m["pooled_q_bh"] < 0.05) & m["direction_consistent"]]
    rep = rep[["genus", "pooled_effect", "pooled_q_bh", "max_abs_dz", "n_cohorts"]].copy()
    rep["s5_direction"] = np.where(rep["pooled_effect"] > 0, "up_in_tumour", "down_in_tumour")
    return rep.sort_values("pooled_effect", ascending=False)


def signature_scores(host: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Mean z-score of log-normalised expression across each pre-specified program."""
    # values are already normalised (not integer counts) -> log1p to stabilise, then z
    lg = np.log1p(host.clip(lower=0))
    z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1).replace(0, np.nan), axis=0)
    scores, used = {}, {}
    for name, genes in SIGNATURES.items():
        present = [g for g in genes if g in z.index]
        used[name] = present
        if present:
            scores[name] = z.loc[present].mean(axis=0)
    return pd.DataFrame(scores), used


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> None:
    summary: dict = {"series": "GSE227919"}

    mic = read_matrix(MIC)
    host = read_matrix(HOST)
    ann = load_annotations()
    summary["input_shapes"] = {"microbe": list(mic.shape), "host": list(host.shape),
                               "annotations": list(ann.shape)}
    summary["class_counts"] = ann["Class"].value_counts().to_dict()

    rel, mic_meta = bacterial_genus_matrix(mic)
    rel.columns = [norm_sample(c) for c in rel.columns]
    summary["microbial"] = mic_meta

    # ---- align samples -----------------------------------------------------
    common = [s for s in rel.columns if s in host.columns and s in ann.index]
    summary["alignment"] = {
        "microbe_samples": int(rel.shape[1]), "host_samples": int(host.shape[1]),
        "annotation_samples": int(ann.shape[0]), "common": len(common),
        "note": "microbe ids use '-' and host/annotation ids use '_'; normalised"}
    if len(common) < 20:
        raise RuntimeError(f"only {len(common)} common samples")
    rel = rel[common]
    host = host[common]
    ann = ann.loc[common]

    # ================= S6a: direction check =================================
    expct = s5_expectations()
    grp_order = ["Control", "HkNR", "Dysplasia", "Cancer"]

    grp_means = (rel.T.groupby(ann["Class"]).mean().T
                 .reindex(columns=[g for g in grp_order if g in ann["Class"].unique()]))
    grp_means.index.name = "genus"
    grp_means.to_csv(TAB / "s6a_group_abundance.csv")

    rows = []
    for _, e in expct.iterrows():
        g = e["genus"]
        if g not in rel.index:
            rows.append({**e.to_dict(), "present_in_rnaseq": False,
                         "rna_cancer_minus_control": np.nan,
                         "rna_direction": "not_detected", "direction_agrees": np.nan,
                         "p_cancer_vs_control": np.nan,
                         "median_rel_abundance": np.nan})
            continue
        absent_frac = float((rel.loc[g] == 0).mean())
        cancer = rel.loc[g, ann["Class"] == "Cancer"]
        control = rel.loc[g, ann["Class"] == "Control"]
        diff = float(cancer.mean() - control.mean())
        rna_dir = "up_in_tumour" if diff > 0 else "down_in_tumour"
        # non-parametric; exploratory only given n=9 cancers
        try:
            p = float(stats.mannwhitneyu(cancer, control, alternative="two-sided").pvalue)
        except ValueError:
            p = np.nan
        rows.append({**e.to_dict(), "present_in_rnaseq": True,
                     "pct_samples_absent": round(100 * absent_frac, 1),
                     "median_rel_abundance": float(rel.loc[g].median()),
                     "rna_cancer_minus_control": diff,
                     "rna_direction": rna_dir,
                     "direction_agrees": bool(rna_dir == e["s5_direction"]),
                     "p_cancer_vs_control": p})
    dir_tab = pd.DataFrame(rows)
    dir_tab.to_csv(TAB / "s6a_direction_check.csv", index=False)

    tested = dir_tab[dir_tab["present_in_rnaseq"].astype(bool)]
    n_agree = int(tested["direction_agrees"].sum())
    n_tested = int(len(tested))
    # sign test against the S5 expectation
    if n_tested:
        sign_p = float(stats.binomtest(n_agree, n_tested, 0.5).pvalue)
    else:
        sign_p = np.nan
    summary["s6a"] = {"genera_from_s5": int(len(dir_tab)),
                      "detected_in_rnaseq": n_tested,
                      "direction_agrees": n_agree,
                      "direction_agreement_pct": round(100 * n_agree / n_tested, 1) if n_tested else None,
                      "sign_test_p": sign_p,
                      "note": "direction check only; effect sizes are not comparable "
                              "between 16S amplicons and RNA-seq-mined microbial reads"}

    # ================= S6b: host coupling ===================================
    sig, used = signature_scores(host)
    sig.index = [norm_sample(i) for i in sig.index]
    sig = sig.loc[common]
    sig.to_csv(TAB / "s6_host_signature_scores.csv")
    summary["s6b_signatures"] = {k: v for k, v in used.items()}

    # couple: genera that replicated in S5 AND are detected here
    focus = [g for g in expct["genus"] if g in rel.index]
    corr_rows = []
    for g in focus:
        for s in sig.columns:
            rho, p = stats.spearmanr(rel.loc[g], sig[s])
            corr_rows.append({"genus": g, "signature": s,
                              "spearman_rho": round(float(rho), 3),
                              "p": float(p)})
    corr = pd.DataFrame(corr_rows)
    corr["q_bh"] = bh(corr["p"].values).round(4)
    corr = corr.sort_values("p")
    corr.to_csv(TAB / "s6b_genus_signature_corr.csv", index=False)

    summary["s6b"] = {
        "genera_tested": len(focus),
        "signatures_tested": int(sig.shape[1]),
        "tests": int(len(corr)),
        "significant_q05": int((corr["q_bh"] < 0.05).sum()),
        "top": corr.head(8)[["genus", "signature", "spearman_rho", "q_bh"]]
                    .to_dict("records"),
        "note": "pre-specified host programs; n is small, q-values are within this "
                "coupling panel only",
    }

    (TAB / "s6_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False, default=str), encoding="utf-8")

    print(json.dumps(summary, indent=1, ensure_ascii=False, default=str)[:5000])
    print("\n--- S6a direction check ---")
    cols = ["genus", "s5_direction", "rna_direction", "direction_agrees",
            "p_cancer_vs_control"]
    print(dir_tab[cols].to_string(index=False))
    print("\n--- S6b top correlations ---")
    print(corr.head(10).to_string(index=False))
    print(f"\nwrote {TAB}/s6_*.csv|json")


if __name__ == "__main__":
    main()
