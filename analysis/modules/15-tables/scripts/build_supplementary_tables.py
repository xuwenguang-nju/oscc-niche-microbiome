#!/usr/bin/env python3
"""Build the supplementary tables.

S1  all 615 genera x cross-cohort meta-analysis
S2  TCMA / TCGA-HNSC sample annotation (the 177 microbiomes and their covariates)
S3  clinical association + survival models in TCGA-HNSC
S4  species resolution and the genus-unassigned family test
S5  platform concordance (RNA-seq-mined TCMA versus amplicon 16S)

Nothing here is hand-typed: each table is a projection of a results/tables file that an
earlier analysis module wrote.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
TAB = ROOT / "results" / "tables"


def s1_all_genera() -> None:
    m = pd.read_csv(TAB / "s5_meta_analysis.csv")
    keep = ["genus", "n_cohorts", "cohorts", "pooled_effect", "pooled_se", "pooled_p",
            "pooled_q_bh", "Q", "I2_pct", "p_heterogeneity", "n_cohorts_positive",
            "n_cohorts_negative", "direction_consistent", "max_abs_dz",
            "n_cohorts_large_dz"]
    t = m[keep].rename(columns={
        "genus": "Genus", "n_cohorts": "Cohorts tested (n)", "cohorts": "Cohorts",
        "pooled_effect": "Pooled log2 fold change", "pooled_se": "SE",
        "pooled_p": "Pooled p", "pooled_q_bh": "FDR q", "Q": "Cochran Q",
        "I2_pct": "Heterogeneity I2 (%)", "p_heterogeneity": "Heterogeneity p",
        "n_cohorts_positive": "Cohorts positive", "n_cohorts_negative": "Cohorts negative",
        "direction_consistent": "Direction consistent",
        "max_abs_dz": "Largest |d_z|", "n_cohorts_large_dz": "Cohorts with |d_z|>=0.8",
    }).sort_values(["FDR q", "Pooled log2 fold change"])
    t.to_csv(TAB / "tableS1_all_genera_meta.csv", index=False, encoding="utf-8-sig")
    print(f"S1 all-genera meta        {t.shape}  "
          f"(q<0.05: {(t['FDR q'] < 0.05).sum()}, tested in >=3 cohorts: "
          f"{(t['Cohorts tested (n)'] >= 3).sum()})")


def s2_tcma_samples() -> None:
    # master.csv carries the covariates but not the tissue type; the TCMA export
    # (s6v2_tcma_samples.csv) carries sample_type. Join on the sample barcode.
    m = pd.read_csv(TAB / "master.csv")
    tc = pd.read_csv(TAB / "s6v2_tcma_samples.csv")
    if "sample_type" not in m.columns:
        m = m.merge(tc[["sample_barcode", "sample_type"]].drop_duplicates("sample_barcode"),
                    on="sample_barcode", how="left")
    keep = {
        "aliquot": "sample_barcode", "case": "case_id", "sample_type": "sample_type",
        "site": "primary_site", "site_group": "site_group",
        "stage": "ajcc_pathologic_stage", "grade": "tumor_grade",
        "hpv_ish": "hpv_ish", "hpv_p16": "hpv_p16",
        "smoking": "tobacco_smoking_status", "pack_years": "pack_years",
        "alcohol": "alcohol_history", "age": "age", "sex": "sex",
        "shannon": "shannon", "richness": "richness",
        "parvimonas_pct": "parvimonas_rel",   # fraction, converted below
        # os_time is GDC days_to_death / days_to_last_follow_up, never converted. It used to be
        # exported as "os_months" without conversion, so the column said months and held days
        # (a maximum of 3,981 read as 332 years). Both units are now exported and labelled.
        "os_days": "os_time", "os_event": "os_event",
    }
    missing = [c for c in keep.values() if c not in m.columns]
    if missing:
        raise SystemExit(f"master.csv missing: {missing}")
    t = m[list(keep.values())].rename(columns={v: k for k, v in keep.items()})
    # parvimonas_rel is a proportion of the community, not a percentage; convert so the
    # column name matches the values
    t["parvimonas_pct"] = (t["parvimonas_pct"] * 100).round(4)
    # months are derived, not delivered; 30.4375 d/month is the same divisor the Kaplan-Meier
    # panels use (365.25 / 12), so the table and the figure cannot disagree
    t["os_months"] = (t["os_days"] / 30.4375).round(1)
    t = t.rename(columns={"parvimonas_pct": "Parvimonas (% of community)"})
    cols = list(t.columns)
    cols.insert(cols.index("os_days") + 1, cols.pop(cols.index("os_months")))
    t = t[cols].sort_values(["case", "aliquot"])
    t.to_csv(TAB / "tableS2_tcma_sample_annotation.csv", index=False,
             encoding="utf-8-sig")
    print(f"S2 TCMA sample annotation {t.shape}  "
          f"(with OS: {t.os_months.notna().sum()}, deaths: {int(t.os_event.fillna(0).sum())})")


def s3_clinical_survival() -> None:
    """Emit one file per model block.

    A single stacked table would be 25 columns wide with each block filling only
    about six of them, so four fifths of the cells would be em dashes. Splitting by
    block keeps every table readable at full width.
    """
    blocks = [
        ("tableS3a_clinical_association.csv", "s11_clinical_assoc.csv", {
            "metric": "Microbiome metric", "variable": "Clinical variable",
            "test": "Test", "n": "n", "n_groups": "Groups (n)",
            "statistic": "Statistic", "p": "p", "q_bh": "FDR q",
            "spearman_rho": "Spearman rho", "spearman_p": "Spearman p",
            "group_summary": "Group summary (median)",
        }),
        ("tableS3b_univariable_cox.csv", "s11_survival_univariable.csv", {
            "metric": "Microbiome metric", "HR_per_SD": "HR per SD",
            "lower95": "95% CI lower", "upper95": "95% CI upper", "p": "p",
            "n": "Patients (n)", "events": "Deaths (n)", "concordance": "C-index",
        }),
        ("tableS3c_multivariable_cox.csv", "s11_survival_multivariable.csv", {
            "metric": "Microbiome metric", "term": "Model term", "HR": "HR",
            "lower95": "95% CI lower", "upper95": "95% CI upper", "p": "p",
        }),
        ("tableS3d_likelihood_ratio_test.csv", "s11_survival_lrt.csv", {
            "metric": "Microbiome metric", "loglik_clin": "log-likelihood (clinical only)",
            "loglik_clin_plus_metric": "log-likelihood (clinical + metric)",
            "LRT_chisq": "LRT chi-square", "LRT_p": "LRT p",
            "AIC_clin": "AIC (clinical only)", "AIC_clin_plus_metric": "AIC (clinical + metric)",
            "n": "Patients (n)", "events": "Deaths (n)",
        }),
        ("tableS3e_stage_coding_sensitivity.csv", "s11_survival_stage_sensitivity.csv", {
            "metric": "Microbiome metric", "term": "Model term", "HR": "HR", "p": "p",
        }),
    ]
    for out, src, ren in blocks:
        d = pd.read_csv(TAB / src)
        d = d.rename(columns=ren)
        d.to_csv(TAB / out, index=False, encoding="utf-8-sig")
        print(f"  {out:42} {d.shape}")
    # the stacked version is superseded by the five blocks above
    stale = TAB / "tableS3_clinical_and_survival.csv"
    if stale.exists():
        stale.unlink()
        print("  removed the superseded stacked tableS3_clinical_and_survival.csv")


def s4_species_and_hidden() -> None:
    """Two small tables shipped together as one workbook-style pair."""
    g = pd.read_csv(TAB / "s8b_species_summary.csv").rename(columns={
        "study": "Cohort", "labels": "Species labels resolved",
        "median_parvimonas_pct_of_bacteria": "Parvimonas (% of bacteria, median)",
        "median_resolved_frac_of_genus": "Resolved fraction of genus (median)",
        "frac_samples_with_target": "Samples with target species",
    })
    g.to_csv(TAB / "tableS4a_species_resolution_summary.csv", index=False,
             encoding="utf-8-sig")

    h = pd.read_csv(TAB / "s8c_hidden_family_tests.csv").rename(columns={
        "study": "Cohort", "comparison": "Comparison",
        "family": "Family", "n_pairs": "Pairs (n)",
        "mean_clr_diff": "Mean CLR difference", "median_clr_diff": "Median CLR difference",
        "frac_up_in_tumour": "Fraction up in tumour",
        "mean_rel_tumour": "Mean rel. abundance tumour",
        "mean_rel_normal": "Mean rel. abundance normal",
        "p": "p", "q": "FDR q",
        "family_total_rel_tumour": "Family total rel. abundance (tumour)",
        "frac_of_family_genusless": "Genus-unassigned share of family",
    })
    h = h[["Cohort", "Comparison", "Family", "Pairs (n)",
           "Mean CLR difference", "Median CLR difference", "Fraction up in tumour",
           "Mean rel. abundance tumour", "Mean rel. abundance normal",
           "Family total rel. abundance (tumour)",
           "Genus-unassigned share of family", "p", "FDR q"]]
    h = h.sort_values("FDR q")
    h.to_csv(TAB / "tableS4b_hidden_family_tests.csv", index=False, encoding="utf-8-sig")
    print(f"S4a species resolution    {g.shape}")
    print(f"S4b genus-unassigned fam. {h.shape}  (q<0.05: {(h['FDR q'] < 0.05).sum()})")


def s5_platform() -> None:
    c = pd.read_csv(TAB / "s11_genus_concordance.csv")
    c = c.rename(columns={
        "cohort": "Cohort", "n_tumour": "Tumour samples (n)",
        "n_genera_shared": "Genera shared (n)",
        "spearman_rho_all": "Spearman rho (all shared genera)",
        "p_all": "p (all shared)",
        "n_genera_prevfiltered": "Genera after prevalence filter (n)",
        "spearman_rho_filtered": "Spearman rho (prevalence-filtered)",
        "p_filtered": "p (prevalence-filtered)",
        "top20_jaccard": "Top-20 Jaccard",
        "top20_overlap_n": "Top-20 overlap (n genera)",
    })
    c.to_csv(TAB / "tableS5a_platform_concordance.csv", index=False, encoding="utf-8-sig")

    d = pd.read_csv(TAB / "s11_genus_concordance_disagreement.csv").rename(columns={
        "genus": "Genus", "cohort": "Cohort",
        "tcma_mean_pct": "TCMA mean (%)", "cohort_mean_pct": "16S mean (%)",
        "tcma_rank": "TCMA rank", "cohort_rank": "16S rank",
        "rank_gap": "Rank gap (TCMA - 16S)"})
    d = d[["Cohort", "Genus", "TCMA mean (%)", "16S mean (%)", "TCMA rank",
           "16S rank", "Rank gap (TCMA - 16S)"]]
    d = d.reindex(d["Rank gap (TCMA - 16S)"].abs().sort_values(ascending=False).index)
    d.to_csv(TAB / "tableS5b_platform_disagreement.csv", index=False,
             encoding="utf-8-sig")
    print(f"S5a platform concordance  {c.shape}")
    print(f"S5b discordant genera     {d.shape}")


def main() -> None:
    s1_all_genera()
    s2_tcma_samples()
    s3_clinical_survival()
    s4_species_and_hidden()
    s5_platform()


if __name__ == "__main__":
    main()
