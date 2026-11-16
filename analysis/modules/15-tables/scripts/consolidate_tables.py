#!/usr/bin/env python3
"""Reduce the supplementary tables from eleven to five.

Decision rule
-------------
* keep a table only if it carries numbers that a reader must look up AND that no
  figure already shows; otherwise the figure stands in for it
* two things cannot be a figure at all - a per-sample annotation and the full
  meta-analysis across every genus - so those stay as tables, shipped as data files
* the five survival blocks collapse into one tidy table, because each block has
  only a handful of rows and they all describe the same three metrics

What changes
------------
S1 + S3(species) -> S3 (volcano figure now carries the selection landscape;
                    s5_meta_analysis.csv ships as a data file)
S2               -> S2 (kept: per-sample annotation, no figure can show it)
S3a-e            -> S4 (one long-format table: metric, model, term, estimate, CI, p, q)
S4a + S5a        -> dropped (fully covered by F4-B/C and F7-A/B)
S4b              -> folded into S4 as three extra rows
S5b              -> dropped (F7-C shows the largest rank gaps; full list ships as data)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
TAB = ROOT / "results" / "tables"


def tidy_survival() -> pd.DataFrame:
    """Collapse the five survival blocks into one long-format table.

    Every block answers the same question about the same three metrics, so a single
    table with a 'model' column reads better than five tables of three rows each.
    """
    rows = []

    for _, r in pd.read_csv(TAB / "s11_survival_univariable.csv").iterrows():
        rows.append({
            "Metric": r.metric.replace("_rel", ""),
            "Model": "univariable Cox",
            "Term": "per SD",
            "Estimate": round(r.HR_per_SD, 3),
            "CI low": round(r.lower95, 3),
            "CI high": round(r.upper95, 3),
            "p": r.p, "FDR q": pd.NA,
            "n": int(r.n), "Note": f"C-index {r.concordance:.3f}",
        })

    for _, r in pd.read_csv(TAB / "s11_survival_multivariable.csv").iterrows():
        term = {"stage_num": "stage (ordinal 1-4)", "age": "age",
                "sexmale": "sex (male vs female)"}.get(r.term, "per SD")
        rows.append({
            "Metric": r.metric.replace("_rel", ""),
            "Model": "multivariable Cox",
            "Term": term,
            "Estimate": round(r.HR, 3),
            "CI low": round(r.lower95, 3),
            "CI high": round(r.upper95, 3),
            "p": r.p, "FDR q": pd.NA,
            "n": pd.NA, "Note": "adjusted for stage, age, sex",
        })

    for _, r in pd.read_csv(TAB / "s11_survival_stage_sensitivity.csv").iterrows():
        term = {"shannon_z": "per SD", "richness_z": "per SD",
                "stage_binlate (III-IV)": "stage (III-IV vs I-II)",
                "age": "age", "sexmale": "sex (male vs female)"}.get(r.term, r.term)
        if not r.term.endswith("_z"):
            continue          # keep only the microbiome terms; covariates are in S4c
        rows.append({
            "Metric": r.metric.replace("_rel", ""),
            "Model": "multivariable Cox, stage binarised",
            "Term": term,
            "Estimate": round(r.HR, 3),
            "CI low": pd.NA, "CI high": pd.NA,
            "p": r.p, "FDR q": pd.NA,
            "n": pd.NA, "Note": "sensitivity to stage coding",
        })

    for _, r in pd.read_csv(TAB / "s11_survival_lrt.csv").iterrows():
        rows.append({
            "Metric": r.metric.replace("_rel", ""),
            "Model": "likelihood-ratio test",
            "Term": f"clinical vs clinical + metric",
            "Estimate": round(r.LRT_chisq, 3),
            "CI low": pd.NA, "CI high": pd.NA,
            "p": r.LRT_p, "FDR q": pd.NA,
            "n": int(r.n),
            "Note": f"AIC {r.AIC_clin:.1f} vs {r.AIC_clin_plus_metric:.1f}",
        })

    for _, r in pd.read_csv(TAB / "s11_clinical_assoc.csv").iterrows():
        rows.append({
            "Metric": r.metric.replace("_rel", ""),
            "Model": "clinical association",
            "Term": f"{r.variable} ({r.test})",
            "Estimate": round(r.statistic, 3) if pd.notna(r.statistic) else pd.NA,
            "CI low": pd.NA, "CI high": pd.NA,
            "p": r.p, "FDR q": r.q_bh,
            # keep the full group summary: truncating here previously cut the last
            # group's median mid-number (e.g. "... 4.0 (n=1) med=") in Table S2
            "n": int(r.n), "Note": str(r.group_summary),
        })

    return pd.DataFrame(rows)


def main() -> None:
    t = tidy_survival().rename(columns={"p": "p", "n": "n"})
    out = TAB / "tableS2_survival_models.csv"
    t.to_csv(out, index=False, encoding="utf-8-sig")
    print(f"S4 survival models  {t.shape}")
    print(t.groupby("Model").size().to_string())

    drop = ["tableS1_all_genera_meta.csv",
            "tableS3a_clinical_association.csv", "tableS3b_univariable_cox.csv",
            "tableS3c_multivariable_cox.csv", "tableS3d_likelihood_ratio_test.csv",
            "tableS3e_stage_coding_sensitivity.csv",
            "tableS4a_species_resolution_summary.csv",
            "tableS4b_hidden_family_tests.csv",
            "tableS5a_platform_concordance.csv",
            "tableS5b_platform_disagreement.csv"]
    for f in drop:
        p = TAB / f
        if p.exists():
            p.unlink()
            print(f"  removed {f}")
    # keep S2 as-is; rename so the numbering follows the new order
    old = TAB / "tableS2_tcma_sample_annotation.csv"
    if old.exists():
        old.rename(TAB / "tableS2_tcma_sample_annotation.csv")
    print("\nremaining table files:")
    for p in sorted(list(TAB.glob("table[0-9]*.csv"))):
        d = pd.read_csv(p)
        print(f"  {p.name:44} {len(d):>4} x {d.shape[1]}")


if __name__ == "__main__":
    main()
