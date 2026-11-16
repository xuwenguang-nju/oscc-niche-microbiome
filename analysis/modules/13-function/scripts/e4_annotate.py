#!/usr/bin/env python3
"""E4: annotate pathway IDs with MetaCyc names and assemble the final E tables.

Why: PWY-6902 and friends are unreadable. PICRUSt2 ships the mapping in its own package, so we
join on it rather than inventing labels. Any ID missing from the map is kept with its raw ID and
flagged, never silently renamed.

Also adds a per-comparison effect-size summary, because with n=36 paired samples a Wilcoxon test
will flag thousands of pathways whose median difference is a few tenths of a percent. Reporting
"277/544 significant" without effect size would overstate the result.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("/home/v3724/oscc-mb/work/picrust/PRJNA866676/out")
RES = Path("/home/v3724/results/tables")
INFO = Path("/home/v3724/oscc-mb/envs/picrust2/lib/python3.12/site-packages/"
            "picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz")


def main() -> None:
    tests = pd.read_csv(RES / "e2_pathway_paired_tests.csv")

    # ---- descriptions ----
    info = pd.read_csv(INFO, sep="\t", header=None, names=["pathway", "description"],
                       on_bad_lines="skip")
    info["pathway"] = info["pathway"].str.strip()
    print(f"MetaCyc map: {len(info)} entries")
    merged = tests.merge(info, on="pathway", how="left")
    n_missing = int(merged.description.isna().sum())
    print(f"tests rows without a description: {n_missing}/{len(merged)} "
          f"({100*n_missing/len(merged):.1f}%)")
    merged["description"] = merged.description.fillna(merged.pathway)
    merged.to_csv(RES / "e4_pathway_tests_annotated.csv", index=False)

    # ---- effect-size-aware summary per comparison ----
    rows = []
    for comp, g in merged.groupby("comparison"):
        sig = g[g.q < 0.05]
        # "meaningful" = significant AND |median difference| >= 0.5% of functional potential
        strong = sig[sig.median_diff.abs() >= 0.005]
        rows.append({
            "comparison": comp, "n_pairs": int(g.n_pairs.iloc[0]),
            "n_pathways": len(g), "n_q_lt_0.05": len(sig),
            "pct_q_lt_0.05": round(100 * len(sig) / len(g), 1),
            "n_absdiff_ge_0.005": len(strong),
            "max_abs_median_diff": round(float(g.median_diff.abs().max()), 4),
            "median_abs_median_diff_of_sig": round(
                float(sig.median_diff.abs().median()), 5) if len(sig) else np.nan,
        })
    S = pd.DataFrame(rows)
    S.to_csv(RES / "e4_comparison_summary.csv", index=False)
    print("\n=== effect-size-aware summary ===")
    print(S.to_string(index=False))

    # ---- top named pathways per comparison ----
    lines = []
    for comp in ["tumor_vs_paracancerous", "tumor_vs_saliva", "tumor_vs_cancer_surface"]:
        g = merged[merged.comparison == comp]
        if not len(g):
            continue
        sig = g[g.q < 0.05].reindex(
            g[g.q < 0.05].median_diff.abs().sort_values(ascending=False).index).head(10)
        lines.append(f"\n### {comp} (top 10 by |median diff|)\n")
        for _, r in sig.iterrows():
            lines.append(f"- `{r.pathway}` {r.description}: "
                         f"median diff {r.median_diff:+.4f}, q={r.q:.2g}\n")
    (RES / "e4_top_pathways.md").write_text(
        "# Top differential KEGG/MetaCyc pathways\n" + "".join(lines), encoding="utf-8")
    print("\n=== top pathways by |median diff| ===")
    print("".join(lines))

    # ---- QC + key numbers into one JSON for the report ----
    qc = json.loads((RES / "e2_qc_nsti.json").read_text())
    perm = json.loads((RES / "e2_permanova.json").read_text())
    dist = json.loads((RES / "e3_summary.json").read_text())
    (RES / "e5_key_numbers.json").write_text(json.dumps({
        "qc_nsti": qc, "permanova": perm, "distance_summary": dist,
        "effect_size_summary": S.to_dict("records"),
    }, indent=1), encoding="utf-8")
    print("\nwrote e4_* and e5_key_numbers.json")


if __name__ == "__main__":
    main()
