#!/usr/bin/env python3
"""Build the source CSVs for the per-rank supplementary table.

Why this exists
---------------
Results 5 reports "1 of 20 phyla, 3 of 30 classes, 7 of 73 orders, 9 of 114 families, 6 of 209 genera
and 2 of 408 species pass FDR". Those six fractions appeared nowhere but in the body text: the only
file carrying them was fw2024_420/rank_level_tests.csv, an 854-row analysis intermediate that a reader
cannot be pointed at and that the submission package does not contain.

This writes the two blocks a reader actually needs:
  (a) one row per taxonomic rank: tests run, tests passing FDR, and the direction of those that pass;
  (b) every one of the 28 tests that pass FDR, with the effect, the fraction of patients agreeing and
      the q value, so the direction quoted for each rank can be checked taxon by taxon.

Rank order follows the hierarchy, not the alphabet, and no number is recomputed here: both blocks are
projections of the existing test output.
"""
from __future__ import annotations

import pathlib

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
TAB = ROOT / "results/tables"
SRC = TAB / "fw2024_420/rank_level_tests.csv"

RANK_ORDER = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def main() -> None:
    rl = pd.read_csv(SRC, encoding="utf-8-sig")
    rl["rank"] = pd.Categorical(rl["rank"], categories=RANK_ORDER, ordered=True)
    rl = rl.sort_values(["rank", "q", "taxon"])

    # (a) per-rank summary
    rows = []
    for rank in RANK_ORDER:
        sub = rl[rl["rank"] == rank]
        sig = sub[sub.q < 0.05]
        rows.append({
            "Rank": rank,
            "Tests": len(sub),
            "FDR < 0.05": len(sig),
            "Lower in tumour": int((sig.direction == "down_in_tumour").sum()),
            "Higher in tumour": int((sig.direction == "up_in_tumour").sum()),
            "Smallest q": f"{sub.q.min():.3g}",
        })
    summary = pd.DataFrame(rows)
    summary.loc[len(summary)] = {
        "Rank": "All ranks",
        "Tests": int(summary["Tests"].sum()),
        "FDR < 0.05": int(summary["FDR < 0.05"].sum()),
        "Lower in tumour": int(summary["Lower in tumour"].sum()),
        "Higher in tumour": int(summary["Higher in tumour"].sum()),
        "Smallest q": f"{rl.q.min():.3g}",
    }
    out_a = TAB / "tableS4_rank_level_summary.csv"
    summary.to_csv(out_a, index=False, encoding="utf-8-sig")

    # (b) the tests that pass, one per row
    sig = rl[rl.q < 0.05].copy()
    out_b = TAB / "tableS4_rank_level_significant.csv"
    pd.DataFrame({
        "Rank": sig["rank"].astype(str).values,
        "Taxon": sig.taxon.values,
        # The column names say "mean" because that is what rank_level_tests.csv carries: the
        # arithmetic mean across the 20 patients. Fig. 5 reports medians of the same quantity, so an
        # unnamed "(%)" column invited the two to be read as the same statistic.
        "Tumour, mean (%)": (sig.mean_T * 100).map(lambda v: f"{v:.3f}").values,
        "Normal, mean (%)": (sig.mean_N * 100).map(lambda v: f"{v:.3f}").values,
        "Difference (points)": (sig.delta * 100).map(lambda v: f"{v:+.2f}").values,
        # frac_up is the share of pairs with tumour > normal, so for a taxon that is lower in
        # tumour the agreeing share is its complement; report the share that matches the direction
        "Patients agreeing": [
            f"{(f if d > 0 else 1 - f) * 100:.0f}%"
            for f, d in zip(sig.frac_up, sig.delta)].__iter__(),
        "Direction": sig.direction.str.replace("_in_tumour", "", regex=False)
                                  .str.replace("down", "lower", regex=False)
                                  .str.replace("up", "higher", regex=False).values,
        "p": sig.p.map(lambda v: f"{v:.3g}").values,
        "q": sig.q.map(lambda v: f"{v:.3g}").values,
    }).to_csv(out_b, index=False, encoding="utf-8-sig")

    print(summary.to_string(index=False))
    print(f"\nwrote {out_a.name} ({len(summary)} rows) and {out_b.name} ({len(sig)} rows)")
    # the six fractions the Results text quotes must come out of block (a)
    quoted = {r["Rank"]: f"{r['FDR < 0.05']} of {r['Tests']}" for _, r in summary.iterrows()}
    print("fractions for Results 5:", quoted)


if __name__ == "__main__":
    main()
