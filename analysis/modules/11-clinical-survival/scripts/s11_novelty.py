#!/usr/bin/env python3
"""A4: how novel are the genera we replicate? (literature-novelty assessment)

Why: one of the retrieved papers (Front Immunol 2023, mesothelioma) made its novelty claim
quantitative -- "107 survival-associated genera, only 1 previously reported for this cancer".
Nine of the 13 full-text papers do some version of this. It is cheap and it converts "we found
X" into "we found X, of which only Y were already known", which is what a reviewer asks.

Method: take the genera that survived our own cross-cohort meta-analysis (S5), then look each
one up in the ORMCKB OSCC curation already in this project
(results/tables/ormckb_oscc_genus_votes.csv and _entries.csv), which counts papers reporting
the genus as up or down in OSCC. Report:
  - whether the genus has any OSCC report
  - how many papers, and in which direction
  - whether the literature direction agrees with ours

Runs locally: both inputs are project files.

Run: python3 s11_novelty.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

def locate_root() -> Path:
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results" / "tables").is_dir() and (cand / "analysis").is_dir():
            return cand
    return here.parents[4]


ROOT = locate_root()
T = ROOT / "results" / "tables"


def main() -> None:
    meta = pd.read_csv(T / "s5_meta_analysis.csv")
    votes = pd.read_csv(T / "ormckb_oscc_genus_votes.csv")
    entries = pd.read_csv(T / "ormckb_oscc_genus_entries.csv")

    # our replicated genera: >=3 cohorts, direction-consistent, FDR<0.05
    ours = meta[(meta.n_cohorts >= 3) & meta.direction_consistent
                & (meta.pooled_q_bh < 0.05)].copy()
    ours["our_direction"] = ours.pooled_effect.map(
        lambda x: "up_in_tumour" if x > 0 else "down_in_tumour")
    print(f"our replicating genera: {len(ours)}")

    v = votes.set_index("Microbe")
    rows = []
    for _, r in ours.iterrows():
        g = r.genus
        if g in v.index:
            vv = v.loc[g]
            n_up = int(vv.n_papers_up)
            n_dn = int(vv.n_papers_down)
            lit = ("up" if n_up > n_dn else ("down" if n_dn > n_up else "conflict"))
            agree = (("up_in_tumour" if lit == "up" else "down_in_tumour") == r.our_direction)
            rows.append({"genus": g, "pooled_effect": round(r.pooled_effect, 3),
                         "pooled_q": r.pooled_q_bh, "our_direction": r.our_direction,
                         "oscc_papers": int(vv.n_papers), "papers_up": n_up,
                         "papers_down": n_dn, "lit_majority": lit,
                         "lit_available": True, "direction_agrees": bool(agree)})
        else:
            rows.append({"genus": g, "pooled_effect": round(r.pooled_effect, 3),
                         "pooled_q": r.pooled_q_bh, "our_direction": r.our_direction,
                         "oscc_papers": 0, "papers_up": 0, "papers_down": 0,
                         "lit_majority": "", "lit_available": False,
                         "direction_agrees": None})

    A = pd.DataFrame(rows).sort_values(["lit_available", "oscc_papers"],
                                       ascending=[True, True])
    A.to_csv(T / "s11_novelty_assessment.csv", index=False)

    n = len(A)
    rep = int(A.lit_available.sum())
    print(f"\n=== novelty: {n - rep}/{n} of our replicating genera have NO prior OSCC report "
          f"in ORMCKB ===")
    print(A.to_string(index=False))

    agree = A[A.lit_available & A.direction_agrees.notna()]
    n_agree = int(agree.direction_agrees.astype(bool).sum())
    print(f"\nwith prior report: {rep} | direction agrees: {n_agree}/{len(agree)}")
    print("  CAVEAT: genera renamed recently (Schaalia split from Actinomyces,")
    print("          Lancefieldella from Atopobium) may be absent from ORMCKB because of")
    print("          nomenclature rather than true novelty. Treat 'no prior report' as")
    print("          'no prior report under this name'.")
    if len(agree):
        print("\n=== where the literature disagrees with us ===")
        bad = agree[agree.direction_agrees.astype(bool) == False]  # noqa: E712
        print(bad[["genus", "our_direction", "lit_majority", "papers_up", "papers_down"]]
              .to_string(index=False) if len(bad) else "  none")

    # context: how well studied is OSCC microbiome overall?
    print(f"\n=== context ===")
    print(f"  ORMCKB OSCC genus-level entries: {len(entries)} rows, "
          f"{entries.Microbe.nunique()} genera, {entries.PMID.nunique()} papers")
    tc = entries.groupby("Microbe").size().sort_values(ascending=False)
    print(f"  most-studied genera: {', '.join(f'{k}({v})' for k, v in tc.head(6).items())}")
    print(f"  median papers per genus: {tc.median():.0f}")


if __name__ == "__main__":
    main()
