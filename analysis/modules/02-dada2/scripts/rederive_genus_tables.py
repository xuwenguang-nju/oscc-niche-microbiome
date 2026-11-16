"""Re-derive genus tables for every cohort with host reads removed.

Why this exists: the original genus collapse kept eukaryotic (host) amplicons in
the denominator. Measured host read share differs enormously across the cohorts
(PRJNA1415359 0.0%, PRJNA597251 0.0%, PRJNA1415268 2.4%, PRJNA666746 49.0%), so
without removal every bacterial genus in the worst-affected cohort is diluted
about twofold relative to the others, which biases cross-cohort effect sizes.

This applies exactly the rule now in assign_taxonomy.R::collapse_to_genus -- keep
only Kingdom in {Bacteria, Archaea} -- to cohorts whose taxonomy is already
computed, so all cohorts are derived by one implementation.

Inputs per cohort  : results/dada2/<study>/asv_table.csv, taxonomy_silva.csv
Outputs per cohort : genus_table_silva.csv (overwritten), genus_table_summary.csv

Run:  python analysis/modules/02-dada2/scripts/rederive_genus_tables.py
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
RES = ROOT / "results" / "dada2"
PROKARYOTA = {"Bacteria", "Archaea"}


def rederive(study: str) -> dict:
    od = RES / study
    seqtab = pd.read_csv(od / "asv_table.csv", index_col=0)
    tax = pd.read_csv(od / "taxonomy_silva.csv", index_col=0)

    # asv_table.csv orientation has differed between hosts: the R pipeline writes
    # ASVs x samples (sequence-named columns) for some cohorts while others are
    # stored samples x ASVs. Detect it instead of assuming, otherwise the run
    # fails with a confusing "ASV count differs" error.
    n_tax = tax.shape[0]
    if seqtab.shape[1] == n_tax:
        pass                                  # samples x ASVs, already correct
    elif seqtab.shape[0] == n_tax:
        seqtab = seqtab.T                     # ASVs x samples -> transpose
    else:
        raise ValueError(
            f"{study}: neither orientation matches the taxonomy "
            f"(table {seqtab.shape}, taxonomy {n_tax} ASVs)")

    # Align ASV columns positionally to the taxonomy rownames. Both derive from the
    # same seqtab column order, so position is the correct correspondence; the
    # stored column labels may be raw sequences rather than ASV ids.
    seqtab.columns = tax.index

    reads_total = int(seqtab.values.sum())
    kingdom = tax["Kingdom"]
    is_pro = kingdom.isin(PROKARYOTA)
    host_reads = int(seqtab.loc[:, ~is_pro.values].values.sum())

    sub = seqtab.loc[:, is_pro.values]
    genus = tax.loc[is_pro.values, "Genus"].fillna("unassigned").replace("", "unassigned")
    genus_tab = sub.T.groupby(genus.values).sum().T
    genus_tab = genus_tab.loc[:, genus_tab.sum().sort_values(ascending=False).index]

    genus_tab.to_csv(od / "genus_table_silva.csv")

    pro_reads = int(genus_tab.values.sum())
    summary = {
        "study": study,
        "samples": int(seqtab.shape[0]),
        "asvs_total": int(seqtab.shape[1]),
        "asvs_prokaryote": int(is_pro.sum()),
        "reads_total": reads_total,
        "reads_prokaryote": pro_reads,
        "host_reads_pct": round(100 * host_reads / reads_total, 1),
        "bacterial_no_genus_reads_pct": round(
            100 * int(genus_tab["unassigned"].sum()) / reads_total, 1)
        if "unassigned" in genus_tab else 0.0,
        "n_genera": int(genus_tab.shape[1]),
    }
    # conservation check: nothing invented or lost inside the prokaryote subset
    assert pro_reads + host_reads == reads_total, f"{study}: reads not conserved"
    return summary


def main() -> None:
    studies = [p.name for p in sorted(RES.iterdir())
               if (p / "asv_table.csv").exists() and (p / "taxonomy_silva.csv").exists()]
    print("cohorts:", studies)
    rows = [rederive(s) for s in studies]
    df = pd.DataFrame(rows)
    out = ROOT / "results" / "tables" / "genus_table_summary.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nwrote {out.relative_to(ROOT)}")
    print(json.dumps(rows, indent=1)[:400])


if __name__ == "__main__":
    main()
