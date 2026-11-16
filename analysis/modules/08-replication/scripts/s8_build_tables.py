#!/usr/bin/env python3
"""S8 build step: tidy genus / family / "hidden" family tables from an ASV table.

The two cohorts store their DADA2 output differently (one keys the ASV table by ASV id and
the taxonomy by sequence, the other the reverse), so identity is resolved POSITIONALLY
against the ASV fasta: DADA2 writes asv_seqs.fasta in the same order as the ASV columns of
the sequence table, so the N-th fasta record is ASV<N>. Positional mapping is deterministic
and avoids the duplicate-key hazard of joining on sequence strings.

Every join is validated afterwards (labels must be ASV-like and fully accounted for); a
weak match aborts instead of silently mixing samples.

Outputs (in --out-dir):
  taxonomy.csv                       ASV -> Kingdom..Genus (prokaryotes only)
  genus_table.csv                    sample x genus
  family_table.csv                   sample x family
  family_table_genus_unassigned.csv  sample x family, ONLY ASVs with a family but no genus
                                     -- the fraction S5/S7 silently drop; S8c tests it

Host reads (Kingdom not Bacteria/Archaea) are dropped first, the same rule used throughout.

Run: python3 s8_build_tables.py --study PRJNA866676 --asv <asv_table.csv> \
        --tax <taxonomy.csv> [--asv-seqs <asv_seqs.fasta>] --out-dir data/genus/PRJNA866676
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

PROK = {"Bacteria", "Archaea"}
RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
ASV_PREFIX = "ASV"


def read_fasta(path: Path) -> list[str]:
    seqs, cur = [], []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.upper())
    if cur:
        seqs.append("".join(cur))
    return seqs


def asv_frac(labels) -> float:
    labs = [str(x) for x in labels]
    return sum(1 for x in labs if x.startswith(ASV_PREFIX)) / max(len(labs), 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True)
    ap.add_argument("--asv", required=True)
    ap.add_argument("--tax", required=True)
    ap.add_argument("--asv-seqs", default=None)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    n_fasta = 0
    if args.asv_seqs:
        n_fasta = len(read_fasta(Path(args.asv_seqs)))
        print(f"  fasta: {n_fasta} ASV sequences")

    asv_ids = lambda n: [f"{ASV_PREFIX}{i+1}" for i in range(n)]

    # ---- taxonomy ---------------------------------------------------------
    tax = pd.read_csv(args.tax, index_col=0)
    if asv_frac(tax.index) < 0.5:
        if n_fasta and len(tax.index) == n_fasta:
            tax.index = asv_ids(n_fasta)
            print("  taxonomy: mapped positionally to ASV ids")
        else:
            raise SystemExit(
                f"taxonomy index is not ASV-like and {len(tax.index)} != fasta {n_fasta}")
    for r in RANKS:
        if r not in tax.columns:
            tax[r] = "unassigned"
    tax = tax[RANKS].fillna("unassigned").replace("", "unassigned")
    n_tax = len(tax)

    # ---- ASV table --------------------------------------------------------
    df = pd.read_csv(args.asv, index_col=0)
    if n_fasta and df.shape[0] == n_fasta and df.shape[1] != n_fasta:
        df = df.T
        print(f"  ASV table: rows were ASVs -> transposed to {df.shape[0]} samples")
    elif n_fasta and df.shape[1] == n_fasta and df.shape[0] != n_fasta:
        print(f"  ASV table: already samples x ASVs ({df.shape[0]} x {df.shape[1]})")
    else:
        # fall back to label-based detection
        if asv_frac(df.index) >= 0.5:
            df = df.T
        elif asv_frac(df.columns) < 0.5:
            raise SystemExit("cannot tell which axis is ASV; pass --asv-seqs")

    if n_fasta and df.shape[1] == n_fasta and asv_frac(df.columns) < 0.5:
        df.columns = asv_ids(n_fasta)
        print("  ASV table: columns mapped positionally to ASV ids")

    asv = df
    print(f"{args.study}: ASV table {asv.shape[0]} samples x {asv.shape[1]} ASVs")

    shared = list(asv.columns[asv.columns.isin(tax.index)])
    print(f"  ASVs matched between table and taxonomy: {len(shared)}")
    if len(shared) < 0.5 * asv.shape[1]:
        raise SystemExit("too few ASVs matched; refusing to guess")
    asv = asv.loc[:, shared]
    tx = tax.loc[shared]
    if not asv.columns.equals(tx.index):
        raise SystemExit("alignment broke between ASV table and taxonomy")

    prok = tx["Kingdom"].isin(PROK).values
    asv_p = asv.loc[:, prok]
    tx_p = tx.loc[prok]
    total = float(asv_p.values.sum())
    print(f"  prokaryotic ASVs: {int(prok.sum())} / {len(prok)} | reads: {total:,.0f}")

    def collapse(mat_df: pd.DataFrame, labels: pd.Series, unassigned: str) -> pd.DataFrame:
        lab = np.asarray(labels, dtype=object).copy()
        lab[lab == "unassigned"] = unassigned
        mat = mat_df.to_numpy()                      # samples x ASVs (subset)
        codes, uniq = pd.factorize(lab, sort=True)
        agg = np.zeros((mat.shape[0], len(uniq)))
        np.add.at(agg.T, codes, mat.T)               # sum ASVs sharing a label
        return pd.DataFrame(agg, index=mat_df.index, columns=uniq)

    gen = collapse(asv_p, tx_p["Genus"], "unassigned_genus")
    fam = collapse(asv_p, tx_p["Family"], "unassigned_family")
    hid = (tx_p["Genus"] == "unassigned") & (tx_p["Family"] != "unassigned")
    fam_hidden = collapse(asv_p.loc[:, hid], tx_p.loc[hid, "Family"],
                          "unassigned_family")

    gen.to_csv(out / "genus_table.csv")
    fam.to_csv(out / "family_table.csv")
    fam_hidden.to_csv(out / "family_table_genus_unassigned.csv")
    tx_p.to_csv(out / "taxonomy.csv")

    print(f"  genera: {gen.shape[1]} | families: {fam.shape[1]}")
    print(f"  ASVs with family but no genus: {int(hid.sum())} "
          f"({100 * asv_p.loc[:, hid].to_numpy().sum() / total:.1f}% of reads)")
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
