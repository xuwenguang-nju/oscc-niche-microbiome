#!/usr/bin/env python3
"""S8b: does the Parvimonas signal come from Parvimonas micra?

Why this matters: S7 showed Parvimonas is the one genus that is enriched in tumour versus
every oral niche. Whether that enrichment is P. micra (the species repeatedly implicated in
colorectal and oral carcinogenesis) or a mixture is a different claim, and only exact-match
species calls can support it.

Resolution is deliberately conservative (see species_only.R): only ASVs equal to a reference
sequence get a species. So every statement below is (a) conditional on the resolved fraction,
which is reported, and (b) accompanied by the unresolved remainder so the reader can see how
much of the genus is still unaccounted for.

Run: python3 s8b_species_micra.py --genus-dir data/genus --out results/tables
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# cohort -> (tumour niche, comparator niche, extra niches of interest)
COHORTS = {
    "PRJNA866676": ("tumor", "paracancerous",
                    ["cancer_surface", "saliva", "tongue_coat",
                     "contralateral_normal_mucosa"]),
    "PRJNA412445": ("tumor", "contralateral_normal_mucosa",
                    ["lesion_side_mucosa", "saliva", "subgingival_plaque"]),
}


def load(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path, index_col=0) if path.is_file() else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus-dir", default="data/genus")
    ap.add_argument("--manifest", default="results/tables/oscc_niche_sample_manifest.csv")
    ap.add_argument("--out", default="results/tables")
    ap.add_argument("--genus", default="Parvimonas")
    ap.add_argument("--species", default="micra")
    ap.add_argument("--target", default=None,
                    help="display label; defaults to '<genus> <species>'")
    args = ap.parse_args()

    # the species table columns hold the bare species epithet (e.g. "micra")
    col = args.species
    label = args.target or f"{args.genus} {args.species}"
    args.target = label
    gd, out = Path(args.genus_dir), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(args.manifest, encoding="utf-8-sig")

    summ, tests, niche_tab = [], [], []
    for study, (tum, comp, others) in COHORTS.items():
        sp = load(gd / study / "species_table_Parvimonas.csv")
        gen = load(gd / study / "genus_table.csv")
        if sp is None or gen is None:
            print(f"[skip] {study}: missing species/genus table")
            continue

        meta = man[man.study == study].drop_duplicates("run").set_index("run")
        common = sp.index.intersection(gen.index).intersection(meta.index)
        if not len(common):
            print(f"[skip] {study}: no sample overlap")
            continue

        sp = sp.loc[common]
        gen_c = gen.loc[common]
        total = gen_c.sum(axis=1)
        gen_col = args.genus
        parv = gen_c[gen_col] if gen_col in gen_c.columns else sp.sum(axis=1)

        # absolute abundance relative to all bacterial reads, and share within the genus
        rel_total = (sp.div(total.replace(0, np.nan), axis=0)).fillna(0.0)
        share_in_genus = (sp.div(parv.replace(0, np.nan), axis=0)).fillna(0.0)

        resolved_cols = [c for c in sp.columns if c != "Parvimonas_unresolved"]
        resolved_reads = sp[resolved_cols].sum(axis=1)
        res_frac = (resolved_reads / parv.replace(0, np.nan)).fillna(0.0)

        summ.append({
            "study": study,
            "n_samples": int(len(common)),
            "labels": ", ".join(sorted(sp.columns)),
            "median_parvimonas_pct_of_bacteria": round(float(
                100 * (parv / total.replace(0, np.nan)).median()), 4),
            "median_resolved_frac_of_genus": round(float(res_frac.median()), 3),
            "frac_samples_with_target": round(float(
                (rel_total.get(col, pd.Series(0, index=sp.index)) > 0).mean()), 3),
        })

        # per-niche medians
        for niche, grp in meta.loc[common].groupby("niche"):
            if col in rel_total.columns:
                niche_tab.append({
                    "study": study, "niche": niche, "n": int(len(grp)),
                    f"{args.target}_pct_of_bacteria": round(float(
                        100 * rel_total.loc[grp.index, col].median()), 4),
                    f"{args.target}_pct_of_parvimonas": round(float(
                        100 * share_in_genus.loc[grp.index, col].median()), 3),
                })

        # paired test inside patient: target abundance tumour vs comparator
        if col in sp.columns:
            pairs = []
            for pid, grp in meta.loc[common].groupby("patient_id"):
                a = grp.index[grp["niche"] == tum]
                b = grp.index[grp["niche"] == comp]
                if len(a) == 1 and len(b) == 1:
                    pairs.append((pid, rel_total.loc[a[0], col],
                                  rel_total.loc[b[0], col]))
            if len(pairs) >= 5:
                a = np.array([p[1] for p in pairs])
                b = np.array([p[2] for p in pairs])
                with np.errstate(divide="ignore"):
                    lr = np.log2((a + 1e-6) / (b + 1e-6))
                try:
                    p_val = stats.wilcoxon(lr).pvalue
                except ValueError:
                    p_val = np.nan
                tests.append({
                    "study": study, "taxon": args.target,
                    "comparison": f"{tum} vs {comp}",
                    "n_pairs": len(pairs),
                    "median_pct_tumour": round(float(100 * np.median(a)), 4),
                    "median_pct_comparator": round(float(100 * np.median(b)), 4),
                    "median_log2fc": round(float(np.median(lr)), 3),
                    "frac_up_in_tumour": round(float((lr > 0).mean()), 3),
                    "p_paired_wilcoxon": p_val,
                })

    if summ:
        S = pd.DataFrame(summ).drop(columns=[c for c in summ[0] if str(c).startswith("n_")])
        S.to_csv(out / "s8b_species_summary.csv", index=False)
        print("=== species-level summary ===")
        print(S.to_string(index=False))
    if niche_tab:
        N = pd.DataFrame(niche_tab)
        N.to_csv(out / "s8b_species_by_niche.csv", index=False)
        print("\n=== P. micra median % of bacterial reads, by niche ===")
        print(N.to_string(index=False))
    if tests:
        T = pd.DataFrame(tests)
        T.to_csv(out / "s8b_species_paired_tests.csv", index=False)
        print("\n=== paired test inside patient ===")
        print(T.to_string(index=False))


if __name__ == "__main__":
    main()
