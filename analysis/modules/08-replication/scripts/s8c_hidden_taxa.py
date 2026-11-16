#!/usr/bin/env python3
"""S8c: are tumour-enriched taxa hiding in the ASVs that never reached genus level?

Motivation: in PRJNA866676, 42% of ASVs receive no genus. Every genus-level analysis in this
project (S5, S7) silently drops them, so "the only tumour-enriched genus is Parvimonas" could
be an artefact of that blind spot. Two numbers have to be separated, and they differ by an
order of magnitude:
    ~42% of ASVs          <- how many taxa are invisible
    ~7% of reads          <- how much of the data that actually is
So the honest test is at the rank where most of these ASVs DO have an assignment: family.

Method: take only ASVs with a curated family but no genus, collapse to family, and test each
family's relative abundance inside the patient (tumour vs that patient's normal tissue) with a
paired Wilcoxon on CLR values, BH-corrected. Families are reported with the fraction of their
reads that came from genus-less ASVs, so a hit can be traced to "unclassified member of a
known family" versus "a well-classified family".

Run: python3 s8c_hidden_taxa.py --genus-dir data/genus --out results/tables
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# Tumour vs the patient's own adjacent mucosa, per cohort.
COMPARISONS = {
    "PRJNA866676": ("tumor", "paracancerous"),
    "PRJNA412445": ("tumor", "contralateral_normal_mucosa"),
}


def clr(mat: np.ndarray, pseudo: float = 1e-6) -> np.ndarray:
    lg = np.log(mat + pseudo)
    return lg - lg.mean(axis=1, keepdims=True)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    r = pv[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    q[ok] = out
    return q


def paired_families(rel: pd.DataFrame, meta: pd.DataFrame, tumour: str,
                    normal: str, min_prev: float = 0.20,
                    min_rel: float = 0.005) -> pd.DataFrame:
    # Two filters, for different reasons: prevalence keeps families that are actually present,
    # and a minimum relative abundance keeps CLR log-ratios from being computed on near-zero
    # values, where the pseudocount dominates and produces spuriously huge differences.
    keep = ((rel > 0).mean(axis=0) >= min_prev) & (rel.median(axis=0) >= min_rel)
    r = rel.loc[:, keep]
    cl = pd.DataFrame(clr(r.values), index=r.index, columns=r.columns)
    m = meta.loc[cl.index]

    rows = []
    for fam in r.columns:
        a, b = [], []
        for pid, grp in m.groupby("patient_id"):
            t = grp.index[grp["niche"] == tumour]
            n = grp.index[grp["niche"] == normal]
            if len(t) == 1 and len(n) == 1:
                a.append(cl.loc[t[0], fam])
                b.append(cl.loc[n[0], fam])
        if len(a) < 5:
            continue
        a, b = np.array(a), np.array(b)
        d = a - b
        try:
            p = stats.wilcoxon(d).pvalue
        except ValueError:
            p = np.nan
        rows.append({
            "family": fam, "n_pairs": len(d),
            "mean_clr_diff": float(d.mean()),
            "median_clr_diff": float(np.median(d)),
            "frac_up_in_tumour": float((d > 0).mean()),
            "mean_rel_tumour": float(r.loc[m.index[m["niche"] == tumour], fam].mean()),
            "mean_rel_normal": float(r.loc[m.index[m["niche"] == normal], fam].mean()),
            "p": p,
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q"] = bh(out["p"].values)
    return out.sort_values("p")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus-dir", default="data/genus")
    ap.add_argument("--manifest", default="results/tables/oscc_niche_sample_manifest.csv")
    ap.add_argument("--out", default="results/tables")
    args = ap.parse_args()

    gd, out = Path(args.genus_dir), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(args.manifest, encoding="utf-8-sig")

    all_res, all_share = [], []
    for study, (tumour, normal) in COMPARISONS.items():
        f_hid = gd / study / "family_table_genus_unassigned.csv"
        f_all = gd / study / "family_table.csv"
        f_gen = gd / study / "genus_table.csv"
        if not all(p.is_file() for p in (f_hid, f_all, f_gen)):
            print(f"[skip] {study}: tables missing")
            continue

        hid = pd.read_csv(f_hid, index_col=0)
        fam = pd.read_csv(f_all, index_col=0)
        gen = pd.read_csv(f_gen, index_col=0)

        meta = man[man.study == study].drop_duplicates("run").set_index("run")
        common = hid.index.intersection(meta.index)
        if not len(common):
            print(f"[skip] {study}: no sample overlap with manifest")
            continue

        # how much of the bacterial signal lives in genus-less ASVs, per sample
        tot = gen.loc[common].sum(axis=1)
        hid_tot = hid.loc[common].sum(axis=1)
        share = (hid_tot / tot.replace(0, np.nan)).dropna()
        all_share.append(pd.DataFrame({
            "study": study, "run": share.index,
            "niche": meta.loc[share.index, "niche"].values,
            "genus_unassigned_read_share": share.values,
        }))
        print(f"{study}: genus-less ASVs carry "
              f"{100*share.median():.1f}% of bacterial reads (median across samples)")

        # relative abundance WITHIN the genus-less fraction
        rel_h = hid.loc[common].div(hid.loc[common].sum(axis=1).replace(0, np.nan), axis=0)
        res_h = paired_families(rel_h.fillna(0.0), meta, tumour, normal)
        if res_h.empty:
            print(f"  no family testable for {study}")
            continue

        # attribute each family to how much of it is genus-less
        rel_all = fam.loc[common].div(fam.loc[common].sum(axis=1).replace(0, np.nan), axis=0)
        res_all = paired_families(rel_all.fillna(0.0), meta, tumour, normal)
        share_map = {}
        if not res_all.empty:
            share_map = dict(zip(res_all.family, res_all.mean_rel_tumour))
        res_h["family_total_rel_tumour"] = res_h.family.map(share_map)
        # Correct share: (reads from genus-less ASVs) / (all reads of that family) in tumour
        # samples. The previous version divided two different normalisations and could
        # exceed 1, which is meaningless as a share.
        t_runs = meta.loc[common].index[meta.loc[common]["niche"] == tumour]
        fam_cols = [f for f in res_h.family if f in hid.columns and f in fam.columns]
        num = hid.loc[t_runs, fam_cols].mean()
        den = fam.loc[t_runs, fam_cols].mean().replace(0, np.nan)
        share = (num / den).rename("frac_of_family_genusless").reset_index()
        share.columns = ["family", "frac_of_family_genusless"]
        res_h = res_h.drop(columns=["frac_of_family_genusless"], errors="ignore")
        res_h = res_h.merge(share, on="family", how="left")
        res_h["study"] = study
        res_h["comparison"] = f"{tumour} vs {normal}"
        all_res.append(res_h)

        print(f"  families tested: {len(res_h)} | q<0.05: {int((res_h.q < 0.05).sum())}")
        top = res_h.head(8)
        print("  top by p (within the genus-less fraction):")
        for _, r in top.iterrows():
            print(f"    {r.family:<28} n={int(r.n_pairs):2d} "
                  f"clrDiff={r.mean_clr_diff:+.2f} q={r.q:.3g} "
                  f"frac_genusless={r.frac_of_family_genusless:.2f}")

    if all_res:
        R = pd.concat(all_res, ignore_index=True)
        R = R.round({"mean_clr_diff": 4, "median_clr_diff": 4, "frac_up_in_tumour": 3,
                     "mean_rel_tumour": 6, "mean_rel_normal": 6, "p": 6, "q": 6,
                     "family_total_rel_tumour": 6, "frac_of_family_genusless": 3})
        R.to_csv(out / "s8c_hidden_family_tests.csv", index=False)
    if all_share:
        S = pd.concat(all_share, ignore_index=True)
        S.to_csv(out / "s8c_genusless_read_share.csv", index=False)
        print("\n=== genus-less read share by cohort x niche (median) ===")
        print(S.groupby(["study", "niche"])["genus_unassigned_read_share"]
              .median().round(4).to_string())


if __name__ == "__main__":
    main()
