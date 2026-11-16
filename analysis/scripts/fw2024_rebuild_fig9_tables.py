#!/usr/bin/env python3
"""Rebuild the per-patient and per-rank tables behind Fig. 9 from the vendor taxonomy tables.

Why this exists: these tables (`paired_genus_tests`, `rank_level_tests`, `species_level_tests`,
`fusobacterium_species_per_sample`, `whole_profile_vs_public`, `consortium_ratio_per_patient`) were
originally produced in an ad-hoc session and had no generator in the repository. When the
tumour/normal labels for patient P16 were corrected on 2026-11-12 (see
`analysis/scripts/fw2024_downsampling.py::apply_sample_corrections`), every one of them became
stale. This script regenerates them from the corrected genus/species counts, using the same tests
the figure scripts expect: paired Wilcoxon per taxon, BH within rank or within table.

Input : data/raw/FW2024-420/results/03.taxonomy.zip  (vendor absolute abundances, all ranks)
        results/tables/s5_meta_analysis.csv          (for the whole-profile comparison)
Output: results/tables/fw2024_420/{paired_genus_tests,rank_level_tests,species_level_tests,
        fusobacterium_species_per_sample,whole_profile_vs_public,consortium_ratio_per_patient}.csv

Run: python analysis/scripts/fw2024_rebuild_fig9_tables.py
"""
from __future__ import annotations

import re
import zipfile
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
ZIP = ROOT / "data/raw/FW2024-420/results/03.taxonomy.zip"
TAB = ROOT / "results/tables/fw2024_420"
RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

pid = lambda s: re.sub(r"^[NTY]", "", s)
grp = lambda s: "T" if s[0] in "TY" else "N"


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def read_rank(rank):
    with zipfile.ZipFile(ZIP) as z:
        b = z.read(f"03.taxonomy/abundance/absolute/all.{rank}.absolute.xls")
    for enc in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            d = pd.read_csv(StringIO(b.decode(enc)), sep="\t")
            break
        except UnicodeDecodeError:
            continue
    samples = [c for c in d.columns if re.match(r"^[NTY]\d+$", str(c))]
    # vendor labels carry a rank prefix ("g__Actinomyces"); strip it so the names join to the
    # public meta-analysis and to the genus lists in fw2024_downsampling.py
    label = d[rank].astype(str).str.replace(r"^[a-z]__", "", regex=True)
    m = d.set_index(label)[samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # relabel the corrected sample name so downstream joins work
    m = m.rename(columns={"YP18": "TP18"})
    # swap P16 as documented in fw2024_downsampling.apply_sample_corrections
    a, b_ = "TP16", "NP16"
    if a in m.columns and b_ in m.columns:
        va, vb = m[a].copy(), m[b_].copy()
        m[a], m[b_] = vb, va
    return m


def paired_table(m, keep_all=True):
    """Paired Wilcoxon per taxon, tumour vs normal, on relative abundance."""
    rel = m.div(m.sum(axis=0), axis=1)
    T = [c for c in rel.columns if grp(c) == "T"]
    N = [c for c in rel.columns if grp(c) == "N"]
    tp = {pid(c): c for c in T}
    np_ = {pid(c): c for c in N}
    pids = sorted(set(tp) & set(np_), key=int)
    rows = []
    for taxon, v in rel.iterrows():
        t = np.array([v[tp[p]] for p in pids], float)
        n = np.array([v[np_[p]] for p in pids], float)
        if np.all(t == 0) and np.all(n == 0):
            continue
        d = t - n
        try:
            p = float(stats.wilcoxon(t, n).pvalue) if np.any(d != 0) else 1.0
        except ValueError:
            p = 1.0
        rows.append(dict(taxon=taxon, n_pairs=len(pids), mean_T=float(t.mean()),
                         mean_N=float(n.mean()), delta=float(d.mean()), p=p,
                         frac_up=float((d > 0).mean()),
                         mean_rel_T=float(t.mean()), mean_rel_N=float(n.mean())))
    out = pd.DataFrame(rows)
    out["q"] = bh(out.p.values)
    out["direction"] = np.where(out.delta > 0, "up_in_tumour", "down_in_tumour")
    return out


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)

    # ---- rank-level tests -------------------------------------------------
    rt_parts = []
    for rk in RANKS:
        m = read_rank(rk)
        t = paired_table(m)
        t.insert(0, "rank", rk)
        rt_parts.append(t)
    RT = pd.concat(rt_parts, ignore_index=True)
    RT = RT[["rank", "taxon", "n_pairs", "mean_T", "mean_N", "delta", "frac_up", "p", "q", "direction"]]
    RT.to_csv(TAB / "rank_level_tests.csv", index=False, encoding="utf-8-sig")
    print(f"rank_level_tests.csv      {len(RT):5} rows  ({dict(RT.groupby('rank').size())})")

    # ---- genus-level paired tests ----------------------------------------
    gm = read_rank("Genus")
    G = paired_table(gm)
    G = G[["taxon", "n_pairs", "mean_T", "mean_N", "delta", "frac_up", "p", "q", "direction"]]
    G = G.rename(columns={"taxon": "genus"})
    G.to_csv(TAB / "paired_genus_tests.csv", index=False, encoding="utf-8-sig")
    print(f"paired_genus_tests.csv    {len(G):5} rows   FDR<0.05: {int((G.q < 0.05).sum())}")

    # ---- species-level tests ---------------------------------------------
    sm = read_rank("Species")
    S = paired_table(sm)
    S = S[["taxon", "n_pairs", "mean_T", "mean_N", "delta", "frac_up", "p", "q", "direction"]]
    S.to_csv(TAB / "species_level_tests.csv", index=False, encoding="utf-8-sig")
    print(f"species_level_tests.csv   {len(S):5} rows")

    # ---- Fusobacterium species per sample (fn = nucleatum, fp = periodonticum) ----
    def species_col(needle):
        hit = [t for t in sm.index if needle in str(t).lower()]
        return hit[0] if hit else None

    fn_t = species_col("nucleatum")
    fp_t = species_col("periodonticum")
    unc = [t for t in sm.index if "unclassified" in str(t).lower() and "fusobacter" in str(t).lower()]
    if fn_t and fp_t:
        rel = sm.div(sm.sum(axis=0), axis=1)
        rows = []
        for c in rel.columns:
            rows.append(dict(fn=float(rel.loc[fn_t, c]) * 100, fp=float(rel.loc[fp_t, c]) * 100,
                             pid=pid(c), grp=grp(c)))
        FS = pd.DataFrame(rows)
        FS.to_csv(TAB / "fusobacterium_species_per_sample.csv", index=False, encoding="utf-8-sig")
        print(f"fusobacterium_species_per_sample.csv {len(FS)} rows  "
              f"({fn_t.split(';')[-1]} / {fp_t.split(';')[-1]})")

    # ---- whole profile vs public meta-analysis ---------------------------
    meta = pd.read_csv(ROOT / "results/tables/s5_meta_analysis.csv")
    rel = gm.div(gm.sum(axis=0), axis=1)
    T = [c for c in rel.columns if grp(c) == "T"]
    N = [c for c in rel.columns if grp(c) == "N"]
    tp = {pid(c): c for c in T}
    np_ = {pid(c): c for c in N}
    pids = sorted(set(tp) & set(np_), key=int)
    rows = []
    for genus in rel.index:
        v = rel.loc[genus]
        t = np.array([v[tp[p]] for p in pids], float)
        n = np.array([v[np_[p]] for p in pids], float)
        # keep every genus shared with the meta-analysis, including those absent from one group;
        # the epsilon in log2FC handles the zero case, and Fig. 9J's prevalence axis separates them
        rows.append(dict(
            genus=genus,
            log2FC=float(np.log2((t.mean() + 1e-9) / (n.mean() + 1e-9))),
            prev_T=float((t > 0).mean()),
            prev_N=float((n > 0).mean()),
        ))
    WP = pd.DataFrame(rows)
    m = meta.set_index("genus")[["pooled_effect"]]
    WP = WP.join(m, on="genus", how="inner").dropna(subset=["pooled_effect"])
    WP.to_csv(TAB / "whole_profile_vs_public.csv", index=False, encoding="utf-8-sig")
    print(f"whole_profile_vs_public.csv  {len(WP)} shared genera")

    # ---- original-set consortium ratio per patient -----------------------
    import importlib.util
    spec = importlib.util.spec_from_file_location("down", ROOT / "analysis/scripts/fw2024_downsampling.py")
    down = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(down)
    A = [g for g in down.ANAEROBES_ORIGINAL if g in rel.index]
    E = [g for g in down.AEROBES_ORIGINAL if g in rel.index]
    d = pd.DataFrame({"A": rel.reindex(A).sum(), "E": rel.reindex(E).sum()})
    d["logratio"] = np.log10((d.A + 1e-9) / (d.E + 1e-9))
    d["grp"] = [grp(c) for c in d.index]
    d["pid"] = [pid(c) for c in d.index]
    d = d[["A", "E", "grp", "pid", "logratio"]].sort_values(["pid", "grp"])
    d.to_csv(TAB / "consortium_ratio_per_patient.csv", encoding="utf-8-sig")
    print(f"consortium_ratio_per_patient.csv {len(d)} rows  "
          f"(original sets: {len(A)} ana / {len(E)} aer)")


if __name__ == "__main__":
    main()
