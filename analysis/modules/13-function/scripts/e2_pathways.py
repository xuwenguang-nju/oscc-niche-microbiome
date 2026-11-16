#!/usr/bin/env python3
"""E: functional potential of the six OSCC niches (PICRUSt2 / KEGG pathways).

Two questions
-------------
1. QC first: how trustworthy is the prediction for THIS community? PICRUSt2's NSTI measures
   how far each ASV sits from a reference genome. Most tumour-microbiome papers report KEGG
   results without ever showing NSTI, which hides exactly the case where the prediction is
   extrapolation. We report per-ASV and per-sample NSTI and the read mass they carry.

2. Does the FUNCTIONAL profile reproduce the taxonomic result? S7 found, at genus level, that
   the intratumoural community is closest to adjacent normal mucosa, and that tumour surface
   is closer to saliva than deep tumour is (backflow gradient). If the same ordering appears in
   KEGG pathway space, the source-attribution result is not an artefact of 16S taxonomy.
   This is an independent axis of evidence on the SAME samples.

Statistical choices
-------------------
- Compositional data: Bray-Curtis on relative abundance for ordination and distance ranking.
- Paired designs use Wilcoxon signed-rank within patient (tumour vs adjacent normal etc.).
- BH correction within each comparison family (implemented here, not relying on scipy version).
- Distance-to-niche ranking mirrors S7: for each patient, tumour profile vs the mean profile of
  that niche computed from OTHER patients (leave-self-out), so a patient cannot match itself.

Run: python3 e2_pathways.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import braycurtis, pdist, squareform

OUT = Path("/home/v3724/oscc-mb/work/picrust/PRJNA866676/out")
RES = Path("/home/v3724/results/tables")
MANIFEST = Path("/home/v3724/oscc-mb/work/s8/results/tables/"
                "oscc_niche_sample_manifest.csv")
STUDY = "PRJNA866676"


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = n - rank
        val = min(prev, p[idx] * n / i)
        q[idx] = val
        prev = val
    return np.clip(q, 0, 1)


def permanova(D: np.ndarray, labels, nperm: int = 999, seed: int = 0):
    lab = np.asarray(labels)
    groups = np.unique(lab)
    n = len(lab)

    def ssw(l):
        s = 0.0
        for g in groups:
            idx = np.where(l == g)[0]
            if len(idx) < 2:
                continue
            sub = D[np.ix_(idx, idx)]
            s += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
        return s

    sw = ssw(lab)
    st = D[np.triu_indices(n, 1)].sum() / n
    F = ((st - sw) / (len(groups) - 1)) / (sw / (n - len(groups)))
    rng = np.random.default_rng(seed)
    cnt = sum(ssw(rng.permutation(lab)) <= sw for _ in range(nperm))
    R2 = (F * (len(groups) - 1)) / (F * (len(groups) - 1) + (n - len(groups)))
    return float(F), float((cnt + 1) / (nperm + 1)), float(R2)


def main() -> None:
    # ---------------- QC ----------------
    nst_s = pd.read_csv(OUT / "KO_metagenome_out/weighted_nsti.tsv.gz",
                        sep="\t", index_col=0).iloc[:, 0]
    mtab = pd.read_csv(OUT / "combined_marker_predicted_and_nsti.tsv.gz",
                       sep="\t", index_col=0)
    ncol = [c for c in mtab.columns if "NSTI" in c.upper()][0]
    nst_a = mtab[ncol]
    seqtab = pd.read_csv(OUT / "KO_metagenome_out/seqtab_norm.tsv.gz",
                         sep="\t", index_col=0)
    tot = seqtab.sum(axis=1).reindex(nst_a.index).fillna(0)
    share = tot / tot.sum() * 100

    qc = {
        "n_asv_input": 7910, "n_asv_placed": int(len(nst_a)),
        "n_asv_lost_in_placement": 7910 - int(len(nst_a)),
        "pct_asv_lost": round(100 * (7910 - len(nst_a)) / 7910, 2),
        "sample_weighted_nsti": {
            "median": round(float(nst_s.median()), 4),
            "p95": round(float(np.percentile(nst_s, 95)), 4),
            "max": round(float(nst_s.max()), 4),
            "n_gt_0.5": int((nst_s > 0.5).sum())},
        "asv_nsti": {
            "median": round(float(nst_a.median()), 4),
            "max": round(float(nst_a.max()), 4),
            "n_gt_0.5": int((nst_a > 0.5).sum()),
            "n_gt_2.0": int((nst_a > 2.0).sum())},
        "read_mass_pct_nsti_gt_0.5": round(float(share[nst_a > 0.5].sum()), 3),
        "read_mass_pct_nsti_gt_1.0": round(float(share[nst_a > 1.0].sum()), 3),
    }
    print("=== QC: NSTI ===")
    print(json.dumps(qc, indent=1))

    # ---------------- data ----------------
    P = pd.read_csv(OUT / "pathways_out/path_abun_unstrat.tsv.gz", sep="\t", index_col=0)
    print(f"\npathway table: {P.shape[0]} KEGG pathways x {P.shape[1]} samples")
    rel = P.div(P.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)

    man = pd.read_csv(MANIFEST, encoding="utf-8-sig")
    man = man[man.study == STUDY].drop_duplicates("run").set_index("run")
    common = [s for s in rel.columns if s in man.index]
    print(f"samples with niche annotation: {len(common)}/{rel.shape[1]}")
    rel = rel[common]
    meta = man.loc[common]
    niches = sorted(meta.niche.unique())
    print("niches:", niches)

    # ---------------- A: PERMANOVA across niches ----------------
    D = bray = squareform(pdist(rel.T.values, metric="braycurtis"))
    F, p, R2 = permanova(bray, meta.niche.values)
    print("\n=== A. PERMANOVA (pathway, Bray-Curtis, niche) ===")
    print(f"  F={F:.2f}  p={p:.4f}  R2={R2:.3f}")

    # ---------------- B: pathway-space niche ranking ----------------
    mat = rel.T.values                      # samples x pathways
    idx_niche = meta.niche.values
    idx_pat = meta.patient_id.values
    prof = {}


    per_patient = []
    for p_ in pd.unique(idx_pat):
        for n_ in niches:
            sel = (idx_pat == p_) & (idx_niche == n_)
            if sel.sum():
                prof[(p_, n_)] = mat[sel].mean(axis=0)

    for p_ in pd.unique(idx_pat):
        if (p_, "tumor") not in prof:
            continue
        srcs = [n for n in niches if n != "tumor"]
        d = {}
        for n_ in srcs:
            others = [prof[(q, n_)] for q in pd.unique(idx_pat)
                      if q != p_ and (q, n_) in prof]
            if len(others) >= 3:
                d[n_] = float(braycurtis(prof[(p_, "tumor")],
                                         np.mean(np.vstack(others), axis=0)))
        if len(d) < 2:
            continue
        order = sorted(d, key=d.get)
        per_patient.append({"patient_id": p_, "n_sources": len(d),
                            "nearest": order[0], "order": ">".join(order),
                            "closest_distance": round(d[order[0]], 4),
                            **{f"d_{k}": round(v, 4) for k, v in d.items()}})
    R = pd.DataFrame(per_patient)
    R.to_csv(RES / "e2_pathway_niche_ranking.csv", index=False)
    print(f"\n=== B. pathway-space niche ranking (n={len(R)} patients) ===")
    if len(R):
        print("  nearest source:")
        print("   " + R.nearest.value_counts().to_string().replace("\n", "\n   "))
        # binomial-ish test: is adjacent normal nearest more often than chance?
        k = int((R.nearest == "paracancerous").sum())
        n_src = int(R.n_sources.median())
        pbin = stats.binomtest(k, len(R), 1 / n_src).pvalue
        print(f"  adjacent normal nearest in {k}/{len(R)} patients "
              f"(chance 1/{n_src}) binom p={pbin:.4g}")

    # ---------------- C: paired pathway differences ----------------
    families = [("tumor", "paracancerous"), ("tumor", "saliva"),
                ("tumor", "contralateral_normal"), ("tumor", "cancer_surface")]
    tests = []
    for a_, b_ in families:
        pairs = []
        for p_ in pd.unique(idx_pat):
            sa = (idx_pat == p_) & (idx_niche == a_)
            sb = (idx_pat == p_) & (idx_niche == b_)
            if sa.sum() and sb.sum():
                pairs.append((mat[sa].mean(axis=0), mat[sb].mean(axis=0)))
        if len(pairs) < 5:
            print(f"\n=== C. {a_} vs {b_}: only {len(pairs)} pairs, skipped ===")
            continue
        A = np.vstack([x[0] for x in pairs])
        B = np.vstack([x[1] for x in pairs])
        rows = []
        for i, pw in enumerate(rel.index):
            diff = A[:, i] - B[:, i]
            if np.allclose(diff, 0):
                rows.append((pw, 0.0, 1.0))
                continue
            try:
                _, pv = stats.wilcoxon(A[:, i], B[:, i])
            except ValueError:
                pv = 1.0
            rows.append((pw, float(np.median(diff)), float(pv)))
        df = pd.DataFrame(rows, columns=["pathway", "median_diff", "p"])
        df["q"] = bh(df.p.values)
        df["comparison"] = f"{a_}_vs_{b_}"
        df["n_pairs"] = len(pairs)
        sig = df[df.q < 0.05]
        print(f"\n=== C. {a_} vs {b_} (n={len(pairs)} pairs) ===")
        print(f"  q<0.05: {len(sig)}/{len(df)} pathways")
        if len(sig):
            top = sig.reindex(sig.median_diff.abs().sort_values(ascending=False).index).head(10)
            print(top.round(4).to_string(index=False))
        tests.append(df)

    if tests:
        allt = pd.concat(tests, ignore_index=True)
        allt.to_csv(RES / "e2_pathway_paired_tests.csv", index=False)
        print(f"\nwrote e2_pathway_paired_tests.csv ({len(allt)} rows)")

    (RES / "e2_qc_nsti.json").write_text(json.dumps(qc, indent=1), encoding="utf-8")
    (RES / "e2_permanova.json").write_text(json.dumps(
        {"niche_permanova_pathway": {"F": round(F, 2), "p": p, "R2": round(R2, 3),
                                     "n_samples": int(len(common))}}, indent=1),
        encoding="utf-8")
    print("\ndone")


if __name__ == "__main__":
    main()
