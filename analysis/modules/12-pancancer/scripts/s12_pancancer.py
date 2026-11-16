#!/usr/bin/env python3
"""C (pan-cancer context): where does the HNSC/oral microbiome sit among other TCGA tumours?

Why this matters for the paper
------------------------------
Our whole argument is that the OSCC intratumoral community is sourced from the oral cavity.
The sharpest external support for "oral" is a contrast: oral/aerodigestive tumours should look
different from gastrointestinal tumours in exactly the taxa that distinguish mouth from gut.
TCMA v2 lets us test that because it ships HNSC, ESCA, STAD, COAD and READ through one pipeline
(same reference, same decontamination), so project-to-project differences are not pipeline
artefacts the way they would be across separately processed 16S studies.

What is computed
----------------
1. Composition: mean relative abundance of the major genera per project.
2. Ordination: PCoA on Bray-Curtis of relative abundance, coloured by project, to locate HNSC.
3. A single discriminating axis: the ratio of "oral" to "gut" indicator genera, per project.
4. Parvimonas (our key genus from S5/S7): mean abundance and prevalence per project.

The oral/gut indicator sets are declared explicitly below rather than learned from the data,
so the test is not circular: taxa are assigned by prior knowledge of body-site distribution.

Run: python3 s12_pancancer.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

def locate_root() -> Path:
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results" / "tables").is_dir() and (cand / "analysis").is_dir():
            return cand
    return here.parents[4]


ROOT = locate_root()
RAW = ROOT / "data" / "raw" / "tcma"
T = ROOT / "results" / "tables"

# Prior-knowledge indicator sets (not derived from these data).
ORAL_INDICATORS = ["Streptococcus", "Veillonella", "Fusobacterium", "Prevotella",
                   "Neisseria", "Rothia", "Actinomyces", "Haemophilus",
                   "Porphyromonas", "Capnocytophaga", "Parvimonas"]
GUT_INDICATORS = ["Bacteroides", "Faecalibacterium", "Roseburia", "Blautia",
                  "Eubacterium", "Ruminococcus", "Alistipes", "Akkermansia",
                  "Barnesiella", "Odoribacter"]

PROJ_ORDER = ["HNSC", "ESCA", "STAD", "COAD", "READ"]


def pcoa(dist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = dist.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (dist ** 2) @ J
    w, v = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1]
    w, v = w[idx], v[:, idx]
    pos = w > 0
    coords = v[:, pos] * np.sqrt(w[pos])
    return coords, w[pos] / w[pos].sum()


def permanova(D: np.ndarray, labels, nperm: int = 999, seed: int = 0):
    """One-way PERMANOVA on a distance matrix (Anderson 2001)."""
    lab = np.asarray(labels)
    groups = np.unique(lab)
    n = len(lab)

    def ss_within(l):
        s = 0.0
        for gg in groups:
            idx = np.where(l == gg)[0]
            if len(idx) < 2:
                continue
            sub = D[np.ix_(idx, idx)]
            s += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
        return s

    sw = ss_within(lab)
    st = D[np.triu_indices(n, 1)].sum() / n
    F = ((st - sw) / (len(groups) - 1)) / (sw / (n - len(groups)))
    rng = np.random.default_rng(seed)
    cnt = sum(ss_within(rng.permutation(lab)) <= sw for _ in range(nperm))
    return float(F), float((cnt + 1) / (nperm + 1))


def main() -> None:
    g = pd.read_csv(RAW / "pancancer.genus.relabund.txt", sep="\t", index_col=0)
    m = pd.read_csv(RAW / "pancancer.sample_metadata.txt", sep="\t", index_col=0)
    m.index.name = "Sample"

    tot = g.sum(axis=0)
    keep = tot[tot > 0].index
    print(f"dropped {len(tot)-len(keep)} samples with zero total abundance")
    g, m = g[keep], m.loc[keep]
    print(f"samples: {len(keep)} | genera: {g.shape[0]}")
    print(m.project.value_counts().to_string())

    # ---- 1. composition per project -------------------------------------
    top = g.mean(axis=1).nlargest(12).index.tolist()
    comp = g.loc[top].T.groupby(m.project).mean().T.reindex(columns=PROJ_ORDER)
    comp.index.name = "genus"
    comp.to_csv(T / "s12_project_composition.csv")

    # ---- 2. ordination --------------------------------------------------
    bc = squareform(pdist(g.T.values, metric="braycurtis"))
    coords, evr = pcoa(bc)
    O = pd.DataFrame(coords[:, :3], index=g.columns,
                     columns=["PCoA1", "PCoA2", "PCoA3"])
    O["project"] = m.project.values
    O["site"] = m.AnatomicSite.values
    O["sample_type"] = m.sample_type.values
    O.to_csv(T / "s12_pcoa_coordinates.csv")
    print(f"\nPCoA variance explained (%): {np.round(evr[:3]*100, 1)}")

    F, p = permanova(bc, m.project.values)
    n, k = len(m), len(PROJ_ORDER)
    R2 = (F * (k - 1)) / (F * (k - 1) + (n - k))
    print(f"PERMANOVA project: F = {F:.2f}, p = {p:.3f}, R2 = {R2:.3f}")

    # ---- 3. oral vs gut indicator axis ----------------------------------
    or_ = [x for x in ORAL_INDICATORS if x in g.index]
    gu = [x for x in GUT_INDICATORS if x in g.index]
    print(f"\noral indicators found: {len(or_)}/{len(ORAL_INDICATORS)} -> {or_}")
    print(f"gut  indicators found: {len(gu)}/{len(GUT_INDICATORS)} -> {gu}")

    rows, per_sample = [], {}
    for proj in PROJ_ORDER:
        sel = m.index[m.project == proj]
        if len(sel) == 0:
            continue
        sub = g[sel]
        o = sub.loc[or_].sum(axis=0)
        u = sub.loc[gu].sum(axis=0)
        per_sample[proj] = (100 * o - 100 * u).values
        rows.append({
            "project": proj, "n": len(sel),
            "oral_pct": 100 * o.mean(), "gut_pct": 100 * u.mean(),
            "oral_minus_gut_pct": float((100 * o - 100 * u).mean()),
            "oral_over_gut": float((o.sum() + 1e-9) / (u.sum() + 1e-9)),
            "median_n_oral_genera": float((sub.loc[or_] > 0).sum(axis=0).median()),
        })
    A = pd.DataFrame(rows)
    A.to_csv(T / "s12_oral_gut_axis.csv", index=False)
    # per-sample values, so the figure can show distributions rather than only means
    per_sample_df = pd.DataFrame([
        {"sample": s_, "project": m.loc[s_, "project"],
         "site": m.loc[s_, "AnatomicSite"],
         "oral_minus_gut_pct": float(100 * g.loc[or_, s_].sum()
                                     - 100 * g.loc[gu, s_].sum())}
        for s_ in g.columns])
    per_sample_df.to_csv(T / "s12_oral_gut_per_sample.csv", index=False)
    print("\n=== oral vs gut indicator axis ===")
    print(A.round(3).to_string(index=False))

    H, pk = stats.kruskal(*[per_sample[p_] for p_ in per_sample])
    print(f"\nKruskal-Wallis oral-minus-gut across projects: H = {H:.1f}, p = {pk:.3e}")

    # ---- 4. focus genera including Parvimonas ---------------------------
    focus = [x for x in ["Parvimonas", "Fusobacterium", "Streptococcus", "Veillonella",
                         "Prevotella", "Actinomyces", "Porphyromonas", "Bacteroides"]
             if x in g.index]
    F_ = []
    for proj in PROJ_ORDER:
        sel = m.index[m.project == proj]
        if len(sel) == 0:
            continue
        sub = g[sel]
        for gn in focus:
            v = sub.loc[gn]
            F_.append({"project": proj, "genus": gn,
                       "mean_pct": float(100 * v.mean()),
                       "prevalence_pct": float(100 * (v > 0).mean())})
    Fd = pd.DataFrame(F_)
    Fd.to_csv(T / "s12_focus_genera_by_project.csv", index=False)
    print("\n=== focus genera by project (mean %) ===")
    print(Fd.pivot(index="genus", columns="project", values="mean_pct")
          .reindex(columns=PROJ_ORDER).round(3).to_string())

    pv = [100 * g.loc["Parvimonas", m.index[m.project == p_]].values for p_ in PROJ_ORDER]
    H2, p2 = stats.kruskal(*pv)
    print(f"\nKruskal-Wallis Parvimonas across projects: H = {H2:.1f}, p = {p2:.3e}")

    # ---- 5. sensitivity: primary tumour only ----------------------------
    # Normal-tissue samples in a project mean would blur the signal, so the main claim is
    # re-checked on Primary Tumor alone (the tumour-microbiome question is about tumours).
    sens = []
    for stype in ["Primary Tumor", "Solid Tissue Normal"]:
        ms = m[m.sample_type == stype]
        for proj in PROJ_ORDER:
            sel = ms.index[ms.project == proj]
            if len(sel) < 3:
                continue
            sub = g[sel]
            o = sub.loc[or_].sum(axis=0)
            u = sub.loc[gu, :].sum(axis=0)
            sens.append({"sample_type": stype, "project": proj, "n": len(sel),
                         "oral_pct": float(100 * o.mean()),
                         "gut_pct": float(100 * u.mean()),
                         "oral_minus_gut_pct": float((100 * o - 100 * u).mean()),
                         "parvimonas_pct": float(100 * sub.loc["Parvimonas"].mean())})
    Sd = pd.DataFrame(sens)
    Sd.to_csv(T / "s12_sample_type_sensitivity.csv", index=False)
    print("\n=== sensitivity: primary tumour vs solid tissue normal ===")
    print(Sd.round(2).to_string(index=False))

    # One-way ANOVA on the axis score itself. PERMANOVA above tests the whole genus composition;
    # the score is a single derived number, so its between-project share of variance needs its own
    # test, and the manuscript quotes both. Computing it here keeps that claim regenerable.
    axis = pd.read_csv(T / "s12_oral_gut_per_sample.csv", encoding="utf-8-sig")
    groups = [gg["oral_minus_gut_pct"].values for _, gg in axis.groupby("project")]
    Fa, pa = stats.f_oneway(*groups)
    grand = axis["oral_minus_gut_pct"].mean()
    ss_b = sum(len(gg) * (gg.mean() - grand) ** 2 for gg in groups)
    ss_t = ((axis["oral_minus_gut_pct"] - grand) ** 2).sum()
    eta2 = ss_b / ss_t
    print(f"\none-way ANOVA on the axis score by project: F = {Fa:.1f}, p = {pa:.2g}, eta2 = {eta2:.3f}")

    (T / "s12_pancancer_summary.json").write_text(json.dumps({
        "n_samples": int(len(m)), "n_genera": int(g.shape[0]),
        "projects": m.project.value_counts().to_dict(),
        "pcoa_var_explained_pct": [round(float(v) * 100, 2) for v in evr[:3]],
        "permanova_project": {"F": round(F, 2), "p": p, "R2": round(R2, 3)},
        "anova_score_by_project": {"F": round(float(Fa), 1), "p": float(pa),
                                   "eta2": round(float(eta2), 3),
                                   "df_between": len(groups) - 1,
                                   "df_within": len(axis) - len(groups)},
        "kruskal_oral_minus_gut": {"H": round(float(H), 1), "p": float(pk)},
        "kruskal_parvimonas": {"H": round(float(H2), 1), "p": float(p2)},
        "oral_indicators_used": or_, "gut_indicators_used": gu,
    }, indent=1), encoding="utf-8")
    print("\nwrote s12_* tables")


if __name__ == "__main__":
    main()
