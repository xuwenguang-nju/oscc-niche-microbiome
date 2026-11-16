"""S5: cross-cohort replication of the tumour-vs-adjacent-normal effect.

The question this answers is not "which genera are significant" but "does the
SMALL effect size replicate". PRJNA1415359 produced zero residual genera with a
maximum |Cohen's dz| of 0.674; a single negative cohort cannot distinguish
"there is no tumour-specific signal" from "this cohort is unusual". Three
independent paired cohorts test which of those holds.

What is computed per cohort
---------------------------
For every genus, the per-patient paired log-ratio
    log((tumour + e)/sum) - log((adjacent normal + e)/sum)
is averaged, giving an effect estimate and its standard error. This is the same
estimator as 04-niche-analysis S4a, restricted to the tumour-vs-normal contrast
because that is the only one all four cohorts share (only PRJNA1415359 has
saliva).

Then across cohorts
-------------------
DerSimonian-Laird random-effects pooling with I^2 and Cochran's Q. Random
effects, not fixed: the cohorts differ in country, platform, read length and
depth, so a fixed-effect model would assert a common true effect the design does
not support.

Readouts
--------
1. pooled effect per genus, with direction consistency across cohorts
2. per-genus replication of LARGE effects (|dz| >= 0.8 in >= 2 cohorts)
3. the per-cohort max |dz|, which is the actual replication claim

Usage
-----
python analysis/modules/05-meta-analysis/scripts/cross_cohort_meta.py \
    --studies PRJNA1415359 PRJNA1415268 PRJNA597251 PRJNA666746
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[4]
DADA2 = ROOT / "results" / "dada2"
MANIFEST = ROOT / "results" / "tables" / "oscc_niche_sample_manifest.csv"
TAB = ROOT / "results" / "tables"
OUT = ROOT / "analysis" / "out"

ALPHA = 0.05
DEPTH_MIN = 1000
PSEUDO = 0.5
LARGE_DZ = 0.8

# Control types differ across cohorts and this matters for interpretation:
#   adjacent_normal_tissue / matched_normal / paracancerous -> peri-tumoural normal
#   contralateral_normal_mucosa                             -> distant normal
# The alias list must contain every name in the manifest, otherwise paired_effect()
# returns an empty frame and the cohort is skipped SILENTLY (this is how PRJNA412445
# was lost from the default run between 2026-10-23 and 2026-10-28).
NORMAL_ALIASES = ("adjacent_normal_tissue", "matched_normal", "paracancerous",
                  "contralateral_normal_mucosa")


def bh_q(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.clip(ranked, 0.0, 1.0)
    return q


def load_cohort(study: str):
    """Return (counts, meta) aligned to the manifest, or None if unavailable."""
    path = DADA2 / study / "genus_table_silva.csv"
    if not path.exists():
        print(f"  {study}: no genus table -> skipped")
        return None
    df = pd.read_csv(path, index_col=0)
    df.index = [str(i).replace("_F_filt.fastq.gz", "").replace("_R_filt.fastq.gz", "")
                for i in df.index]
    man = pd.read_csv(MANIFEST)
    man = man[man.study == study].drop_duplicates("run").set_index("run")
    keys = pd.Index([k.split("_")[0] for k in df.index])
    keep = keys.isin(man.index)
    df, keys = df.loc[keep], keys[keep]
    df.index = keys
    df = df[~df.index.duplicated(keep="first")]
    meta = man.reindex(df.index)[["niche", "patient_id"]]
    counts = df.drop(columns=[c for c in ("unassigned",) if c in df.columns])
    counts = counts.loc[counts.sum(axis=1) >= DEPTH_MIN]
    meta = meta.loc[counts.index]
    return counts, meta


def paired_effect(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Per-genus mean paired log-ratio (tumour - adjacent normal) and its SE."""
    tumour = meta[meta.niche == "tumor"]
    normal_niche = next((n for n in NORMAL_ALIASES if (meta.niche == n).any()), None)
    if normal_niche is None or tumour.empty:
        return pd.DataFrame()
    normals = meta[meta.niche == normal_niche]

    L = []
    for pid, trow in tumour.groupby("patient_id"):
        nrow = normals[normals.patient_id == pid]
        if nrow.empty:
            continue
        t = counts.loc[trow.index[0]].astype(float)
        n = counts.loc[nrow.index[0]].astype(float)
        if min(t.sum(), n.sum()) < DEPTH_MIN:
            continue
        pt = (t + PSEUDO) / (t.sum() + PSEUDO * len(t))
        pn = (n + PSEUDO) / (n.sum() + PSEUDO * len(n))
        L.append(np.log(pt.values) - np.log(pn.values))

    if len(L) < 3:
        return pd.DataFrame()
    M = np.vstack(L)
    mu = M.mean(axis=0)
    sd = M.std(axis=0, ddof=1)
    n = M.shape[0]
    se = sd / np.sqrt(n)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = mu / se
    p = 2 * stats.t.sf(np.abs(tstat), df=n - 1)
    p[~np.isfinite(p)] = 1.0
    dz = np.where(sd > 0, mu / sd, 0.0)

    return pd.DataFrame({
        "genus": counts.columns.values,
        "n_patients": n,
        "effect": mu,
        "se": se,
        "cohens_dz": dz,
        "p_value": p,
        "q_value": bh_q(p),
        "prevalence": (counts > 0).mean(axis=0).values,
        "normal_niche_used": normal_niche,
    })


def dersimonian_laird(effects: np.ndarray, variances: np.ndarray) -> dict:
    k = len(effects)
    if k < 2:
        return {"pooled": float(effects[0]) if k else np.nan,
                "se": float(np.sqrt(variances[0])) if k else np.nan,
                "k": k, "Q": np.nan, "I2": np.nan, "tau2": 0.0,
                "p_het": np.nan, "p_pooled": np.nan}
    w = 1.0 / variances
    fixed = float(np.sum(w * effects) / np.sum(w))
    Q = float(np.sum(w * (effects - fixed) ** 2))
    dfree = k - 1
    p_het = float(stats.chi2.sf(Q, dfree))
    C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - dfree) / C) if C > 0 else 0.0
    I2 = max(0.0, (Q - dfree) / Q) * 100 if Q > 0 else 0.0
    wr = 1.0 / (variances + tau2)
    pooled = float(np.sum(wr * effects) / np.sum(wr))
    se = float(np.sqrt(1.0 / np.sum(wr)))
    z = pooled / se if se > 0 else 0.0
    return {"pooled": pooled, "se": se, "k": k, "Q": Q, "I2": I2, "tau2": tau2,
            "p_het": p_het, "p_pooled": float(2 * stats.norm.sf(abs(z)))}


def main() -> None:
    ap = argparse.ArgumentParser()
    # Default = the four cohorts whose normal is peri-tumoural, so that the pooled
    # effect compares like with like. PRJNA412445 uses contralateral (distant)
    # normal and is therefore analysed as its own stratum, together with the
    # prospective cohort: see s5_reciprocal_shift_per_cohort.csv.
    ap.add_argument("--studies", nargs="+", default=[
        "PRJNA1415359", "PRJNA1415268", "PRJNA597251", "PRJNA666746"])
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    outdir = (Path(args.outdir).resolve() if args.outdir else TAB)
    outdir.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    per_cohort: dict[str, pd.DataFrame] = {}
    print("per-cohort paired effects (tumour vs adjacent normal):")
    for s in args.studies:
        got = load_cohort(s)
        if got is None:
            continue
        eff = paired_effect(*got)
        if eff.empty:
            print(f"  {s}: no usable tumour/normal pairs -> skipped")
            continue
        per_cohort[s] = eff
        eff.to_csv(outdir / f"s5_effect_{s}.csv", index=False)
        big = int((eff.cohens_dz.abs() >= LARGE_DZ).sum())
        print(f"  {s}: {eff.n_patients.iloc[0]} pairs x {len(eff)} genera | "
              f"max|dz|={eff.cohens_dz.abs().max():.3f} | |dz|>={LARGE_DZ}: {big} | "
              f"q<{ALPHA}: {int((eff.q_value < ALPHA).sum())}")

    if len(per_cohort) < 2:
        print("\nfewer than two cohorts available -- cannot pool. "
              "Process the remaining cohorts first.")
        return

    all_genera = sorted(set().union(*[set(e.genus) for e in per_cohort.values()]))
    rows, dz_rows = [], []
    for g in all_genera:
        eff, var, cohorts, tstats = [], [], [], []
        for s, e in per_cohort.items():
            r = e[e.genus == g]
            if r.empty:
                continue
            r = r.iloc[0]
            if not np.isfinite(r.se) or r.se <= 0:
                continue
            eff.append(r.effect)
            var.append(r.se ** 2)
            tstats.append(r.cohens_dz)
            cohorts.append(s)
        if not eff:
            continue
        eff_a, var_a = np.array(eff), np.array(var)
        pooled = dersimonian_laird(eff_a, var_a)
        n_pos = int((eff_a > 0).sum())
        n_neg = int((eff_a < 0).sum())
        rows.append({
            "genus": g, "n_cohorts": pooled["k"], "cohorts": ";".join(cohorts),
            "pooled_effect": pooled["pooled"], "pooled_se": pooled["se"],
            "pooled_p": pooled["p_pooled"], "Q": pooled["Q"],
            "I2_pct": pooled["I2"], "p_heterogeneity": pooled["p_het"],
            "n_cohorts_positive": n_pos, "n_cohorts_negative": n_neg,
            "direction_consistent": bool(max(n_pos, n_neg) == pooled["k"]),
            "max_abs_dz": float(np.max(np.abs(tstats))),
            "n_cohorts_large_dz": int(np.sum(np.abs(tstats) >= LARGE_DZ)),
        })
        for s, e_val, t_val in zip(cohorts, eff_a, tstats):
            dz_rows.append({"genus": g, "study": s, "effect": e_val,
                            "cohens_dz": t_val})

    meta = pd.DataFrame(rows)
    if meta.empty:
        print("no genus testable in >=2 cohorts")
        return
    meta["pooled_q_bh"] = bh_q(meta.pooled_p.fillna(1.0).values)
    meta = meta.sort_values("pooled_effect")
    meta.to_csv(outdir / "s5_meta_analysis.csv", index=False)
    pd.DataFrame(dz_rows).to_csv(outdir / "s5_effect_by_cohort.csv", index=False)

    k_max = int(meta.n_cohorts.max())
    summary = {
        "cohorts_used": list(per_cohort),
        "n_genera_pooled": int(len(meta)),
        "n_genera_in_all_cohorts": int((meta.n_cohorts == k_max).sum()),
        "large_effect_replicated": {
            "definition": f"|Cohen's dz| >= {LARGE_DZ} in >= 2 cohorts",
            "n_genera": int((meta.n_cohorts_large_dz >= 2).sum()),
            "genera": meta.loc[meta.n_cohorts_large_dz >= 2, "genus"].tolist(),
        },
        "max_abs_dz_by_cohort": {
            s: round(float(e.cohens_dz.abs().max()), 4) for s, e in per_cohort.items()},
        "n_large_dz_by_cohort": {
            s: int((e.cohens_dz.abs() >= LARGE_DZ).sum()) for s, e in per_cohort.items()},
        "direction_consistent_genera": int(meta.direction_consistent.sum()),
        "pooled_significant_bh": int((meta.pooled_q_bh < ALPHA).sum()),
        "median_I2_pct": round(float(meta.I2_pct.median()), 1),
    }
    (OUT / "s5_replication_summary.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"\n=== S5 replication summary (up to {k_max} cohorts) ===")
    print(f"genera pooled: {len(meta)}  (testable in all {k_max}: "
          f"{summary['n_genera_in_all_cohorts']})")
    print("\nmax |Cohen's dz| per cohort:")
    for s, v in summary["max_abs_dz_by_cohort"].items():
        print(f"  {s:<16} max|dz| = {v:.3f}   |dz|>={LARGE_DZ}: "
              f"{summary['n_large_dz_by_cohort'][s]}")
    print(f"\nlarge effects replicated in >=2 cohorts: "
          f"{summary['large_effect_replicated']['n_genera']} "
          f"{summary['large_effect_replicated']['genera'][:10]}")
    print(f"direction-consistent genera: {summary['direction_consistent_genera']}/{len(meta)}")
    print(f"pooled significant after BH: {summary['pooled_significant_bh']}")
    print(f"median I2 (heterogeneity): {summary['median_I2_pct']}%")
    print(f"\nwrote {outdir.relative_to(ROOT)}/s5_meta_analysis.csv")
    print(f"wrote {OUT.relative_to(ROOT)}/s5_replication_summary.json")


if __name__ == "__main__":
    main()
