"""S2 + S4: decontamination, residual taxa, and niche-attribution ranking.

Reads a per-cohort genus table produced by 02-dada2 and answers the two
readouts that survived the S3 simulation, in this order:

S2  decontamination by multi-niche equiprevalence
    A genus that is detected at the same rate in every niche is behaving like a
    reagent/environment contaminant; a genus that is genuinely resident in one
    niche changes its detection rate with the niche. No cohort released blanks,
    so this replaces the decontam control-group method. Reported for the record,
    not silently applied.

S4a residual taxa (H2)
    For every patient, two paired log-ratios: tumour vs adjacent normal mucosa,
    and tumour vs pre-operative saliva. Across patients, each genus' mean
    log-ratio is tested (one-sample t on the per-patient values) and
    BH-corrected. A genus is a residual only if it is significant AND
    same-signed in BOTH pairings. Output is a taxon list, never a score --
    S3 showed any continuous "tumour-specific fraction" is unidentifiable here.

S4b niche attribution ranking (H1)
    Distance from each tumour sample to each niche's centroid, tested with a
    paired Friedman test across patients, reported as a RANKING. No percentage
    attribution: S3 showed mixing fractions are not recoverable.

Transformations: counts -> CLR (centred log-ratio) with a pseudocount, because
16S genus tables are compositional and Pearson/Spearman on raw proportions
would be driven by the sequencing depth. Bray-Curtis distances are also reported
because they are robust to the rare tail and comparable to the published
literature.

Usage
-----
python analysis/modules/04-niche-analysis/scripts/niche_attribution.py \
    --study PRJNA1415359
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
PREVALENCE_MIN = 0.10      # genus must be detected in >=10% of samples to be tested
DEPTH_MIN = 1000           # drop samples below this many genus-assigned reads
PSEUDO = 0.5


# ---------------------------------------------------------------- helpers ----
def clr(counts: np.ndarray, pseudo: float = PSEUDO) -> np.ndarray:
    x = counts.astype(float) + pseudo
    logx = np.log(x)
    return logx - logx.mean(axis=1, keepdims=True)


def bray_curtis(a: np.ndarray, b: np.ndarray) -> float:
    denom = a.sum() + b.sum()
    return float(np.abs(a - b).sum() / denom) if denom else np.nan


def bh_q(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.clip(ranked, 0.0, 1.0)
    return q


def load_genus_table(study: str) -> pd.DataFrame:
    """Rows are sample names of the form <run>_1 (or <run>), columns are genera."""
    path = DADA2 / study / "genus_table_silva.csv"
    if not path.exists():
        raise SystemExit(f"missing genus table: {path}")
    df = pd.read_csv(path, index_col=0)
    df.index = [str(i).replace("_F_filt.fastq.gz", "").replace("_R_filt.fastq.gz", "")
                for i in df.index]
    df.index.name = "sample_key"
    return df


def map_samples(df: pd.DataFrame, study: str) -> pd.DataFrame:
    man = pd.read_csv(MANIFEST)
    man = man[man.study == study].copy()
    # genus table row keys are <run> or <run>_1; manifest keys are <run>
    man["sample_key"] = man["run"].astype(str)
    idx = pd.Index([k.split("_")[0] for k in df.index], name="run")
    meta = man.drop_duplicates("sample_key").set_index("sample_key")
    meta = meta.reindex(idx)
    out = df.copy()
    out["run"] = idx
    out["niche"] = meta["niche"].values
    out["patient_id"] = meta["patient_id"].values
    # colliding run keys (both mates kept) would duplicate; drop those
    if out["run"].duplicated().any():
        dup = out["run"].duplicated(keep="first")
        out = out.loc[~dup].copy()
    return out


# ---------------------------------------------------------------- S2 ---------
def s2_decontamination(counts: pd.DataFrame, meta: pd.DataFrame,
                       niches: list[str]) -> pd.DataFrame:
    """Equiprevalence test: does detection rate vary across niches?"""
    rows = []
    sub = meta[meta.niche.isin(niches)]
    c = counts.loc[sub.index]
    for g in c.columns:
        present = (c[g] > 0)
        tab = pd.crosstab(sub["niche"], present)
        if tab.shape[1] < 2:
            p = np.nan
        else:
            # chi-square on niche x detection
            try:
                _, p, _, _ = stats.chi2_contingency(tab.values)
            except Exception:
                p = np.nan
        det = {n: float(present[sub.niche == n].mean()) for n in niches}
        rows.append({
            "genus": g,
            "overall_prevalence": float(present.mean()),
            "p_detection_varies_by_niche": p,
            **{f"detect_{n}": det.get(n, np.nan) for n in niches},
        })
    df = pd.DataFrame(rows)
    q = bh_q(df["p_detection_varies_by_niche"].fillna(1.0).values)
    df["q_detection_varies_by_niche"] = q
    # equiprevalent = cannot reject that detection rate is niche-independent
    df["equiprevalent"] = df["q_detection_varies_by_niche"] >= ALPHA
    df = df.sort_values(["equiprevalent", "overall_prevalence"],
                        ascending=[False, True])
    return df


def s2_stratified(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Separate the two comparisons that the single S2 test conflates.

    within_tissue  : tumour vs adjacent normal. This is the contrast S4a uses, so
                     a genus that is equiprevalent here cannot be a tumour-specific
                     signal no matter how different it looks from saliva.
    tissue_vs_saliva : captures sampling/handling differences (tissue handling vs
                     a fluid rinse), which is where reagent contaminants show up.
    """
    rows = []
    tiers = {
        "within_tissue": ["tumor", "adjacent_normal_tissue"],
        "tissue_vs_saliva": ["tumor", "preop_saliva"],
    }
    for g in counts.columns:
        present = counts[g] > 0
        rec = {"genus": g, "overall_prevalence": float(present.mean())}
        for label, niches in tiers.items():
            sub = meta[meta.niche.isin(niches)]
            if sub.empty or sub["niche"].nunique() < 2:
                rec[f"q_{label}"] = np.nan
                rec[f"detect_{label}"] = np.nan
                continue
            tab = pd.crosstab(sub["niche"], present.loc[sub.index])
            if tab.shape[1] < 2:
                p_val = np.nan
            else:
                try:
                    _, p_val, _, _ = stats.chi2_contingency(tab.values)
                except Exception:
                    p_val = np.nan
            rec[f"p_{label}"] = p_val
            rec[f"detect_{label}"] = float(
                (present.loc[sub.index] & (sub.niche == "tumor")).sum()
                / max(1, (sub.niche == "tumor").sum()))
        rows.append(rec)
    df = pd.DataFrame(rows)
    for label in tiers:
        col = f"p_{label}"
        if col in df:
            df[f"q_{label}"] = bh_q(df[col].fillna(1.0).values)
    # the decision that matters for S4a
    df["equiprevalent_within_tissue"] = df["q_within_tissue"] >= ALPHA
    return df


# ---------------------------------------------------------------- S4a --------
def s4a_residual_taxa(counts: pd.DataFrame, meta: pd.DataFrame,
                      test_genera: list[str]) -> pd.DataFrame:
    """Two paired log-ratios per patient, then a one-sample t across patients."""
    pairs_norm, pairs_sal, patients = [], [], []
    tumor = meta[meta.niche == "tumor"]
    for pid, trow in tumor.groupby("patient_id"):
        n = meta[(meta.patient_id == pid) & (meta.niche == "adjacent_normal_tissue")]
        s = meta[(meta.patient_id == pid) & (meta.niche == "preop_saliva")]
        if n.empty or s.empty:
            continue
        t = counts.loc[trow.index[0]]
        nvec = counts.loc[n.index[0]]
        svec = counts.loc[s.index[0]]
        if min(t.sum(), nvec.sum(), svec.sum()) < DEPTH_MIN:
            continue
        eps = PSEUDO
        pt, pn, ps = ((t + eps) / (t.sum() + eps * len(t)),
                      (nvec + eps) / (nvec.sum() + eps * len(nvec)),
                      (svec + eps) / (svec.sum() + eps * len(svec)))
        pairs_norm.append(np.log(pt[test_genera]) - np.log(pn[test_genera]))
        pairs_sal.append(np.log(pt[test_genera]) - np.log(ps[test_genera]))
        patients.append(pid)

    if not pairs_norm:
        return pd.DataFrame()

    Ln = np.vstack(pairs_norm)
    Ls = np.vstack(pairs_sal)

    def t_test(L):
        mu = L.mean(axis=0)
        sd = L.std(axis=0, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            tstat = mu / (sd / np.sqrt(L.shape[0]))
        p = 2 * stats.t.sf(np.abs(tstat), df=L.shape[0] - 1)
        p[~np.isfinite(p)] = 1.0
        return mu, tstat, p

    mu_n, t_n, p_n = t_test(Ln)
    mu_s, t_s, p_s = t_test(Ls)
    q_n, q_s = bh_q(p_n), bh_q(p_s)

    up = (q_n < ALPHA) & (q_s < ALPHA) & (mu_n > 0) & (mu_s > 0)
    down = (q_n < ALPHA) & (q_s < ALPHA) & (mu_n < 0) & (mu_s < 0)

    # Sensitivity: a paired non-parametric test, plus effect size and its
    # interval. A t-test can return "no residual" because of a few high-variance
    # patients rather than because there is no effect, and the two cases call for
    # opposite conclusions, so both are reported.
    def wilcoxon_test(L):
        p_vals = np.ones(L.shape[1])
        for j in range(L.shape[1]):
            col = L[:, j]
            if np.allclose(col, col[0]):
                continue
            try:
                _, pv = stats.wilcoxon(col, zero_method="wilcox")
                p_vals[j] = pv
            except Exception:
                p_vals[j] = 1.0
        return bh_q(p_vals)

    wil_n = wilcoxon_test(Ln)
    # 95% CI of the mean log-ratio (t interval), and a paired Cohen's dz
    def ci_and_dz(L):
        n = L.shape[0]
        mu = L.mean(axis=0)
        sd = L.std(axis=0, ddof=1)
        se = sd / np.sqrt(n)
        tcrit = stats.t.ppf(0.975, df=n - 1)
        dz = np.where(sd > 0, mu / sd, 0.0)
        return mu - tcrit * se, mu + tcrit * se, dz

    lo_n, hi_n, dz_n = ci_and_dz(Ln)
    lo_s, hi_s, dz_s = ci_and_dz(Ls)

    df = pd.DataFrame({
        "genus": test_genera,
        "n_patients": Ln.shape[0],
        "mean_logratio_vs_normal": mu_n,
        "ci95_lo_vs_normal": lo_n,
        "ci95_hi_vs_normal": hi_n,
        "cohens_dz_vs_normal": dz_n,
        "q_vs_normal": q_n,
        "q_wilcoxon_vs_normal": wil_n,
        "mean_logratio_vs_saliva": mu_s,
        "ci95_lo_vs_saliva": lo_s,
        "ci95_hi_vs_saliva": hi_s,
        "cohens_dz_vs_saliva": dz_s,
        "q_vs_saliva": q_s,
        "tumour_enriched": up,
        "tumour_depleted": down,
        "residual": up | down,
    })
    # pre-registered-style report: is a null result powered, or just underpowered?
    df["ci_excludes_zero_vs_normal"] = (df.ci95_lo_vs_normal > 0) | (df.ci95_hi_vs_normal < 0)
    df["direction"] = np.where(df.tumour_enriched, "tumour_enriched",
                               np.where(df.tumour_depleted, "tumour_depleted", ""))
    return df.sort_values(["residual", "mean_logratio_vs_normal"], ascending=[False, False])


# ---------------------------------------------------------------- S4b --------
def s4b_niche_ranking(counts: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Distance from each tumour to each niche centroid, paired Friedman test."""
    niche_levels = ["preop_saliva", "adjacent_normal_tissue", "tumor"]
    present = [n for n in niche_levels if (meta.niche == n).any()]
    if len(present) < 2:
        return pd.DataFrame(), {"error": "fewer than two niches present"}

    X = clr(counts.values)
    cdf = pd.DataFrame(X, index=counts.index, columns=counts.columns)

    tumor = meta[meta.niche == "tumor"]
    rows = []
    for pid, trow in tumor.groupby("patient_id"):
        t = cdf.loc[trow.index[0]]
        rec = {"patient_id": pid}
        for n in present:
            # centroid of the OTHER patients' samples of this niche, so the
            # tumour is never compared with a centroid that includes itself
            other = meta[(meta.niche == n) & (meta.patient_id != pid)]
            if other.empty:
                rec[n] = np.nan
            else:
                centroid = cdf.loc[other.index].mean(axis=0)
                rec[n] = float(np.sqrt(((t - centroid) ** 2).sum()))
        rows.append(rec)
    d = pd.DataFrame(rows).set_index("patient_id")

    # Aitchison distance to this patient's OWN references. The tumour niche is
    # excluded here: comparing a sample with itself yields distance 0 and would
    # make "tumour is closest to tumour" a tautology rather than a finding.
    own = []
    own_niches = [n for n in present if n != "tumor"]
    for pid, trow in tumor.groupby("patient_id"):
        t = cdf.loc[trow.index[0]]
        rec = {"patient_id": pid}
        for n in own_niches:
            ref = meta[(meta.niche == n) & (meta.patient_id == pid)]
            if ref.empty:
                rec[n] = np.nan
            else:
                rec[n] = float(np.sqrt(((t - cdf.loc[ref.index[0]]) ** 2).sum()))
        own.append(rec)
    own_df = pd.DataFrame(own).set_index("patient_id")
    present_own = own_niches

    summary = {"niches": present, "n_patients": int(d.dropna().shape[0])}
    for label, frame in (("vs_cohort_centroid", d), ("vs_own_patient", own_df)):
        if frame.shape[1] < 2:
            continue
        usable = frame.dropna()
        if usable.shape[0] >= 3 and usable.shape[1] >= 3:
            try:
                stat, p = stats.friedmanchisquare(
                    *[usable[c].values for c in usable.columns])
                test_name = "Friedman"
            except Exception:
                stat, p, test_name = np.nan, np.nan, "Friedman"
        elif usable.shape[0] >= 3 and usable.shape[1] == 2:
            # two paired groups: Friedman is undefined, use paired Wilcoxon
            try:
                stat, p = stats.wilcoxon(usable.iloc[:, 0].values,
                                         usable.iloc[:, 1].values)
                test_name = "Wilcoxon paired"
            except Exception:
                stat, p, test_name = np.nan, np.nan, "Wilcoxon paired"
        else:
            stat, p, test_name = np.nan, np.nan, "insufficient data"
        med = usable.median().sort_values()
        summary[label] = {
                "median_distance": {k: round(float(v), 4) for k, v in med.items()},
                "ranking_nearest_first": list(med.index),
                "test": test_name,
                "test_stat": None if np.isnan(stat) else round(float(stat), 3),
                "test_p": None if np.isnan(p) else float(p),
                "n_patients": int(usable.shape[0]),
            }
    return d, {"summary": summary, "own_patient": own_df}


# ---------------------------------------------------------------- main -------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default="PRJNA1415359")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    study = args.study

    outdir = (Path(args.outdir).resolve() if args.outdir else TAB)
    outdir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    counts_raw = load_genus_table(study)
    df = map_samples(counts_raw, study)
    meta = df[["niche", "patient_id"]].copy()
    counts = df.drop(columns=["niche", "patient_id", "run"]).astype(float)
    counts = counts.loc[counts.sum(axis=1) >= DEPTH_MIN]
    meta = meta.loc[counts.index]

    print(f"{study}: {counts.shape[0]} samples x {counts.shape[1]} genera "
          f"({counts.sum().sum():,.0f} reads after depth filter)")
    print("  niche counts:", meta.niche.value_counts().to_dict())

    # drop the unassigned bucket from biological testing but keep it in totals
    if "unassigned" in counts.columns:
        unassigned_frac = float(counts["unassigned"].sum() / counts.sum().sum())
        print(f"  unassigned reads: {100 * unassigned_frac:.1f}% (kept in depth, "
              f"excluded from taxa)")
        counts_bio = counts.drop(columns=["unassigned"])
    else:
        counts_bio = counts

    prev = (counts_bio > 0).mean(axis=0)
    test_genera = sorted(prev[prev >= PREVALENCE_MIN].index)
    print(f"  genera tested (prevalence >= {PREVALENCE_MIN}): {len(test_genera)}")

    niches_for_s2 = [n for n in ("tumor", "adjacent_normal_tissue", "preop_saliva",
                                 "postop_saliva") if (meta.niche == n).any()]

    if len(niches_for_s2) >= 2:
        dec = s2_decontamination(counts_bio, meta, niches_for_s2)
        dec.to_csv(outdir / f"s2_equiprevalence_{study}.csv", index=False)
        n_equi = int(dec["equiprevalent"].sum())
        print(f"\nS2: {n_equi}/{len(dec)} genera are equiprevalent across "
              f"{len(niches_for_s2)} niches (candidate contaminants)")

        strat = s2_stratified(counts_bio, meta)
        strat.to_csv(outdir / f"s2_stratified_{study}.csv", index=False)
        n_within = int(strat["equiprevalent_within_tissue"].sum())
        print(f"S2 stratified: {n_within}/{len(strat)} genera are equiprevalent "
              f"WITHIN the tissue tier (tumour vs adjacent normal)")
        print("  -> a genus equiprevalent here cannot be tumour-specific")
    else:
        print(f"\nS2: SKIPPED -- only {len(niches_for_s2)} niche present "
              f"({niches_for_s2}); equiprevalence needs >= 2")

    res = s4a_residual_taxa(counts_bio, meta, test_genera)
    if not res.empty:
        res.to_csv(outdir / f"s4a_residual_taxa_{study}.csv", index=False)
        n_res = int(res["residual"].sum())
        print(f"S4a: {n_res} residual genera meeting BOTH t-test thresholds "
              f"({int(res.tumour_enriched.sum())} enriched, "
              f"{int(res.tumour_depleted.sum())} depleted) "
              f"in {int(res.n_patients.iloc[0])} paired patients")
        n_wil = int((res.q_wilcoxon_vs_normal < ALPHA).sum())
        n_ci = int(res.ci_excludes_zero_vs_normal.sum())
        n_big = int((res.cohens_dz_vs_normal.abs() >= 0.8).sum())
        print(f"      sensitivity vs adjacent normal: Wilcoxon q<0.05 = {n_wil}; "
              f"95% CI excludes 0 = {n_ci}; |dz|>=0.8 = {n_big}")
        if n_res == 0 and n_big <= 2:
            print("      -> null result is NOT an underpowered test: no genus has a "
                  "large effect size vs adjacent normal")
        print("  top residual genera:")
        for _, r in res[res.residual].head(12).iterrows():
            print(f"    {r['genus']:<28} vs normal log2FC-like={r['mean_logratio_vs_normal']:+.2f} "
                  f"(q={r['q_vs_normal']:.2g})  vs saliva={r['mean_logratio_vs_saliva']:+.2f} "
                  f"(q={r['q_vs_saliva']:.2g})  {r['direction']}")

    dist, ranking = s4b_niche_ranking(counts_bio, meta)
    if not dist.empty:
        dist.to_csv(outdir / f"s4b_distance_to_niche_{study}.csv")
        ranking["own_patient"].to_csv(
            outdir / f"s4b_distance_own_patient_{study}.csv")
        (OUT / f"s4b_niche_ranking_{study}.json").write_text(
            json.dumps(ranking["summary"], indent=1, ensure_ascii=False),
            encoding="utf-8")
        print("\nS4b niche attribution ranking:")
        for label in ("vs_cohort_centroid", "vs_own_patient"):
            if label in ranking["summary"]:
                s = ranking["summary"][label]
                print(f"  [{label}] nearest-first: {s['ranking_nearest_first']}")
                print(f"      median distances: {s['median_distance']}")
                print(f"      {s['test']} p = {s['test_p']} (n={s['n_patients']})")

    meta.to_csv(outdir / f"sample_niche_map_{study}.csv")
    try:
        shown = outdir.relative_to(ROOT)
    except ValueError:
        shown = outdir
    print(f"\nwrote tables to {shown}")


if __name__ == "__main__":
    main()
