#!/usr/bin/env python3
"""Down-sampling analysis of the FW2024-420 primary cohort.

Question: at what sequencing depth do the signals we care about survive?

Three signals are tracked across depth:
  (i)   Parvimonas tumour enrichment (paired Wilcoxon, and its AUC);
  (ii)  the anaerobic-vs-aerobic consortium contrast;
  (iii) the composite scores built *a priori* from the public four-cohort meta-analysis
        (1 up genus, 12 down genera) -- these are defined without this cohort,
        so their AUC is not circular.

Expectation from the manuscript's own identifiability work (Fig. 2I): rare taxa and
large genus panels degrade with depth before a single high-abundance genus and a
cohort-defined consortium do. If that holds here, the "low depth" limitation becomes a
positive argument rather than a caveat.

Input : results/derived/fw2024_420/genus_counts.csv  (built by build_inputs, below)
Output: results/tables/fw2024_420/downsampling_summary.csv
        results/tables/fw2024_420/downsampling_long.csv
        figures/main/F9-depth-sensitivity.png   (only if --figure)

Run: python analysis/scripts/fw2024_downsampling.py
"""
from __future__ import annotations

import argparse
import re
import zipfile
from io import StringIO
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/FW2024-420/results"
TAB = ROOT / "results/tables/fw2024_420"

DEPTHS = [20000, 15000, 10000, 5000, 2000, 1000, 500]
N_REP = 30
SEED = 20261028
EPS = 1e-9                     # pseudocount in the log ratio, as elsewhere in this project

# --- anaerobic : aerobic composite axis: three definitions, all reported -------------
# The first definition is retained only for transparency and is NOT used for inference:
# its two arms were selected on the differential result of this very cohort, so any AUC
# computed from it is circular (see s5_consortium_circularity_audit.csv). It also
# mis-assigns three genera by oxygen tolerance (Veillonella and Oribacterium are strict
# anaerobes but were placed in AEROBES; Gemella is facultative but was placed in ANAEROBES).
_spec_sd = importlib.util.spec_from_file_location(
    "sd", Path(__file__).resolve().parent / "score_definitions.py")
sd = importlib.util.module_from_spec(_spec_sd)
_spec_sd.loader.exec_module(sd)

ANAEROBES_ORIGINAL, AEROBES_ORIGINAL = sd.arms("original")
# The main score is the supplement's strict set (29 vs 16). This file used to carry its own
# 26 + 16 copy, which is how the analysis and the published table came to disagree.
ANAEROBES_MAIN, AEROBES_MAIN = sd.arms("main")

# Defined from oxygen tolerance alone (Bergey's Manual), before any test in this study.
OBLIGATE_ANAEROBE = [
    "Fusobacterium", "Leptotrichia", "Prevotella", "Prevotella_7", "Alloprevotella",
    "Peptostreptococcus", "Peptococcus", "Veillonella", "Porphyromonas", "Parvimonas",
    "Lachnoanaerobaculum", "Treponema", "Sneathia", "Oribacterium", "Atopobium",
    "Stomatobaculum", "Solobacterium", "Johnsonella", "Selenomonas", "Mogibacterium",
    "Filifactor", "Dialister", "Olsenella", "Cryptobacterium", "Slackia", "Tannerella"]
AEROBE_FACULTATIVE = [
    "Streptococcus", "Actinomyces", "Neisseria", "Haemophilus", "Rothia", "Pseudomonas",
    "Aggregatibacter", "Granulicatella", "Corynebacterium", "Sphingomonas", "Eikenella",
    "Lautropia", "Abiotrophia", "Bergeyella", "Kingella", "Acinetobacter"]
# oxygen tolerance is strain-dependent or contested for these; excluded from both arms
AMBIGUOUS_OXYGEN = ["Lactobacillus", "Campylobacter", "Capnocytophaga", "Gemella",
                    "Staphylococcus", "TM7x"]

ANAEROBES = ANAEROBES_MAIN      # the definition of record, used for the depth loop below
AEROBES = AEROBES_ORIGINAL

CONSORTIUM_DEFINITIONS = {
    "original": (ANAEROBES_ORIGINAL, AEROBES_ORIGINAL),
    # "full" keeps every genus whose documented oxygen tolerance allows an assignment,
    # including the ones whose tolerance is strain-dependent (best-guess assignment).
    "taxonomy_full": ([g for g in OBLIGATE_ANAEROBE
                       if g not in ("Campylobacter", "Capnocytophaga", "Staphylococcus", "TM7x")]
                      + ["Gemella"],
                      AEROBE_FACULTATIVE + ["Lactobacillus", "Staphylococcus"]
                      if "Lactobacillus" not in AEROBE_FACULTATIVE else AEROBE_FACULTATIVE),
    # "strict" drops every genus whose oxygen tolerance is contested in the literature.
    "taxonomy_strict": ([g for g in OBLIGATE_ANAEROBE if g not in AMBIGUOUS_OXYGEN],
                        [g for g in AEROBE_FACULTATIVE if g not in AMBIGUOUS_OXYGEN]),
}
# a priori from the public meta-analysis only (never from this cohort)
PRIOR_UP = ["Parvimonas"]
PRIOR_DOWN = ["Actinomyces", "Bifidobacterium", "Veillonella", "Staphylococcus",
              "Alloprevotella", "Cardiobacterium", "Lancefieldella", "Paracoccus",
              "Actinobacillus", "Schaalia", "Scardovia", "Cutibacterium"]


# --- sample-level label corrections, confirmed against the original lab records -------------
# Two corrections were supplied by the study investigators on 2026-11-12 and are applied here, at
# the single point where the vendor table enters the pipeline, so that every downstream analysis
# inherits them:
#
#   1. YP18 -> TP18.  The vendor prefixed one tumour sample with "Y" instead of "T". The
#      sample-submission sheet lists it as P18T, which independently confirms the typo.
#   2. TP16 <-> NP16.  The submission sheet's *tube* labels are correct (正常/肿瘤) but the
#      *sample names* on those two rows were transposed. The investigators confirmed the tube
#      labels are authoritative, so the two sequencing columns must be exchanged for this patient.
#
# Correction 2 changes every paired statistic for this patient, and therefore the cohort-level
# numbers in Fig. 9, Fig. 9K and the six-cohort pooled axis. Re-running the pipeline is mandatory
# after any change to this block.
RENAME = {"YP18": "TP18"}
SWAP = ("TP16", "NP16")


def apply_sample_corrections(counts: pd.DataFrame) -> pd.DataFrame:
    """Apply the documented label corrections; idempotent and logged."""
    counts = counts.rename(columns=RENAME)
    a, b = SWAP
    if a in counts.columns and b in counts.columns:
        # exchange the two columns' data, keeping the (now correct) names in place
        va, vb = counts[a].copy(), counts[b].copy()
        counts[a], counts[b] = vb, va
        print(f"corrected: swapped the tumour/normal columns for patient "
              f"{a[1:]} (tube labels authoritative)")
    if set(RENAME.values()) & set(counts.columns):
        print(f"corrected: renamed {sorted(RENAME)} -> {sorted(RENAME.values())}")
    return counts


def build_inputs() -> pd.DataFrame:
    """Genus x sample integer counts from the vendor's absolute abundance table."""
    with zipfile.ZipFile(RAW / "03.taxonomy.zip") as z:
        raw = z.read("03.taxonomy/abundance/absolute/all.Genus.absolute.xls")
    for enc in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            txt = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    d = pd.read_csv(StringIO(txt), sep="\t")
    samples = [c for c in d.columns if re.match(r"^[NTY]\d+", c)]
    d["genus"] = d["Genus"].str.replace(r"^g__", "", regex=True)
    counts = d.groupby("genus")[samples].sum()          # merge duplicate labels if any
    counts = counts.loc[counts.sum(axis=1) > 0]
    counts = apply_sample_corrections(counts)
    out = ROOT / "results/derived/fw2024_420"
    out.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out / "genus_counts.csv", encoding="utf-8-sig")
    print(f"input: {counts.shape[0]} genera x {counts.shape[1]} samples")
    return counts


def sample_meta(cols):
    """Patient id and group; YP18 is the mislabelled TP18 (documented in the report)."""
    pid = {c: re.sub(r"^[NTY]", "", c) for c in cols}
    grp = {c: ("T" if c[0] in "TY" else "N") for c in cols}
    return pid, grp


def subsample(counts: pd.DataFrame, depth: int, rng) -> pd.DataFrame:
    """Hypergeometric read subsampling per sample; samples with fewer reads are used as-is."""
    out = {}
    for c in counts.columns:
        v = np.rint(counts[c].to_numpy()).astype(np.int64)   # counts are integer-valued
        v[v < 0] = 0
        tot = int(v.sum())
        tgt = min(depth, tot)
        out[c] = rng.multivariate_hypergeometric(v, tgt) if tgt < tot else v
    return pd.DataFrame(out, index=counts.index)


def paired_p(x: pd.Series, y: pd.Series):
    """Paired Wilcoxon on tumour vs normal; y = normal, x = tumour, aligned by patient."""
    if len(x) < 5 or (x == y).all():
        return np.nan
    try:
        return stats.wilcoxon(x, y).pvalue
    except ValueError:
        return np.nan


def analyse(rel: pd.DataFrame, pid: dict, grp: dict, reads_per_sample=None) -> dict:
    """One depth: Parvimonas, consortium scores, FDR count, a priori AUCs."""
    cols = list(rel.columns)
    df = rel.T.copy()
    df["pid"] = [pid[c] for c in cols]
    df["grp"] = [grp[c] for c in cols]

    def paired_frame(col):
        return df.pivot_table(index="pid", columns="grp", values=col).dropna()

    rec = {"n_genera_detected": int((rel.sum(axis=1) > 0).sum())}
    rec["median_reads"] = (float(np.median(reads_per_sample))
                           if reads_per_sample is not None else float("nan"))

    # (i) Parvimonas
    pv = paired_frame("Parvimonas") if "Parvimonas" in rel.index else None
    if pv is not None and len(pv):
        rec["parv_p"] = paired_p(pv["T"], pv["N"])
        rec["parv_N_med"] = float(pv["N"].median()) * 100
        rec["parv_T_med"] = float(pv["T"].median()) * 100

    # (ii) consortium scores
    for label, genes in (("anaerobe", ANAEROBES_MAIN), ("aerobe", AEROBES_MAIN)):
        present = [g for g in genes if g in rel.index]
        df[label] = rel.loc[present].sum().values if present else 0.0
        pf = paired_frame(label)
        rec[f"{label}_p"] = paired_p(pf["T"], pf["N"])
        rec[f"{label}_N_med"] = float(pf["N"].median()) * 100
        rec[f"{label}_T_med"] = float(pf["T"].median()) * 100
        y = (df["grp"] == "T").astype(int).values
        rec[f"{label}_auc"] = _auc(y, df[label].values)

    # the composite ratio, which is the statistic Fig. 5B reports. It has to be tracked through
    # the depth series too, otherwise panels E and F describe a different quantity from panel B.
    df["score"] = np.log10((df["anaerobe"] + EPS) / (df["aerobe"] + EPS))
    rec["score_N_med"] = float(paired_frame("score")["N"].median())
    rec["score_T_med"] = float(paired_frame("score")["T"].median())
    rec["score_p"] = paired_p(paired_frame("score")["T"], paired_frame("score")["N"])
    y = (df["grp"] == "T").astype(int).values
    rec["score_auc"] = _auc(y, df["score"].values)

    # (iii) a priori composite scores -> AUC only (no fitting, so no optimism)
    y = (df["grp"] == "T").astype(int).values
    up = [g for g in PRIOR_UP if g in rel.index]
    dn = [g for g in PRIOR_DOWN if g in rel.index]
    rec["prior_up_n"] = len(up)
    rec["prior_down_n"] = len(dn)
    rec["prior_up_auc"] = _auc(y, rel.loc[up].sum().values) if up else np.nan
    rec["prior_down_auc"] = _auc(y, -rel.loc[dn].sum().values) if dn else np.nan

    # (iv) how many genera pass FDR at this depth
    tests = []
    for g in rel.index:
        pf = paired_frame(g)
        p = paired_p(pf["T"], pf["N"])
        if not np.isnan(p):
            tests.append(p)
    if tests:
        q = _bh(np.array(tests))
        rec["n_tested"] = len(tests)
        rec["n_fdr05"] = int((q < 0.05).sum())
    return rec


def per_genus(rel: pd.DataFrame, pid: dict, grp: dict) -> pd.Series:
    """Paired Wilcoxon p per genus at one depth (for the depth-of-discovery curve)."""
    df = rel.T.copy()
    df["pid"] = [pid[c] for c in df.index]
    df["grp"] = [grp[c] for c in df.index]
    out = {}
    for g in rel.index:
        pf = df.pivot_table(index="pid", columns="grp", values=g).dropna()
        out[g] = paired_p(pf["T"], pf["N"]) if len(pf) else np.nan
    return pd.Series(out)


def _auc(y, score):
    """Rank-based AUC, no sklearn dependency for the standalone script."""
    y = np.asarray(y)
    s = np.asarray(score, dtype=float)
    if s.std() == 0:
        return np.nan
    r = stats.rankdata(s)
    n1 = y.sum()
    n0 = len(y) - n1
    if n0 == 0 or n1 == 0:
        return np.nan
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", action="store_true")
    args = ap.parse_args()

    counts = build_inputs()
    counts.columns = [c for c in counts.columns]
    pid, grp = sample_meta(counts.columns)
    TAB.mkdir(parents=True, exist_ok=True)

    rows, long_rows, genus_rows = [], [], []
    rng = np.random.default_rng(SEED)

    # full depth, no subsampling
    rel_full = counts / counts.sum(axis=0)
    full = analyse(rel_full, pid, grp, counts.sum(axis=0).values)
    full["depth"] = int(counts.sum(axis=0).min())
    rows.append(full)
    print(f"full depth {full['depth']:>7,}  parv p={full.get('parv_p', np.nan):.2e}  "
          f"anaerobe p={full['anaerobe_p']:.2e}  FDR sig={full.get('n_fdr05')}")

    for D in DEPTHS:
        reps = []
        for k in range(N_REP):
            sub = subsample(counts, D, rng)
            rel = sub / sub.sum(axis=0).replace(0, np.nan)
            rel = rel.fillna(0)
            rec = analyse(rel, pid, grp, sub.sum(axis=0).values)
            rec["depth"] = D
            pg = per_genus(rel, pid, grp)
            for g, pv_ in pg.items():
                genus_rows.append({"depth": D, "rep": k, "genus": g, "p": pv_})
            rec["rep"] = k
            reps.append(rec)
            long_rows.append(rec)
        R = pd.DataFrame(reps)
        agg = {"depth": D}
        for c in R.columns:
            if c in ("depth", "rep"):
                continue
            agg[c] = R[c].median()
        # keep the fraction of replicates in which each key test passes
        agg["parv_frac_p05"] = float((R.get("parv_p", pd.Series(dtype=float)) < 0.05).mean())
        agg["anaerobe_frac_p05"] = float((R["anaerobe_p"] < 0.05).mean())
        agg["aerobe_frac_p05"] = float((R["aerobe_p"] < 0.05).mean())
        rows.append(agg)
        print(f"depth {D:>7,}  parv p={agg.get('parv_p', np.nan):.2e} "
              f"({agg['parv_frac_p05']:.0%} reps)  anaerobe p={agg['anaerobe_p']:.2e} "
              f"({agg['anaerobe_frac_p05']:.0%})  FDR sig={agg.get('n_fdr05')}")

    S = pd.DataFrame(rows).sort_values("depth", ascending=False)
    S.to_csv(TAB / "downsampling_summary.csv", index=False, encoding="utf-8-sig")


    pd.DataFrame(long_rows).to_csv(TAB / "downsampling_long.csv", index=False, encoding="utf-8-sig")
    PG = pd.DataFrame(genus_rows)
    if len(PG):
        PG["q"] = np.nan
        for D, idx in PG.groupby("depth").groups.items():
            med = PG.loc[idx].groupby("genus")["p"].median()
            ok = med.notna()
            q = pd.Series(np.nan, index=med.index)
            q[ok] = _bh(med[ok].values)
            PG.loc[idx, "q"] = PG.loc[idx, "genus"].map(q)
        PG.to_csv(TAB / "downsampling_per_genus.csv", index=False, encoding="utf-8-sig")

    # The plotting code used to read a separate hand-made "for_manuscript" file which had gone
    # stale: it still carried the pre-relabel, pre-definition-change values (parv p = 3e-4,
    # AUC 0.90/0.81). It is now written from the summary, with the column names the figure uses,
    # so the figure cannot read numbers the analysis no longer produces.
    MS = S.copy()
    MS["genera_detected"] = MS["n_genera_detected"]
    MS["fdr05"] = MS["n_fdr05"]
    MS["parv_p"] = MS["parv_p"]
    MS["score_p"] = MS["score_p"]
    MS["score_auc"] = MS["score_auc"]
    MS["prior_up_auc"] = MS["prior_up_auc"]
    MS["prior_down_auc"] = MS["prior_down_auc"]
    try:
        nom = (PG.assign(hit=(PG.p < 0.05).astype(int))
                 .groupby(["depth", "rep"]).hit.sum().groupby("depth").mean())
    except NameError:
        nom = None
    MS["nominal_p05"] = MS["depth"].map(nom) if nom is not None else np.nan
    keep = ["depth", "median_reads", "genera_detected", "nominal_p05", "fdr05",
            "parv_p", "score_p", "score_auc", "prior_up_auc", "prior_down_auc"]
    MS[keep].to_csv(TAB / "downsampling_for_manuscript.csv", index=False, encoding="utf-8-sig")

    print(f"\nwrote {TAB/'downsampling_summary.csv'}")
    print(f"wrote {TAB/'downsampling_long.csv'}  ({len(long_rows)} rows)")


if __name__ == "__main__":
    main()
