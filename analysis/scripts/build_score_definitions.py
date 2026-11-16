#!/usr/bin/env python3
"""Build the authoritative score-definition table and recompute the main composite score.

Why this exists
---------------
The oxygen-tolerance axis existed in three inconsistent forms: a hard-coded list in the plotting
code (26 + 16), a table in the supplement (29 + 17 / 29 + 16), and an undocumented 8 + 6 set
derived from the cohort's own differential result. The two that were supposed to be the same set
disagreed, and the one used for the main figure was not the one published in the supplement.

This script makes the supplement's taxonomy sets the definition of record (author's decision:
main score = 29 anaerobic + 16 aerobic/facultative, Gemella excluded), writes every set out
member by member so no set is undocumented, and recomputes the composite score from the
definition rather than from a list embedded in a plotting script.

Definitions written to results/tables/score_definition_members.csv
    main     29 obligate anaerobes vs 16 aerobic/facultative (Gemella excluded)  <- used for inference
    full     30 vs 17, Gemella counted as aerobic/facultative
    eHOMD    the original set restricted to genera the eHOMD annotation covers (8 vs 3)
    original 8 vs 6, chosen from this cohort's own differential result (transparency only)

Score: log10[(sum of relative abundance in the anaerobic arm + 1e-9) /
             (sum in the aerobic arm + 1e-9)], per sample.

Outputs
    results/tables/score_definition_members.csv
    results/tables/tableS3_genus_set_definitions.csv   (regenerated, four sets explicit)
    results/tables/s5_score_inhouse.csv                (per-patient, per definition)
    results/tables/s5_score_per_cohort.csv             (per cohort, per definition)
    results/tables/s5_score_pooled.csv                 (random-effects pooled, per definition)
"""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/tables"
COUNTS = ROOT / "results/derived/fw2024_420/genus_counts.csv"
EPS = 1e-9

_spec = importlib.util.spec_from_file_location("ccm", ROOT / "analysis/modules/05-meta-analysis/"
                                                      "scripts/cross_cohort_meta.py")
ccm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ccm)

pid = lambda s: re.sub(r"^[NTY]", "", s)                      # noqa: E731
grp = lambda s: "T" if s[0] in "TY" else "N"                  # noqa: E731


_spec_defs = importlib.util.spec_from_file_location(
    "sd", ROOT / "analysis/scripts/score_definitions.py")
sd = importlib.util.module_from_spec(_spec_defs)
_spec_defs.loader.exec_module(sd)

MAIN_ANA, MAIN_AER = sd.arms("main")
EXCLUDED = sd.EXCLUDED
ORIGINAL_ANA, ORIGINAL_AER = sd.arms("original")
EHOMD_ANA, EHOMD_AER = sd.arms("eHOMD")


def members() -> pd.DataFrame:
    """One row per genus with its oxygen class and its four set memberships."""
    t3 = pd.read_csv(TAB / "tableS3_genus_set_definitions.csv", encoding="utf-8-sig")
    t3 = t3.rename(columns={c: c.split(" (")[0] for c in t3.columns})
    ehomd_ana, ehomd_aer = EHOMD_ANA, EHOMD_AER
    aer_full = MAIN_AER + EXCLUDED

    def arm(g, ana, aer):
        return "anaerobe" if g in ana else "aerobe" if g in aer else ""

    rows = []
    for g in sorted(t3.Genus):
        rows.append({
            "Genus": g,
            "Oxygen class": t3.set_index("Genus").loc[g, "Oxygen class"],
            "main": arm(g, MAIN_ANA, MAIN_AER),
            "full": arm(g, MAIN_ANA, aer_full),
            "eHOMD": arm(g, ehomd_ana, ehomd_aer),
            "original": arm(g, ORIGINAL_ANA, ORIGINAL_AER),
        })
    return pd.DataFrame(rows)


def inhouse(defs: pd.DataFrame, tag: str) -> tuple[pd.DataFrame, int, int]:
    """Per-patient score in the hospital cohort for one definition."""
    c = pd.read_csv(COUNTS, index_col=0, encoding="utf-8-sig")
    rel = c.div(c.sum(axis=0), axis=1)
    a = [g for g in defs.loc[defs[tag] == "anaerobe", "Genus"] if g in rel.index]
    e = [g for g in defs.loc[defs[tag] == "aerobe", "Genus"] if g in rel.index]
    d = pd.DataFrame({"ana": rel.reindex(a).sum(), "aer": rel.reindex(e).sum()})
    d["logratio"] = np.log10((d.ana + EPS) / (d.aer + EPS))
    d["grp"] = [grp(s) for s in d.index]
    d["pid"] = [pid(s) for s in d.index]
    p = d.pivot_table(index="pid", columns="grp", values=["ana", "aer", "logratio"]).dropna()
    return p, len(a), len(e)


def cohort(study: str, defs: pd.DataFrame, tag: str):
    g = ccm.load_cohort(study)
    if g is None:
        return None
    c, m = g
    eff = ccm.paired_effect(c, m)
    if eff.empty:
        return None
    nn = eff.normal_niche_used.iloc[0]
    R = c.div(c.sum(axis=1), axis=0)
    a = [x for x in defs.loc[defs[tag] == "anaerobe", "Genus"] if x in R.columns]
    e = [x for x in defs.loc[defs[tag] == "aerobe", "Genus"] if x in R.columns]
    d = pd.DataFrame({"ana": R[a].sum(axis=1), "aer": R[e].sum(axis=1)})
    d["logratio"] = np.log10((d.ana + EPS) / (d.aer + EPS))
    d["niche"] = m.reindex(d.index)["niche"].values
    d["pid"] = m.reindex(d.index)["patient_id"].values
    p = d[d.niche.isin(["tumor", nn])].pivot_table(
        index="pid", columns="niche", values="logratio").dropna()
    if len(p) < 5:
        return None
    dif = (p["tumor"] - p[nn]).values
    return dict(n=len(p), effect=float(dif.mean()),
                se=float(dif.std(ddof=1) / np.sqrt(len(dif))),
                p=float(stats.wilcoxon(dif).pvalue))


def sensitivity_table(pooled: pd.DataFrame, per: pd.DataFrame) -> pd.DataFrame:
    """Pooled estimate per definition with the three interval types and leave-one-out range.

    The intervals quoted in Results 5 come from here, so the table is generated rather than typed:
    z from the DerSimonian-Laird standard error, t(k-2) and HKSJ from the same statistics, and the
    leave-one-cohort-out range by refitting with each cohort removed.
    """
    from scipy import stats as st
    rows = []
    for tag, r in pooled.set_index("definition").iterrows():
        T = per[per.definition == tag]
        k = len(T)
        se = r["se"]
        tcrit = st.t.ppf(0.975, k - 2)
        # HKSJ: variance scaled by the mean squared standardised residual
        w = 1 / T.se.values ** 2
        mu = (w * T.effect.values).sum() / w.sum()
        q_hk = ((T.effect.values - mu) ** 2 / T.se.values ** 2).sum() / (k - 1)
        se_hksj = np.sqrt(q_hk) * np.sqrt(1 / w.sum()) if q_hk > 0 else se
        loo = []
        for drop in T.study:
            T2 = T[T.study != drop]
            w2 = 1 / T2.se.values ** 2
            loo.append((w2 * T2.effect.values).sum() / w2.sum())
        rows.append({
            "definition": tag, "k": k, "pooled": r["pooled"], "se": se,
            "I2_pct": r["I2"], "tau2": np.nan,
            "lo_z": r["pooled"] - 1.96 * se, "hi_z": r["pooled"] + 1.96 * se,
            "lo_t": r["pooled"] - tcrit * se, "hi_t": r["pooled"] + tcrit * se,
            "lo_hksj": r["pooled"] - tcrit * se_hksj, "hi_hksj": r["pooled"] + tcrit * se_hksj,
            "loo_pooled_min": min(loo), "loo_pooled_max": max(loo),
            "loo_all_positive": all(v > 0 for v in loo),
            "loo_z_all_exclude0": all(v - 1.96 * se > 0 for v in loo),
            "drop_study_cohort": np.nan, "drop_ci_low": np.nan, "drop_ci_high": np.nan,
        })
    return pd.DataFrame(rows)


def main() -> None:
    defs = members()
    defs.to_csv(TAB / "score_definition_members.csv", index=False, encoding="utf-8-sig")

    # the supplement table now carries all four sets, so no set is undocumented
    t3 = pd.DataFrame({
        "Genus": defs.Genus,
        "Oxygen class": defs["Oxygen class"],
        "Set 1: taxonomy, full (29 vs 17)": defs["full"].replace("", "not assigned"),
        "Set 2: taxonomy, strict (29 vs 16)": defs["main"].replace("", "excluded"),
        "Set 3: eHOMD annotation": defs["eHOMD"].replace("", "not annotated"),
        "Set 4: original, study-cohort derived": defs["original"].replace("", "not assigned"),
    })
    t3.to_csv(TAB / "tableS3_genus_set_definitions.csv", index=False, encoding="utf-8-sig")

    inh, per, pooled = [], [], []
    for tag in ["main", "full", "eHOMD", "original"]:
        p, na, ne = inhouse(defs, tag)
        lr = p[("logratio", "T")] - p[("logratio", "N")]
        y = np.r_[np.zeros(len(p)), np.ones(len(p))]
        s = np.r_[p[("logratio", "N")].values, p[("logratio", "T")].values]
        inh.append(dict(definition=tag, n_ana=na, n_aer=ne, n=len(p),
                        ana_N=p[("ana", "N")].median() * 100,
                        ana_T=p[("ana", "T")].median() * 100,
                        aer_N=p[("aer", "N")].median() * 100,
                        aer_T=p[("aer", "T")].median() * 100,
                        lr_up=int((lr > 0).sum()), p_lr=float(stats.wilcoxon(lr).pvalue),
                        auc=float(roc_auc_score(y, s))))
        rows = []
        for study in ["PRJNA1415359", "PRJNA1415268", "PRJNA597251", "PRJNA666746", "PRJNA412445"]:
            r = cohort(study, defs, tag)
            if r:
                rows.append(dict(definition=tag, study=study, **r))
        rows.append(dict(definition=tag, study="hospital cohort", n=len(lr),
                         effect=float(lr.mean()),
                         se=float(lr.std(ddof=1) / np.sqrt(len(lr))),
                         p=float(stats.wilcoxon(lr).pvalue)))
        per.extend(rows)
        T = pd.DataFrame(rows)
        pl = ccm.dersimonian_laird(T.effect.values, T.se.values ** 2)
        pooled.append(dict(definition=tag, n_ana=na, n_aer=ne, k=len(T),
                           n_pos=int((T.effect > 0).sum()), pooled=pl["pooled"],
                           se=pl["se"], lo=pl["pooled"] - 1.96 * pl["se"],
                           hi=pl["pooled"] + 1.96 * pl["se"], I2=pl["I2"], p=pl["p_pooled"]))
        print(f"{tag:9} arms {na}+{ne}  k={len(T)}  positive={int((T.effect>0).sum())}  "
              f"pooled={pl['pooled']:+.3f} I2={pl['I2']:.0f}%")

    pd.DataFrame(inh).to_csv(TAB / "s5_score_inhouse.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(per).to_csv(TAB / "s5_score_per_cohort.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(pooled).to_csv(TAB / "s5_score_pooled.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(per).pipe(lambda d: d[d.definition == "main"]).to_csv(
        TAB / "s5_reciprocal_shift_per_cohort.csv", index=False, encoding="utf-8-sig")
    sensitivity_table(pd.DataFrame(pooled), pd.DataFrame(per)).to_csv(
        TAB / "s5_score_sensitivity.csv", index=False, encoding="utf-8-sig")
    print("\nwrote score_definition_members.csv, tableS3_genus_set_definitions.csv, "
          "s5_score_{inhouse,per_cohort,pooled}.csv")


if __name__ == "__main__":
    main()
