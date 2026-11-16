#!/usr/bin/env python3
"""D: is the TCMA profile quantitatively concordant with amplicon 16S, genus by genus?

Why this exists
---------------
S6 concluded that the tumour-microbiome signal does not reproduce on the TCMA platform, but
that conclusion rested on a DIRECTION test (do the ten S5 genera point the same way). Direction
agreement is a weak claim. Genome Medicine 2023 (PMID 37131219) did this properly: they had
RNA-seq, 16S and FISH on overlapping samples and compared them genus by genus with Spearman,
after prevalence-filtering because zero-inflated rare taxa distort rank correlation.

We do not have matched samples between TCMA and our 16S cohorts, so the strongest available
test is cross-sectional: does the two platforms rank genera by abundance the same way? If the
rank concordance is high, the platforms agree about the community and the S6 discrepancy is
about specific taxa; if low, the platforms are not interchangeable at all and any comparison is
uninterpretable.

Reported per cohort:
  - Spearman rho across intersecting genera, unfiltered and prevalence-filtered
  - Jaccard overlap of the top-20 genera
  - the genera that disagree most, so the failure mode is visible

Run: python3 s11_genus_concordance.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

WORK = Path("/home/v3724/oscc-mb/work")
TCMA = WORK / "clin" / "s6v2_tcma_genus.csv"
MANIFEST = WORK / "s8" / "results" / "tables" / "oscc_niche_sample_manifest.csv"
COHORTS = {
    "PRJNA866676": WORK / "s8/data/genus/PRJNA866676/genus_table.csv",
    "PRJNA412445": WORK / "s8/data/genus/PRJNA412445/genus_table.csv",
    "PRJNA1415359": WORK / "s8/data/genus/PRJNA1415359/genus_table.csv",
}
TUMOUR = {"tumor"}


def relabund(df: pd.DataFrame) -> pd.DataFrame:
    return df.div(df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)


def main() -> None:
    tcma = pd.read_csv(TCMA, index_col=0)
    tcma_r = relabund(tcma)
    tcma_mean = tcma_r.mean(axis=0)
    print(f"TCMA: {tcma.shape[0]} samples x {tcma.shape[1]} genera "
          f"(all primary tumour aliquots)")

    man = pd.read_csv(MANIFEST, encoding="utf-8-sig")

    rows, disc_rows = [], []
    for name, path in COHORTS.items():
        if not path.is_file():
            print(f"[skip] {name}: {path} missing")
            continue
        g = pd.read_csv(path, index_col=0)
        # normalise run ids (SRRxxxx_1_F_filt -> SRRxxxx)
        import re
        g.index = [re.match(r"(SRR\d+|ERR\d+|DRR\d+)", str(i)).group(1)
                   if re.match(r"(SRR\d+|ERR\d+|DRR\d+)", str(i)) else str(i)
                   for i in g.index]
        m = man[man.study == name].drop_duplicates("run").set_index("run")
        tum = [i for i in g.index if i in m.index and m.loc[i, "niche"] in TUMOUR]
        if len(tum) < 5:
            print(f"[skip] {name}: only {len(tum)} tumour samples")
            continue
        gt = relabund(g.loc[tum])
        gt_mean = gt.mean(axis=0)
        print(f"\n{name}: {len(tum)} tumour samples x {g.shape[1]} genera")

        shared = tcma_mean.index.intersection(gt_mean.index)
        shared = [s for s in shared if s != "unassigned" and not s.startswith("unassigned")]
        a = tcma_mean[shared].values
        b = gt_mean[shared].values

        # unfiltered
        rho, p = stats.spearmanr(a, b)
        # prevalence filter (Genome Med lesson): keep taxa seen in >=20% of samples on BOTH sides
        prev_t = (tcma_r[shared] > 0).mean(axis=0).values
        prev_g = (gt[shared] > 0).mean(axis=0).values
        keep = (prev_t >= 0.20) & (prev_g >= 0.20)
        if keep.sum() >= 5:
            rho_f, p_f = stats.spearmanr(a[keep], b[keep])
        else:
            rho_f, p_f = np.nan, np.nan

        top_t = set(tcma_mean[shared].nlargest(20).index)
        top_g = set(gt_mean[shared].nlargest(20).index)
        jac = len(top_t & top_g) / len(top_t | top_g)

        rows.append({
            "cohort": name, "n_tumour": len(tum),
            "n_genera_shared": len(shared),
            "spearman_rho_all": round(float(rho), 3), "p_all": float(p),
            "n_genera_prevfiltered": int(keep.sum()),
            "spearman_rho_filtered": round(float(rho_f), 3) if rho_f == rho_f else np.nan,
            "p_filtered": float(p_f) if p_f == p_f else np.nan,
            "top20_jaccard": round(jac, 3),
            "top20_overlap_n": len(top_t & top_g),
        })

        # rank disagreement among prevalence-filtered taxa
        d = pd.DataFrame({"genus": np.array(shared)[keep],
                          "tcma_mean_pct": a[keep] * 100,
                          "cohort_mean_pct": b[keep] * 100})
        d["tcma_rank"] = d.tcma_mean_pct.rank(ascending=False)
        d["cohort_rank"] = d.cohort_mean_pct.rank(ascending=False)
        d["rank_gap"] = (d.tcma_rank - d.cohort_rank).abs()
        d["cohort"] = name
        disc_rows.append(d.sort_values("rank_gap", ascending=False).head(8))

        print(f"  shared genera: {len(shared)} | rho(all)={rho:.3f} p={p:.2e}")
        print(f"  prevalence-filtered: {int(keep.sum())} genera | rho={rho_f:.3f} "
              f"p={p_f:.2e}" if rho_f == rho_f else "  filtered: too few genera")
        print(f"  top-20 Jaccard: {jac:.3f} ({len(top_t & top_g)}/20 shared)")

    R = pd.DataFrame(rows)
    R.to_csv("s11_genus_concordance.csv", index=False)

    if disc_rows:
        D = pd.concat(disc_rows, ignore_index=True)
        D["tcma_mean_pct"] = D.tcma_mean_pct.round(3)
        D["cohort_mean_pct"] = D.cohort_mean_pct.round(3)
        D.to_csv("s11_genus_concordance_disagreement.csv", index=False)
        print("\n=== genera whose rank disagrees most (prevalence-filtered) ===")
        print(D[["cohort", "genus", "tcma_mean_pct", "cohort_mean_pct",
                 "tcma_rank", "cohort_rank", "rank_gap"]].to_string(index=False))

    print("\n=== summary ===")
    print(R.to_string(index=False))
    print("\nInterpretation guide:")
    print("  rho > 0.6  -> platforms rank the community similarly; S6 discrepancy is taxon-specific")
    print("  rho < 0.3  -> platforms are not interchangeable even at genus level")


if __name__ == "__main__":
    main()
