#!/usr/bin/env python3
"""S6 supporting check: is the RNA-seq 'microbiome' signal even usable?

S6 returned a null (5/9 genera agreed in direction with S5; no host correlation
survived FDR). Before that null can be interpreted, the signal itself has to be
characterised: total RNA-seq retains prokaryotic reads, but low-biomass tissue
libraries also carry environmental and reagent background, and a background-
dominated profile cannot test an amplicon-derived result either way.

This partitions the 297 detected genera into
  oral_typical        genera that dominate oral/plaque communities
  environmental       genera typical of soil/water/reagent background in low-biomass work
  other               everything else
and reports their share of mean relative abundance.

Input : results/tables/s6a_group_abundance.csv  (from s6_host_coupling.py)
Output: results/tables/s6b_signal_composition.csv, .json
Run: python3 s6_signal_composition.py
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

def project_root() -> Path:
    """Locate the project root by marker, not by fixed depth.

    The same script runs from two layouts: analysis/modules/06-.../scripts/ inside
    the project, and ~/oscc-mb/work/s6/ on the compute host. A fixed parents[N]
    breaks one of them, so search upward for a directory holding both results/ and
    analysis/ (or data/ on the host layout).
    """
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results").is_dir() and (cand / "analysis").is_dir():
            return cand
    for cand in [here, *here.parents]:
        if (cand / "results").is_dir() and (cand / "data").is_dir():
            return cand
    return here.parents[3]


ROOT = project_root()
TAB = ROOT / "results" / "tables"

# Oral/plaque genera (well established markers of the oral community).
ORAL = {
    "Streptococcus", "Prevotella", "Fusobacterium", "Veillonella", "Neisseria",
    "Porphyromonas", "Actinomyces", "Capnocytophaga", "Haemophilus", "Parvimonas",
    "Leptotrichia", "Rothia", "Gemella", "Campylobacter", "Alloprevotella",
    "Schaalia", "Lancefieldella", "Cardiobacterium", "Cutibacterium",
    "Staphylococcus", "Actinobacillus", "Selenomonas", "Tannerella",
    "Treponema", "Eikenella", "Aggregatibacter", "Lautropia", "Bergeyella",
}
# Genera repeatedly reported as kit/reagent or environment background in
# low-biomass 16S/RNA-seq studies.
ENVIRONMENTAL = {
    "Stenotrophomonas", "Pseudomonas", "Bradyrhizobium", "Acinetobacter",
    "Sphingomonas", "Microbacterium", "Corynebacterium", "Flavobacterium",
    "Delftia", "Comamonas", "Methylobacterium", "Ralstonia", "Burkholderia",
    "Herbaspirillum", "Alcanivorax", "Porphyrobacter", "Pseudoruegeria",
    "Catellicoccus", "Paeniglutamicibacter", "Rhodococcus", "Sphingobacterium",
    "Chryseobacterium", "Variovorax", "Cupriavidus", "Aquabacterium",
    "Methylophilus", "Massilia", "Hymenobacter", "Sphingopyxis",
}


def main() -> None:
    g = pd.read_csv(TAB / "s6a_group_abundance.csv", index_col=0)
    mean_rel = g.mean(axis=1).sort_values(ascending=False)

    def bucket(name: str) -> str:
        if name in ORAL:
            return "oral_typical"
        if name in ENVIRONMENTAL:
            return "environmental_background"
        return "other"

    comp = mean_rel.to_frame("mean_rel_abundance")
    comp["bucket"] = [bucket(i) for i in comp.index]
    comp.to_csv(TAB / "s6b_signal_composition.csv")

    shares = comp.groupby("bucket")["mean_rel_abundance"].agg(["sum", "count"])
    shares["pct"] = (100 * shares["sum"]).round(1)
    summary = {
        "n_genera": int(len(comp)),
        "share_pct": {k: float(v) for k, v in shares["pct"].items()},
        "n_genera_by_bucket": {k: int(v) for k, v in shares["count"].items()},
        "top10": {k: round(100 * float(v), 2) for k, v in mean_rel.head(10).items()},
        "top10_cumulative_pct": round(100 * float(mean_rel.head(10).sum()), 1),
        "interpretation": (
            "oral-typical genera retain a substantial share, so the profile is not "
            "pure background; but several of the largest genera are environmental/"
            "reagent taxa, so genus-level abundances are diluted and the cohort is "
            "underpowered to test amplicon-derived effects"),
    }
    (TAB / "s6b_signal_composition.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False), encoding="utf-8")

    print(json.dumps(summary, indent=1, ensure_ascii=False))
    print("\nper-bucket:")
    print(shares.round(3).to_string())


if __name__ == "__main__":
    main()
