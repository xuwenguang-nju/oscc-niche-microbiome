#!/usr/bin/env python3
"""The oxygen-tolerance score definitions, in one place.

The axis was previously defined three times: a 26 + 16 list inside the downsampling script, a
29 + 17 / 29 + 16 table in the supplement, and an undocumented 8 + 6 set derived from the cohort's
own differential result. The supplement and the analysis therefore disagreed about which genera
count, and the main figure used a set the paper did not publish.

Everything that needs a definition imports it from here. Table S3 is generated from these lists
(`build_score_definitions.py`), so table and analysis cannot drift apart again.

Author's decision: the main score is the 29 obligate anaerobes against the 16 aerobic/facultative
genera, which is the supplement's "taxonomy, strict" set, with Gemella excluded because sources
classify it inconsistently.
"""
from __future__ import annotations

MAIN_ANA = [
    "Alloprevotella", "Atopobium", "Bulleidia", "Cryptobacterium", "Dialister", "Filifactor",
    "Fretibacterium", "Fusobacterium", "Johnsonella", "Lachnoanaerobaculum", "Leptotrichia",
    "Mogibacterium", "Olsenella", "Oribacterium", "Parvimonas", "Peptococcus",
    "Peptostreptococcus", "Porphyromonas", "Prevotella", "Prevotella_7", "Pyramidobacter",
    "Selenomonas", "Slackia", "Sneathia", "Solobacterium", "Stomatobaculum", "Tannerella",
    "Treponema", "Veillonella",
]

MAIN_AER = [
    "Abiotrophia", "Acinetobacter", "Actinomyces", "Aggregatibacter", "Bergeyella",
    "Corynebacterium", "Eikenella", "Granulicatella", "Haemophilus", "Kingella", "Lautropia",
    "Neisseria", "Pseudomonas", "Rothia", "Sphingomonas", "Streptococcus",
]

# sources disagree, so the main score excludes it and the "full" variant counts it as aerobic
EXCLUDED = ["Gemella"]
FULL_AER = MAIN_AER + EXCLUDED

# Transparency only: derived from this cohort's own differential result, so any AUC computed from
# it is circular. Reported as a sensitivity analysis, never as the main estimate.
ORIGINAL_ANA = ["Parvimonas", "Peptostreptococcus", "Peptococcus", "Gemella", "Porphyromonas",
                "Fusobacterium", "Treponema", "Filifactor"]
ORIGINAL_AER = ["Actinomyces", "Neisseria", "Veillonella", "Rothia", "Oribacterium", "Lautropia"]

# eHOMD is not an independent classification: it keeps the genera the annotation happens to cover
EHOMD_ANA = ["Bulleidia", "Parvimonas", "Prevotella", "Pyramidobacter", "Slackia", "Solobacterium",
             "Treponema", "Veillonella"]
EHOMD_AER = ["Abiotrophia", "Actinomyces", "Streptococcus"]

DEFINITIONS = {
    "main": (MAIN_ANA, MAIN_AER),
    "full": (MAIN_ANA, FULL_AER),
    "eHOMD": (EHOMD_ANA, EHOMD_AER),
    "original": (ORIGINAL_ANA, ORIGINAL_AER),
}


def arms(tag: str = "main") -> tuple[list[str], list[str]]:
    """(anaerobic, aerobic/facultative) genus lists for one definition."""
    return DEFINITIONS[tag]
