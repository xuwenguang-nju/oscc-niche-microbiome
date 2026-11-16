#!/usr/bin/env python3
"""Apply the panel-reference corrections to results/reports/figure-captions-jom.md.

Every replacement below fixes a caption sentence that does not describe the panel the plotting
script actually draws. Ground truth for each was read from the script that renders the panel
(`analysis/scripts/audit_caption_panels.py` prints the letters and titles) and from the result
table behind the panel. Each `old` string must occur exactly once, so a silent no-op is
impossible; if any assertion fails nothing is written.
"""
import pathlib

P = pathlib.Path("results/reports/figure-captions-jom.md")
text = P.read_text(encoding="utf-8")
orig = text

EDITS = [
    # ---------------------------------------------------------------- Figure 1
    # D is a per-cohort settings table (platform, primer trim, median reads), not a per-sample
    # depth plot; E counts niches per patient for the six public cohorts only; F is a horizontal
    # bar chart of samples analysed, not a scatter of niches against cohort size.
    ("""which it is resolved. (**D**) Sequencing depth per sample for the seven datasets. (**E**) Patients by
number of niches sampled; the public cohorts contribute one niche per patient, the hospital sample
two. (**F**) Niches sampled per patient against cohort size, with the two cohorts shared with that
pooled analysis marked.""",
     """which it is resolved. (**D**) Per-cohort technical settings: platform, whether the 5' primer was
trimmed, and median reads per sample. (**E**) Patients by number of niches sampled, per cohort.
(**F**) Samples analysed per cohort, with tissue types per patient inside each bar and the two cohorts
that also appear in that pooled analysis marked."""),

    ("""*Data:* 898 archived runs from six ENA/SRA projects, of which 822 samples from 245 patients entered
analysis after removing 76 runs from PRJNA1415359 that had no tissue; plus 40 hospital samples from 20
patients (20 tumour, 20 distant normal). (**A**) 782 public samples across 245 patients. (**C**) No
question is testable in fewer than three cohorts. (**D**) Median bacterial reads per sample 0.56–17.3 million across the seven datasets.
(**F**) Six niches per patient in PRJNA866676, two in the hospital sample, one in every other cohort.""",
     """*Data:* 898 archived runs from six ENA/SRA projects, of which 782 entered analysis after 116 runs from
PRJNA1415359 were removed (donors with no tissue, or with no tumour sample); plus 40 hospital samples
from 20 patients (20 tumour, 20 distant normal). (**A**) 782 public samples across 245 patients.
(**C**) Testable in 1 to 5 cohorts and confirmed in 5 and 4 for the two best-supported conclusions,
in 1 for the saliva gradient, and in none for the RNA-seq check. (**D**) Median reads per sample
12k–83k across the six public cohorts; the hospital libraries were delivered as genus-level counts
and have no equivalent setting. (**F**) At most six tissue types per patient in PRJNA866676, five in
PRJNA412445, four in PRJNA1415359 and two elsewhere; 674 of the 822 samples are not in that pool."""),

    # ---------------------------------------------------------------- Figure 2
    # D plots the SD of the estimate (one patient vs a cohort of 20); no interval-coverage panel
    # exists in this figure. F is the detection floor by depth; G the excess over the null by
    # source distinctness; H the sensitivity/specificity of calling a taxon present.
    ("""(**D**) Interval coverage by design. (**E**) Sensitivity to 25% reference leakage. (**F**) Median
excess over the null against the true fraction. (**G**) Power against the true fraction. (**H**)
Coverage against the true fraction. (**I**) Genera detected at a true fraction of 0.10, by sequencing
depth.""",
     """(**D**) Standard deviation of the estimate, one patient against a cohort of 20. (**E**) Sensitivity to
25% reference leakage. (**F**) Smallest detectable tumour-specific fraction, by sequencing depth.
(**G**) Excess over the null against the true fraction, by source distinctness. (**H**) Sensitivity and
specificity for calling a taxon present, against the true fraction. (**I**) Genera detected at a true
fraction of 0.10, by sequencing depth."""),

    ("""(**D**) Coverage is nominal or better for cohort-level fits and poor for
single-patient fits; the standard deviation of the estimate falls 2.3–5.3-fold from one patient to a
cohort of 20. (**E**) ≥0.97 sensitivity at every true fraction. (**I**) 9, 10, 11, 11 and 14 genera.""",
     """(**D**) The standard deviation of the estimate falls 2.3–5.3-fold from one
patient to a cohort of 20, and 2.7–4.1-fold at a true fraction of 0.10. (**E**) ≥0.97 sensitivity at
every true fraction. (**F**) 10% at every depth up to 130k reads and 5% at 968k. (**G**) At a true
fraction of 0.30 the excess runs from −0.03 with distinct sources to +0.26 with near-identical ones.
(**H**) Specificity is 1.0 at every true fraction and sensitivity 1.0 wherever novel taxa are
present, so the taxon call is not what fails. (**I**) 9, 10, 11, 11 and 14 genera."""),

    # ---------------------------------------------------------------- Figure 3
    # B pairs each patient's deep-tumour-to-saliva distance against the cancer surface's.
    ("""(**B**) Distance from tumour to
cancer surface against distance from tumour to saliva.""",
     """(**B**) Distance from deep
tumour to saliva against distance from the cancer surface to saliva, paired by patient."""),

    # ---------------------------------------------------------------- Figure 4
    # D is the genus-UNASSIGNED share (the blind spot), not the assigned share.
    ("""(**D**) Cumulative
share of reads carried by genus-assigned ASVs, by cohort.""",
     """(**D**) Cumulative
distribution of the share of reads in genus-unassigned ASVs, per sample, by cohort."""),

    ("""(**B**) 0.188% in tumour, 0.051% paracancerous, 0.024% cancer surface, 0.043% saliva, 0.000% tongue
coat.""",
     """(**B**) 0.188% in tumour, 0.051% paracancerous, 0.024% cancer surface, 0.043% saliva, 0.003%
contralateral, 0.000% tongue coat."""),

    # Table S1 gives the pooled log-ratios of these thirteen genera: -1.06 to +0.76.
    ("""(**L**) Pooled mean log-ratio +0.11 to +0.69.""",
     """(**L**) Pooled mean log-ratio −1.06 (Actinomyces) to +0.76 (Parvimonas)."""),

    # ---------------------------------------------------------------- Figure 5
    # B plots the strict oxygen-tolerance sets (15/20 up, AUC 0.78); 46.1 -> 63.5 and 37.5 -> 21.4
    # are that definition's percentages, not a mixture of definitions.
    ("""(**B**) Anaerobes 43.0% → 63.5%, aerobes 42.6% → 21.4%.""",
     """(**B**) On the strict oxygen-tolerance sets that this panel plots, anaerobes 46.1% → 63.5% of reads
and aerobes 37.5% → 21.4%."""),

    # J is the whole-profile scatter against the public pooled effect.
    ("""(**J**) +1.21 here against +0.11 to +0.69 in the meta-analysis.""",
     """(**J**) ρ = 0.21 across the 120 shared genera, rising to 0.76 among the 41 prevalent in at least half
the samples on both sides; *Parvimonas* itself is +1.21 here against a pooled +0.76 (0.46–1.05) in the
meta-analysis."""),

    # L: the pooled values are the four in s5_consortium_four_definitions_pooled.csv.
    ("""(**L**) strict +0.29 (0.14–0.44), taxonomy +0.32, taxonomy strict +0.22, eHOMD excludes zero only under a z interval.""",
     """(**L**) full taxonomy +0.27 (0.14–0.41), strict +0.29 (0.14–0.44), eHOMD +0.46 (0.12–0.79), original
sets +0.59 (0.18–0.99) and shown only for transparency; the intervals use a normal approximation, and
eHOMD no longer excludes zero under t(k − 2)."""),
]

for old, new in EDITS:
    n = text.count(old)
    assert n == 1, f"expected 1 occurrence, found {n}:\n{old[:120]}"
    text = text.replace(old, new)

assert text != orig
P.write_text(text, encoding="utf-8", newline="\n")
print(f"applied {len(EDITS)} caption corrections -> {P}")
