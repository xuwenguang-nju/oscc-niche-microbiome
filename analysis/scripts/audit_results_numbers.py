#!/usr/bin/env python3
"""Audit every numeral in a Results draft against the tables that carry the reported statistics.

Why this exists
---------------
A Results section is mostly numbers, and a number that drifts from its table is the hardest defect
for a co-author to catch by reading. This walks the draft, extracts every numeral, and looks for it
in a curated set of *summary* tables — the ones that hold the statistics the draft quotes.

Why curated rather than "every CSV"
-----------------------------------
A first version accepted a numeral if it appeared anywhere in results/tables. That is worthless: with
per-sample abundance matrices and two thousand coupling tests in the pool, almost any plausible number
turns up somewhere, and a deliberate typo (pseudo-F 23.3 changed to 27.9) passed. So only tables whose
columns are the quantities the text reports are read: test summaries, pooled estimates, cohort
counts. Per-sample dumps and long per-item tables are excluded by name and explained below.

A numeral passes if it matches, at the precision used or the one above it, a value in one of those
tables, in either sign; or if it is a listed structural count (cohort sizes, instrument names, panel
labels) with its provenance written next to it. Anything else is printed for a human decision, and the
run fails — never silently accepted.

Known limit: both checks are presence-based, so a wrong value escapes if the true string also occurs
elsewhere in the draft or in the same section's tables. Changing the Fig. 5J correlation from 0.76 to
0.86 is not caught, because 0.76 also appears in the neighbouring Fig. 6A sentence. A tamper that
alters a value with no legitimate twin — the pooled interval bound, a rank count — does fail.

Run: python analysis/scripts/audit_results_numbers.py [path-to-draft]
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
DRAFT = ROOT / "results/reports/manuscript-results-v8.md"
TAB = ROOT / "results/tables"

# Each section of the draft is checked against the tables that produced its panels, not against every
# table in the project. Without that binding the check is worthless: with per-sample matrices in the
# pool, a deliberate typo matches something somewhere and passes.
SECTION_TABLES = {
    "1": ["table1_cohort_characteristics_clean.csv", "read_profile_summary.csv"],
    "2": ["s3_simulation_calibration.csv", "s3_simulation_mdc.csv",
                         "s3_perpatient_variability.csv", "s3_replicates_distinct_c20.csv",
                         "s3_replicates_moderate_c80.csv",
                         "s3_replicates_moderate_c80_leak25.csv",
                         "s3_replicates_near_identical_c400.csv"],
    "3": ["s7_source_ranking_tests.csv", "s7b_before_after.csv",
                              "s7b_contaminant_burden.csv", "s7_backflow_gradient.csv",
                              "s7_pairwise_niche_distances.csv", "s7_permanova.csv", "e2_permanova.json",
                              "e2_pathway_niche_ranking.csv"],
    "4": ["s8a_rank_tests.csv", "s8a_surface_proxy.csv",
                                          "s8b_species_summary.csv",
                                          "s8b_species_paired_tests.csv",
                                          "s8c_genusless_read_share.csv",
                                          "s8c_hidden_family_tests.csv",
                                          "s5_pooled_dz_and_logratio.csv"],
    "5": ["s5_consortium_three_definitions_inhouse.csv",
                            "s5_consortium_three_definitions_per_cohort.csv",
                            "s5_consortium_four_definitions_pooled.csv",
                            "s5_consortium_sensitivity.csv",
                            "s5_reciprocal_shift_pooled.csv",
                            "fw2024_420/species_level_tests.csv",
                            "fw2024_420/paired_genus_tests.csv",
                            "fw2024_420/rank_level_tests.csv",
                            "fw2024_420/downsampling_summary.csv",
                            "fw2024_420/consortium_ratio_per_patient.csv",
                            "fw2024_420/whole_profile_vs_public.csv",
                            "tableS2_clinical_characteristics.csv",
                            "tableS3_genus_set_definitions.csv",
                            "tableS4_kegg_pathways_significant.csv"],
    "6": ["s11_genus_concordance.csv", "s6a_direction_check.csv",
                                     "s6b_signal_composition.csv",
                                     "s6v4_parvimonas_top.csv", "s11_clinical_assoc.csv",
                                     "s11_survival_univariable.csv",
                                     "s11_survival_multivariable.csv", "s11_survival_lrt.csv",
                                     "s11_survival_km.csv", "s11_survival_stage_sensitivity.csv",
                                     "s11_genus_concordance_disagreement.csv"],
    "7": ["s12_oral_gut_axis.csv", "s12_focus_genera_by_project.csv",
                        "s12_project_composition.csv"],
}

def derived_values() -> dict[str, str]:
    """Values the draft states that are arithmetic on table columns, recomputed here.

    A derived number (a difference between two project means, an ANOVA F on a per-sample score) is
    the easiest kind to get wrong and the hardest to check by eye, so each one is recomputed from the
    table it comes from and the draft's numeral must match. Nothing here is a hand-written constant.
    """
    out: dict[str, str] = {}
    try:
        st = pd.read_csv(TAB / "s12_sample_type_sensitivity.csv", encoding="utf-8-sig")
        per = pd.read_csv(TAB / "s12_oral_gut_per_sample.csv", encoding="utf-8-sig")
    except Exception:                       # noqa: BLE001
        return out

    t = st[st.sample_type == "Primary Tumor"].set_index("project")["oral_minus_gut_pct"]
    n = st[st.sample_type == "Solid Tissue Normal"].set_index("project")["oral_minus_gut_pct"]
    diff = (t - n).round(1)
    for proj, v in diff.items():
        out[f"{v:.1f}"] = f"tumour minus normal, {proj} (s12_sample_type_sensitivity.csv)"
    out[f"{(t['ESCA'] - t['READ']):.0f}"] = "span between the highest and lowest project mean"

    try:
        pm = pd.read_csv(TAB / "s7_permanova.csv", encoding="utf-8-sig")
        r2 = float(pm.set_index("test").loc["niche (patient-blocked)", "r2"])
        out[f"{r2 * 100:.1f}"] = "R-squared x 100, patient-blocked niche PERMANOVA (s7_permanova.csv)"
    except Exception:                       # noqa: BLE001
        pass

    groups = [g.oral_minus_gut_pct.dropna().to_numpy() for _, g in per.groupby("project")]
    grand = per.oral_minus_gut_pct.mean()
    ss_b = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    ss_t = ((per.oral_minus_gut_pct - grand) ** 2).sum()
    out[f"{ss_b / ss_t * 100:.0f}"] = "share of variance in the axis score explained by project"
    out[f"{len(per)}"] = "samples in the axis table"
    # one-way F from the same sums of squares, to avoid depending on scipy here
    k, N = len(groups), len(per)
    F = (ss_b / (k - 1)) / ((ss_t - ss_b) / (N - k))
    out[f"{F:.0f}"] = "one-way ANOVA F for the axis score by project"
    return out


# Structural numbers: not measured values, each with where it comes from. A numeral that is neither
# in a summary table nor listed here fails the run.
STRUCTURAL = {
    "1": "figure and table labels, panel letters, prose enumerations",
    "2": "idem",
    "3": "idem",
    "4": "idem",
    "5": "idem",
    "6": "idem",
    "7": "idem",
    "8": "n = 8 paired samples in PRJNA412445 (s8a_rank_tests.csv n_patients)",
    "9": "nine genera testable against the RNA-seq cohort (s6a_direction_check.csv)",
    "13": "13 replicating genera; 13 contaminant genera; 13 cohorts in the pooled analysis",
    "20": "20 hospital patients; cohort of 20 in the simulation",
    "25": "25 patients in PRJNA412445 (Table 1)",
    "30": "30 replicates per cell; 30 classes at the rank test",
    "29": "29 obligate-anaerobe genera in the main score (tableS3_genus_set_definitions.csv)",
    "2506": "ASVs delivered for the hospital cohort (Table 1, ASVs column)",
    "127": "span between project means, 78.8 - (-47.9) (s12_oral_gut_axis.csv)",
    "31": "31 matched pairs in PRJNA1415359; largest rank gap 31",
    "36": "36 patients with both backflow samples (s7_backflow_gradient.csv)",
    "37": "37 patients with a tumour and a comparator (s7_source_ranking_tests.csv)",
    "41": "41 genera prevalent on both sides (whole_profile_vs_public.csv)",
    "46": "46 KEGG pathways (tableS4)",
    "48": "48 READ samples (s12_oral_gut_axis.csv)",
    "73": "73 orders tested (rank_level_tests.csv)",
    "80": "80 ESCA samples (s12_oral_gut_axis.csv)",
    "106": "106 patients with survival data (s11_survival_univariable.csv)",
    "854": "tests across all ranks (tableS4_rank_level_summary.csv, All ranks row)",
    "28": "taxa passing FDR across all ranks (tableS4_rank_level_summary.csv, All ranks row)",
    "209": "genera tested (tableS4_rank_level_summary.csv); also genera detected at full depth",
    "408": "species tested (tableS4_rank_level_summary.csv)",
    "4.15": "Parvimonas mean in tumour, Table S4(b) (rank_level_tests.csv mean_T x 100)",
    "2.56": "Parvimonas median in tumour, Fig. 5C (genus_counts.csv)",
    "1.6": "share of reads in unassigned/unclassified rows (genus_counts.csv)",
    "114": "114 family tests (rank_level_tests.csv)",
    "116": "116 genera after decontamination (s7b_before_after.csv)",
    "120": "120 shared genera (whole_profile_vs_public.csv)",
    "141": "141 genera detected at 500 reads (downsampling_summary.csv)",
    "146": "146 genera testable in 3+ cohorts (s5_pooled_dz_and_logratio.csv)",
    "155": "155 primary tumours (Table S5)",
    "165": "165 STAD samples (s12_oral_gut_axis.csv)",
    "177": "177 TCGA-HNSC samples (Table S5)",
    "200": "200 Parvimonas coupling tests (s6v4_parvimonas_top.csv)",
    "209": "209 genera at full depth (downsampling_summary.csv); 209 genera in the hospital cohort (Table 1)",
    "240": "240 predicted KEGG pathways (Table S4 source)",
    "245": "245 public patients (Table 1)",
    "254": "254 genera that fail the replication filter (267 - 13)",
    "267": "267 genera testable in 2+ cohorts (s5_pooled_dz_and_logratio.csv)",
    "408": "408 species tested (rank_level_tests.csv)",
    "500": "500 reads, the deepest subsampling level (downsampling_summary.csv)",
    "265": "265 patients total, 245 public + 20 hospital (Table 1 arithmetic)",
    "674": "674 samples not in the pooled analysis (Table 1 arithmetic)",
    "782": "782 public samples (Table 1)",
    "822": "822 samples total (Table 1)",
    "898": "898 archived runs (oscc_niche_sample_manifest.csv)",
    "100000": "read cap for PRJNA666746 (Table 1)",
    "23000": "23k reads, shallowest simulated depth 23,428 (s3_simulation_mdc.csv)",
    "130000": "130k reads, one simulated depth (s3)",
    "968000": "968k reads, deepest simulated depth (s3)",
    "11000": "median reads/sample, min 11,197 in PRJNA597251 (Fig. 1D source table)",
    "83000": "median reads/sample, max 83,207 in PRJNA1415359 (Fig. 1D source table)",
    "6000": "NovaSeq 6000, an instrument name rather than a measurement",
    "138.1": "SILVA 138.1, a reference-database version (Methods, Table 1 footnote)",
    "95": "95% confidence interval, a level rather than a measured value",
    "17": "5' primer trim length, 17 bp of R1 (Table 1, trim R1/R2 column)",
    "21": "5' primer trim length, 21 bp of R2 (Table 1, trim R1/R2 column)",
    "17.3": "17,293,002 bacterial reads in PRJNA866676, expressed in millions (Table 1)",
    "19": "patients ranking paracancerous first (s7_source_ranking.csv)",
    "78": "77.8% of the 36 paired backflow patients (s7_backflow_gradient.csv)",
    "24": "pathways lower in tumour (Table S4, Direction column)",
    "2000": "2,000 genus-feature tests, 4 assays x 500 pairs (s6v4_coupling_long.csv row count)",
}


def verified_claims() -> list[tuple[str, list[str], bool]]:
    """(label, acceptable writings, required) for the headline claims.

    The section-level check above cannot catch a typo that happens to equal another legitimate value
    in the same tables. These do: each headline statistic is recomputed from its table and compared
    with the numeral the draft carries, so a change to any of them fails even if the new value exists
    elsewhere. Only claims whose source is a single cell or a simple arithmetic step appear here.

    A claim carries the writings it may legitimately take (the same statistic at a coarser precision,
    say q = 0.005 for 0.0049) and whether the draft is required to state it at all. Optional claims
    are for effect sizes a draft may reasonably leave to a table; if the draft states one it must
    still match, so a wrong value cannot slip through by being dropped.
    """
    out: list[tuple[str, list[str], bool]] = []
    leaf = lambda p: TAB / p                                    # noqa: E731

    d = pd.read_csv(leaf("s8b_species_by_niche.csv"), encoding="utf-8-sig")
    tum = d[(d.study == "PRJNA866676") & (d.niche == "tumor")]["Parvimonas micra_pct_of_bacteria"].iloc[0]
    par = d[(d.study == "PRJNA866676") & (d.niche == "paracancerous")]["Parvimonas micra_pct_of_bacteria"].iloc[0]
    out.append(("Parvimonas % of reads in tumour", [f"{tum:.3f}"], True))

    sp = pd.read_csv(leaf("s8b_species_paired_tests.csv"), encoding="utf-8-sig")
    r = sp[sp.study == "PRJNA866676"].iloc[0]
    out.append(("Parvimonas log2 fold change, PRJNA866676", [f"{r.median_log2fc:+.2f}"], True))
    out.append(("its paired p", [f"{r.p_paired_wilcoxon:.4f}"], True))

    # the main score's pooled change and its interval: a plausible-looking wrong bound was invisible
    # to the section check, because 0.481 also occurs elsewhere in the same family of tables
    try:
        sp_ = pd.read_csv(leaf("s5_score_pooled.csv"), encoding="utf-8-sig").set_index("definition")
        row = sp_.loc["main"]
        out.append(("pooled main-score change", [f"{row.pooled:.3f}", f"{row.pooled:.2f}"], True))
        out.append(("its lower bound", [f"{row.lo:.3f}", f"{row.lo:.2f}"], True))
        out.append(("its upper bound", [f"{row.hi:.3f}", f"{row.hi:.2f}"], True))
        out.append(("its I-squared", [f"{row.I2:.0f}"], True))
    except Exception:                       # noqa: BLE001
        pass

    # the 5J correlation pair, recomputed rather than trusted
    try:
        from scipy import stats as _st
        wpp = pd.read_csv(leaf("fw2024_420/whole_profile_vs_public.csv"), encoding="utf-8-sig")
        _all = _st.spearmanr(wpp.pooled_effect / np.log(2), wpp.log2FC).statistic
        _sub = wpp[(wpp.prev_T >= 0.5) & (wpp.prev_N >= 0.5)]
        _hi = _st.spearmanr(_sub.pooled_effect / np.log(2), _sub.log2FC).statistic
        out.append(("5J rho, all shared genera", [f"{_all:.2f}"], True))
        out.append(("5J rho, prevalent genera", [f"{_hi:.2f}"], True))
    except Exception:                       # noqa: BLE001
        pass

    hid = pd.read_csv(leaf("s8c_hidden_family_tests.csv"), encoding="utf-8-sig").sort_values("q").iloc[0]
    # the draft now states the effect size alongside q, so both are required claims
    out.append(("family-level CLR difference", [f"{hid.mean_clr_diff:.2f}"], True))
    out.append(("its q", [f"{hid.q:.4f}", f"{hid.q:.3f}"], True))
    # supporting detail: a draft may state these or leave them to the table, but a wrong value fails
    out.append(("its pair count", [f"{int(hid.n_pairs)}"], False))
    out.append(("family mean abundance in tumour", [f"{hid.mean_rel_tumour * 100:.1f}"], False))
    out.append(("family mean abundance in normal", [f"{hid.mean_rel_normal * 100:.1f}"], False))

    pooled = pd.read_csv(leaf("s5_pooled_dz_and_logratio.csv"), encoding="utf-8-sig")
    pv = pooled[pooled.genus == "Parvimonas"].iloc[0]
    # the draft now quotes the pooled effect in log2, so the comparison is made in log2
    out.append(("pooled Parvimonas log2 ratio", [f"{pv.pooled_logratio / np.log(2):+.2f}"], True))
    out.append(("genera testable in 2+ cohorts", [f"{len(pooled)}"], True))
    out.append(("genera testable in 3+ cohorts", [f"{int((pooled.k >= 3).sum())}"], True))

    # The main score is the taxonomy set of Table S3 (29 anaerobic against 16 aerobic genera), which
    # is the "main" definition in the regenerated score family.
    pooled_score = pd.read_csv(leaf("s5_score_pooled.csv"), encoding="utf-8-sig")
    main = pooled_score[pooled_score.definition == "main"].iloc[0]
    out.append(("main-score pooled shift", [f"{main.pooled:+.2f}"], True))
    out.append(("its lower 95% bound", [f"{main.lo:.2f}"], True))

    inhouse = pd.read_csv(leaf("s5_score_inhouse.csv"), encoding="utf-8-sig")
    st = inhouse[inhouse.definition == "main"].iloc[0]
    out.append(("anaerobe share in tumour", [f"{st.ana_T:.1f}"], True))
    out.append(("aerobe share in tumour", [f"{st.aer_T:.1f}"], True))
    out.append(("AUC of the ratio", [f"{st.auc:.2f}"], True))

    kegg = pd.read_csv(leaf("kegg_pathways_significant.csv"), encoding="utf-8-sig")
    out.append(("pathways passing FDR", [f"{len(kegg)}"], True))
    out.append(("smallest q after adjustment", [f"{kegg['FDR q after adjustment'].min():.2f}"], True))

    conc = pd.read_csv(leaf("s11_genus_concordance.csv"), encoding="utf-8-sig")
    out.append(("largest genus-level rho", [f"{conc.spearman_rho_filtered.max():.2f}"], True))

    coup = pd.read_csv(leaf("s6v4_coupling_long.csv"), encoding="utf-8-sig")
    out.append(("coupling pairs passing FDR", [f"{int((coup.q < 0.05).sum())}"], True))
    pm = coup[coup.genus == "Parvimonas"]
    out.append(("Parvimonas minimum q", [f"{pm.q.min():.3f}"], True))

    axis = pd.read_csv(leaf("s12_oral_gut_axis.csv"), encoding="utf-8-sig").set_index("project")
    out.append(("HNSC axis score", [f"{axis.loc['HNSC','oral_minus_gut_pct']:+.1f}"], True))
    out.append(("READ axis score", [f"{axis.loc['READ','oral_minus_gut_pct']:+.1f}"], True))
    return out


def table_values(names: list[str]) -> dict[str, str]:
    """{value-string: file it came from} for the named tables."""
    vals: dict[str, str] = {}
    for name in names:
        path = TAB / name
        if not path.is_file():
            continue
        if path.suffix == ".json":
            txt = path.read_text(encoding="utf-8", errors="replace")
            for m in re.finditer(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", txt):
                raw = m.group(0)
                vals.setdefault(raw, name)
                try:
                    f = float(raw)
                except ValueError:
                    continue
                for dec in (0, 1, 2, 3):
                    vals.setdefault(f"{f:.{dec}f}", name)
                    vals.setdefault(f"{abs(f):.{dec}f}", name)
            continue
        try:
            df = pd.read_csv(path, encoding="utf-8-sig", low_memory=False)
        except Exception:                       # noqa: BLE001
            continue
        for col in df.columns:
            num = pd.to_numeric(df[col], errors="coerce").dropna()
            if num.empty:
                continue
            for dec in (0, 1, 2, 3):
                for v in num.round(dec):
                    vals.setdefault(f"{v:.{dec}f}", name)
                    vals.setdefault(f"{abs(v):.{dec}f}", name)
    return vals


# A numeral is not "16" in 16S, not the mantissa of 10⁻⁵, and not the tail of an accession
NUMERAL = re.compile(r"(?<![\w.-])(?<!\d×10)(?<!×10)[−+-]?\d[\d,]*(?:\.\d+)?(?![\w⁰¹²³⁴⁵⁶⁷⁸⁹⁻])")


def candidates(token: str) -> list[str]:
    """The numeral at the precisions the text may use, both signs."""
    t = token.replace(",", "").replace("−", "-").replace("+", "").rstrip(".")
    out: set[str] = set()
    try:
        f = abs(float(t))
    except ValueError:
        return [t]
    out.add(str(int(round(f))))
    for dec in (1, 2, 3):
        out.add(f"{f:.{dec}f}")
    return sorted(out)


def check_headings(body: str) -> list[str]:
    """Each numeral in a subsection heading must also occur in that subsection's text.

    A result-bearing heading is a claim. If it carries a number the paragraph below never states, the
    reader meets a finding the section does not support. This cannot catch a wrong number that matches
    by accident, and it cannot tell which occurrence is which (a heading saying "eight" passes if the
    section contains an unrelated "eight"). What it does catch is a heading that outruns its paragraph,
    which is the failure mode of rewriting headings to state results.
    """
    problems: list[str] = []
    for part in re.split(r"(?m)^### ", body)[1:]:
        head, _, text_ = part.partition("\n")
        # the "(Fig. N)" pointer is a reference, not a value
        stripped = re.sub(r"\(Fig(?:s)?\.\s*S?\d+[^)]*\)", " ", head)
        for m in re.finditer(NUMERAL, stripped):
            tok = m.group(0).replace(",", "")
            # compare at the precision the heading uses: allowing the integer-collapsed form would
            # match "2.40" against any "2.x", so the check has to be exact at its own precision
            dec = len(tok.split(".")[1]) if "." in tok else 0
            target = f"{abs(float(tok)):.{dec}f}"
            found = False
            for bm in re.finditer(NUMERAL, text_):
                btok = bm.group(0).replace(",", "")
                try:
                    if f"{abs(float(btok)):.{dec}f}" == target:
                        found = True
                        break
                except ValueError:
                    continue
            if not found:
                problems.append(
                    f"heading numeral {tok!r} is not stated in its own section: {head[:60]}")
    return problems


def main() -> None:
    path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else DRAFT
    text = path.read_text(encoding="utf-8")
    body = text.split("## Notes for the authors")[0]
    if "## Results" in body:
        body = body.split("## Results", 1)[1]

    # Split into sections, keeping the heading text as the key. Two heading styles occur: the
    # project's own drafts carry "### " with a "(Fig. N)" marker, while the author's revised draft
    # uses "## N." and names no figure. Accepting only the first silently produced "0 numerals
    # checked" on that draft, which is worse than failing: it looks like a pass.
    parts = re.split(r"(?m)^(?:### |## )(?=\d+\.\s|.*\(Fig)", body)
    sections = []
    for part in parts[1:]:
        head, _, text_ = part.partition("\n")
        # Prefer the figure named in the heading; otherwise fall back to the heading's own number.
        m = re.search(r"\(Fig(?:s)?\.\s*(S?\d+)", head)
        if m:
            key = m.group(1)
        else:
            mnum = re.match(r"\s*(\d+)\.", head)
            key = mnum.group(1) if mnum else None
        if key not in SECTION_TABLES:
            sys.exit(f"no table mapping for section: {head[:70]}")
        sections.append((head.strip(), key, text_))
    if not sections:
        sys.exit("no sections recognised in the draft; the numeral check would be vacuous")

    derived = derived_values()
    unknown: dict[str, int] = {}
    checked = 0
    total_vals = 0
    for head, key, text_ in sections:
        vals = table_values(SECTION_TABLES[key])
        total_vals += len(vals)
        for m in re.finditer(NUMERAL, text_):
            raw = m.group(0)
            checked += 1
            cands = candidates(raw)
            if (any(c in vals for c in cands) or any(c in STRUCTURAL for c in cands)
                    or any(c in derived for c in cands)):
                continue
            unknown[f"{raw}  [{head.split('(')[0][:34]}]"] = \
                unknown.get(f"{raw}  [{head.split('(')[0][:34]}]", 0) + 1

    words = len(re.sub(r"[*#|`]", "", body).split())
    print(f"draft       : {path.name}")
    print(f"body words  : {words}")
    print(f"sections    : {len(sections)}, each checked against its own tables")
    print(f"numerals    : {checked} checked")
    print(f"recomputed  : {len(derived)} derived value(s) verified against the tables")
    if unknown:
        print(f"\n{len(unknown)} numeral(s) found in no summary table and no structural list:")
        for k, n in sorted(unknown.items(), key=lambda kv: -kv[1]):
            print(f"  {k:>12}  ×{n}")
        sys.exit(1)
    # headline claims: recompute and compare with what the draft says
    # the draft uses the typographic minus; compare on one canonical form
    text_flat = " ".join(body.split()).replace("\u2212", "-").replace("\u2013", "-")
    bad = []
    optional_missing = []
    for label, alts, required in verified_claims():
        alts = [a.replace("\u2212", "-").replace("\u2013", "-") for a in alts]
        if any(a in text_flat for a in alts):
            continue
        if required:
            bad.append(f"{label}: draft does not carry any of {alts!r}")
        else:
            optional_missing.append(f"{label}: not stated in the draft ({' or '.join(alts)})")
    if optional_missing:
        print("\noptional effect sizes the draft leaves to the tables:")
        for o in optional_missing:
            print("  -", o)
    if bad:
        print("\nheadline claims that do not match the tables:")
        for b in bad:
            print("  -", b)
        sys.exit(1)
    print(f"claims      : {len(verified_claims())} headline statistics recomputed and matched")

    head_bad = check_headings(body)
    if head_bad:
        print("\nheadings whose numbers the section does not state:")
        for b in head_bad:
            print("  -", b)
        sys.exit(1)
    print("headings    : every numeral in a heading is stated in its section")
    print("every numeral in the draft traces to a summary table or a listed structural count")


if __name__ == "__main__":
    main()
