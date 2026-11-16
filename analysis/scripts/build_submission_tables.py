#!/usr/bin/env python3
"""Build submission-ready table files (.docx) from the result CSVs.

One .docx per table, each with a caption block above the table and footnotes below, a repeating
header row, and a fixed layout so the columns keep their widths. Written for the Journal of Oral
Microbiology submission package, whose instructions ask for tables that are editable and
interpretable without the main text.

Run: python analysis/scripts/build_submission_tables.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results" / "tables"
OUT = ROOT / "results" / "reports" / "tables_word"

def _tidy_axis_sensitivity(df: pd.DataFrame) -> pd.DataFrame:
    """Turn the 19 raw sensitivity fields into the six columns a reader needs.

    The three interval types are the point of the table (z is optimistic at k = 6, t(k-2) is the
    honest one, HKSJ the alternative), so they are kept side by side rather than spread over six
    cryptic columns.
    """
    f = lambda lo, hi: f"{lo:+.3f} to {hi:+.3f}"          # noqa: E731
    out = pd.DataFrame({
        "Definition": df["definition"].str.replace("_", " ", regex=False),
        "k": df["k"],
        "Pooled \u0394log\u2081\u2080": df["pooled"].map(lambda v: f"{v:+.3f}"),
        "SE": df["se"].map(lambda v: f"{v:.3f}"),
        "I\u00b2 (%)": df["I2_pct"].map(lambda v: f"{v:.0f}"),
        "95% CI (z)": [f(a, b) for a, b in zip(df["lo_z"], df["hi_z"])],
        "95% CI t(k\u22122)": [f(a, b) for a, b in zip(df["lo_t"], df["hi_t"])],
        "95% CI HKSJ": [f(a, b) for a, b in zip(df["lo_hksj"], df["hi_hksj"])],
        "LOO range": [f(a, b) for a, b in zip(df["loo_pooled_min"], df["loo_pooled_max"])],
    })
    # the drop-hospital-cohort interval is only present in the older table; keep the column when it
    # exists so the table does not silently lose a sensitivity analysis
    if "drop_ci_low" in df and df["drop_ci_low"].notna().any():
        out["Drop hospital cohort"] = [f(a, b) for a, b in
                                       zip(df["drop_ci_low"], df["drop_ci_high"])]
    else:
        out["LOO all positive"] = df["loo_all_positive"].map(lambda v: "yes" if v else "no")
    return out



SPECS = [
    # ---- main tables -------------------------------------------------------------------
    dict(
        name="Table1",
        src="table1_cohort_characteristics_clean.csv",
        caption="**Table 1. The seven datasets: sampling, sequencing and processing.**",
        font=8.0,
        widths=[2.2, 2.0, 1.0, 1.0, 2.4, 1.1, 1.1, 1.2, 1.2, 1.3, 0.9, 1.1,
                1.1, 1.1, 1.8, 1.5],
        footnotes=[
            "Rows 1–6 are public deposits re-processed from raw reads in this study; row 7 is the hospital "
            "cohort described in Methods (independent hospital-based convenience sample, not a registered study).",
            "DADA2 was run with pool = FALSE for every cohort. trim and trunc give the R1/R2 cycle values; "
            "two cohorts carried a 17/21 bp 5\u2032 primer residue that was trimmed first.",
            "PRJNA666746 was down-sampled to 100,000 reads per sample (read cap) before denoising; no other "
            "cohort required capping.",
            "The hospital cohort was processed by the sequencing provider with DADA2 and Bayesian taxonomic "
            "assignment against SILVA 138.1, whereas rows 1-6 were re-processed here against SILVA 138.2; "
            "per-sample settings are therefore not applicable.",
            "What the provider delivered for row 7, and which result each part supports: an ASV table with "
            "its sequences (2,506 ASVs, the ASVs column), absolute and relative abundance tables at all seven "
            "ranks from kingdom to species (the Genera column, and the species and phylum results in Fig. 5H "
            "and 5I), alpha and beta diversity tables, and PICRUSt2 predictions of KEGG, MetaCyc and EC "
            "abundances, which underlie the pathway results in Fig. 5M and 5N and Table S5. Fig. S1 is a "
            "different cohort (PRJNA866676) and does not draw on these tables. No raw reads were transferred.",
            "Genus-unassigned gives the median percentage of reads per sample not assigned to a genus; it is "
            "not available for the hospital cohort because the delivered tables are already collapsed by rank.",
        ],
    ),
    dict(
        name="TableS2",
        src="tableS2_clinical_characteristics.csv",
        caption="**Table S2. Clinical characteristics of the 20 hospital-cohort patients.**",
        font=9.0,
        widths=[1.7, 1.7, 1.6, 3.0, 3.0, 1.5],
        footnotes=[
            "Each patient contributed one tumour and one distant normal mucosal sample, collected at the same "
            "operation; the paired analysis uses all 20 pairs.",
            "Patients are identified by study codes. Direct identifiers are held by the hospital and were not "
            "transferred for analysis.",
            "All 20 tumours were tongue squamous cell carcinoma. T, N and M are reproduced exactly as they "
            "stand in the hospital record, which lists a single TNM field without an edition and without a "
            "clinical or pathological prefix; grade is the reported histological grade.",
            "Two things the record does not settle, and which are therefore not claimed here. (i) Edition: one "
            "patient (P16) is recorded as N3b, a subcategory that exists only from the AJCC/UICC 8th edition "
            "— the 7th subdivides N2 but leaves N3 undivided — so the table is consistent with the 8th "
            "edition but this is an inference from a single entry, not a stated fact. (ii) cTNM or pTNM: the "
            "T4a, N2b and N3b subcategories are used in both clinical and pathological staging, and the "
            "record carries no prefix, so the pathological status is unknown. Both must be confirmed against "
            "the original charts before submission; neither affects any test in this paper, which uses "
            "stage only as an ordinal covariate in the TCGA analyses and not for this cohort.",
            "Patients were not accrued consecutively; see Methods for the sampling frame and its limitations.",
        ],
    ),
    # ---- supplementary tables ----------------------------------------------------------
    dict(
        name="TableS1",
        src="tableS1_replicating_genera.csv",
        caption="**Table S1. The thirteen genera that replicate across three or more cohorts.**",
        font=8.5,
        widths=[3.2, 1.8, 1.9, 1.1, 3.0, 1.5, 1.2, 1.7, 2.0, 2.0, 2.0],
        footnotes=[
            "Ordered by pooled effect. The pooled effect is a log\u2082 ratio of relative abundances, "
            "tumour over normal, from a DerSimonian\u2013Laird random-effects model over the cohorts in "
            "which the genus was testable; the 95% CI is the normal approximation. The same unit is used "
            "on the axis of Fig. 4L.",
            "Largest |d_z| is the biggest standardised effect the genus reached in any single cohort. It "
            "is not evidence about whether that cohort drives the pooled estimate; that question is "
            "answered by the leave-one-cohort-out analysis, not by this column.",
            "Literature columns come from the curated ORMCKB survey described in Methods and are counted "
            "per genus, not per paper direction. \u201cNone reported\u201d marks genera absent from the prior "
            "OSCC literature, for which agreement is not applicable rather than negative; \u201cconflict\u201d "
            "marks genera whose prior reports disagree with each other, and agreement there is not scored.",
        ],
    ),
    dict(
        name="TableS3",
        src="tableS3_genus_set_definitions.csv",
        caption="**Table S3. Genus-set definitions used in the oxygen-tolerance sensitivity analysis.**",
        font=7.5,
        widths=[3.0, 2.0, 3.0, 3.0, 2.6, 3.4],
        footnotes=[
            "Set 2 (29 obligate anaerobes against 16 aerobic or facultative genera) is the main score, used "
            "for every inference in the paper and for panels B, E, F, K and L of Fig. 5. Set 1 adds Gemella "
            "to the aerobic arm (29 against 17). Set 3 restricts the original sets to the genera the eHOMD "
            "annotation covers (8 against 3). Set 4 was chosen from this cohort's own differential result; "
            "it is reported for transparency only and is never used for inference.",
            "The two taxonomy sets (1 and 2) were assigned from documented oxygen tolerance alone, without "
            "reference to any test on these data, and are reproduced in full here so neither is implicit. "
            "Set 3 inherits that assignment and only restricts it. Set 4 is the exception: it was chosen from "
            "this cohort's own differential result, which is why it is not used for inference, and its "
            "membership is recorded here rather than being described as independent. Fig. 5L compares the "
            "four pooled estimates and part (b) tabulates them.",
            "\u201cexcluded\u201d marks Gemella, whose oxygen tolerance is classified inconsistently across "
            "sources; \u201cnot annotated\u201d means the genus is absent from the eHOMD annotation, not that "
            "it is aerotolerant (eHOMD covers 11 of the 46 genera); \u201cnot assigned\u201d means the genus "
            "is in neither arm of that set.",
        ],
        # part b: the pooled axis estimate under each independent definition, so that the intervals
        # quoted in Results 5 (t(k-2), leave-one-cohort-out, drop the hospital cohort) are tabulated.
        extra=dict(
            src="s5_score_sensitivity.csv",
            caption="Table S3 (b). Pooled anaerobe : aerobe shift under each definition.",
            font=7.5,
            widths=[3.0, 0.7, 2.2, 1.2, 1.4, 2.9, 2.9, 2.9, 2.6, 2.2],
            transform=_tidy_axis_sensitivity,
            footnotes=[
                "DerSimonian\u2013Laird random effects over the six cohorts; \u0394log\u2081\u2080 of the "
                "anaerobe : aerobe ratio, tumour minus normal. All six cohort-level estimates are positive "
                "under every definition.",
                "With k = 6 the normal-approximation interval is optimistic, so t(k \u2212 2) and "
                "Hartung\u2013Knapp\u2013Sidik\u2013Jonkman intervals are given alongside it; under t(k \u2212 2) "
                "the eHOMD definition no longer excludes zero, which is why Results 5 claims consistency of "
                "direction rather than a pooled magnitude.",
                "The leave-one-cohort-out range gives the spread of the pooled estimate when any single "
                "cohort is removed; under the main definition every such range stays positive and its "
                "normal-approximation interval excludes zero. A separate drop-the-hospital-cohort fit is "
                "not tabulated because that cohort is one of the six and its removal is one of the "
                "leave-one-out cases already shown.",
            ],
        ),
    ),
    dict(
        name="TableS6",
        src="tcma_sample_annotation.csv",
        caption="**Table S6. Per-sample annotation of the TCGA-HNSC cohort (177 samples).**",
        font=7.5,
        widths=[1.6, 1.2, 1.7, 3.0, 1.1, 1.0, 0.8, 0.9, 0.9, 2.4, 1.1, 1.0, 0.8, 0.8,
                1.4, 1.1, 1.9, 1.2, 1.2, 1.0],
        footnotes=[
            "One row per aliquot. Diversity metrics and Parvimonas abundance are the values used for the "
            "survival analyses.",
            "os_days is GDC days_to_death for a death and days_to_last_follow_up otherwise, and is the "
            "quantity every survival model uses. os_months is that value divided by 30.4375 (365.25/12) for "
            "display and is not an input to any model. os_event is 1 for death and 0 for a censored "
            "patient. An earlier version of this column was labelled os_months but contained days, so "
            "values such as 3,981 were not months.",
            "The table holds 155 primary tumours and 22 solid tissue normal samples, and no metastatic "
            "sample; the paired comparison uses primary tumour against solid tissue normal within a case.",
            "Missing values are shown blank.",
        ],
    ),
    dict(
        # (a) one row per taxonomic rank, because Results 5 quotes six fractions that otherwise
        # live only in an 854-row analysis intermediate
        name="TableS4",
        src="tableS4_rank_level_summary.csv",
        caption="**Table S4. Tests at every taxonomic rank in the hospital cohort, and the ones that pass FDR.**",
        font=8.0,
        widths=[3.2, 2.2, 2.4, 2.8, 2.8, 2.4],
        footnotes=[
            "Paired Wilcoxon tests of tumour against distant normal mucosa in the 20 hospital patients, "
            "run independently at each rank on the provider's relative-abundance tables; q is "
            "Benjamini\u2013Hochberg within that rank. Part (b) lists every test that passes.",
            "The per-rank counts are what Results 5 reports (1 of 20 phyla, 3 of 30 classes, 7 of 73 "
            "orders, 9 of 114 families, 6 of 209 genera, 2 of 408 species). They are not a correction "
            "across ranks: each rank is a separate family of tests, and the claim is only that a "
            "difference is detectable at every rank, not that any single taxon replicates.",
            "Abundances at each rank are percentages of that rank's own total (all rows of the provider's "
            "table for that rank, unassigned and unclassified included), so the percentages sum to 100 at "
            "every rank; a test at one rank is not independent of the test on the same reads at another.",
        ],
        extra=dict(
            src="tableS4_rank_level_significant.csv",
            caption="Table S4 (b). The 28 taxa whose tumour\u2013normal difference survives FDR, by rank.",
            font=8.0,
            widths=[2.0, 4.4, 2.0, 2.0, 2.6, 2.6, 2.0, 1.6, 1.6],
            footnotes=[
                "Ordered by rank, then by q within rank. The two abundance columns are the arithmetic "
                "mean across the 20 paired patients of each sample's percentage, not a median, and "
                "\u201cDifference\u201d is tumour mean minus normal mean. \u201cPatients agreeing\u201d is the "
                "share of the 20 pairs whose sign matches the direction in that row, so it counts "
                "1 \u2212 frac_up for a taxon that is lower in tumour.",
                "Denominator: at every rank, each taxon is divided by the total of all rows at that same "
                "rank in the provider's absolute-abundance table, so the percentages at a rank sum to 100 "
                "independently of the other ranks. Rows for unassigned and unclassified taxa are retained "
                "\u2014 41 of the 209 genus-level rows, carrying 1.6% of reads \u2014 because dropping them "
                "would silently renormalise every other value.",
                "These means are not interchangeable with the values in Fig. 5: that figure reports "
                "medians of the same per-sample percentages, which is why *Parvimonas* appears here at "
                "4.153% and 0.549% but at 2.56% and 0.39% there. Both describe the same 20 patients and the "
                "same denominator; the difference is the summary statistic, and it is large because this "
                "genus is absent from several samples and has a long right tail.",
                "Actinobacteriota is the one phylum-level difference, and it is the shift drawn in "
                "Fig. 5I. Each rank is corrected separately within its own family of tests, so a q at one "
                "rank is not comparable with a q at another; the text therefore reports the per-rank counts "
                "instead of nominating a single winning taxon.",
            ],
        ),
    ),
    dict(
        name="TableS5",
        src="kegg_pathways_significant.csv",
        caption="**Table S5. Predicted KEGG pathways that differ between tumour and normal mucosa.**",
        font=7.5,
        widths=[1.5, 6.2, 2.4, 1.5, 1.5, 1.6, 1.3, 1.3, 1.7, 1.9],
        footnotes=[
            "All 46 of the 240 PICRUSt2-predicted pathways that pass FDR in the hospital cohort, split by "
            "direction and ordered by q. PICRUSt2 infers pathway abundance from 16S taxonomy, so these are "
            "predictions, not measurements.",
            "Several maps carry host-disease names (for example Prostate cancer, Bladder cancer, Systemic lupus "
            "erythematosus) because bacterial orthologues populate those KEGG maps; the names do not indicate "
            "that a human disease process is present.",
            "The right-hand column gives FDR q after adjusting each pathway for the anaerobe-to-aerobe axis "
            "score. No pathway survives, which is the result reported in Fig. 5N.",
            "Names were resolved from KEGG on 2026-11-15. One identifier (ko00471) no longer resolves and is "
            "left as its accession.",
        ],
    ),
]

def _shade(cell):
    tcPr = cell._tc.get_or_add_tcPr()
    sh = OxmlElement("w:shd")
    sh.set(qn("w:val"), "clear")
    sh.set(qn("w:color"), "auto")
    sh.set(qn("w:fill"), "EDF1F7")
    tcPr.append(sh)


def _repeat_header(row):
    trPr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    trPr.append(el)


def _fixed_layout(table):
    table.autofit = False
    lay = OxmlElement("w:tblLayout")
    lay.set(qn("w:type"), "fixed")
    table._tbl.tblPr.append(lay)


def add_table(doc, df, font_pt=8.0, widths=None):
    ncol = len(df.columns)
    t = doc.add_table(rows=0, cols=ncol)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    _fixed_layout(t)
    if widths is None or len(widths) != ncol:
        usable, w0 = 25.5, 3.0
        widths = [w0] + [(usable - w0) / max(1, ncol - 1)] * (ncol - 1)
    hdr = t.add_row()
    for j, col in enumerate(df.columns):
        c = hdr.cells[j]
        c.text = ""
        p = c.paragraphs[0]
        p.paragraph_format.space_after = Pt(0)
        r = p.add_run(str(col))
        r.bold = True
        r.font.size = Pt(font_pt)
        _shade(c)
    _repeat_header(hdr)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, v in enumerate(row):
            cells[j].text = ""
            p = cells[j].paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing = 1.0
            txt = "" if pd.isna(v) else str(v)
            r = p.add_run(txt)
            r.font.size = Pt(font_pt)
    for row in t.rows:
        for c, w in zip(row.cells, widths):
            c.width = Cm(w)
    return t


def build(spec):
    df = pd.read_csv(TAB / spec["src"], encoding="utf-8-sig")
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(10.5)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = Cm(29.7), Cm(21)
    for a in ("left_margin", "right_margin"):
        setattr(sec, a, Cm(2.0))
    sec.top_margin = sec.bottom_margin = Cm(1.8)

    cap = doc.add_paragraph()
    cap.paragraph_format.space_after = Pt(8)
    r = cap.add_run(spec["caption"].replace("**", ""))
    r.bold = True
    r.font.size = Pt(11)

    add_table(doc, df, font_pt=spec["font"], widths=spec["widths"])

    doc.add_paragraph()
    for fn in spec["footnotes"]:
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(2)
        rr = p.add_run(fn)
        rr.font.size = Pt(8.5)
        rr.font.color.rgb = RGBColor(0x33, 0x33, 0x33)

    # optional second block in the same file (Table S3 carries the genus assignments as part a
    # and the pooled sensitivity of the axis under those definitions as part b)
    extra = spec.get("extra")
    if extra:
        doc.add_paragraph()
        cap2 = doc.add_paragraph()
        cap2.paragraph_format.space_after = Pt(6)
        r2 = cap2.add_run(extra["caption"])
        r2.bold = True
        r2.font.size = Pt(10)
        df2 = pd.read_csv(TAB / extra["src"], encoding="utf-8-sig")
        if "transform" in extra:
            df2 = extra["transform"](df2)
        add_table(doc, df2, font_pt=extra.get("font", spec["font"]), widths=extra["widths"])
        doc.add_paragraph()
        for fn in extra["footnotes"]:
            p = doc.add_paragraph()
            p.paragraph_format.space_after = Pt(2)
            rr = p.add_run(fn)
            rr.font.size = Pt(8.5)
            rr.font.color.rgb = RGBColor(0x33, 0x33, 0x33)

    OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT / f"{spec['name']}.docx"
    doc.save(dst)
    return dst, df.shape


def main():
    for spec in SPECS:
        dst, shape = build(spec)
        print(f"{dst.name:14} {shape[0]:3} rows x {shape[1]:2} cols -> {dst.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
