#!/usr/bin/env python3
"""Assemble every manuscript table into one Word file for review.

Layout choices:
  * landscape throughout - most of these tables are 10-25 columns wide and would be
    unreadable in portrait
  * 8 pt body, 7 pt for the very wide supplementary tables, header row bold + shaded
  * header row repeats on every page so a table can be read after a page break
  * missing values print as an em dash rather than "nan"
  * each table carries its own caption and its own footnote listing the abbreviations
    used in that table, so the file can be read without the manuscript next to it
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[4]
TAB = ROOT / "results" / "tables"
OUT = ROOT / "results" / "reports" / "tables-review.docx"

# (file, table number, caption, footnote, body font size)
SPEC = [
    ("table1_cohort_characteristics.csv", "Table 1",
     "Characteristics of the six public oral squamous cell carcinoma 16S cohorts.",
     "Runs, patients and niche counts come from the sample manifest; platform, "
     "trimLeft, truncLen, maxEE, pool and the read cap come from cohort_params.yaml; "
     "median read length and the 5-prime primer fraction are measured from the raw reads "
     "and are given per read, because R1 and R2 differ (for example 85 and 94 % in "
     "PRJNA1415359). PRJNA666746 carries a partial 806R match that its submitter did "
     "not trim, so it is labelled rather than reported as clean. ASV, genus and read "
     "counts come from the genus tables. Where a cell holds several settings it is one "
     "merged field, and the units are in the header: 'Platform / layout' pairs the "
     "sequencer with single- or paired-end, and 'DADA2 settings' lists trim, trunc, "
     "maxEE, pool and the per-sample read cap in that order. The number in brackets "
     "under Region is how many runs actually carried a country value in the ENA "
     "metadata; two cohorts leave it unset, so their country comes from the cohort "
     "notes instead. ASV, amplicon sequence variant; maxEE, maximum expected "
     "errors; pool, DADA2 pooling mode.", 6.8),

    ("table2_replicating_genera.csv", "Table 2",
     "Genera that replicate across cohorts: tested in at least three cohorts with "
     "FDR q < 0.05 and a consistent direction of effect.",
     "Pooled log2 fold change is tumour versus normal from a random-effects model on "
     "Cohen's d_z. Fifteen genera reach q < 0.05; two are excluded here - Collinsella "
     "was testable in only two cohorts and Olsenella has an inconsistent direction. "
     "Supplementary Figure S3 shows all of them, and the full statistics for all 615 "
     "genera ship as data (s5_meta_analysis.csv), so the exclusion can be audited. "
     "|d_z| >= 0.8 is the conventional large-effect threshold; no genus reaches it. "
     "FDR, false discovery rate; I2, proportion of variance from heterogeneity.", 7.5),

    ("tableS1_tcma_sample_annotation.csv", "Supplementary Table S1",
     "Sample-level annotation for the 177 TCGA-HNSC microbiomes in The Cancer "
     "Microbiome Atlas.",
     "Tissue type comes from the TCMA export, covariates from the GDC API. HPV status "
     "is available for a subset only: ISH for 53 samples and p16 for a different subset, "
     "which is why the two are separate columns. OS time is available for 156 of 177 "
     "samples (134 patients after de-duplication, 73 deaths). alpha-diversity is "
     "computed on the genus table and converted to percentages, because the source "
     "column is a proportion of the community; Parvimonas is given as a percentage of "
     "the bacterial community. OS, overall survival; ISH, in situ hybridisation; "
     "HPV, human papillomavirus.", 5.8),

    ("tableS2_survival_models.csv", "Supplementary Table S2",
     "Clinical associations and overall-survival models for the three microbiome "
     "metrics, in one long-format table.",
     "Five model families share the same three metrics, so they are stacked in one "
     "table with a 'Model' column rather than split across five tables. Stage is "
     "entered ordinally (1-4) in the primary model: as a factor it hits complete "
     "separation at stage I, where no deaths occur, so the hazard ratio is not "
     "estimable. One block repeats the model with stage binarised as a sensitivity "
     "analysis. No clinical association survives FDR (smallest q = 0.22) and no "
     "microbiome metric improves the clinical model. HR, hazard ratio; "
     "LRT, likelihood-ratio test; AIC, Akaike information criterion.", 6.2),
]


def shade(cell, hex_fill: str) -> None:
    tcPr = cell._tc.get_or_add_tcPr()
    el = OxmlElement("w:shd")
    el.set(qn("w:val"), "clear")
    el.set(qn("w:fill"), hex_fill)
    tcPr.append(el)


def repeat_header(row) -> None:
    trPr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    trPr.append(el)


def fmt(v) -> str:
    """Format a cell: 3 significant decimals, scientific for very small p/q values."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "\u2014"
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, (int, float)):
        v = float(v)
        if v == int(v) and abs(v) < 1e15:   # counts stored as floats
            return f"{int(v)}"
        a = abs(v)
        if a < 1e-3:                 # p and q values
            return f"{v:.2e}"
        if a < 1:
            return f"{v:.3f}"
        if a < 100:
            return f"{v:.2f}"
        return f"{v:.0f}"
    s = str(v).strip()
    if s.lower() in ("nan", "none", ""):
        return "\u2014"
    return _wrappable(s)


def _wrappable(text: str) -> str:
    """Insert break opportunities into long runs joined by separators.

    Cell values such as "PRJNA1415359;PRJNA1415268;PRJNA597251;PRJNA666746" are a single
    49-character token, so Word cannot wrap them and instead splits an accession number
    down the middle. A space after each separator gives Word a place to break; the CSV
    is untouched, so downstream users still see the original string.
    """
    # keep explicit line breaks; only collapse the doubles the separator pass creates
    text = text.replace("\n", "\x00")
    for sep in (";", "|", ","):
        text = text.replace(sep, sep + " ")
    parts = [" ".join(p.split()) for p in text.split("\x00")]
    return "\n".join(parts)


def add_table(doc, df: pd.DataFrame, num: str, caption: str, footnote: str,
              body_pt: float) -> None:
    h = doc.add_paragraph()
    r = h.add_run(f"{num}. {caption}")
    r.bold = True
    r.font.size = Pt(9)
    h.paragraph_format.space_before = Pt(14)
    h.paragraph_format.space_after = Pt(4)
    h.paragraph_format.keep_with_next = True

    tbl = doc.add_table(rows=1, cols=df.shape[1])
    tbl.style = "Table Grid"
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    tbl.autofit = False
    tbl.allow_autofit = False

    sec = doc.sections[0]
    usable_cm = sec.page_width.cm - sec.left_margin.cm - sec.right_margin.cm
    cell_pad = 0.06                     # cm each side; the Word default is 0.19
    tblPr = tbl._tbl.tblPr
    mar = OxmlElement("w:tblCellMar")
    for side in ("left", "right"):
        el = OxmlElement(f"w:{side}")
        el.set(qn("w:w"), str(int(cell_pad / 2.54 * 1440)))
        el.set(qn("w:type"), "dxa")
        mar.append(el)
    tblPr.append(mar)
    widths = column_widths(df, usable_cm, body_pt)
    for j, w in enumerate(widths):
        tbl.columns[j].width = Cm(w)

    hdr = tbl.rows[0]
    repeat_header(hdr)
    for j, col in enumerate(df.columns):
        cell = hdr.cells[j]
        cell.width = Cm(widths[j])
        cell.text = ""
        p = cell.paragraphs[0]
        p.paragraph_format.space_after = Pt(0)
        p.paragraph_format.line_spacing = 1.0
        run = p.add_run(header_text(col))
        run.bold = True
        run.font.size = Pt(body_pt)
        shade(cell, "E8EEF4")

    for _, row in df.iterrows():
        cells = tbl.add_row().cells
        # a row must not break across pages: a split cell is unreadable and leaves a
        # repeated header with a single orphaned fragment under it
        no_split(tbl.rows[-1])
        for j, col in enumerate(df.columns):
            cells[j].width = Cm(widths[j])
            cells[j].text = ""
            p = cells[j].paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing = 1.0
            run = p.add_run(fmt(row[col]))
            run.font.size = Pt(body_pt)

    f = doc.add_paragraph()
    fr = f.add_run(footnote)
    fr.italic = True
    fr.font.size = Pt(7)
    fr.font.color.rgb = RGBColor(0x44, 0x44, 0x44)
    f.paragraph_format.space_before = Pt(3)
    f.paragraph_format.space_after = Pt(10)


def header_text(col) -> str:
    """Make a column name wrappable without changing its meaning.

    Machine names such as site_group or hpv_p16 have no spaces, so Word has no place
    to break them and instead splits mid-token ("site_grou / p"). Replacing the
    separators with spaces lets the header wrap at a sensible point; the CSV keeps the
    original column name, so nothing downstream changes.
    """
    t = str(col).replace('_', ' ').replace('.', ' ')
    return ' '.join(t.split())


def _advance_cm(text: str, pt: float, bold: bool = False) -> float:
    """Measured width of a string in cm at the given size.

    Counting characters is not accurate enough: "Bifidobacterium" is 15 characters but
    needs 2.27 cm at 7.5 pt, and a column sized by character count wraps its last letter
    onto a second line. Matplotlib ships a font that reports real advance widths, so the
    measurement is taken from the glyph outlines instead.
    """
    try:
        from matplotlib.font_manager import FontProperties
        from matplotlib.textpath import TextPath
        fp = FontProperties(family="DejaVu Sans", size=pt,
                            weight="bold" if bold else "normal")
        return TextPath((0, 0), text, prop=fp).get_extents().width / 72 * 2.54
    except Exception:
        return 0.20 * len(text) * pt / 7.0


def column_widths(df: pd.DataFrame, usable_cm: float, body_pt: float) -> list[float]:
    """Divide the page width so that no header word or value is split mid-token.

    Word breaks inside a word when the cell is narrower than that word, which is how
    "ASVs" became "ASV / s" and "PRJNA1415359" became "PRJNA1415353 / 59". Each column
    therefore needs a floor equal to its longest *unbreakable* token: a numeric value is
    unbreakable as a whole, while free text may wrap at its longest word.

    Widths are built as floor + share of the leftover space, so the total is exactly the
    usable width and no column ends up below the width its content demands.
    """
    pad = 2 * 0.06 + 0.08          # cell margins (set below) plus a hair of slack
    floors: list[float] = []
    weight: list[float] = []

    for col in df.columns:
        vals = [fmt(v) for v in df[col].tolist()]
        values = [v for v in vals if v and v != "\u2014"]

        numeric = bool(values) and all(
            len(v) <= 12 and re.fullmatch(r"-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?", v)
            for v in values)
        if numeric:
            token = max(values, key=len)     # a number cannot wrap at all
        else:
            words = [w for v in vals for w in str(v).replace("\n", " ").split()]
            token = max(words, key=len) if words else ""

        need_val = _advance_cm(token, body_pt, bold=False)
        # a header wraps at its spaces, so only its longest WORD must fit; measuring the
        # whole header inflated the floor by 3-4x for labels like
        # "Largest |d_z| in any cohort" and starved every other column
        hdr_words = header_text(col).split() or [""]
        need_hdr = max(_advance_cm(w, body_pt, bold=True) for w in hdr_words)
        floors.append(max(need_val, need_hdr) * 1.02 + pad)

        # extra space goes where content is long, but capped so one verbose column
        # cannot swallow the table
        weight.append(max(1.0, min(max((len(v) for v in vals), default=1), 12)))

    base = float(sum(floors))
    if base >= usable_cm:
        # the table genuinely cannot honour every floor at this font size; scale
        # proportionally and let the text columns wrap
        k = usable_cm / base
        return [f * k for f in floors]

    extra = usable_cm - base
    wsum = float(sum(weight))
    return [f + extra * w / wsum for f, w in zip(floors, weight)]


def no_split(row) -> None:
    trPr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:cantSplit")
    trPr.append(el)



def main() -> None:
    doc = Document()

    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = Cm(29.7), Cm(21.0)
    for m in ("left_margin", "right_margin"):
        setattr(sec, m, Cm(1.4))
    sec.top_margin = sec.bottom_margin = Cm(1.6)

    style = doc.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(9)

    t = doc.add_paragraph()
    tr = t.add_run("Tables for the intratumoral microbiome study")
    tr.bold = True
    tr.font.size = Pt(15)

    sub = doc.add_paragraph()
    sr = sub.add_run(
        "Two main tables and two supplementary tables. Every value was written by a "
        "script from the analysis outputs in results/tables; nothing is hand-typed. "
        "Most analyses are already shown as figure panels, so only the material a "
        "reader must look up is tabulated. Anything that could be shown as a figure "
        "is a figure (see the supplementary figure series) and is not repeated here.")
    sr.font.size = Pt(8.5)
    sr.italic = True
    sr.font.color.rgb = RGBColor(0x33, 0x33, 0x33)

    index = doc.add_paragraph()
    ir = index.add_run("Contents")
    ir.bold = True
    ir.font.size = Pt(10)
    index.paragraph_format.space_before = Pt(10)

    for f, num, cap, _fn, _pt in SPEC:
        p = doc.add_paragraph(style="List Bullet")
        run = p.add_run(f"{num}  ")
        run.bold = True
        run.font.size = Pt(8.5)
        df = pd.read_csv(TAB / f)
        tail = f"  ({len(df)} rows x {df.shape[1]} columns)"
        run2 = p.add_run(cap + tail)
        run2.font.size = Pt(8.5)

    usable = sec.page_width.cm - sec.left_margin.cm - sec.right_margin.cm
    for f, num, cap, fn, pt in SPEC:
        df = pd.read_csv(TAB / f)
        add_table(doc, df, num, cap, fn, pt)

    # page numbers, so a printed copy can be referenced
    footer_p = sec.footer.paragraphs[0]
    footer_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = footer_p.add_run()
    run.font.size = Pt(7.5)
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    run._r.addnext(fld)

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}  "
          f"({OUT.stat().st_size/1024:.0f} KB)")
    for f, num, cap, _fn, _pt in SPEC:
        df = pd.read_csv(TAB / f)
        print(f"  {num:22} {len(df):>4} x {df.shape[1]:>2} rows x cols")


if __name__ == "__main__":
    main()