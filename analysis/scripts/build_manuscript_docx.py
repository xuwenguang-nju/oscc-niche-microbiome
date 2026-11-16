#!/usr/bin/env python3
"""Render a Markdown manuscript to .docx for review.

Used for `results/reports/manuscript-v2-assembled.md` -> `.docx`. Handles the subset of
Markdown these drafts use: ATX headings, bold/italic/code inline spans, sub/superscript
characters, pipe tables with a header separator, blockquotes, and dash bullets.

Run: python analysis/scripts/build_manuscript_docx.py <in.md> <out.docx>
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

SUP = {"\u2070": "0", "\u00b9": "1", "\u00b2": "2", "\u00b3": "3", "\u2074": "4",
       "\u2075": "5", "\u2076": "6", "\u2077": "7", "\u2078": "8", "\u2079": "9",
       "\u207a": "+", "\u207b": "\u2212", "\u207f": "n"}
SUB = {"\u2080": "0", "\u2081": "1", "\u2082": "2", "\u2083": "3", "\u2084": "4",
       "\u2085": "5", "\u2086": "6", "\u2087": "7", "\u2088": "8", "\u2089": "9"}

H_RE = re.compile(r"^(#{1,6})\s+(.*)$")
TABLE_SEP = re.compile(r"^\|[-: |]+\|$")
COMB = re.compile(r"\*\*(?P<bold>.+?)\*\*|`(?P<code>[^`]+)`"
                   r"|<sub>(?P<sub>[^<]+)</sub>|<sup>(?P<sup>[^<]+)</sup>"
                   r"|(?<![\*\w])\*(?P<ital>[^*\n]+?)\*(?!\*)")


def _plain(par, text, bold, ital, mono, sub=False, sup=False):
    r = par.add_run(text)
    r.bold, r.italic = bold, ital
    if mono:
        r.font.name = "Consolas"
        r.font.size = Pt(9)
    if sub:
        r.font.subscript = True
    if sup:
        r.font.superscript = True
    return r


def add_runs(par, text, bold=False, ital=False, mono=False):
    """Emit inline spans, converting U+207x / U+208x runs into real super/subscripts."""
    pos = 0
    for m in COMB.finditer(text):
        if m.start() > pos:
            _emit(par, text[pos:m.start()], bold, ital, mono)
        if m.group("bold"):
            add_runs(par, m.group("bold"), True, ital, mono)
        elif m.group("code"):
            _plain(par, m.group("code"), bold, ital, True)
        elif m.group("sub"):
            _plain(par, m.group("sub"), bold, ital, mono, sub=True)
        elif m.group("sup"):
            _plain(par, m.group("sup"), bold, ital, mono, sup=True)
        elif m.group("ital"):
            add_runs(par, m.group("ital"), bold, True, mono)
        pos = m.end()
    if pos < len(text):
        _emit(par, text[pos:], bold, ital, mono)


def _emit(par, text, bold, ital, mono):
    buf, mode = "", None
    for ch in text:
        new = "sup" if ch in SUP else ("sub" if ch in SUB else None)
        if new != mode:
            if buf:
                _plain(par, buf, bold, ital, mono, sub=(mode == "sub"), sup=(mode == "sup"))
            buf, mode = "", new
        buf += SUP.get(ch, SUB.get(ch, ch))
    if buf:
        _plain(par, buf, bold, ital, mono, sub=(mode == "sub"), sup=(mode == "sup"))


def style_doc(d):
    s = d.styles["Normal"]
    s.font.name = "Times New Roman"
    s.font.size = Pt(11)
    s.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    pf = s.paragraph_format
    pf.line_spacing = 1.5
    pf.space_after = Pt(6)
    for name, size in (("Heading 1", 16), ("Heading 2", 13.5), ("Heading 3", 11.5)):
        st = d.styles[name]
        st.font.name = "Times New Roman"
        st.font.size = Pt(size)
        st.font.color.rgb = RGBColor(0, 0, 0)
        st.font.bold = True
    sec = d.sections[0]
    sec.page_width, sec.page_height = Cm(21), Cm(29.7)
    for a in ("left_margin", "right_margin"):
        setattr(sec, a, Cm(2.3))
    sec.top_margin = sec.bottom_margin = Cm(2.0)


def set_widths(t, widths):
    t.autofit = False
    tblPr = t._tbl.tblPr
    lay = OxmlElement("w:tblLayout")
    lay.set(qn("w:type"), "fixed")
    tblPr.append(lay)
    for row in t.rows:
        for c, w in zip(row.cells, widths):
            c.width = Cm(w)


def shade(cell):
    tcPr = cell._tc.get_or_add_tcPr()
    sh = OxmlElement("w:shd")
    sh.set(qn("w:val"), "clear")
    sh.set(qn("w:color"), "auto")
    sh.set(qn("w:fill"), "EDF1F7")
    tcPr.append(sh)


def build(md_path, out_path, table_widths=None):
    lines = Path(md_path).read_text(encoding="utf-8").splitlines()
    d = Document()
    style_doc(d)
    i, tno = 0, 0
    while i < len(lines):
        raw = lines[i]
        s = raw.strip()
        if not s:
            i += 1
            continue
        if s.startswith("<!--"):
            i += 1
            continue
        m = H_RE.match(s)
        if m:
            lvl = min(len(m.group(1)), 3)
            p = d.add_heading("", level=lvl)
            add_runs(p, m.group(2))
            i += 1
            continue
        if s.startswith("|") and i + 1 < len(lines) and TABLE_SEP.match(lines[i + 1].strip()):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                if not TABLE_SEP.match(lines[i].strip()):
                    rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            ncol = max(len(r) for r in rows)
            t = d.add_table(rows=0, cols=ncol)
            t.style = "Table Grid"
            t.alignment = WD_TABLE_ALIGNMENT.CENTER
            for k, r in enumerate(rows):
                cells = t.add_row().cells
                for j in range(ncol):
                    cells[j].text = ""
                    para = cells[j].paragraphs[0]
                    para.paragraph_format.space_after = Pt(1)
                    para.paragraph_format.line_spacing = 1.0
                    add_runs(para, r[j] if j < len(r) else "", bold=(k == 0))
                    if k == 0:
                        shade(cells[j])
                    for run in para.runs:
                        run.font.size = Pt(8.5)
            t.rows[0]._tr.get_or_add_trPr().append(OxmlElement("w:tblHeader"))
            if table_widths and tno < len(table_widths):
                set_widths(t, table_widths[tno])
            tno += 1
            d.add_paragraph()
            continue
        if s.startswith(">"):
            buf = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                buf.append(lines[i].strip().lstrip("> ").strip())
                i += 1
            p = d.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.6)
            p.paragraph_format.space_before = Pt(6)
            add_runs(p, " ".join(x for x in buf if x))
            for run in p.runs:
                run.font.size = Pt(9.5)
                run.font.color.rgb = RGBColor(0x44, 0x44, 0x44)
                run.italic = True
            continue
        if re.match(r"^([-*]|\d+\.)\s", s):
            p = d.add_paragraph(style="List Bullet" if s[:1] in "-*" else None)
            add_runs(p, re.sub(r"^([-*]|\d+\.)\s+", "", s))
            i += 1
            continue
        buf = [s]
        i += 1
        while (i < len(lines) and lines[i].strip()
               and not H_RE.match(lines[i].strip())
               and not lines[i].strip().startswith(("|", ">", "<!--"))
               and not re.match(r"^([-*]|\d+\.)\s", lines[i].strip())):
            buf.append(lines[i].strip())
            i += 1
        p = d.add_paragraph()
        add_runs(p, " ".join(buf))
    d.save(out_path)
    return len(d.paragraphs), len(d.tables)


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    widths = [[1.0, 4.6, 3.6, 2.9], [1.2, 5.0, 3.2, 4.0]] if "--refs" in sys.argv else None
    print(build(src, dst, widths))
