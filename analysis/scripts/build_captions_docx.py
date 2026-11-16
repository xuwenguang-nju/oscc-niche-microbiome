#!/usr/bin/env python3
"""Rebuild results/reports/figure-captions-jom.docx from the .md, and rewrap the .md paragraphs.

The .docx is a rendering of the .md, not an independent document: the figure-title sentence is
bold, the panel letters stay bold, and the *Data:* paragraph is italic. Building it here keeps the
two in step, so an edited caption cannot ship as Word with the old text.

One paragraph in the .md = one paragraph in Word: consecutive non-blank lines are joined, so the
hard wraps used for reading the markdown do not become separate Word paragraphs (they did on the
first attempt, and the bold spans then broke at the wrap).
"""
import pathlib
import re
import textwrap

from docx import Document
from docx.shared import Pt

MD = pathlib.Path("results/reports/figure-captions-jom.md")
DOCX = MD.with_suffix(".docx")
WIDTH = 100

TOKEN = re.compile(r"\*\*(.+?)\*\*|\*([^*\n]+)\*")


def rewrap(md: str) -> str:
    """Re-wrap prose blocks at WIDTH; leave headings and rules alone."""
    out, buf = [], []

    def flush():
        if buf:
            joined = " ".join(x.strip() for x in buf)
            out.extend(textwrap.wrap(joined, width=WIDTH, break_long_words=False,
                                     break_on_hyphens=False) or [""])
            buf.clear()

    for line in md.split("\n"):
        s = line.strip()
        if not s or s.startswith("#") or s == "---":
            flush()
            out.append(line)
        else:
            buf.append(s)
    flush()
    return "\n".join(out) + "\n"


def blocks(md: str):
    """Yield ('h1'|'h2'|'p', text) — one entry per markdown block."""
    buf = []
    for line in md.split("\n"):
        s = line.strip()
        if not s or s == "---":
            if buf:
                yield "p", " ".join(buf)
                buf = []
            continue
        if s.startswith("## "):
            if buf:
                yield "p", " ".join(buf)
                buf = []
            yield "h2", s[3:]
        elif s.startswith("# "):
            if buf:
                yield "p", " ".join(buf)
                buf = []
            yield "h1", s[2:]
        else:
            buf.append(s)
    if buf:
        yield "p", " ".join(buf)


def add_runs(par, text: str) -> None:
    pos = 0
    for m in TOKEN.finditer(text):
        if m.start() > pos:
            par.add_run(text[pos:m.start()])
        if m.group(1) is not None:
            par.add_run(m.group(1)).bold = True
        else:
            par.add_run(m.group(2)).italic = True
        pos = m.end()
    if pos < len(text):
        par.add_run(text[pos:])


def main() -> None:
    md = rewrap(MD.read_text(encoding="utf-8"))
    MD.write_text(md, encoding="utf-8", newline="\n")

    doc = Document()
    doc.styles["Normal"].font.name = "Calibri"
    doc.styles["Normal"].font.size = Pt(10.5)

    paras = list(blocks(md))
    n_caps = sum(1 for kind, text in paras
                 if kind == "p" and text.lstrip("*").startswith("Figure "))
    assert n_caps == 11, f"expected 11 captions, found {n_caps}"

    for kind, text in paras:
        if kind == "h1":
            doc.add_heading(text, level=1)
        elif kind == "h2":
            doc.add_heading(text, level=2)
        else:
            add_runs(doc.add_paragraph(), text)

    doc.save(DOCX)
    print(f"wrote {DOCX} — {n_caps} captions")


if __name__ == "__main__":
    main()
