#!/usr/bin/env python3
"""Render a Results draft to .docx for circulation, leaving the markdown untouched.

Why a separate renderer
-----------------------
The caption builder (build_revised_kit_captions.py) renders captions from markdown without
rewrapping the source, and this reuses its block/italic parser for the same reason: the author's
line breaks are deliberate and a renderer must not rewrite them.

What it does that a plain paste would not
-----------------------------------------
- Headings become real Word headings, so the document has a navigable outline.
- Italic taxon names (*Parvimonas micra*), which are meaningful in this manuscript, are preserved
  rather than left as asterisks.
- First-line indent and the journal's line spacing are set so the draft reads like a manuscript
  page instead of a wall of text.
- A footer records the source file and its checksum, so a circulated copy can be traced back to the
  exact draft it came from.

Usage
-----
    python analysis/scripts/build_results_docx.py [draft.md] [-o out.docx]

Default source is the latest revised draft; default output sits beside it.
"""
from __future__ import annotations

import hashlib
import importlib.util
import pathlib
import sys

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Cm, Pt

ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_SRC = (ROOT / "submission" / "drafting-kit-revised" / "drafts"
               / "Results_JOM_v4_2026-11-16.md")
BUILDER = ROOT / "analysis" / "scripts" / "build_captions_docx.py"


def load_parser():
    """Reuse the caption builder's markdown block/italic parser."""
    spec = importlib.util.spec_from_file_location("capmod", BUILDER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def add_runs_with_code(mod, par, text: str) -> None:
    """Emit a paragraph, giving inline code monospace and bold/italic their meaning.

    The caption builder's add_runs knows only ** and *. A Methods draft is full of parameter names in
    backticks (pool = FALSE, /cases, survival), and passing them straight through printed the
    backticks and let the asterisk in a glob such as analysis/modules/*/README.md start a spurious
    italic run. Segmenting on backticks first fixes both, and identifiers read as code in Word.
    """
    for i, chunk in enumerate(text.split("`")):
        if not chunk:
            continue
        if i % 2 == 1:                                  # inside a code span
            run = par.add_run(chunk)
            run.font.name = "Consolas"
            run.font.size = Pt(10)
        else:
            mod.add_runs(par, chunk)


def parse_blocks(text: str):
    """Yield ('h1'|'h2'|'h3'|'p', text), one entry per markdown block.

    The caption builder handles only '#' and '##'. A Methods draft uses '###' for its subsections, and
    reusing that parser silently turned nine headings into body paragraphs, which is why this walker
    handles three levels instead.
    """
    buf: list[str] = []
    for line in text.split("\n"):
        s = line.strip()
        if not s or s == "---":
            if buf:
                yield "p", " ".join(buf)
                buf = []
            continue
        level = 0
        for ch in s:
            if ch == "#":
                level += 1
            else:
                break
        if level and s[level:level + 1] == " " and level <= 3:
            if buf:
                yield "p", " ".join(buf)
                buf = []
            yield f"h{level}", s[level:].strip()
        else:
            buf.append(s)
    if buf:
        yield "p", " ".join(buf)


def main() -> None:
    args = sys.argv[1:]
    out = None
    if "-o" in args:
        k = args.index("-o")
        out = pathlib.Path(args[k + 1])
        del args[k:k + 2]
    src = pathlib.Path(args[0]) if args else DEFAULT_SRC
    if not src.is_file():
        sys.exit(f"source not found: {src}")
    out = out or src.with_suffix(".docx")

    mod = load_parser()
    text = src.read_text(encoding="utf-8")          # read only; the markdown is never rewritten
    parablock = list(parse_blocks(text))

    doc = Document()
    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(12)

    for kind, body in parablock:
        if kind in ("h1", "h2", "h3"):
            h = doc.add_heading(body, level=int(kind[1]))
            if kind == "h1":
                h.alignment = WD_ALIGN_PARAGRAPH.LEFT
            continue
        p = doc.add_paragraph()
        p.paragraph_format.first_line_indent = Cm(0.75)
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.line_spacing = 2.0
        add_runs_with_code(mod, p, body)

    sha = hashlib.sha256(src.read_bytes()).hexdigest()[:12]
    foot = doc.sections[0].footer.paragraphs[0]
    foot.text = (f"Rendered from {src.name} (sha256 {sha}) — markdown source unchanged; "
                 f"this copy is for reading, the .md remains the editable original.")
    foot.runs[0].font.size = Pt(7)

    doc.save(str(out))
    n_h2 = sum(1 for k, _ in parablock if k in ("h2", "h3"))
    n_p = sum(1 for k, _ in parablock if k == "p")
    try:
        shown = out.resolve().relative_to(ROOT)
    except ValueError:
        shown = out
    print(f"wrote {shown} — {n_h2} sections, {n_p} paragraphs, source {src.name} sha {sha}")


if __name__ == "__main__":
    main()
