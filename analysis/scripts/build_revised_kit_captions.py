#!/usr/bin/env python3
"""Render the revised kit's figure-captions.md to .docx without touching the .md.

Why this exists rather than reusing build_captions_docx.py
---------------------------------------------------------
That script rewraps its markdown in place, which is fine for the file this project generates but
not for the author's revised caption file: reflowing their deliberate line breaks would be an
unrequested edit. This renders the same Word structure (bold title sentence, bold panel letters,
italics) from the markdown exactly as it stands.

Usage
-----
    python analysis/scripts/build_revised_kit_captions.py [kit_dir]

Defaults to submission/drafting-kit-revised.
"""
import importlib.util
import pathlib
import sys

from docx import Document
from docx.shared import Pt

ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_KIT = ROOT / "submission" / "drafting-kit-revised"
BUILDER = ROOT / "analysis" / "scripts" / "build_captions_docx.py"


def load_builder():
    spec = importlib.util.spec_from_file_location("capmod", BUILDER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    kit = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_KIT
    md_path = kit / "figure-captions.md"
    docx_path = kit / "figure-captions.docx"

    mod = load_builder()
    md = md_path.read_text(encoding="utf-8")          # read only; never rewritten
    paras = list(mod.blocks(md))
    n_caps = sum(1 for kind, text in paras
                 if kind == "p" and text.lstrip("*").startswith("Figure "))
    assert n_caps == 11, f"expected 11 captions in {md_path.name}, found {n_caps}"

    doc = Document()
    doc.styles["Normal"].font.name = "Calibri"
    doc.styles["Normal"].font.size = Pt(10.5)
    for kind, text in paras:
        if kind == "h1":
            doc.add_heading(text, level=1)
        elif kind == "h2":
            doc.add_heading(text, level=2)
        else:
            mod.add_runs(doc.add_paragraph(), text)
    doc.save(str(docx_path))
    print(f"wrote {docx_path.relative_to(ROOT)} — {n_caps} captions (source .md untouched)")


if __name__ == "__main__":
    main()
