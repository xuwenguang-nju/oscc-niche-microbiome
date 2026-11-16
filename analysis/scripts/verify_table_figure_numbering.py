#!/usr/bin/env python3
"""Check that tables and figures are numbered by first citation, and that the files agree.

The rule this enforces is the journal's: every table and figure is numbered in the order of its
first mention in the main text, main items first and then the supplement in its own sequence. A
table that is only mentioned in a "*Tables:*" footer does not count as cited — the reader has to
meet it in a sentence.

Sources of truth
----------------
citation order : results/reports/manuscript-methods-v1.md then manuscript-results-v7.md
                 (the section files, in manuscript order; the assembled .docx predates them)
table files    : submission/jom-2026/04-tables/*.csv and *.docx
figure files   : results/tables/figure_manifest.csv

Run: python analysis/scripts/verify_table_figure_numbering.py
"""
from __future__ import annotations

import csv
import pathlib
import re
import sys

from docx import Document

ROOT = pathlib.Path(__file__).resolve().parents[2]
# the manuscript sections in document order: Methods, then the Results being submitted. A new
# Results draft is picked up by name, so the check always runs against the file in play.
_RESULTS = next((ROOT / "results/reports" / f"manuscript-results-v{n}.md"
                 for n in range(20, 7, -1)
                 if (ROOT / "results/reports" / f"manuscript-results-v{n}.md").is_file()),
                ROOT / "results/reports/manuscript-results-v7.md")
TEXT = [ROOT / "results/reports/manuscript-methods-v1.md", _RESULTS]
TABDIR = ROOT / "submission/jom-2026/04-tables"
MANIFEST = ROOT / "results/tables/figure_manifest.csv"

fails: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        fails.append(msg)


def citations() -> list[tuple[str, int]]:
    """[(token, line_number)] in document order, body text only."""
    pat = re.compile(r"Table\s+(S?)(\d+)|Fig(?:s)?\.?\s*(S?)(\d+)")
    out = []
    for path in TEXT:
        if not path.is_file():
            raise SystemExit(f"missing text source: {path}")
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            s = line.strip()
            if s.startswith("*Tables:*") or s.startswith(">") or s.startswith("#"):
                continue                      # footers, header notes, section headings
            for m in pat.finditer(line):
                if m.group(2):
                    out.append((f"Table {m.group(1)}{m.group(2)}", i))
                else:
                    out.append((f"Fig. {m.group(3)}{m.group(4)}", i))
    return out


def first_mentions(cites, prefix: str) -> list[str]:
    seen: list[str] = []
    for tok, _ in cites:
        if tok.startswith(prefix) and tok not in seen:
            seen.append(tok)
    return seen


def main() -> None:
    cites = citations()

    # ---------------------------------------------------------------- tables
    tables = first_mentions(cites, "Table")
    main_tables = [t for t in tables if not t.startswith("Table S")]
    supp_tables = [t for t in tables if t.startswith("Table S")]
    check(main_tables == ["Table 1"],
          f"main tables cited: {main_tables}; exactly one main table is allowed")
    check(supp_tables == [f"Table S{i}" for i in range(1, len(supp_tables) + 1)],
          f"supplementary tables are not in ascending first-citation order: {supp_tables}")

    files = sorted(p.stem for p in TABDIR.glob("*.csv"))
    check(files == sorted(t.replace(" ", "") for t in tables),
          f"table files on disk {files} do not match the tables cited {tables}")
    check(not (TABDIR / "Table2.csv").exists(),
          "04-tables/Table2.csv still exists: the second main table was supposed to move to the "
          "supplement")

    # each .docx caption must open with the table's own number
    for stem in files:
        docx = TABDIR / f"{stem}.docx"
        check(docx.is_file(), f"{stem}.docx missing")
        if not docx.is_file():
            continue
        first = next((p.text.strip() for p in Document(docx).paragraphs if p.text.strip()), "")
        want = stem.replace("TableS", "Table S").replace("Table1", "Table 1")
        check(first.startswith(want + "."),
              f"{stem}.docx caption starts with {first[:40]!r}, expected it to start with {want!r}")

    # ---------------------------------------------------------------- figures
    with open(MANIFEST, encoding="utf-8") as fh:
        manifest = {r["figure"].split("-")[0].replace("Fig", ""): r
                    for r in csv.DictReader(fh)}
    main_figs = first_mentions(cites, "Fig. ")
    main_seq = [f.replace("Fig. ", "") for f in main_figs if not f.startswith("Fig. S")]
    supp_seq = [f.replace("Fig. S", "") for f in main_figs if f.startswith("Fig. S")]

    check(main_seq == [str(i) for i in range(1, len(main_seq) + 1)],
          f"main figures are not cited in ascending order: {main_seq}")
    check(supp_seq == [str(i) for i in range(1, len(supp_seq) + 1)],
          f"supplementary figures are not cited in ascending order: {supp_seq}")

    on_disk_main = [k for k, r in manifest.items() if r["set"] == "main"]
    on_disk_supp = [k.replace("S", "") for k, r in manifest.items() if r["set"] == "supplementary"]
    check(sorted(main_seq) == sorted(on_disk_main),
          f"main figures cited {main_seq} do not match the manifest {on_disk_main}")
    check(sorted(supp_seq) == sorted(on_disk_supp),
          f"supplementary figures cited {supp_seq} do not match the manifest {on_disk_supp}")

    # ---------------------------------------------------------------- report
    print(f"sources       : {TEXT[0].name}, {TEXT[1].name}")
    print("table order   :", " -> ".join(tables))
    print("main figures  :", " -> ".join(main_seq))
    print("supp figures  :", " -> ".join("S" + s for s in supp_seq))
    print(f"table files   : {', '.join(files)}")

    if fails:
        print("\nFAILURES")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print("\nNUMBERING OK — first-citation order holds for every table and figure")


if __name__ == "__main__":
    main()
