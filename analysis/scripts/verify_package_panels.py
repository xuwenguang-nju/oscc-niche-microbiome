#!/usr/bin/env python3
"""Verify the figure/panel cross-references in the submission documents against the deliverables.

What it checks
--------------
1. `figure_manifest.csv` matches the files on disk: every row's file exists, its dpi is 600, its
   width is inside the 182 mm two-column limit, and its panel_letters string has the recorded count.
2. Every "Fig. N<letter>" / "Figure N" reference in the submission checklist names a figure that
   exists and a panel letter that figure actually has.
3. Every panel letter a caption describes exists in that figure, and every figure has 11 captions
   of two paragraphs (body + *Data:*), inside the house limits.
4. The caption .docx contains the same 11 captions as the markdown, with no Markdown residue.

It reads the delivered files, so it fails when a document is edited without the figures.
"""
from __future__ import annotations

import csv
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "results/tables/figure_manifest.csv"
CAPTIONS = ROOT / "results/reports/figure-captions-jom.md"
CHECKLIST = ROOT / "submission/jom-2026/08-checklist/submission-checklist.md"
WIDTH_LIMIT_MM = 182.0

fails: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        fails.append(msg)


def load_manifest():
    with open(MANIFEST, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    figs = {}
    for r in rows:
        # "Fig4-replication..." -> "4", "FigS1-niche..." -> "S1"
        figs[r["figure"].split("-")[0].replace("Fig", "")] = r
    return rows, figs


def check_manifest(rows, figs) -> None:
    check(len(rows) == 11, f"manifest has {len(rows)} rows, expected 11")
    main = [r for r in rows if r["set"] == "main"]
    supp = [r for r in rows if r["set"] == "supplementary"]
    check(len(main) == 7, f"{len(main)} main figures, manuscript uses 7")
    check(len(supp) == 4, f"{len(supp)} supplementary figures, manuscript uses 4")
    for r in rows:
        path = ROOT / "figures/manuscript-v6" / r["set"] / r["figure"]
        check(path.is_file(), f"manifest row has no file: {r['figure']}")
        check(r["dpi"] == "600", f"{r['figure']}: dpi {r['dpi']}, journal colour rule is 300+ and "
                                 f"this set is delivered at 600")
        check(float(r["width_mm"]) <= WIDTH_LIMIT_MM,
              f"{r['figure']}: {r['width_mm']} mm wide, limit {WIDTH_LIMIT_MM}")
        letters = r["panel_letters"].split(",")
        check(len(letters) == int(r["panels"]),
              f"{r['figure']}: {len(letters)} letters for {r['panels']} panels")
        # the letters must start at A and be consecutive, or a figure's lettering has drifted
        check(letters == [chr(65 + i) for i in range(len(letters))],
              f"{r['figure']}: panel letters not consecutive from A: {r['panel_letters']}")


def figure_key(token: str) -> str:
    """Manifest keys are already '1'..'7' and 'S1'..'S4'."""
    return token.upper() if token.upper().startswith("S") else token


def check_references(figs) -> None:
    text = CHECKLIST.read_text(encoding="utf-8")
    # "Fig. 4L", "Fig. 1A", "Figs. 6 and S3", "Fig. S1" — the letter is optional
    n = 0
    for m in re.finditer(r"Fig(?:s)?\.?\s+(S?\d+)([A-N])?\b", text):
        fig, letter = m.group(1), m.group(2)
        key = figure_key(fig)
        check(key in figs, f"checklist cites Fig. {fig}, which is not in the manifest")
        if key in figs and letter:
            check(letter in figs[key]["panel_letters"].split(","),
                  f"checklist cites Fig. {fig}{letter}, but that figure has "
                  f"{figs[key]['panels']} panels ({figs[key]['panel_letters']})")
        n += 1
    print(f"checklist: {n} figure references checked")


def check_captions(figs) -> None:
    blocks = [b.strip() for b in CAPTIONS.read_text(encoding="utf-8").split("\n\n")]
    caps: dict[str, list[str]] = {}
    cur = None
    for b in blocks:
        flat = " ".join(b.split())
        if flat.lstrip("*").startswith("Figure "):
            cur = flat.split(".")[0].replace("*", "").strip().split()[1]
            caps[cur] = [flat]
            # every panel letter the caption names must exist in that figure
            key = figure_key(cur)
            check(key in figs, f"caption for Figure {cur} has no manifest row")
            if key in figs:
                have = set(figs[key]["panel_letters"].split(","))
                named = set(re.findall(r"\(\*\*([A-N])\*\*\)", flat))
                check(named <= have,
                      f"Figure {cur} caption names panels {sorted(named - have)} that the figure "
                      f"does not have")
                check(named == have or len(have) > 14,
                      f"Figure {cur} caption covers {len(named)} of {len(have)} panels")
        elif flat.startswith("*Data:") and cur:
            caps[cur].append(flat)

    check(len(caps) == 11, f"{len(caps)} captions, expected 11")
    for k, v in caps.items():
        check(len(v) == 2, f"Figure {k}: {len(v)} paragraphs, expected body + Data")
        title = v[0].split(".**", 1)[0].replace("**", "")
        n_title = len(" ".join(title.split()[2:]).split())
        n_body = len(re.sub(r"[*]", "", " ".join(v)).split())
        check(n_title <= 20, f"Figure {k}: title {n_title} words, house limit 20")
        # This file is the master, in which the Data: note is a separate paragraph. The submission
        # text merges it into the prose and is held to 300 words there (340 for the fourteen-panel
        # Figure 5); build_submission_captions.py enforces that. Here the cap only catches runaways.
        check(n_body <= 340, f"Figure {k}: {n_body} words, over the 340 master cap")
    print(f"captions: {len(caps)} checked, all panel letters present in their figures")


def check_docx(caps_count: int = 11) -> None:
    from docx import Document

    d = Document(CAPTIONS.with_suffix(".docx"))
    # headings are style Heading 2, so counting style Normal finds only the caption paragraphs
    paras = [p.text for p in d.paragraphs if p.style.name == "Normal"]
    n = sum(1 for t in paras if t.startswith("Figure "))
    check(n == caps_count, f"caption .docx has {n} captions, markdown has {caps_count}")
    residue = [t for t in paras if "**" in t]
    check(not residue, f"caption .docx carries Markdown residue in {len(residue)} paragraphs")
    print(f"caption .docx: {n} captions, no Markdown residue")


TABDIR = ROOT / "submission/jom-2026/04-tables"


def check_tables() -> None:
    """The table inventory must be 1 main + N supplementary, numbered without gaps.

    Added when the rank-level table (Table S4) was inserted ahead of KEGG, which pushed the KEGG
    table to S5 and the TCGA annotation to S6. Without this, a renamed file and a stale inventory
    entry can disagree while every other check still passes.
    """
    present = sorted(p.stem for p in TABDIR.glob("*.docx"))
    mains = [s for s in present if s.startswith("Table") and not s.startswith("TableS")]
    supp = sorted((s for s in present if s.startswith("TableS")),
                  key=lambda s: int(s.replace("TableS", "")))
    check(mains == ["Table1"], f"main tables on disk: {mains}, expected ['Table1']")
    check([int(s.replace("TableS", "")) for s in supp] == list(range(1, len(supp) + 1)),
          f"supplementary table numbers are not 1..N without gaps: {supp}")
    # each docx must carry its own number in the caption, and open with that number
    from docx import Document

    for stem in present:
        d = Document(TABDIR / f"{stem}.docx")
        cap = next((p.text for p in d.paragraphs if p.text.strip().startswith("Table")), "")
        expect = f"Table {stem.replace('Table', '')}."   # TableS3 -> "Table S3."
        check(cap.startswith(expect), f"{stem}.docx caption starts {cap[:40]!r}, expected {expect!r}")
        check((TABDIR / f"{stem}.csv").is_file(), f"{stem}.docx has no source CSV beside it")
    print(f"tables: {len(mains)} main + {len(supp)} supplementary, captions and CSV sources present")


def main() -> None:
    rows, figs = load_manifest()
    check_manifest(rows, figs)
    check_references(figs)
    check_captions(figs)
    check_docx()
    check_tables()

    if fails:
        print("\nFAILURES")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
