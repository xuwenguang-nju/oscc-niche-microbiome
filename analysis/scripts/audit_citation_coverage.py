#!/usr/bin/env python3
"""Report which figures, supplementary figures and tables each section cites, and which are orphaned.

Why this exists
---------------
"Are all the figures and tables cited?" cannot be answered from one file: the manuscript is written as
four section files plus a figure-caption file, and a figure may be cited in Results while a table is
cited only in Methods. This walks every section file, records every Fig/Table token with the section it
came from, and prints a coverage matrix. It also reports panel-level coverage for the main figures,
because a figure whose panel N is never mentioned is a figure a reader stops reading.

A token counts only when it appears in body text: a "*Tables:*" footer, a "> " header note and a
"## Notes for the authors" block are excluded, since the journal's numbering rule is about the first
citation a reader meets in a sentence.

Run: python analysis/scripts/audit_citation_coverage.py
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
SEC = ROOT / "results/reports"

SECTIONS = [
    ("Introduction", "manuscript-introduction-v1.md"),
    ("Methods", "manuscript-methods-v1.md"),
    ("Results", "manuscript-results-v8.md"),
    ("Discussion", "manuscript-discussion-v1.md"),
]
CAPTIONS = ROOT / "results/reports/figure-captions-submission.md"

# Fig. 2A / Fig. 2A, 2B / Figs. 3-4 / Fig. S1 / Figure 6
def body_of(path: pathlib.Path) -> list[str]:
    """Lines that count as body text: no headers, no footers, no author notes."""
    out, skip = [], False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("## Notes for the authors"):
            skip = True
        if skip:
            continue
        s = line.strip()
        if not s or s.startswith(">") or s.startswith("*Tables:*") or s.startswith("*Figures:*"):
            continue
        out.append(line)
    return out


FIG = re.compile(r"Fig(?:ure)?s?\.?\s*(S?\d+)")
TAB = re.compile(r"Table\s+(S?\d+)")
ELEMENT = re.compile(r"^\s*(?:S?(\d+))?\s*([A-N])?\s*$")


def _run_after(line: str, start: int, first: str) -> tuple[list[tuple[str, str]], int]:
    """Parse a citation run such as "6G, 6H, 6I, 6J" or "3A, 3I" or "S1".

    Returns ([(figure, panel or "")], end index). A bare letter belongs to the figure named just
    before it; a bare number starts a new figure; "and" and en/em dashes separate elements, so
    "Figs. 3-4" and "Figs. 6 and 7" both parse.
    """
    out: list[tuple[str, str]] = []
    i = start
    cur = first
    after_dash = False
    while i < len(line):
        if line[i] in ",; \t":
            if line[i] in ",;":
                after_dash = False
            i += 1
            continue
        if line[i] in "\u2013\u2014-":
            after_dash = True
            i += 1
            continue
        # a range such as "7G-7I" names every letter between its ends
        if i < len(line) and line[i] in ",;":
            pass
        m = re.match(r"(?:and\s+)?(S?)(\d+)?([A-N])?", line[i:])
        if not m or (m.group(0).strip() in ("", "and")):
            break
        tok = m.group(0).strip()
        digits, letter = m.group(2), m.group(3)
        if digits:
            cur = ("S" if m.group(1) else "") + digits
            if after_dash and out and not m.group(3):
                # "Figs. 3-4": the end of a figure range; the start was recorded already
                pass
            if (letter and after_dash and out and out[-1][0] == cur and out[-1][1]):
                for code in range(ord(out[-1][1]) + 1, ord(letter)):
                    out.append((cur, chr(code)))
            out.append((cur, letter or ""))
        elif letter and cur:
            # extend a dash range: "G-7I" fills H, "7G-7I" fills H as well
            if after_dash and out and out[-1][0] == cur and out[-1][1]:
                lo, hi = out[-1][1], letter
                for code in range(ord(lo) + 1, ord(hi)):
                    out.append((cur, chr(code)))
            out.append((cur, letter))
        else:
            break
        after_dash = False
        i += m.end()
        # stop when the citation run is clearly over (a word character follows without a separator)
        if i < len(line) and line[i] not in ",; \u2013\u2014-" and line[i:i+4] != " and":
            break
    return out, i


def scan(path: pathlib.Path):
    """({fig: set(panels)}, {fig: occurrences}, {table: occurrences}) for one section file."""
    figs: dict[str, set[str]] = {}
    fig_n: dict[str, int] = {}
    tabs: dict[str, int] = {}
    for line in body_of(path):
        for m in FIG.finditer(line):
            fig = ("S" if m.group(1).startswith("S") else "") + m.group(1).lstrip("S")
            figs.setdefault(fig, set())
            fig_n[fig] = fig_n.get(fig, 0) + 1
            run, _ = _run_after(line, m.end(1), fig)
            for f, letter in run:
                figs.setdefault(f, set())
                if letter:
                    figs[f].add(letter)
        for m in TAB.finditer(line):
            t = ("S" if m.group(1).startswith("S") else "") + m.group(1).lstrip("S")
            tabs[t] = tabs.get(t, 0) + 1
    return figs, fig_n, tabs


def main() -> None:
    results = {name: scan(SEC / f) for name, f in SECTIONS}

    fig_ids = [str(i) for i in range(1, 8)] + [f"S{i}" for i in range(1, 5)]
    tab_ids = ["1"] + [f"S{i}" for i in range(1, 6)]

    print("Figure citations by section. Each cell is occurrences (distinct panels named).\n")
    header = f"{'figure':9}" + "".join(f"{n:>20}" for n, _ in SECTIONS) + "   panels named ever"
    print(header)
    panels_all: dict[str, set[str]] = {}
    for f in fig_ids:
        cells = []
        for name, _ in SECTIONS:
            figs, fig_n, _ = results[name]
            if f in fig_n:
                cells.append(f"{fig_n[f]} ({len(figs[f])})")
                panels_all.setdefault(f, set()).update(figs[f])
            else:
                cells.append("-")
                panels_all.setdefault(f, set())
        print(f"{('Fig. ' + f):9}" + "".join(f"{c:>20}" for c in cells) +
              f"   {''.join(sorted(panels_all[f])) or '-'}")

    print("\nTable citations by section (occurrences)\n")
    print(f"{'table':9}" + "".join(f"{n:>14}" for n, _ in SECTIONS))
    for t in tab_ids:
        cells = []
        for name, _ in SECTIONS:
            _, _, tabs = results[name]
            cells.append(str(tabs[t]) if t in tabs else "-")
        print(f"{('Table ' + t):9}" + "".join(f"{c:>14}" for c in cells))

    captioned = set()
    if CAPTIONS.is_file():
        for m in re.finditer(r"\*\*Figure (S?\d+)\.", CAPTIONS.read_text(encoding="utf-8")):
            captioned.add(m.group(1))
    print(f"\ncaptions present for: "
          f"{', '.join(sorted(captioned, key=lambda x: (x[0] == 'S', int(x.lstrip('S')))))}")

    cited_figs = {f for name, _ in SECTIONS for f in results[name][1]}
    cited_tabs = {t for name, _ in SECTIONS for t in results[name][2]}
    print("\nfigures never cited in a sentence:", ", ".join(f for f in fig_ids if f not in cited_figs) or "none")
    print("tables never cited in a sentence: ", ", ".join(t for t in tab_ids if t not in cited_tabs) or "none")

    # panel coverage: which panel letters are never mentioned in any section
    LIMITS = {"1": "F", "2": "I", "3": "I", "4": "M", "5": "N", "6": "L", "7": "I",
              "S1": "C", "S2": "H", "S3": "K", "S4": "G"}
    missing_panels = []
    print("\npanel coverage (letters never mentioned in a sentence):")
    for f in fig_ids:
        allp = {chr(65 + i) for i in range(ord(LIMITS[f]) - 64)}
        named = panels_all.get(f, set())
        missing = sorted(allp - named)
        print(f"  Fig. {f:3} {len(named):2}/{len(allp):2} panels named" +
              (f"; never named: {''.join(missing)}" if missing else ""))
        if missing and not f.startswith("S"):
            missing_panels.append((f, "".join(missing)))

    print("\nmain-figure panels with no callout: " +
          (", ".join(f"Fig. {f}{p}" for f, p in missing_panels) or "none"))


if __name__ == "__main__":
    main()
