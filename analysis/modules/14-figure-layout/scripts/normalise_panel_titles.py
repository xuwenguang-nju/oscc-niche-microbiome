#!/usr/bin/env python3
"""Normalise panel titles across every figure script to one convention.

The target convention (set by the user):
  1. no figure-level title inside the image -- the title belongs in the figure legend
  2. every panel may keep its title, but it must be CENTRED
  3. panel identifiers are uppercase letters: A, B, C, D ...

Before this pass the scripts used three incompatible styles, and three figures carried the
panel letter twice (a standalone bold letter drawn by `panel_letter()` plus a "(a)" prefix in
the title). This pass removes `suptitle`, folds the letter into a single centred title, and
drops the duplicated standalone letters.

Why a script and not hand edits: 13 files, ~40 title strings, three different prefix forms
("a  ", "b  ", "(a) ", f-strings, and a variable holding the title). A paren-balanced scan is
safer than 40 regexes, and it is reproducible if a figure is added later.

Run: python3 normalise_panel_titles.py [--dry-run]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import project_root  # noqa: E402

SCRIPTS = [
    "analysis/modules/03-simulation/scripts/plot_simulation_calibration.py",
    "analysis/modules/04-niche-analysis/scripts/plot_niche_attribution.py",
    "analysis/modules/05-meta-analysis/scripts/plot_forest.py",
    "analysis/modules/06-host-coupling/scripts/plot_s6.py",
    "analysis/modules/06-host-coupling/scripts/s6v2_plot.py",
    "analysis/modules/06-host-coupling/scripts/s6v3_plot.py",
    "analysis/modules/06-host-coupling/scripts/s6v4_plot.py",
    "analysis/modules/06-host-coupling/scripts/s6v5_plot.py",
    "analysis/modules/07-six-niche/scripts/s7_plot.py",
    "analysis/modules/08-replication/scripts/s8_plot.py",
    "analysis/modules/09-figures/scripts/fig1_study_design.py",
    "analysis/modules/11-clinical-survival/scripts/s11_plot.py",
    "analysis/modules/12-pancancer/scripts/s12_plot.py",
    "analysis/modules/13-function/scripts/e5_plot.py",
]

# a panel prefix inside a title string: optional "(", one letter a-j, optional ")", whitespace
PREFIX = re.compile(r'^(\(?)([a-j])(\)?)(\s{1,3})')


def _strip_call(text: str, needle: str) -> tuple[str, int]:
    """Delete every `needle(...)` statement, keeping paren balance. Returns (text, count)."""
    out, i, removed = [], 0, 0
    while True:
        j = text.find(needle, i)
        if j < 0:
            out.append(text[i:])
            break
        # find the opening paren of the call
        k = text.find("(", j + len(needle) - 1)
        if k < 0:
            out.append(text[i:])
            break
        depth, m = 0, k
        while m < len(text):
            ch = text[m]
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    break
            m += 1
        # swallow the rest of the logical line (trailing newline / line continuation)
        end = m + 1
        while end < len(text) and text[end] in " \t":
            end += 1
        if end < len(text) and text[end] == "\n":
            end += 1
        # drop the indentation that preceded it, if the line becomes empty
        start = j
        while start > 0 and text[start - 1] in " \t":
            start -= 1
        out.append(text[i:start])
        i = end
        removed += 1
    return "".join(out), removed


def _fix_titles(text: str) -> tuple[str, int]:
    """Uppercase panel letters and centre every set_title() call."""
    out, i, changed = [], 0, 0
    while True:
        j = text.find(".set_title(", i)
        if j < 0:
            out.append(text[i:])
            break
        k = j + len(".set_title(")
        depth, m = 1, k
        while m < len(text):
            ch = text[m]
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    break
            m += 1
        call = text[j:m + 1]

        # loc="left" -> default centre. keep any other kwargs untouched.
        new_call = call.replace(', loc="left"', "").replace('loc="left", ', "")
        new_call = new_call.replace('loc="left"', "")

        # first argument, if a literal/f-string, may carry a panel prefix
        arg_match = re.match(r'\.set_title\(\s*(f?)"([^"]*)"', new_call)
        if arg_match:
            is_f, body = arg_match.group(1), arg_match.group(2)
            pref = PREFIX.match(body)
            if pref and not body.startswith(("$", "\\")):
                letter = pref.group(2).upper()
                rest = body[pref.end():]
                new_body = f"{letter}  {rest}"
                new_call = new_call.replace(f'{is_f}"{body}"', f'{is_f}"{new_body}"', 1)
                changed += 1
        out.append(text[i:j])
        out.append(new_call)
        i = m + 1
    return "".join(out), changed


def _drop_panel_letter_loops(text: str) -> tuple[str, int]:
    """Remove `for ... in zip(...): panel_letter(...)` blocks (the letter is now in the title)."""
    pat = re.compile(
        r"\n[ \t]*for [A-Za-z_]+, [A-Za-z_]+ in zip\([^\n]*\):\n"
        r"[ \t]*panel_letter\([^\n]*\)\n")
    text, n = pat.subn("\n", text)
    single = re.compile(r"\n[ \t]*panel_letter\([^\n]*\)\n")
    text, n2 = single.subn("\n", text)
    return text, n + n2


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    root = project_root(Path(__file__).resolve().parent)

    for rel in SCRIPTS:
        p = root / rel
        if not p.is_file():
            print(f"  [skip] missing {rel}")
            continue
        src = p.read_text(encoding="utf-8")
        text, n_sup = _strip_call(src, ".suptitle(")
        text, n_loop = _drop_panel_letter_loops(text)
        text, n_ttl = _fix_titles(text)
        if text == src:
            print(f"  [same] {p.name}")
            continue
        if not a.dry_run:
            p.write_text(text, encoding="utf-8")
        print(f"  [edit] {p.name:36} suptitle-{n_sup}  letter-loop-{n_loop}  retitled-{n_ttl}")


if __name__ == "__main__":
    main()
