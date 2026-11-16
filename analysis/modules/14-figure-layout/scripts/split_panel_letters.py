#!/usr/bin/env python3
"""Move panel letters OUT of the title and into an independent top-left letter.

Target convention (user's instruction):
  - the letter stands alone at the panel's top-left corner, bold
  - the panel title text is centred, WITHOUT the letter in it
  - S2 (the 2x5 small-multiple figure) keeps its current merged A..J titles

Why ast and not regex: the previous regex pass over the same files silently corrupted
`figsize=`, `fig_*()` and `figures/` (see output/panel-convention-change.md). Parsing the
module and splicing only the recorded spans cannot touch unrelated code.

Each `ax.set_title("A  rest")` becomes
    ax.set_title("rest")
    panel_letter(ax, "A", case="upper", dy=1.02)

Run: python3 split_panel_letters.py [--dry-run]
"""
from __future__ import annotations

import argparse
import ast
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figlib import project_root  # noqa: E402

SCRIPTS = [
    "03-simulation/scripts/plot_simulation_calibration.py",
    "04-niche-analysis/scripts/plot_niche_attribution.py",
    "05-meta-analysis/scripts/plot_forest.py",
    "06-host-coupling/scripts/plot_s6.py",
    "06-host-coupling/scripts/s6v2_plot.py",
    "06-host-coupling/scripts/s6v3_plot.py",
    "06-host-coupling/scripts/s6v4_plot.py",
    "06-host-coupling/scripts/s6v5_plot.py",
    "07-six-niche/scripts/s7_plot.py",
    "08-replication/scripts/s8_plot.py",
    "09-figures/scripts/fig1_study_design.py",
    "11-clinical-survival/scripts/s11_plot.py",
    "12-pancancer/scripts/s12_plot.py",
    "13-function/scripts/e5_plot.py",
]

PREFIX = re.compile(r"^([A-J])[ ]{1,3}")
LETTER_CALL = 'panel_letter({ax}, "{L}", case="upper", dx=-0.10, dy=1.02)'


def _ensure_import(src: str) -> str:
    """Make panel_letter importable without disturbing the existing import style.

    NOTE: check for the *import*, not for any mention. After the splice the file always
    contains `panel_letter(ax, ...)`, so testing `"panel_letter" in src` here silently
    skipped the import and produced NameError at render time.
    """
    if re.search(r"^from runtime import .*\bpanel_letter\b", src, re.M):
        return src
    m = re.search(r"^from runtime import (.+)$", src, re.M)
    if m:
        names = m.group(1).split("#")[0].strip()
        if "panel_letter" not in names:
            new = f"from runtime import {names}, panel_letter"
            src = src[:m.start()] + new + src[m.end():]
        return src
    # exec()-style load of runtime.py puts every function in module globals: nothing to do
    return src


def _split_title(raw: str):
    """Strip the 'A  ' prefix from the ORIGINAL source of a title argument.

    Working on the raw slice (rather than unparsing the AST) keeps f-strings, quoting and
    line continuations byte-identical, and only the three prefix characters move.

    Returns (letter, remaining_source) or (None, None).
    """
    m = re.match(r'^(f?)(["\'])([A-J])[ ]{1,3}', raw)
    if not m:
        return None, None
    fprefix, quote, letter = m.group(1), m.group(2), m.group(3)
    return letter, f"{fprefix}{quote}" + raw[m.end():]


def transform(path: Path, dry: bool) -> tuple[int, list[str]]:
    src = path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    lines = src.splitlines(keepends=True)

    # absolute offsets for (lineno, col_offset) -> index in src
    starts = [0]
    for ln in lines:
        starts.append(starts[-1] + len(ln))

    def off(lineno: int, col: int) -> int:
        return starts[lineno - 1] + col

    # S2 is figure 2 in s7_plot.py: everything after this line keeps its merged titles
    skip_from = None
    if path.name == "s7_plot.py":
        for i, ln in enumerate(lines, 1):
            if ln.lstrip().startswith("fig2 = plt.figure("):
                skip_from = i
                break

    edits = []          # (start, end, replacement)
    letters = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "set_title" and node.args):
            continue
        if skip_from is not None and node.lineno >= skip_from:
            continue
        ax = node.func.value
        ax_name = ax.id if isinstance(ax, ast.Name) else None
        if ax_name is None:
            continue
        a0 = node.args[0]
        raw_first = src[off(a0.lineno, a0.col_offset):off(a0.end_lineno, a0.end_col_offset)]
        letter, new_text = _split_title(raw_first)
        if letter is None:
            continue

        stmt_start = off(node.lineno, node.col_offset)
        stmt_end = off(node.end_lineno, node.end_col_offset)
        title_only = src[stmt_start:stmt_end]

        # rebuild the set_title call with the new first argument
        new_call = title_only.replace(raw_first, new_text, 1)

        indent = " " * node.col_offset
        letter_call = LETTER_CALL.format(ax=ax_name, L=letter)
        edits.append((stmt_start, stmt_end,
                      new_call + "\n" + indent + letter_call))
        letters.append(f"{letter}->{ax_name}")

    if not edits:
        return 0, letters

    out = src
    for start, end, repl in sorted(edits, reverse=True):
        out = out[:start] + repl + out[end:]

    out = _ensure_import(out)
    ast.parse(out)          # fail loudly rather than write a broken file
    if not dry:
        path.write_text(out, encoding="utf-8")
    return len(edits), letters


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    root = project_root(Path(__file__).resolve().parent)

    for rel in SCRIPTS:
        p = root / "analysis/modules" / rel
        if not p.is_file():
            print(f"  [skip] {rel}")
            continue
        n, letters = transform(p, a.dry_run)
        tag = "dry" if a.dry_run else "edit"
        print(f"  [{tag}] {p.name:34} panels={n:2}  {', '.join(letters)}")


if __name__ == "__main__":
    main()
