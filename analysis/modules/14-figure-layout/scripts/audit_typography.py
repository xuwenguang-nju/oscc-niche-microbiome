#!/usr/bin/env python3
"""Per-script typography audit using the AST (regex misses multi-line calls).

For each supplementary-figure script, list every set_title / set_xlabel / set_ylabel
call with the fontsize it resolves to, flagging any role that carries more than one
size, plus any text drawn below the axes (which stretches the canvas and leaves dead
space next to sibling panels that have no such overhang).
"""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]

SCRIPTS = {
    "S1": "13-function/scripts/e5_plot.py",
    "S2": "07-six-niche/scripts/s7b_plot_genera_by_niche.py",
    "S3": "04-niche-analysis/scripts/plot_niche_attribution.py",
    "S4": "06-host-coupling/scripts/s6v2_plot.py",
    "S5": "06-host-coupling/scripts/s6v3_plot.py",
    "S6": "06-host-coupling/scripts/s6v5_plot.py",
    "S7": "06-host-coupling/scripts/plot_s6.py",
}

ROLES = {"set_title": "title", "set_xlabel": "xlabel", "set_ylabel": "ylabel"}


def base_size(tree: ast.AST) -> float | None:
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "apply_figure_style"):
            for kw in node.keywords:
                if kw.arg == "sizes" and isinstance(kw.value, ast.Tuple):
                    elts = kw.value.elts
                    if elts and isinstance(elts[0], ast.Constant):
                        return float(elts[0].value)
    return None


def kw(call: ast.Call, name: str):
    for k in call.keywords:
        if k.arg == name:
            return k.value
    return None


def const(v):
    return v.value if isinstance(v, ast.Constant) else None


def audit(tag: str, rel: str) -> None:
    path = ROOT / "analysis" / "modules" / rel
    if not path.is_file():
        print(f"{tag}  MISSING {rel}")
        return
    tree = ast.parse(path.read_text(encoding="utf-8"))
    base = base_size(tree)

    per_role: dict[str, list[tuple[int, float | None]]] = {"title": [], "xlabel": [], "ylabel": []}
    below: list[str] = []

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        fn = node.func.attr
        if fn in ROLES:
            fs = kw(node, "fontsize")
            size = float(const(fs)) if const(fs) is not None else base
            per_role[ROLES[fn]].append((node.lineno, size))
        elif fn in ("legend", "text", "annotate"):
            bbox = kw(node, "bbox_to_anchor")
            transform = kw(node, "transform")
            y = kw(node, "y")
            # a legend anchored at negative axes-fraction y, or text placed at
            # negative axes coords, is drawn outside the axes
            if isinstance(bbox, ast.Tuple) and bbox.elts:
                last = const(bbox.elts[-1])
                if isinstance(last, float) and last < 0:
                    below.append(f"L{node.lineno} {fn} bbox y={last}")
            if (isinstance(transform, ast.Attribute) and transform.attr == "transAxes"
                    and isinstance(y, ast.UnaryOp) and isinstance(const(y.operand), float)
                    and const(y.operand) < 0):
                below.append(f"L{node.lineno} {fn} y={const(y.operand)}")

    print(f"\n{tag}  {rel}   (base={base})")
    bad = False
    for role, rows in per_role.items():
        if not rows:
            continue
        sizes = {s for _, s in rows}
        tag2 = "OK " if len(sizes) == 1 else "!! "
        if len(sizes) > 1:
            bad = True
        detail = ", ".join(f"L{l}={s}" for l, s in rows)
        print(f"   {tag2}{role:8} sizes={sorted(sizes)}  [{detail}]")
    mark = "!! " if below else "OK "
    print(f"   {mark}below-axes: {below or 'none'}")


def main() -> None:
    for tag, rel in SCRIPTS.items():
        audit(tag, rel)


if __name__ == "__main__":
    main()
