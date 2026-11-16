#!/usr/bin/env python3
"""Regenerate the submission package inventory from disk.

The inventory is derived data: it lists what the package actually contains, with sizes and
timestamps. It had been hand-maintained and still listed `04-tables/Table2.csv`, a file that no
longer exists after Table 2 was folded into the supplement (Table 2 -> Table S1, and the clinical
table that had been Table 2 became Table S2).

Run: python analysis/scripts/refresh_package_inventory.py
"""
from __future__ import annotations

import csv
import datetime as dt
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
PKG = ROOT / "submission" / "jom-2026"
TARGETS = [PKG / "package-inventory.csv",
           PKG / "07-source-and-reproducibility" / "package-inventory.csv"]
# the inventory must not list itself or the transient Office lock files
SKIP = {"package-inventory.csv"}
SKIP_SUFFIX = {"~$"}


def main() -> None:
    rows = []
    for f in sorted(PKG.rglob("*")):
        if not f.is_file() or f.name in SKIP or f.name.startswith("~$"):
            continue
        rel = f.relative_to(PKG).as_posix()
        rows.append({
            "path": rel,
            "bytes": f.stat().st_size,
            "modified": dt.datetime.fromtimestamp(f.stat().st_mtime).strftime("%Y-%m-%d %H:%M"),
        })

    for out in TARGETS:
        with open(out, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["path", "bytes", "modified"])
            w.writeheader()
            w.writerows(rows)

    tables = [r["path"] for r in rows if r["path"].startswith("04-tables/")]
    print(f"{len(rows)} files inventoried -> {TARGETS[0].relative_to(ROOT)} (+ package mirror)")
    print("04-tables:", ", ".join(p.split('/')[-1] for p in tables))


if __name__ == "__main__":
    main()
