#!/usr/bin/env python3
"""Build Table S2 from the hospital's original clinical record.

Why this exists
---------------
`results/tables/tableS2_clinical_characteristics.csv` had no generator: it was transcribed by hand
from `data/raw/患者信息.xlsx`. The transcription itself was correct — all 20 rows match the source
field by field — but the column header read "pTNM (AJCC)", and neither that edition nor the
pathological status is stated anywhere in the source. A hand-written header is exactly where such an
overstatement survives review, so the table is now derived from the record and the header says only
what the record says.

What the source states, and what it does not
--------------------------------------------
Columns: 住院号 (admission number), 姓名 (name), 性别, 年龄, 部位 (subsite), TNM分期, 组织学分级.
There is no edition field and no c/p prefix, so "AJCC 8th, pathological" cannot be asserted here.

One internal clue does constrain the edition: patient P16 is recorded as T3N3bM0, and N3b is a
subcategory that exists only from the AJCC/UICC 8th edition — the 7th subdivides N2 but leaves N3
undivided. So the record is consistent with the 8th edition (and cannot be a faithful 7th-edition
transcription), which is reported as an inference in the table footnote and flagged for the authors
to confirm, not adopted as fact.

Output: results/tables/tableS2_clinical_characteristics.csv
Run:    python analysis/scripts/build_table_s2.py
"""
from __future__ import annotations

import pathlib
import re

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / "data/raw/患者信息.xlsx"
OUT = ROOT / "results/tables/tableS2_clinical_characteristics.csv"

SEX = {"男": "male", "女": "female"}
SITE = {"左": "Left", "右": "Right"}


def main() -> None:
    src = pd.read_excel(SRC)
    rows = []
    for i, r in src.iterrows():
        tnm = re.sub(r"\s+", "", str(r["TNM分期"]))
        m = re.match(r"T(\w+?)N(\w+?)M(\d)", tnm)
        if not m:
            raise SystemExit(f"unparseable TNM on row {i + 2}: {r['TNM分期']!r}")
        t_, n_, m_ = m.groups()
        rows.append({
            "Patient": f"P{i + 1:02d}",
            "Sex": SEX[str(r["性别"]).strip()],
            "Age (y)": int(r["年龄"]),
            "Subsite": f"{SITE[str(r['部位']).strip()[0]]} tongue",
            # the record gives T, N and M without a c/p prefix, so they are reported as recorded
            "TNM stage (as recorded)": f"T{t_} N{n_} M{m_}",
            "Grade": str(r["组织学分级"]).strip(),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False, encoding="utf-8-sig")

    print(out.to_string(index=False))
    col = out["TNM stage (as recorded)"]
    n_cat = col.str.extract(r"N(\w+?) M")[0]
    t_cat = col.str.extract(r"^T(\w+?) ")[0]
    m_cat = col.str.extract(r"M(\d)$")[0]
    print(f"\nwrote {OUT.name}: {len(out)} patients")
    print(f"T categories: {sorted(t_cat.unique())}")
    print(f"N categories: {sorted(n_cat.unique())}  (N3b in {int((n_cat == '3b').sum())} patient)")
    print(f"M categories: {sorted(m_cat.unique())}")


if __name__ == "__main__":
    main()
