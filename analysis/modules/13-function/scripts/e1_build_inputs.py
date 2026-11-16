#!/usr/bin/env python3
"""Build PICRUSt2 inputs (fasta + BIOM) from a DADA2 ASV table, with a documented filter.

Why a filter at all
-------------------
PICRUSt2's cost is dominated by phylogenetic placement (EPA-NG), which scales with the number
of *unique* sequences. In our tumour data 27,465 of 35,375 ASVs (77.6%) are singletons, i.e.
observed in exactly one sample, and together they carry 2.78% of reads. Singletons in amplicon
data are overwhelmingly residual sequencing error; placing them costs hours and changes
community function by a fraction of a percent. We therefore drop them and state the filter
explicitly rather than silently.

Retained: prevalence >= 2 samples -> 7,910 ASVs, 97.22% of reads (asserted below).

Orientation
-----------
The DADA2 CSV is ASV x sample (rows = sequences, columns = runs). BIOM v2.1 requires
rows = observations, columns = samples, so no transpose is needed, but IDs must be assigned
positionally (DADA2's fasta and table share row order).

Run:
  python3 e1_build_inputs.py --study PRJNA866676 --min-prev 2
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

WORK = Path("/home/v3724/oscc-mb/work")
DADA = {
    "PRJNA866676": WORK / "n866676/results/dada2/PRJNA866676",
    "PRJNA412445": WORK / "n412445/results/dada2/PRJNA412445",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True, choices=sorted(DADA))
    ap.add_argument("--min-prev", type=int, default=2)
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()

    d = DADA[a.study]
    tab = pd.read_csv(d / "asv_table.csv", index_col=0)
    print(f"{a.study}: ASV table {tab.shape[0]} x {tab.shape[1]}")

    if not all(set(str(i).upper()) <= set("ACGTN") for i in tab.index[:50]):
        raise SystemExit("row index does not look like DNA sequences; aborting")

    tot = tab.sum(axis=1)
    prev = (tab > 0).sum(axis=1)
    keep = prev >= a.min_prev
    reads_kept = 100 * tot[keep].sum() / tot.sum()
    print(f"filter prev >= {a.min_prev}: keep {int(keep.sum())} ASVs "
          f"({100*keep.sum()/len(tab):.1f}%) carrying {reads_kept:.2f}% of reads")
    if reads_kept < 90:
        raise SystemExit(f"filter keeps only {reads_kept:.1f}% of reads; refusing")

    sub = tab.loc[keep]

    out = Path(a.outdir) if a.outdir else (WORK / "picrust" / a.study)
    out.mkdir(parents=True, exist_ok=True)

    fa = out / f"{a.study}.prev{a.min_prev}.fasta"
    with open(fa, "w", encoding="utf-8") as fh:
        for i, seq in enumerate(sub.index, 1):
            fh.write(f">ASV{i}\n{seq}\n")
    print(f"wrote {fa.name} ({fa.stat().st_size/1e6:.1f} MB, ASV1..ASV{len(sub)})")

    samples = list(sub.columns)
    obs_ids = [f"ASV{i}" for i in range(1, len(sub) + 1)]
    data, vals = [], sub.values
    for r in range(vals.shape[0]):
        row = vals[r]
        for c in np.nonzero(row)[0]:
            v = int(row[c])
            if v:
                data.append([r, int(c), v])
    biom = {
        "id": None,
        "format": "Biological Observation Matrix 2.1.0",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "wisp-science e1_build_inputs.py",
        "date": pd.Timestamp.now().isoformat(),
        # biom-format indexes col['metadata'] unconditionally in Table.from_json,
        # so both axes must carry an explicit (null) metadata key.
        "rows": [{"id": o, "metadata": None} for o in obs_ids],
        "columns": [{"id": s, "metadata": None} for s in samples],
        "matrix_type": "sparse",
        "matrix_element_type": "int",
        "shape": [len(obs_ids), len(samples)],
        "data": data,
    }
    bm = out / f"{a.study}.prev{a.min_prev}.biom"
    bm.write_text(json.dumps(biom), encoding="utf-8")
    print(f"wrote {bm.name} ({bm.stat().st_size/1e6:.1f} MB, {len(data)} non-zero entries)")

    meta = {"study": a.study, "min_prev": a.min_prev,
            "n_asv_total": int(len(tab)), "n_asv_kept": int(keep.sum()),
            "pct_asv_kept": round(100 * keep.sum() / len(tab), 2),
            "pct_reads_kept": round(float(reads_kept), 2),
            "n_samples": int(len(samples)),
            "total_reads_kept": int(tot[keep].sum()),
            "fasta": fa.name, "biom": bm.name}
    (out / f"{a.study}.input_meta.json").write_text(json.dumps(meta, indent=1),
                                                   encoding="utf-8")
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
