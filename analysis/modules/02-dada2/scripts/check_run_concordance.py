"""Concordance between two independent DADA2 runs of the same cohort.

Why this exists: PRJNA666746 was processed twice with the same pipeline and the
same target subsampling depth but on different hardware -- once locally (Windows,
DADA2 multithreading unavailable, 12.8 h at 1-2 cores) and once on the Linux node
(32 threads, 3.7 h). Both draw an independent random subsample, so they are not
expected to be bit-identical; measuring how far apart they land is a direct
reproducibility check on the step that produced the fourth cohort of S5.

Inputs
  results/dada2/PRJNA666746/asv_table.csv        + taxonomy_silva.csv   (server run)
  results/dada2/PRJNA666746_localrun/asv_table.csv + taxonomy.csv       (local run)

Output
  results/tables/prjna666746_run_concordance.csv
  results/tables/prjna666746_run_concordance_summary.json

Run:  python analysis/modules/02-dada2/scripts/check_run_concordance.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
RES = ROOT / "results" / "dada2"
PROKARYOTA = {"Bacteria", "Archaea"}
# local run kept the "_filt" suffix in sample names; the server run did not
_SUFFIX = re.compile(r"\.rawds_[FR]_filt$")


def genus_matrix(seqtab_csv: Path, tax_csv: Path) -> pd.DataFrame:
    st = pd.read_csv(seqtab_csv, index_col=0)
    tx = pd.read_csv(tax_csv, index_col=0)
    n = tx.shape[0]
    if st.shape[1] != n and st.shape[0] == n:
        st = st.T                                  # ASVs x samples -> samples x ASVs
    if st.shape[1] != n:
        raise ValueError(f"{seqtab_csv}: cannot align {st.shape} to {n} taxonomy rows")
    st.columns = tx.index
    st.index = [_SUFFIX.sub("", str(i)) for i in st.index]
    keep = tx["Kingdom"].isin(PROKARYOTA).values
    genus = tx.loc[keep, "Genus"].fillna("unassigned").replace("", "unassigned").values
    return st.loc[:, keep].T.groupby(genus).sum().T


def main() -> None:
    a = genus_matrix(RES / "PRJNA666746" / "asv_table.csv",
                     RES / "PRJNA666746" / "taxonomy_silva.csv")
    b = genus_matrix(RES / "PRJNA666746_localrun" / "asv_table.csv",
                     RES / "PRJNA666746_localrun" / "taxonomy.csv")

    samples = a.index.intersection(b.index)
    genera = a.columns.intersection(b.columns)
    pa = a.loc[samples, genera]; pa = pa.div(pa.sum(axis=1), axis=0)
    pb = b.loc[samples, genera]; pb = pb.div(pb.sum(axis=1), axis=0)

    bc = 0.5 * np.abs(pa.values - pb.values).sum(axis=1)
    rho = np.array([pa.iloc[i].corr(pb.iloc[i], method="spearman")
                    for i in range(len(samples))])

    per = pd.DataFrame({"sample": samples,
                        "bray_curtis": bc.round(5),
                        "spearman_rho": rho.round(5)})
    out_tab = ROOT / "results" / "tables" / "prjna666746_run_concordance.csv"
    out_tab.parent.mkdir(parents=True, exist_ok=True)
    per.to_csv(out_tab, index=False)

    summary = {
        "cohort": "PRJNA666746",
        "samples_compared": int(len(samples)),
        "genera_compared": int(len(genera)),
        "server_run": {"asvs": 36477, "reads": 2582449, "platform": "Linux, 32 threads",
                       "wall_hours": 3.7},
        "local_run": {"asvs": 36428, "reads": 2572088, "platform": "Windows, 1-2 cores",
                      "wall_hours": 12.8},
        "bray_curtis": {"median": round(float(np.median(bc)), 4),
                        "p95": round(float(np.percentile(bc, 95)), 4),
                        "max": round(float(bc.max()), 4)},
        "spearman_rho": {"median": round(float(np.median(rho)), 4),
                         "min": round(float(rho.min()), 4)},
    }
    out_json = ROOT / "results" / "tables" / "prjna666746_run_concordance_summary.json"
    out_json.write_text(json.dumps(summary, indent=1), encoding="utf-8")

    print(json.dumps(summary, indent=1))
    print(f"\nwrote {out_tab.relative_to(ROOT)}")
    print(f"wrote {out_json.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
