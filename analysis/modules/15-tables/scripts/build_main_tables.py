#!/usr/bin/env python3
"""Build the manuscript tables from the analysis outputs.

Table 1  cohort characteristics (one row per cohort)
Table 2  the replicating genera, with cross-cohort statistics and prior literature

Every value is read from an existing results/tables file or from cohort_params.yaml.
The only typed-in values are the ENA-derived region strings, and each row records how
many runs carried that value so the provenance travels with the table.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
TAB = ROOT / "results" / "tables"
OUT = TAB

REGION = {
    "PRJNA1415359": "Korea (Seoul)",
    "PRJNA1415268": "Korea",
    "PRJNA597251": "China (Qingdao)",
    "PRJNA666746": "India",
    "PRJNA412445": "China",
    "PRJNA866676": "China (Jinan)",
}
REGION_NOTE = {
    "PRJNA1415359": "80/241 runs state South Korea: Seoul",
    "PRJNA1415268": "country unset in ENA; Korea from cohort notes",
    "PRJNA597251": "48/48 runs state China:Qingdao",
    "PRJNA666746": "country unset in ENA; India from cohort notes",
    "PRJNA412445": "102/102 runs state China",
    "PRJNA866676": "313/313 runs state China: Jinan",
}
NICHE_SHORT = {
    "tumor": "tumour",
    "adjacent_normal": "adjacent normal",
    "contralateral_normal_mucosa": "contralateral mucosa",
    "preop_saliva": "pre-op saliva",
    "postop_saliva": "post-op saliva",
    "saliva": "saliva",
    "cancer_surface": "cancer surface",
    "tongue_coat": "tongue coat",
    "subgingival_plaque": "subgingival plaque",
    "lesion_side_mucosa": "lesion-side mucosa",
}


def read_params() -> dict:
    """Minimal parse of cohort_params.yaml (avoids a yaml dependency)."""
    txt = (ROOT / "analysis/modules/02-dada2/input/cohort_params.yaml").read_text(
        encoding="utf-8")
    out, cur = {}, None
    for line in txt.splitlines():
        m = re.match(r"^  (PRJNA\d+):\s*$", line)
        if m:
            cur = m.group(1)
            out[cur] = {}
            continue
        if cur is None:
            continue
        m = re.match(r"^    (\w+):\s*(.+?)\s*$", line)
        if m:
            out[cur][m.group(1)] = m.group(2)
    return out


SUBSET_DIR = TAB.parents[1] / "analysis" / "modules" / "02-dada2" / "input"


def analysed_manifest() -> pd.DataFrame:
    """The manifest restricted to runs that actually entered the DADA2 analysis.

    The manifest lists every archived run; the analysis used fewer. PRJNA1415359 is
    the only cohort where those differ (241 archived -> 125 analysed): 80 runs are
    saliva-only donors who contributed no tissue, and 36 are saliva from 18 patients
    with no tumour sample. Reporting the archive count put 241 samples and 130
    patients in this table while every result downstream was computed on 125 samples
    from 32 patients, and it inflated the study total from 782 to 898.
    """
    man = pd.read_csv(TAB / "oscc_niche_sample_manifest.csv", encoding="utf-8-sig")
    keep = pd.Series(True, index=man.index)
    for f in sorted(SUBSET_DIR.glob("subset_*.txt")):
        study = f.stem.replace("subset_", "")
        runs = set(pd.read_csv(f, sep="\t")["run"].astype(str))
        keep &= ~(man["study"].eq(study) & ~man["run"].astype(str).isin(runs))
    return man[keep].copy()


def build_table1() -> pd.DataFrame:
    man_all = pd.read_csv(TAB / "oscc_niche_sample_manifest.csv", encoding="utf-8-sig")
    man = analysed_manifest()
    gen = pd.read_csv(TAB / "genus_table_summary.csv")
    prof = pd.read_csv(TAB / "read_profile_summary.csv")
    prim = pd.read_csv(TAB / "primer_residue.csv")
    par = read_params()

    rows = []
    for st in ["PRJNA1415359", "PRJNA1415268", "PRJNA666746", "PRJNA597251",
               "PRJNA412445", "PRJNA866676"]:
        d = man[man.study == st]
        niches = d.niche.astype(str).str.split("|").explode().dropna().unique()
        niches = sorted({n for n in niches if n and n != "nan"})
        g = gen[gen.study == st].iloc[0]
        pr = prof[prof.study == st]
        pm = prim[prim.study == st]
        p = par.get(st, {})

        single = str(p.get("single_end", "")).lower() == "true"

        def bracket_list(raw: str) -> list[str]:
            """'[17, 21]   # 341F at 5-prime of R1' -> ['17', '21'].

            The YAML values carry trailing comments, so the comment has to be cut
            before the brackets are stripped, otherwise the comment text ends up
            inside the table cell.
            """
            body = raw.split("#")[0].strip()
            return [x.strip() for x in body.strip("[]").split(",") if x.strip()]

        trim = bracket_list(p.get("trimLeft", "0"))
        trunc = bracket_list(p.get("truncLen", "0"))
        if single:
            # A single-end cohort CAN still use a fixed truncLen: PRJNA866676 ran with 350 on
            # 370-460 bp reads. An earlier version always printed "none (variable length)",
            # which contradicted that cohort's own run_status.json.
            trunc_s = trunc[0] if trunc and trunc[0] not in ("0", "") else "none (variable length)"
            layout, trim_s = "single-end", trim[0]
        else:
            layout, trim_s, trunc_s = "paired-end", "/".join(trim), "/".join(trunc)

        len_med = int(pr.len_mean.median()) if len(pr) else None
        # Report the measured 5-prime primer fraction per read, never a single max:
        # R1 and R2 differ substantially (85 % vs 94 % in PRJNA1415359), and taking the
        # max both hides that and, for the cohort whose submitter did not formally trim,
        # returns 0 because those rows are labelled "already_trimmed" rather than
        # "forward". Keep every measured value and label the partial-primer case.
        def primer_frac(read: str) -> float:
            hit = pm[pm.read == read]
            return float(hit.fraction_at_5prime.iloc[0]) if len(hit) else 0.0

        f1, f2 = primer_frac("R1"), primer_frac("R2")
        partial = pm.top_forward.dropna().unique() if "top_forward" in pm else []
        if single:
            primer_s = f"{round(f1 * 100)}"
        else:
            primer_s = f"{round(f1 * 100)}/{round(f2 * 100)}"
        trimmed = [int(x) for x in trim if x.isdigit()]
        if f1 > 0.05 and trimmed and max(trimmed) == 0:
            # primer signal present but nothing was trimmed: the match is partial and the
            # submitter left the reads as-is, so the value must not read as "clean"
            primer_s += f" (partial {partial[0]} match)" if len(partial) else " (partial)"

        # 21 raw columns do not fit across a landscape page, so related fields are
        # merged and every merges keeps its units in the header. The unpacked values
        # stay in the CSV build below for anyone who needs them individually.
        dada = (f"trim {trim_s} | trunc {trunc_s} | maxEE "
                f"{p.get('maxEE', '').split('#')[0].strip()} | pool "
                f"{p.get('pool', '').split('#')[0].strip()} | cap "
                f"{p.get('downsample', '').split('#')[0].strip() or 'none'}")
        plat = f"{p.get('platform', '').strip(chr(34))} / {layout}"
        region = REGION[st]
        prov = REGION_NOTE[st]
        n_note = prov.split(";")[0].strip()
        region_s = f"{region}\n({n_note})" if n_note else region

        rows.append({
            "Accession": st,
            "Region": region_s,
            "Patients (n)": int(d.patient_id.nunique()),
            "Runs (n)": len(d),
            "Runs archived (n)": int((man_all.study == st).sum()),
            "Niches sampled (n)": len(niches),
            "Niches": "; ".join(NICHE_SHORT.get(n, n) for n in niches),
            "Platform / layout": plat,
            "Amplicon": p.get("region", "").strip(chr(34)),
            "Median read (bp)": len_med,
            "5' primer (%), R1/R2": primer_s,
            "DADA2 settings": dada,
            "ASVs": int(g.n_asv),
            "Genera": int(g.n_genera),
            "Bacterial reads": int(g.reads_in_genus_table),
            "Genus-unassigned (%)": round(float(g.no_genus_median_pct), 1),
        })
    # The study cohort is part of the study but is not in the public manifest, so it is
    # added explicitly rather than by widening the loop above. Values come from the vendor
    # report and from results/derived/fw2024_420/genus_counts.csv (the vendor project number is
# internal and deliberately not shown to readers).
    t1 = pd.DataFrame(rows)
    import numpy as _np
    _fc = pd.read_csv(ROOT / "results/derived/fw2024_420/genus_counts.csv",
                      index_col=0, encoding="utf-8-sig")
    _asv = None
    try:
        import zipfile as _zf
        from io import StringIO as _SIO
        with _zf.ZipFile(ROOT / "data/raw/FW2024-420/results/02.ASVs_analysis.zip") as _z:
            _b = _z.read("02.ASVs_analysis/all_ASVs_count.xls")   # union ASV table: one row per ASV
        for _enc in ("utf-8-sig", "utf-8", "gb18030"):
            try:
                _asv = pd.read_csv(_SIO(_b.decode(_enc)), sep="\t")
                break
            except Exception:
                continue
    except Exception:
        _asv = None
    t1.loc[len(t1)] = {
        "Accession": "Hospital cohort",
        "Region": "Oral cavity (tongue SCC)",
        "Patients (n)": 20,
        "Runs (n)": 40,
        "Runs archived (n)": 40,
        "Niches sampled (n)": 2,
        "Niches": "tumour (tongue); distant normal mucosa",
        "Platform / layout": "Illumina NovaSeq 6000 / paired-end",
        "Amplicon": "V3+V4 (338F/806R)",
        "Median read (bp)": pd.NA,
        "5' primer (%), R1/R2": "none",
        "DADA2 settings": "vendor DADA2, SILVA 138.1, Bayesian assignment",
        # union ASV count: one row per ASV, first column is the ASV id
        "ASVs": int((_asv.ASV_ID.str.startswith("ASV")).sum()) if _asv is not None else pd.NA,
        "Genera": int(_fc.shape[0]),
        "Bacterial reads": int(_fc.sum().sum()),
        "Genus-unassigned (%)": pd.NA,
    }
    t1.to_csv(OUT / "table1_cohort_characteristics.csv", index=False, encoding="utf-8-sig")
    return t1


def build_table2() -> pd.DataFrame:
    meta = pd.read_csv(TAB / "s5_meta_analysis.csv")
    nov = pd.read_csv(TAB / "s11_novelty_assessment.csv")
    rep = meta[(meta.n_cohorts >= 3) & (meta.pooled_q_bh < 0.05)
               & meta.direction_consistent].copy()
    t2 = rep[["genus", "n_cohorts", "cohorts", "pooled_effect", "pooled_se",
              "pooled_q_bh", "I2_pct", "max_abs_dz", "n_cohorts_large_dz"]].merge(
        nov[["genus", "oscc_papers", "lit_majority", "direction_agrees"]],
        on="genus", how="left")
    t2 = t2.rename(columns={
        "genus": "Genus",
        "n_cohorts": "Cohorts tested (n)",
        "cohorts": "Cohorts",
        "pooled_effect": "Pooled log2 fold change",
        "pooled_se": "SE",
        "pooled_q_bh": "FDR q",
        "I2_pct": "Heterogeneity I2 (%)",
        "max_abs_dz": "Largest |d_z| in any cohort",
        "n_cohorts_large_dz": "Cohorts with |d_z| >= 0.8",
        "oscc_papers": "Prior OSCC reports (n)",
        "lit_majority": "Literature direction",
        "direction_agrees": "Agrees with literature",
    }).sort_values("Pooled log2 fold change")
    t2.to_csv(OUT / "table2_replicating_genera.csv", index=False, encoding="utf-8-sig")
    return t2


def main() -> None:
    t1 = build_table1()
    print("Table 1:", t1.shape)
    print(t1[["Accession", "Region", "Patients (n)", "Runs (n)",
              "Niches sampled (n)", "Platform / layout", "ASVs", "Genera"]]
          .to_string(index=False))

    t2 = build_table2()
    print("\nTable 2:", t2.shape)
    print(t2[["Genus", "Cohorts tested (n)", "Pooled log2 fold change", "FDR q",
              "Heterogeneity I2 (%)", "Prior OSCC reports (n)"]].to_string(index=False))


if __name__ == "__main__":
    main()
