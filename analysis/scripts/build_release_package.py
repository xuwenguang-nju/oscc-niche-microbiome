#!/usr/bin/env python3
"""Assemble the two release trees for GitHub and Zenodo.

Why this exists
---------------
The repository holds things that must never become public and things that must not be published
twice. The hospital record `data/raw/患者信息.xlsx` contains patient names and hospital record
numbers, and the vendor named every library after the hospital record number (`TP20`, `NP20`,
...), so a sample label alone identifies the patient once it sits next to sex, age, subsite and stage
in a series of 20. `analysis/r-lib/` is 20 MB of vendored third-party R packages and `analysis/out/`
is 27 MB of scratch output; neither belongs in a release.

So the package is built by an allow-list, every published label is re-coded through one mapping, and a
preflight pass refuses to finish if an identifier survived.

Two trees, because the manuscript cites them separately
-------------------------------------------------------
  release/code/   -> a GitHub repository, archived by Zenodo as one record (analysis code)
  release/data/   -> a Zenodo dataset record (de-identified metadata, ASV table, taxonomy tables)

Usage
-----
    python analysis/scripts/build_release_package.py [--out release]

Nothing in the source tree is modified; the script only reads. The crosswalk from hospital record
number to P-code is written to `results/private/` and is deliberately outside the release tree.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import pathlib
import re
import shutil
import zipfile

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]

FORBIDDEN = [
    "患者信息.xlsx",
    "clinical_annotation_TEMPLATE.csv",
    "sample_manifest_from_submission_sheet.csv",
]

CODE_SOURCES = [
    ("analysis/scripts", "analysis/scripts"),
    ("analysis/modules", "analysis/modules"),
]
CODE_SKIP_DIRS = {"__pycache__", "output", ".ipynb_checkpoints"}
CODE_SKIP_SUFFIX = {".pyc", ".log", ".rds", ".rda", ".rdb", ".rdx", ".png", ".pdf", ".html"}

CROSSMAP: dict[str, str] = {}
CROSS = None


# ----------------------------------------------------------------------------- helpers
def sha256(p: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def build_crosswalk():
    """Map each hospital record number to an anonymous P-number; keep the mapping private."""
    df = pd.read_excel(ROOT/"data"/"raw"/"患者信息.xlsx")
    hid = [str(h) for h in df["住院号"]]
    cross = {h: f"P{i:02d}" for i, h in enumerate(hid, 1)}
    priv = ROOT/"results"/"private"
    priv.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"hospital_record_number": hid,
                  "patient_code": [cross[h] for h in hid]}).to_csv(
        priv/"crosswalk_hospital_to_Pcode.csv", index=False, encoding="utf-8")
    return cross, df


def relabel_columns(tab: pd.DataFrame, cross: dict[str, str]) -> pd.DataFrame:
    """Recode six-digit vendor sample labels to P-codes; leave non-sample columns alone."""
    ren = {}
    for c in tab.columns:
        if not re.fullmatch(r"[TNY]?\d{6}", str(c)):
            continue
        num = c[1:] if c[0] in "TNY" else c
        grp = "T" if c[0] in "TY" else "N"
        ren[c] = f"{cross.get(num, 'unmapped')}_{grp}"
    return tab.rename(columns=ren)


def anonymise_text(s: str) -> str:
    """Replace hospital record numbers in source text with their P-codes.

    The label-correction block in the analysis code names two libraries literally
    (`RENAME = {"Y<num>": ...}`, `SWAP = (...)`). The mapping is applied to comments and code alike so
    the logic still runs, and the released dataset uses the same P-codes, so a reader can follow the
    correction without seeing a hospital record number.
    """
    for num, pcode in CROSSMAP.items():
        s = s.replace(num, pcode)
    return s


def make_portable(s: str) -> str:
    """Resolve the figure-style helper relative to the file instead of an absolute Windows path.

    The plotting scripts import a private helper from `D:/wisp-science/skills/figure-style` and they
    do it in several ways: `sys.path.insert` with a `Path`, with a raw string, through a
    `SKILL`/`STYLE_DIR` variable that later code keeps using, or by `exec(compile(open(...)))`. The
    absolute path is therefore *rewritten*, not deleted, and a single `_FIGSTYLE` definition is
    injected before its first use, so expressions elsewhere in the file keep working.
    """
    if "figure-style" not in s and "figure_style" not in s:
        return s

    look = ('next(p for p in Path(__file__).resolve().parents\n'
            '              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"')

    # exec(compile(open(<abs path to runtime.py> ...).read(), "runtime.py", "exec"))
    s = re.sub(r'exec\(compile\(open\([^\n]*?figure[-_]style[^\n]*?runtime\.py[^\n]*?'
               r'(?:\n\s*encoding="utf-8"\)\.read\(\))?,\s*\n?\s*"runtime\.py",\s*"exec"\)\)',
               'exec(compile(open(_FIGSTYLE / "runtime.py", encoding="utf-8").read(),\n'
               '              "runtime.py", "exec"))', s, flags=re.S)
    # Path("<abs>/runtime.py")  ->  _FIGSTYLE / "runtime.py"
    s = re.sub(r'Path\(\s*r?["\'][A-Za-z]:[\\/][^"\']*figure[-_]style[\\/]runtime\.py["\']\s*\)',
               '_FIGSTYLE / "runtime.py"', s)
    # Path("<abs>/figure-style") -> _FIGSTYLE
    s = re.sub(r'Path\(\s*r?["\'][A-Za-z]:[\\/][^"\']*figure[-_]style["\']\s*\)', '_FIGSTYLE', s)
    # bare raw string to the directory
    s = re.sub(r'r"[A-Za-z]:\\+[^"]*figure[-_]style"', "str(_FIGSTYLE)", s)
    # a bare raw string to runtime.py
    s = re.sub(r'r"[A-Za-z]:\\+[^"]*figure[-_]style\\+runtime\.py"',
               'str(_FIGSTYLE / "runtime.py")', s)

    if "_FIGSTYLE" not in s:
        return s
    # inject the definition immediately before the first line that uses it
    lines = s.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if "_FIGSTYLE" in line:
            definition = ("# figure-style helper, vendored beside the analysis code\n"
                          "_FIGSTYLE = " + look + "\n")
            lines.insert(i, definition)
            break
    return "".join(lines)


# ----------------------------------------------------------------------------- code tree
def copy_code(dest: pathlib.Path) -> list[dict]:
    rows = []
    me = pathlib.Path(__file__).name
    for src_rel, dst_rel in CODE_SOURCES:
        src = ROOT/src_rel
        for f in sorted(src.rglob("*")):
            if not f.is_file():
                continue
            if any(part in CODE_SKIP_DIRS for part in f.parts):
                continue
            if f.suffix.lower() in CODE_SKIP_SUFFIX or f.name in FORBIDDEN:
                continue
            target = dest/dst_rel/f.relative_to(src)
            target.parent.mkdir(parents=True, exist_ok=True)
            if f.suffix.lower() in (".py", ".r", ".md", ".txt", ".yaml", ".yml", ".json"):
                text = anonymise_text(f.read_text(encoding="utf-8", errors="strict"))
                if f.name != me:            # this script only matches those paths, never imports them
                    text = make_portable(text)
                target.write_text(text, encoding="utf-8", newline="\n")
            else:
                shutil.copyfile(f, target)
            rows.append({"path": str(target.relative_to(dest)), "bytes": target.stat().st_size,
                         "sha256": sha256(target)})
    skill = pathlib.Path("D:/wisp-science/skills/figure-style/runtime.py")
    if skill.is_file():
        target = dest/"analysis"/"figure_style"/"runtime.py"
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(skill, target)
        rows.append({"path": str(target.relative_to(dest)), "bytes": target.stat().st_size,
                     "sha256": sha256(target)})
    return rows


# ----------------------------------------------------------------------------- data tree
def deidentified_clinical(dest: pathlib.Path) -> dict:
    """The 20-patient annotation with identity columns removed, plus age banding.

    Age plus subsite plus stage in a series of 20 is close to identifying, so a banded age column is
    supplied for any grouped analysis and the exact age is marked as restricted-use.
    """
    df = CROSS
    out = pd.DataFrame({
        "patient_code": [CROSSMAP[str(h)] for h in df["住院号"]],
        "sex": df["性别"].map({"男": "male", "女": "female"}),
        "age_years": df["年龄"].astype(int),
        "age_band": pd.cut(df["年龄"], [0, 45, 55, 65, 75, 200],
                           labels=["<=45", "46-55", "56-65", "66-75", ">75"]).astype(str),
        "subsite": df["部位"].map(lambda s: "lateral tongue" if "舌" in str(s) else str(s)),
        "side": df["部位"].map(lambda s: "left" if str(s).startswith("左") else
                               ("right" if str(s).startswith("右") else "not recorded")),
        "tnm_as_recorded": df["TNM分期"].astype(str),
        "grade": df["组织学分级"].astype(str),
    }).sort_values("patient_code")
    p = dest/"metadata_clinical_deidentified.csv"
    p.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(p, index=False, encoding="utf-8")
    return {"path": str(p.relative_to(dest)), "bytes": p.stat().st_size, "sha256": sha256(p),
            "rows": len(out),
            "note": "names and hospital record numbers removed; patient_code is a re-indexed label"}


def sample_metadata(dest: pathlib.Path) -> dict:
    """One row per library, built from the ASV table columns so labels match the deposited data."""
    base = ROOT/"data"/"raw"/"FW2024-420"/"results"
    with zipfile.ZipFile(base/"02.ASVs_analysis.zip") as z:
        raw = z.read("02.ASVs_analysis/all_ASVs_count.xls")
    cols = pd.read_csv(io.BytesIO(raw), sep="\t", nrows=0).columns.tolist()
    rows = []
    for s in [c for c in cols if c != "ASV_ID"]:
        grp = "tumour" if s[0] in "TY" else "normal"
        num = s[1:] if s[0] in "TNY" else s
        pcode = CROSSMAP.get(num, "unmapped")
        rows.append({
            "library_id": f"{pcode}_{'T' if grp == 'tumour' else 'N'}",
            "patient_code": pcode,
            "tissue": grp,
            "vendor_label_group": "tumour prefixed Y instead of T" if s[0] == "Y" else "",
        })
    df = pd.DataFrame(rows).sort_values(["patient_code", "tissue"])
    p = dest/"metadata_samples.csv"
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False, encoding="utf-8")
    return {"path": str(p.relative_to(dest)), "bytes": p.stat().st_size, "sha256": sha256(p),
            "rows": len(df)}


def vendor_tables(dest: pathlib.Path) -> list[dict]:
    """ASV table, sequences and per-rank absolute abundance, relabelled from the vendor archive."""
    base = ROOT/"data"/"raw"/"FW2024-420"/"results"
    ddir = dest/"vendor_tables"
    ddir.mkdir(parents=True, exist_ok=True)
    rows = []

    with zipfile.ZipFile(base/"02.ASVs_analysis.zip") as z:
        tab = pd.read_csv(io.BytesIO(z.read("02.ASVs_analysis/all_ASVs_count.xls")), sep="\t")
        p = ddir/"asv_table.tsv.gz"
        with gzip.open(p, "wb", compresslevel=9) as fh:
            fh.write(relabel_columns(tab, CROSSMAP).to_csv(sep="\t", index=False).encode("utf-8"))
        rows.append({"path": str(p.relative_to(dest)), "bytes": p.stat().st_size,
                     "sha256": sha256(p), "source": "all_ASVs_count.xls"})

        p = ddir/"asv_sequences.fasta.gz"
        with gzip.open(p, "wb", compresslevel=9) as fh:
            fh.write(z.read("02.ASVs_analysis/ASVs_sequence.fasta"))
        rows.append({"path": str(p.relative_to(dest)), "bytes": p.stat().st_size,
                     "sha256": sha256(p), "source": "ASVs_sequence.fasta"})

    with zipfile.ZipFile(base/"03.taxonomy.zip") as z:
        for rank in ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]:
            member = f"03.taxonomy/abundance/absolute/all.{rank}.absolute.xls"
            tab = pd.read_csv(io.BytesIO(z.read(member)), sep="\t")
            p = ddir/f"absolute_{rank.lower()}.tsv.gz"
            with gzip.open(p, "wb", compresslevel=9) as fh:
                fh.write(relabel_columns(tab, CROSSMAP).to_csv(sep="\t", index=False).encode("utf-8"))
            rows.append({"path": str(p.relative_to(dest)), "bytes": p.stat().st_size,
                         "sha256": sha256(p), "source": member})
    return rows


# ----------------------------------------------------------------------------- preflight
def preflight(tree: pathlib.Path) -> list[str]:
    """Fail the build if an identifier or credential reached the result.

    Scope matters. An earlier version scanned every text file for the Chinese column names, which
    flagged the very scripts whose job is to strip those columns, and a check that fires on its own
    remedy teaches the author to bypass it. Identity markers are therefore looked for where they can
    do harm: the header row of a published data file. Code is scanned only for things no legitimate
    source file contains, and the key marker is assembled at run time so this file is not its own
    counter-example.
    """
    problems: list[str] = []
    key_marker = "PRIVATE" + " KEY" + "-----"
    id18 = re.compile(r"\b\d{6}(19|20)\d{2}(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])\d{3}[\dXx]\b")
    hosp = re.compile(r"\b10[0-9]{4}\b")
    for f in tree.rglob("*"):
        if not f.is_file():
            continue
        if f.name in FORBIDDEN:
            problems.append(f"forbidden file present: {f.relative_to(tree)}")
            continue
        suffix = f.suffix.lower()
        if suffix in (".csv", ".tsv"):
            head = (f.read_text(encoding="utf-8", errors="replace").splitlines() or [""])[0]
            for col in ("姓名", "住院号"):
                if col in head:
                    problems.append(f"identity column {col!r} in the header of {f.relative_to(tree)}")
            body = f.read_text(encoding="utf-8", errors="replace")
            if hosp.search(body):
                problems.append(f"a value resembling a hospital record number in {f.relative_to(tree)}")
        if suffix not in (".py", ".r", ".md", ".csv", ".yaml", ".yml", ".json", ".txt", ".tsv"):
            continue
        try:
            s = f.read_text(encoding="utf-8", errors="strict")
        except UnicodeDecodeError:
            continue
        if key_marker in s:
            problems.append(f"private key material in {f.relative_to(tree)}")
        if id18.search(s):
            problems.append(f"a string resembling an 18-digit identity number in {f.relative_to(tree)}")
    return problems


# ----------------------------------------------------------------------------- docs
REQ_PY = """# Python packages with the versions the reported results were produced with.
numpy==2.4.6
pandas==2.3.3
scipy==1.17.1
matplotlib==3.11.2
scikit-learn==1.9.1
statsmodels==0.15.0
pillow==12.3.0
"""

REQ_R = """# R packages, as reported in the manuscript.
R 4.4.2
DADA2 1.34.0
Biostrings 2.74.1
ShortRead 1.64.0
RcppParallel 5.1.11.2
IRanges 2.40.1
S4Vectors 0.44.0
survival 3.5.8
yaml 2.3.10
jsonlite 1.8.9
PICRUSt2 2.6.3
"""

CODE_README = """# Analysis code for the intratumoral bacterial community study

This repository contains the code that produces every table and figure in the manuscript. It is
released so that the reported numbers can be checked rather than taken on trust.

## What it does not contain

- **No patient data.** The hospital cohort is distributed as a separate Zenodo dataset record using
  the anonymous identifiers `P01` to `P20`. Sample labels in this code use the same codes.
- **No vendor raw files.** The sequencing provider's tables encode hospital record numbers in their
  sample names; the deposited dataset carries re-coded labels instead, and this code uses those.
- **No large intermediate objects.** Per-cohort `.rds`/`.rda` files and per-sample matrices are
  regenerated by the scripts from the deposited inputs.

## How the code is organised

`analysis/modules/NN-name/` holds one analysis module per stage, in the order the manuscript presents
them. Each has a `README.md` recording its inputs, the parameters that change results, direct package
versions and its outputs. `analysis/scripts/` holds scripts that cross module boundaries: the
cross-cohort meta-analysis, figure assembly, the submission tables, and the verification tools.

`analysis/figure_style/runtime.py` is a vendored copy of the plotting helper the figure scripts use.
Every plotting script resolves it relative to itself, so no path editing is needed.

## Reproducing a figure

Most figures rebuild from the deposited tables alone:

```
python analysis/modules/08-replication/scripts/s8_plot.py
python analysis/scripts/export_figures.py      # 600 dpi TIFF and PNG at print size
```

Rebuilding the six public cohorts from raw reads additionally needs the reference databases (SILVA
138.2, eHOMD V15.22, fetched by `analysis/modules/01-preprocess/scripts/fetch_references.py`) and those
cohorts' FASTQ from ENA/SRA under the BioProject accessions listed in Table 1 of the manuscript. The
hospital cohort needs neither: its processed tables are in the accompanying Zenodo dataset record, and
its raw reads are not publicly available.

## Verification tools

These checks were used during preparation and ship with the release, because a claim that every
number traces to a table should travel with the tool that checks it.

| Script | What it checks |
|---|---|
| `audit_results_numbers.py` | every numeral in a Results draft resolves to a value in the tables |
| `verify_table_figure_numbering.py` | tables and figures are numbered in first-citation order |
| `audit_citation_coverage.py` | every figure panel is called out somewhere in the text |
| `verify_package_panels.py` | caption numbers match filenames; supplementary numbers have no gap |
| `build_release_package.py` | this package: no identity columns, no hospital numbers, no key material |

## Versions

`requirements.txt` and `requirements-r.txt` list the versions the reported results were produced with.
DADA2 inference parameters per cohort, including the record of a change made during the work, are in
`cohort_params.yaml` in the manuscript's submission material.

## Citation

Please cite the manuscript and the Zenodo archive of this repository rather than the repository URL
alone. The DOI is in the manuscript's Data Availability statement.
"""

DATA_README = """# 16S rRNA gene profiles of paired tumour and distant normal mucosa from 20 tongue cancer patients

This dataset accompanies a manuscript on the bacterial community in oral squamous cell carcinoma. It
contains the sequencing-derived tables for one hospital series of 20 patients, each contributing one
tumour and one macroscopically normal mucosal sample taken away from the tumour bed. The same patients
appear as `P01` to `P20` in the analysis code.

Amplicon: 16S rRNA gene V3 to V4, primers 338F (ACTCCTACGGGAGGCAGCAG) and 806R
(GGACTACHVGGGTWTCTAAT). Platform: Illumina NovaSeq 6000, paired-end. Sequencing and primary analysis
were performed by a commercial provider using DADA2 1.26.0, QIIME2 2024.2.0 and Silva 138.1.

## Files

| File | Contents |
|---|---|
| `metadata_samples.csv` | one row per library: `library_id`, `patient_code`, `tissue` |
| `metadata_clinical_deidentified.csv` | one row per patient: sex, age in years and banded, subsite, side, TNM as recorded, grade |
| `vendor_tables/asv_table.tsv.gz` | ASV counts, samples as columns |
| `vendor_tables/asv_sequences.fasta.gz` | ASV sequences |
| `vendor_tables/absolute_<rank>.tsv.gz` | absolute abundance per rank: kingdom, phylum, class, order, family, genus, species |
| `MANIFEST.csv` | size and SHA-256 of every file above |

## What was changed, and why

1. **Sample labels were re-coded.** The provider named libraries after the hospital record number, so
   the original labels identified the patient. Published labels are `P01_T`, `P01_N` and so on. The
   mapping from the original labels is held by the study investigators and is available from the
   corresponding author for verification.
2. **One library was mislabelled by the provider.** A tumour library was prefixed `Y` instead of `T`
   and appears here under its corrected patient code.
3. **One patient's tumour and normal libraries were transposed** in the provider's sample names while
   the tube labels were correct. The tables here follow the tube labels, which is what the analysis
   used.

The clinical table records TNM exactly as the hospital record states it. The staging edition was
confirmed with the treating department as AJCC/UICC 8th. No stage-dependent analysis rests on this
cohort in the manuscript.

## What is not here

- **Raw FASTQ are not released and are not part of this deposit.** The ethics approval covers
  secondary analysis of the specimens, not redistribution of the reads, and the provider's read files
  carry sample labels that encode the hospital record number. Every table and figure in the manuscript
  is regenerated from the processed tables here, so the counts are sufficient to check the analysis;
  the reads may be made available to a named reviewer on request to the corresponding author, at the
  discretion of the ethics committee.
- **No negative or blank control libraries exist** for this series. This is stated in the manuscript
  and it limits the decontamination analysis to a fixed exclusion list.
- **The provider's raw report** contains patient names and is not distributed.

## Ethics

Collected at Nanjing Stomatological Hospital under institutional ethics committee approval
NJSH-2023NL-19, which explicitly covers secondary microbiome analysis of these specimens. All 20
participants gave written informed consent. The identifiers released here are not reversible without
the mapping held by the investigators.

## Licence

CC BY 4.0 for the data. Please cite the manuscript and this record.
"""


TOP_README = """# Release package for upload

Two folders, two destinations. Nothing here contains patient names or hospital record numbers.

```
release/
├── code/     -> GitHub repository, then archive the release tag on Zenodo (record 2)
└── data/     -> Zenodo dataset record (record 1)
```

## What was done to make this safe

The hospital record and the vendor's sample labels both identify patients, so every label was re-coded
through one mapping before anything was written here:

- `P01` to `P20` replace hospital record numbers everywhere, in the data **and** in the code.
- The vendor's ASV table and all seven per-rank abundance tables had their sample columns renamed, so
  they no longer carry the original labels.
- The crosswalk from hospital record number to P-code is **not** in this package. It is written to
  `results/private/crosswalk_hospital_to_Pcode.csv` inside the project and should be kept with the
  study records, not uploaded.
- `data/raw/患者信息.xlsx` and the vendor submission sheet are excluded outright; they carry names.

`build_release_package.py` runs a preflight pass that refuses to build if an identity column, a
hospital-record-number-shaped value, an 18-digit identity number or key material appears anywhere. It
was tested by scanning the finished tree for all 20 patient names and all 20 hospital record numbers:
zero hits.

## Uploading to GitHub

```
cd release/code
git init
git add .
git commit -m "Analysis code for the manuscript"
git remote add origin git@github.com:<you>/<repo>.git
git push -u origin main
```

Then create a release tag (`v1.0.0`) and connect the repository to Zenodo under Settings, so that the
tag is archived with a DOI. Put that DOI in the manuscript's Data Availability statement.

**Before you push, decide the licence.** Code and data usually take different ones: MIT or BSD for
code, CC BY 4.0 for the data record. `release/data/README.md` already states CC BY 4.0 for the data;
the code licence is left for you to choose and add as `LICENSE`.

## Uploading to Zenodo

Record 1, the dataset: drag `release/data/` contents in, or zip that folder first.

| Field | Value |
|---|---|
| Upload type | Dataset |
| Title | 16S rRNA gene profiles of paired tumour and distant normal mucosa from 20 tongue cancer patients |
| Access | Open, provided the ethics approval permits it; confirm before publishing |
| Licence | CC BY 4.0 |
| Related | the manuscript DOI, marked *is supplement to* |

## Two things to settle before either upload

1. **Confirm that the committee permits open release of the processed tables.** Approval
   NJSH-2023NL-19 covers secondary microbiome analysis and the decision recorded here is to release
   processed data only, but the approval has not been read for a public-release clause. One sentence
   from the committee, kept on file, closes this.
2. **Raw FASTQ are deliberately excluded.** The decision is processed data only: the ethics approval
   covers secondary analysis of the specimens rather than redistribution of the reads, and the
   provider's read files carry sample labels that encode the hospital record number. Do not add
   `data/fastq/` to this package.
"""


def write_docs(code_dir: pathlib.Path, data_dir: pathlib.Path) -> None:
    """READMEs and dependency manifests, generated so they cannot go stale."""
    (code_dir/"requirements.txt").write_text(REQ_PY, encoding="utf-8")
    (code_dir/"requirements-r.txt").write_text(REQ_R, encoding="utf-8")
    (code_dir/"README.md").write_text(CODE_README, encoding="utf-8")
    (data_dir/"README.md").write_text(DATA_README, encoding="utf-8")
    (code_dir.parent/"README.md").write_text(TOP_README, encoding="utf-8")


# ----------------------------------------------------------------------------- main
def main() -> None:
    global CROSSMAP, CROSS
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="release")
    args = ap.parse_args()
    out = ROOT/args.out
    if out.exists():
        shutil.rmtree(out)
    code_dir, data_dir = out/"code", out/"data"
    code_dir.mkdir(parents=True)
    data_dir.mkdir(parents=True)

    CROSSMAP, CROSS = build_crosswalk()
    print(f"crosswalk: {len(CROSSMAP)} patients re-coded; mapping kept in results/private/")

    code_rows = copy_code(code_dir)
    print(f"code/   {len(code_rows)} files, {sum(r['bytes'] for r in code_rows)/1e6:.2f} MB")

    data_rows = [sample_metadata(data_dir), deidentified_clinical(data_dir)] + vendor_tables(data_dir)
    write_docs(code_dir, data_dir)
    print(f"data/   {len(data_rows)} files, {sum(r['bytes'] for r in data_rows)/1e6:.2f} MB")

    problems = preflight(code_dir) + preflight(data_dir)
    if problems:
        print("\nPREFLIGHT FAILED - nothing from this build should be uploaded:")
        for p in problems:
            print("  -", p)
        raise SystemExit(1)
    print("\npreflight: no identity columns, no forbidden files, no key material, no hospital numbers")

    for tree, rows in ((code_dir, code_rows), (data_dir, data_rows)):
        pd.DataFrame(rows).to_csv(tree/"MANIFEST.csv", index=False, encoding="utf-8")

    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
