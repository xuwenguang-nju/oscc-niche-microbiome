#!/usr/bin/env python3
"""Assemble one folder with the current figures and table files, for drafting a new manuscript.

Contents
--------
figures/main/              7 main figures, PNG (reading) + TIFF (production)
figures/supplementary/     4 supplementary figures, PNG + TIFF
figures/captions/          the caption list that goes with them, .md + .docx
tables/                    one Word file per table
README.md                  what is here, and which numbers belong to which item
inventory.csv              every file with size, sha256, and measured geometry

Files are copied, never moved: the working set in `figures/manuscript-v6/`, the submission
package's own table copies and its figure folders stay exactly as they are.

The script measures what it copies (dpi, pixel size, physical width) and stops if a figure is
missing a TIFF or carries no dpi tag, so the kit cannot silently contain a superseded PNG.

Run: python analysis/scripts/build_drafting_kit.py
"""
from __future__ import annotations

import csv
import hashlib
import pathlib
import shutil

from docx import Document
from PIL import Image

ROOT = pathlib.Path(__file__).resolve().parents[2]
FIG = ROOT / "figures/manuscript-v6"
KIT = ROOT / "submission/drafting-kit"

MAIN = [
    ("Fig1-design-and-cohort-overlap", "Study design, cohort coverage, overlap with a pooled meta-analysis", 6),
    ("Fig2-identifiability", "Which quantities the design can and cannot resolve", 9),
    ("Fig3-six-niche-attribution", "Six-niche source attribution (PRJNA866676)", 9),
    ("Fig4-replication-and-meta-analysis", "Replication, species resolution, pooled meta-analysis", 13),
    ("Fig5-hospital-paired-cohort", "The independent hospital paired cohort", 14),
    ("Fig6-platform-concordance", "TCGA genome-derived profiles vs 16S", 12),
    ("Fig7-pan-cancer", "Pan-cancer context: the signal follows the site", 9),
]

SUPP = [
    ("FigS1-niche-functional", "Functional prediction and source distance in the six niches", 3),
    ("FigS2-cross-cohort-meta", "Detail of the cross-cohort meta-analysis", 8),
    ("FigS3-clinical-negative", "Clinical and survival analysis, all negative", 11),
    ("FigS4-tcma-host-coupling", "Direction, site and platform detail behind Fig. 6", 7),
]

# table file -> (number, what it holds, first cited in the current Results)
TABLES = [
    ("Table1", "Table 1", "The seven datasets: sampling, sequencing and processing", "Results 1"),
    ("TableS1", "Table S1", "The 13 genera replicating across 3+ cohorts, effects in log\u2082",
     "Results 4"),
    ("TableS2", "Table S2", "Clinical characteristics of the 20 hospital patients", "Results 5"),
    ("TableS3", "Table S3", "(a) genus assignments under 4 oxygen-tolerance definitions; "
                             "(b) pooled axis estimate under each definition", "Results 5"),
    ("TableS4", "Table S4", "(a) every taxonomic rank tested: 854 tests, 28 passing FDR; "
                             "(b) those 28 with effect and direction", "Results 5"),
    ("TableS5", "Table S5", "The 46 PICRUSt2 KEGG pathways differing between tumour and normal",
     "Results 5"),
    ("TableS6", "Table S6", "Per-sample annotation of the 177 TCGA-HNSC samples", "Results 6"),
]

# the submission-ready caption set (style B, with the verified clarifications), written by
# analysis/scripts/build_submission_captions.py. Copied to the kit root under a neutral name.
CAPTIONS = [(ROOT / "results/reports/figure-captions-submission.md", "figure-captions.md"),
            (ROOT / "results/reports/figure-captions-submission.docx", "figure-captions.docx")]

TABDIR = ROOT / "submission/jom-2026/04-tables"
MAX_WIDTH_MM = 182.0
rows: list[dict] = []


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def record(src: pathlib.Path, dst: pathlib.Path, note: str, kind: str) -> None:
    im = Image.open(src) if src.suffix.lower() in (".png", ".tif", ".tiff") else None
    if im is not None:
        dpi = float(im.info.get("dpi", (0, 0))[0] or 0)
        if not dpi:
            raise SystemExit(f"{src.name} carries no dpi tag, so its physical size is unknown")
        rows.append({
            "kind": kind, "file": dst.relative_to(KIT).as_posix(), "note": note,
            "format": im.format, "mode": im.mode, "dpi": f"{dpi:.0f}",
            "px": f"{im.width}x{im.height}",
            "width_mm": f"{im.width / dpi * 25.4:.1f}",
            "height_mm": f"{im.height / dpi * 25.4:.1f}",
            "bytes": dst.stat().st_size, "sha256": sha256(dst),
        })
        if im.format in ("TIFF", "PNG") and im.width / dpi * 25.4 > MAX_WIDTH_MM + 0.5:
            raise SystemExit(f"{src.name} is {im.width / dpi * 25.4:.1f} mm wide, over {MAX_WIDTH_MM}")
    else:
        rows.append({
            "kind": kind, "file": dst.relative_to(KIT).as_posix(), "note": note,
            "format": dst.suffix.lstrip("."), "mode": "", "dpi": "", "px": "",
            "width_mm": "", "height_mm": "",
            "bytes": dst.stat().st_size, "sha256": sha256(dst),
        })


def copy_figures(folder: str, items) -> None:
    for stem, _title, _panels in items:
        out = KIT / "figures" / folder
        out.mkdir(parents=True, exist_ok=True)
        for ext in (".png", ".tif"):
            src = FIG / folder / f"{stem}{ext}"
            if not src.is_file():
                raise SystemExit(f"missing figure: {src}")
            dst = out / f"{stem}{ext}"
            shutil.copyfile(src, dst)
            record(src, dst, f"{folder}: {_title}; panels A-{chr(64 + _panels)}", "figure")


def copy_tables() -> None:
    out = KIT / "tables"
    out.mkdir(parents=True, exist_ok=True)
    for stem, number, what, cited in TABLES:
        src = TABDIR / f"{stem}.docx"
        if not src.is_file():
            raise SystemExit(f"missing table: {src}")
        dst = out / f"{stem}.docx"
        shutil.copyfile(src, dst)
        # the caption paragraph must name the table these files sit under
        cap = next((p.text.strip() for p in Document(dst).paragraphs if p.text.strip()), "")
        if not cap.startswith(f"{number}."):
            raise SystemExit(f"{stem}.docx opens with {cap[:50]!r}, expected '{number}.'")
        record(src, dst, f"{number}: {what}; first cited in {cited}", "table")


def copy_captions() -> None:
    for src, name in CAPTIONS:
        if not src.is_file():
            raise SystemExit(f"missing caption file: {src}")
        dst = KIT / name
        shutil.copyfile(src, dst)
        record(src, dst, f"the 11 captions, in upload order (Fig. 1-7, S1-S4), from {src.name}",
               "caption")


def write_readme() -> None:
    main_lines = "\n".join(
        f"| {i} | `{stem}.tif` / `.png` | {t} | A–{chr(64 + n)} |"
        for i, (stem, t, n) in enumerate(MAIN, 1))
    supp_lines = "\n".join(
        f"| S{i} | `{stem}.tif` / `.png` | {t} | A–{chr(64 + n)} |"
        for i, (stem, t, n) in enumerate(SUPP, 1))
    tab_lines = "\n".join(
        f"| **{number}** | `{stem}.docx` | {what} | {cited} |"
        for stem, number, what, cited in TABLES)

    # Heights are measured from the delivered files, not asserted: it was a prose claim ("281 mm")
    # that went stale the moment Fig. S1 was re-laid out.
    geo = []
    for stem, *_ in MAIN + SUPP:
        im = Image.open(KIT / "figures" / (
            "main" if any(stem == m[0] for m in MAIN) else "supplementary") / f"{stem}.png")
        dpi = float(im.info["dpi"][0])
        geo.append((stem, im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4))
    tallest = max(geo, key=lambda g: g[2])
    widest = max(geo, key=lambda g: g[1])
    geo_line = (f"Measured from the files: widths {min(g[1] for g in geo):.0f}–"
                f"{widest[1]:.0f} mm (two-column limit 182 mm, widest is {widest[0]}); "
                f"heights {min(g[2] for g in geo):.0f}–{tallest[2]:.0f} mm "
                f"(tallest is {tallest[0]}). An A4 text block is about 240 mm, so "
                + ("every figure fits." if tallest[2] <= 240
                   else f"{sum(1 for g in geo if g[2] > 240)} figure(s) need splitting."))

    (KIT / "README.md").write_text(f"""# Drafting kit — figures and tables for the new manuscript

Built 2026-11-15 from the working set, copied (not moved), so `figures/manuscript-v6/` and
`submission/jom-2026/` are untouched. Rebuild with
`python analysis/scripts/build_drafting_kit.py`.

**Main figures: 7 · Supplementary figures: 4 · Main table: 1 · Supplementary tables: 6**

## Figures

`figures/main/` — 7 main figures. Every one exists as a **TIFF** (600 dpi, LZW, RGB) for production
and a same-sized **PNG** for reading and for pasting into a draft.

| No. | File | What it shows | Panels |
|---|---|---|---|
{main_lines}

`figures/supplementary/` — 4 supplementary figures, same two formats.

| No. | File | What it shows | Panels |
|---|---|---|---|
{supp_lines}

{geo_line}

`figure-captions.md` / `.docx` (kit root) — the 11 captions in upload order, one paragraph each:
the panel description and that panel's statistics in the same sentence, which is the form journals
print. Written against the rendered panels, panel by panel, and then checked twice — every panel
letter against the figure, every numeral against the previous caption file plus an explicit list of
verified additions. If a panel changes, its caption has to change with it; the build
(`analysis/scripts/build_submission_captions.py`) is what enforces that.

## Tables

`tables/` — one Word file per table, self-contained: caption above the table, footnotes below, a
repeating header row, and a fixed column layout. All landscape A4.

| No. | File | What it holds | First cited in |
|---|---|---|---|
{tab_lines}

Numbering is by first citation in the text, and the supplement follows the same rule:
Table 1 → S1 → S2 → S3 → S4 → S5 → S6. Table 1 is the only main table.

`inventory.csv` lists every file here with its size, sha256 and measured geometry, so you can tell
a stale copy from the current one without opening anything.

## What changed after the supervisor's review

The kit was rebuilt once the fourteen review points were worked through. What a drafter has to know:

- **One score, one definition.** The main oxygen-tolerance score is the *taxonomy* set of Table S3 —
  **29 obligate-anaerobe genera against 16 aerobic or facultative**, *Gemella* excluded — and it is
  the same set in Fig. 5B, 5E, 5F, 5K and 5L. It was three different sets behind one figure before.
  Pooled shift on this score is +0.29 (95% CI 0.15 to 0.44), I² 25%.
- **Fig. 5E and 5F** plot the composite score's paired p and AUC, not the anaerobe share and not the
  twelve down-moving genera, so the depth series tests the quantity panel B reports.
- **Fig. 5J** has log₂ on both axes, so its identity line is between comparable units; it was ln
  against log₂. **Fig. 5D** y-axis is "% of samples".
- **Fig. 4L, Table S1 and Methods** are in **log₂**. The pooled effect was always a natural-log ratio
  mislabelled "log2FC"; it was converted, not re-estimated.
- **Fig. 7** panel titles no longer overstate, and the span between project means is 127 points.
- **Fig. 1C** separates the genus-level replication (*Parvimonas*, 4 cohorts) from the species-level
  result (*P. micra*, 1 cohort).
- **Table S3** carries all four definitions genus by genus; **Table 1** records what the provider
  actually delivered (ASV table with sequences, all seven ranks, PICRUSt2).
- **Fig. S1 was re-laid out.** Three stacked panels made it 281 mm tall, above an A4 text block;
  A and B now share a row and C keeps full width, giving 168 mm at the same 171 mm width and the
  same font sizes.
- **A second review round changed seven figures and every caption.** Fig. 4B is drawn and now
  labelled as the species (*P. micra*, denominator all bacterial reads) rather than the genus;
  Fig. 4A carries one chance-rank line per cohort, since the two have 2 and 4 candidate sources;
  Fig. 5H's left axis limit no longer clips the single phylum-level FDR hit; Fig. 6C prints the full
  accessions; Fig. 6K is on a linear axis, because a log axis silently dropped the five genera with a
  zero median, *Parvimonas* among them; Fig. S2A and S2C are converted to log₂ so they match
  Table S1; Fig. S3G's legend says "adjusted" because the model holds stage, age and sex. Table S1's
  prior-report column is now generated from the ORMCKB survey instead of hand-copied.
- **A third pass settled the remaining taxonomy and sample-set questions.** Fig. 4J's axis now names
  the species (*P. micra*), matching its input table, rather than the genus. Fig. 6K's x axis says
  RNA-seq cohort, not TCMA: both of its axes come from GSE227919, so the panel is not a
  cross-substrate comparison. Fig. 6's caption now states that E and K are RNA-seq panels. Fig. 4's
  caption gives the sample size per panel rather than one figure-wide n, because its panels come from
  three cohorts. Fig. S4D is described as mixing site-group and tongue-subset strata, which is what
  it draws, and no longer as "the same estimates" as Fig. S3J. Fig. 2A's caption reports the 0.15
  crossing as specific to the plotted scenario, with the range across the four scenarios, instead of
  a general detectability threshold. The survival sample sets (134/73 for Kaplan–Meier and
  univariable, 106/54 for the adjusted models, 28 patients lost to missing stage) are now stated in
  Methods and summarised in the Fig. S3 caption.
- **Two tables are new.** `TableS4.docx` holds every taxonomic rank tested (854 tests, the 28 that
  pass FDR, and each with its effect and direction) — the six per-rank fractions in Results 5 were
  checkable nowhere else. `TableS5.docx` was the KEGG table and moved up from S4 so numbering still
  follows first citation; `TableS6.docx` is the TCGA annotation, formerly S5.

## One thing to decide when you start drafting

1. **File names.** They are kept as they come out of the pipeline (`Fig4-replication-and-meta-analysis`
   and so on) so any file still traces back to `figure_manifest.csv` and to the script that drew it.
   Renaming to `Fig1`…`Fig7` is safe if you prefer it; nothing downstream reads these names.

## What is deliberately not here

- Source data and the analysis tables behind every panel (`results/tables/`) — those are for the
  deposition record, not for the draft.
- The assembled manuscript and its section files (`results/reports/`), which are what you are about to
  replace.
- The previous nine-figure set still sitting in `submission/jom-2026/02-figures-main/`. It is
  superseded; do not draft against it.
""", encoding="utf-8", newline="\n")


def main() -> None:
    if KIT.exists():
        shutil.rmtree(KIT)
    KIT.mkdir(parents=True)

    copy_figures("main", MAIN)
    copy_figures("supplementary", SUPP)
    copy_captions()
    copy_tables()

    with open(KIT / "inventory.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    write_readme()

    figs = [r for r in rows if r["kind"] == "figure"]
    tabs = [r for r in rows if r["kind"] == "table"]
    caps = [r for r in rows if r["kind"] == "caption"]
    print(f"figures {len(figs)} ({sum(r['format'] == 'TIFF' for r in figs)} TIFF, "
          f"{sum(r['format'] == 'PNG' for r in figs)} PNG) | tables {len(tabs)} | captions {len(caps)}")
    print(f"total {sum(r['bytes'] for r in rows) / 1e6:.1f} MB -> {KIT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
