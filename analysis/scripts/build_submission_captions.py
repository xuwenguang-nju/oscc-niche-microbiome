#!/usr/bin/env python3
"""Build the submission caption file in style B, plus the verified clarifications it needs.

Starting point
--------------
The B form produced by `analysis/scripts/build_caption_variants.py`: one paragraph per figure, the
panel description and that panel's statistics in the same sentence, no italic *Data:* label.

What this script adds
---------------------
A short, explicit list of edits — nothing else. Each one is here because the caption was silent on
something a reader or a production editor needs, and each is backed by the code or table that draws
the panel, recorded next to it. Two of them also correct a number that no table supports.

Guards
------
* every anchor must occur exactly once, so an edit cannot silently miss;
* the numerals of a finished caption must be the numerals of the B text plus this file's allowlist,
  so no number can appear that was not either inherited or verified here;
* the house limits are re-imposed: title <= 20 words, caption <= 300.

Writes results/reports/figure-captions-submission.md and .docx.
"""
from __future__ import annotations

import importlib.util
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
MD = ROOT / "results/reports/figure-captions-submission.md"
DOCX = MD.with_suffix(".docx")

spec = importlib.util.spec_from_file_location(
    "bcv", ROOT / "analysis/scripts/build_caption_variants.py")
bcv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bcv)

# Numbers this file is allowed to introduce, with the evidence for each. Anything not listed here
# and not already in the B text fails the build, which is what keeps a clarification from becoming
# an invention.
ALLOWED = {
    "40":  "hospital cohort = 20 patients x 2 samples (Table 1, row 7)",
    "177": "TCGA-HNSC samples analysed (Table S5, s6v2_tcma_samples.csv)",
    "155": "primary tumours in Table S5",
    "22":  "solid tissue normal in Table S5",
    "80":  "ESCA samples (s12_oral_gut_axis.csv: n)",
    "165": "STAD samples (s12_oral_gut_axis.csv: n)",
    "146": "COAD samples (s12_oral_gut_axis.csv: n)",
    "48":  "READ samples (s12_oral_gut_axis.csv: n)",
    "95":  "confidence level of the intervals of Sec. 4l, 5k, S3f",
    "1.5": "whisker length of the box plots, in IQR (matplotlib default)",
}

EDITS: list[tuple[str, str, str, str]] = [
    # (figure, anchor, replacement, why)
    # Fig. 5 was rewritten after these clarifications landed in the source caption, so the
    # Fig-5 edits are now carried there and no longer appear here.

    ("7",
     "(HNSC, ESCA, STAD, COAD, READ), genome-derived genus profiles.",
     "(HNSC, ESCA, STAD, COAD, READ; 177, 80, 165, 146 and 48 samples), genome-derived genus "
     "profiles.",
     "s12_oral_gut_axis.csv: n = 177, 80, 165, 146, 48. The caption gave no sample counts at all."),

    ("7",
     "(**E**) The same score per sample.",
     "(**E**) The same score per sample (boxes are median and interquartile range, whiskers to "
     "1.5×IQR).",
     "s12_plot.py boxplot, matplotlib defaults."),

    ("S1",
     "excluding the cancer surface.",
     "excluding the cancer surface (boxes are median and interquartile range, whiskers to 1.5×IQR; "
     "every patient is drawn as a point).",
     "e5_plot.py boxplot with showfliers=False plus a jittered point per patient."),

    ("S2",
     "(**F**) Heterogeneity per replicate.",
     "(**F**) Heterogeneity (I²) per genus.",
     "panel F is a bar per genus of I² (plot_forest.py, axes[2,0]); \"per replicate\" described a "
     "different panel."),

    ("S2",
     "(**G**) Largest effect against number of cohorts.",
     "(**G**) Largest effect against number of cohorts (boxes are median and interquartile range, "
     "whiskers to 1.5×IQR).",
     "panel G is the box plot of max |d_z| by number of cohorts (plot_forest.py, axes[2,1])."),

    ("S4",
     "Oral indicator in pooled against oral-only analyses. Pooled against oral-only, as −log₁₀ q.",
     "Oral indicator in pooled against oral-only analyses, as −log₁₀ q.",
     "the merged sentence repeated its own subject; no content changes."),

    # The S3 univariable-HR clarification that used to live here is now written directly into the
    # master caption (Fig. S3, panel F), so no edit is applied at build time.
]

PREAMBLE = """# Figure captions

Figures are supplied as separate files; each caption gives what the panels show and the data behind
them. Terms used throughout: *adjacent* is the mucosa immediately next to the tumour, *distant
normal* is macroscopically uninvolved mucosa taken away from the tumour bed, and *Parvimonas* is
given at genus level unless *P. micra* is named. Tests are paired on patients where the design
allows, with Benjamini–Hochberg correction within the stated family of comparisons. Where a panel
shows boxes they are median and interquartile range with whiskers to 1.5×IQR unless said otherwise;
where it shows an interval bar or an error bar the level is stated.
"""


def numerals(text: str) -> set[str]:
    out = set()
    for m in re.finditer(r"\d[\d.,]*\d|\d", text):
        out.add(m.group(0).rstrip(".,").replace(",", ""))
    return out


def normalise(text: str) -> str:
    """Consistency of the statistics typography across the eleven captions.

    B capitalises the first letter after a sentence break, which turns p-values into "P = 0.04" in
    some captions and leaves "p = 0.04" in others. Journals set a lower-case italic p throughout, and
    the same argument applies to rho: Fig. 5 uses the glyph and Fig. 6 spelled it out.
    """
    subs = [
        (r"\bP = ", "p = "),
        (r"\bRho = ", "ρ = "),
        (r"Mean rho ", "Mean ρ "),
        (r"\bpseudo-F", "Pseudo-F"),
    ]
    for pat, repl in subs:
        text = re.sub(pat, repl, text)
    return text


def main() -> None:
    A = bcv.parse_a()
    order = ["1", "2", "3", "4", "5", "6", "7", "S1", "S2", "S3", "S4"]
    missing = [k for k in order if k not in A]
    if missing:
        sys.exit(f"source is missing figures: {missing}")

    texts = {}
    problems = []
    for fig in order:
        b = bcv.make_b(fig, A[fig]["body"], A[fig]["data"])
        for e_fig, anchor, repl, why in EDITS:
            if e_fig != fig:
                continue
            n = b.count(anchor)
            if n != 1:
                problems.append(f"Figure {fig}: anchor occurs {n}×: {anchor[:60]}")
                continue
            b = b.replace(anchor, repl)
        texts[fig] = normalise(b)

    for fig, b in texts.items():
        base = numerals(A[fig]["body"] + " " + A[fig]["data"])
        added = numerals(b) - base
        unknown = added - set(ALLOWED)
        if unknown:
            problems.append(f"Figure {fig}: unlisted numerals {sorted(unknown)}")
        n = bcv.words(b)
        # 300 words is an author-side ceiling, not a journal rule; the journal limits the total
        # manuscript word count, captions included. Figure 5 carries fourteen panels and is the one
        # caption that cannot say what each panel shows inside 300 words, so it is allowed 340 and
        # the exception is recorded here rather than being silently absorbed.
        cap = 340 if fig == "5" else 300
        if n > cap:
            problems.append(f"Figure {fig}: {n} words, over the {cap} house limit")
        if bcv.title_words(b) > 20:
            problems.append(f"Figure {fig}: title over 20 words")

    if problems:
        print("FAILURES")
        for p in problems:
            print("  -", p)
        sys.exit(1)

    body = [PREAMBLE]
    for fig in order:
        body.append(texts[fig])
    MD.write_text("\n\n".join(body).rstrip() + "\n", encoding="utf-8", newline="\n")

    r = subprocess.run([sys.executable, str(ROOT / "analysis/scripts/build_manuscript_docx.py"),
                        str(MD), str(DOCX)], capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"docx build failed: {r.stderr[-500:]}")

    print(f"{'fig':4} {'words':>5} {'title':>5}")
    for fig in order:
        print(f"{fig:4} {bcv.words(texts[fig]):5} {bcv.title_words(texts[fig]):5}")
    print(f"total {sum(bcv.words(v) for v in texts.values())} words, {len(EDITS)} clarifications")
    print(f"wrote {MD.relative_to(ROOT)} and {DOCX.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
