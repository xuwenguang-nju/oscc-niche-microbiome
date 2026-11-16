#!/usr/bin/env python3
"""Bring the two captions that the panel corrections pushed over the 300-word house limit back
under it, and fix the Fig. 5 panel list (E and F describe the wrong panels; panel L has four rows).

House limits: title sentence <= 20 words, whole caption <= 300 words including the Data: note.
Both are the author's instruction, not the journal's.

The edits are applied to the *flattened* caption block and the block is re-wrapped afterwards, so a
hard wrap inside a phrase cannot make a replacement miss.
"""
import pathlib
import re
import textwrap

P = pathlib.Path("results/reports/figure-captions-jom.md")
text = P.read_text(encoding="utf-8")
MISSES: list[tuple[str, str]] = []


def block(fig: str) -> tuple[int, int]:
    """Span of the caption paragraph plus its Data: paragraph for one figure."""
    start = text.index(f"**Figure {fig}.")
    end = text.index("\n\n---", start)
    return start, end


def edit(fig: str, pairs: list[tuple[str, str]]) -> None:
    global text
    s, e = block(fig)
    flat = " ".join(text[s:e].split())
    for old, new in pairs:
        if old in flat:
            flat = flat.replace(old, new, 1)
        elif new in flat:
            pass                       # already applied -- the script is re-runnable
        elif any(later_new in flat for _, later_new in pairs):
            # superseded by a later entry in the same list (two passes over the same phrase)
            MISSES.append((fig, old))
        else:
            raise AssertionError(f"Figure {fig}: neither form present:\n  {old[:140]}")
    # The body sentence and the *Data:* note are two paragraphs in the .md. Flattening them to
    # search must not merge them, so the paragraph break is restored before wrapping.
    body, sep, data = flat.partition("*Data:*")
    wrapped = "\n".join(textwrap.wrap(body.strip(), width=100, break_long_words=False,
                                      break_on_hyphens=False))
    if sep:
        wrapped += "\n\n*Data:* " + "\n".join(textwrap.wrap(data.strip(), width=100,
                                                            break_long_words=False,
                                                            break_on_hyphens=False))
    text = text[:s] + wrapped + text[e:]


edit("2", [
    ("across four leakage scenarios, with the estimate",
     "across four scenarios, with the estimate"),
    ("— the objective closest to zero is χ², the one with the largest bias is L1.",
     "— χ² is closest to zero, L1 furthest."),
    ("(**D**) The standard deviation of the estimate falls",
     "(**D**) The SD of the estimate falls"),
    ("(**H**) Specificity is 1.0 at every true fraction and sensitivity 1.0 wherever novel taxa are present,",
     "(**H**) Specificity is 1.0 throughout and sensitivity 1.0 wherever novel taxa are present,"),
])

edit("5", [
    # E is the paired p-value curve and F the AUC; D also shows the threshold counts.
    ("(**D**) Genera detected as reads are subsampled. (**E**) Genera passing FDR as reads are "
     "subsampled. (**F**) *Parvimonas* and the anaerobe score as reads are subsampled.",
     "(**D**) Genera detected, and genera passing nominal or FDR thresholds, as reads are subsampled. "
     "(**E**) Paired Wilcoxon p for *Parvimonas* and for the anaerobe score over the same depth "
     "series. (**F**) AUC of the anaerobe score and of *Parvimonas* alone over that series."),
    ("reads, against their abundance.", "reads, against abundance."),
    ("(**J**) Agreement between this cohort and the public meta-analysis.",
     "(**J**) Agreement with the public meta-analysis."),
    ("Pooled estimates under three genus-set definitions.",
     "Pooled estimates under four genus-set definitions."),
    ("(**N**) The same pathways after adjusting for composition.",
     "(**N**) The same after adjusting for composition."),
    ("*Data:* 20 paired hospital samples,", "*Data:* 20 paired samples,"),
    ("(**B**) On the strict oxygen-tolerance sets that this panel plots, anaerobes 46.1% → 63.5% of "
     "reads and aerobes 37.5% → 21.4%.",
     "(**B**) Strict oxygen-tolerance sets: anaerobes 46.1% → 63.5% of reads, aerobes 37.5% → 21.4%."),
    ("(**C**) *F. nucleatum* 2.66% → 6.69% of reads, 15 of 20 patients higher;",
     "(**C**) *F. nucleatum* 2.66% → 6.69%, 15 of 20 higher;"),
    ("(**D**, **E**) Depth series from 500 reads to full.",
     "(**D**, **E**, **F**) Depth series from 500 reads to full."),
    ("(**G**) 21 genera retained and 6 lost at 500 reads, median abundance 1.43% against 0.09%.",
     "(**G**) 21 genera retained, 6 lost, median abundance 1.43% against 0.09%."),
    ("(**H**) 12 of 120 family tests pass FDR. (**I**) 1 of 20 phyla.",
     "(**H**) 12 of 120 family tests. (**I**) 1 of 20 phyla."),
    ("ρ = 0.21 across the 120 shared genera, rising to 0.76 among the 41 prevalent in at least half "
     "the samples on both sides; *Parvimonas* itself is +1.21 here against a pooled +0.76 "
     "(0.46–1.05) in the meta-analysis.",
     "ρ = 0.21 across 120 shared genera and 0.76 among the 41 prevalent on both sides; *Parvimonas* "
     "+1.21 here against a pooled +0.76 (0.46–1.05)."),
    ("original sets +0.59 (0.18–0.99) and shown only for transparency; the intervals use a normal "
     "approximation, and eHOMD",
     "original +0.59 shown only for transparency; normal-approximation intervals, and eHOMD"),
    # second pass: the first round left the caption at 305 words
    ("(**A**) Paired ordination of the 20 patients; each arrow joins a patient's normal mucosa to "
     "that patient's tumour.",
     "(**A**) Paired ordination of the 20 patients; each arrow joins normal mucosa to that tumour."),
    ("(**G**) Genera lost and retained at 500 reads, against abundance.",
     "(**G**) Genera lost and retained at 500 reads."),
    ("(**D**, **E**, **F**) Depth series", "(**D**–**F**) Depth series"),
    ("(**M**) 46 of 240 KEGG pathways pass FDR.", "(**M**) 46 of 240 KEGG pathways."),
])

P.write_text(text, encoding="utf-8", newline="\n")

# ---- verify the house limits on the edited file
blocks = [b.strip() for b in P.read_text(encoding="utf-8").split("\n\n")]
caps = {}
cur = None
for b in blocks:
    flat = " ".join(b.split())
    if flat.lstrip("*").startswith("Figure "):
        cur = flat.split(".")[0].lstrip("*")
        caps[cur] = [flat]
    elif flat.startswith("*Data:") and cur:
        caps[cur].append(flat)

ok = True
for k, v in caps.items():
    title = v[0].split(".**", 1)[0].replace("**", "")
    body = re.sub(r"[*]", "", " ".join(v))
    n_body = len(body.split())
    n_title = len(" ".join(title.split()[2:]).split())          # drop the "Figure N." label
    good = n_title <= 20 and n_body <= 300
    ok &= good
    print(f"{k:22} title={n_title:3}w caption={n_body:4}w {'OK' if good else 'OVER'}")
print("all within limits:", ok)
for fig, old in MISSES:
    print(f"note: Figure {fig} entry superseded by a later edit: {old[:70]}")
assert ok
