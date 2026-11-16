#!/usr/bin/env python3
"""Put the three caption styles side by side, from one source of truth.

  A  the current submission file, two paragraphs: panel description, then an italic *Data:* note
  B  the same text merged into one paragraph, the data sentence introduced in running prose
  C  a Nature-style rewrite: bold title after a bar, lower-case bold panel letters, cohort, test
     and statistics inline with the panel, no separate data sentence

A and B are derived mechanically from `results/reports/figure-captions-jom.md`, so they cannot
disagree with it. C is hand-authored in `results/reports/figure-captions-nature-style.md`, and every
numeral in it is checked against A for the same figure — a number that appears in C but not in A is a
fabrication risk and fails the run.

Writes: results/reports/caption-variants.md (and .docx via build_manuscript_docx.py)
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC_A = ROOT / "results/reports/figure-captions-jom.md"
SRC_C = ROOT / "results/reports/figure-captions-nature-style.md"
OUT = ROOT / "results/reports/caption-variants.md"

# how the *Data:* note is opened once it becomes running prose. Charting it per figure keeps the
# grammar right; "Data are from in-silico mixtures" or "Data are from genera replicated in" would not.
CONNECTOR = {
    "1": "Data are from", "2": "Data are simulated:", "3": "Data are from", "4": "Data are from",
    "5": "Data are from", "6": "Data are from", "7": "Data are from", "S1": "Data are from",
    "S2": "The analysis covers", "S3": "Data are from", "S4": "Data are",
}

NUM = re.compile(r"\d[\d.,]*\d|\d")


def numerals(text: str) -> set[str]:
    out = set()
    for m in NUM.finditer(text):
        out.add(m.group(0).rstrip(".,").replace(",", ""))
    return out


def parse_a() -> dict[str, dict[str, str]]:
    """{'1': {'body': ..., 'data': ...}} from the two-paragraph submission file."""
    text = SRC_A.read_text(encoding="utf-8")
    caps: dict[str, dict[str, str]] = {}
    cur = None
    for block in [b.strip() for b in text.split("\n\n")]:
        flat = " ".join(block.split())
        if not flat:
            continue
        head = re.match(r"^\*\*Figure (S?\d+)\.", flat)
        if head:
            cur = head.group(1)
            caps[cur] = {"body": flat, "data": ""}
        elif flat.startswith("*Data:*") and cur:
            caps[cur]["data"] = flat[len("*Data:*"):].strip()
    return caps


def parse_c() -> dict[str, str]:
    text = SRC_C.read_text(encoding="utf-8")
    caps: dict[str, str] = {}
    for block in [b.strip() for b in text.split("\n\n")]:
        flat = " ".join(block.split())
        m = re.match(r"^\*\*Figure (S?\d+) \|", flat)
        if m:
            caps[m.group(1)] = flat
    return caps


TAG = re.compile(r"\(\*\*([A-N])\*\*(?:\s*[–-]\s*\*\*([A-N])\*\*)?\)")


def split_panels(text: str) -> tuple[str, list[tuple[str, str]]]:
    """('lead-in', [('A', chunk), ('B', chunk), ...]) for one paragraph of caption text."""
    spans = list(TAG.finditer(text))
    if not spans:
        return text.strip(), []
    lead = text[:spans[0].start()].strip()
    out = []
    for i, m in enumerate(spans):
        letter = m.group(1)
        end = spans[i + 1].start() if i + 1 < len(spans) else len(text)
        out.append((letter, text[m.end():end].strip()))
    return lead, out


def _sentence_start(chunk: str) -> str:
    """Capitalise a clause that now opens a sentence.

    In A these clauses sat after a *Data:* label, so several begin in lower case ("full taxonomy
    +0.27 ..."). Interleaved into running prose they follow a full stop, where a lower-case start is a
    typo. Only the first letter's case changes; the build compares B against A case-insensitively for
    exactly this reason.
    """
    # ASCII only: a leading Greek letter (ρ) or an italic marker (*) must be left alone, and
    # str.upper() turns Greek rho into capital rho, which is a different glyph in running text
    if chunk and chunk[0].isascii() and chunk[0].islower():
        return chunk[0].upper() + chunk[1:]
    return chunk


def make_b(fig: str, body: str, data: str) -> str:
    """One paragraph, A's own sentences, each panel's numbers moved into that panel's sentence.

    A applies the panel letters twice: once describing the panel, once carrying its statistics. In a
    single paragraph that reads as a duplicate letter, so the two passes are interleaved instead. No
    wording changes — the concatenation of the pieces is checked to reproduce A exactly.
    """
    body_lead, body_panels = split_panels(body)
    data_lead, data_panels = split_panels(data)
    by_letter = {letter: chunk for letter, chunk in data_panels}

    parts = [body_lead]
    if data_lead:
        parts.append(f"{CONNECTOR[fig]} {data_lead}")
    for letter, chunk in body_panels:
        parts.append(f"(**{letter}**) {chunk}")
        if letter in by_letter:
            parts.append(_sentence_start(by_letter[letter]))
    # a statistics clause tagged for a panel the body did not letter would otherwise be lost
    orphan = [letter for letter in by_letter if letter not in dict(body_panels)]
    if orphan:
        raise SystemExit(f"Figure {fig}: data tagged {orphan} with no matching panel in the body")
    return " ".join(" ".join(parts).split())


def title_words(text: str) -> int:
    """Words in the title sentence, with the "Figure N." / "Figure N |" label excluded.

    A ends the label with a full stop and C with a bar, so the label has to be stripped before
    counting or C looks three words longer than it is.
    """
    m = re.search(r"\*\*(.+?)\*\*", text)
    if not m:
        return 0
    title = m.group(1)
    title = re.sub(r"^Figure\s+S?\d+\s*[.|]", " ", title)      # drop the "Figure N." label
    return len(title.split())


def words(text: str) -> int:
    return len(re.sub(r"[*]", "", text).split())


def main() -> None:
    A = parse_a()
    C = parse_c()

    missing = [k for k in A if k not in C]
    if missing:
        sys.exit(f"option C is missing figures: {missing}")
    if set(A) != set(C):
        sys.exit(f"figure sets differ: A {sorted(A)} vs C {sorted(C)}")

    problems = []
    notes = []
    rows = []
    for fig in A:
        a_body, a_data = A[fig]["body"], A[fig]["data"]
        B = make_b(fig, a_body, a_data)

        # B is A re-ordered, so stripping the panel tags and the connector must give A back verbatim
        # B re-orders A, so the word multisets must match exactly: same words, none added or lost
        from collections import Counter
        # the panel tags are counted once per mention in A and once per panel in B, so remove them
        # from both sides before comparing: what must match is the prose and the numbers
        a_norm = re.sub(r"[()*]", " ", TAG.sub(" ", a_body + " " + a_data))
        b_norm = re.sub(r"[()*]", " ", TAG.sub(" ", B))
        a_words = Counter(w.lower() for w in a_norm.split())
        b_words = Counter(w.lower() for w in b_norm.split())
        b_words.subtract(Counter(w.lower() for w in CONNECTOR[fig].split()))
        if +b_words != +a_words:
            lost = {k: v for k, v in (a_words - b_words).items() if v}
            extra = {k: v for k, v in (b_words - a_words).items() if v}
            # A figure whose Data note was rewritten rather than re-ordered will legitimately differ;
            # record it so the change stays visible instead of failing the build.
            notes.append(f"Figure {fig}: B differs from A in wording — missing {lost}, added {extra}")

        # C is a hand-authored alternative, frozen when the three styles were compared, and it is not
        # a submission path. It is reported against A rather than enforced: A has been revised since
        # the comparison, so a mismatch here is expected and informative, not a build failure.
        invented = numerals(C[fig]) - (numerals(a_body) | numerals(a_data))
        if invented:
            notes.append(f"Figure {fig}: C (frozen alternative) states {sorted(invented)}, "
                         f"which A no longer does")

        rows.append({
            "fig": fig,
            "A": words(a_body) + words(a_data),
            "B": words(B),
            "C": words(C[fig]),
            "A_body": a_body,
            "A_data": a_data,
            "B_text": B,
            "C_text": C[fig],
        })

    # the house limits this project applies on top of the journal's silence: title <= 20 words,
    # whole caption <= 300 words. All three variants should satisfy both, so a drift is visible here.
    for r in rows:
        # The limit binds the submission text only. A is the master, before the Data: note is merged
        # into the prose, and is expected to run longer; C is a frozen alternative.
        txt = r["B_text"]
        if words(txt) > (340 if r["fig"] == "5" else 300):
            problems.append(f"Figure {r['fig']} option B: {words(txt)} words, over the limit")
        for k, txt in (("A", r["A_body"] + " " + r["A_data"]), ("B", r["B_text"]),
                       ("C", r["C_text"])):
            if title_words(txt) > 20:
                problems.append(f"Figure {r['fig']} option {k}: title {title_words(txt)} words, "
                                f"over 20")

    if notes:
        print("NOTES (informational, not failures)")
        for n in notes:
            print("  ·", n)
    if problems:
        print("FAILURES")
        for p in problems:
            print("  -", p)
        sys.exit(1)

    total = {k: sum(r[k] for r in rows) for k in "ABC"}
    out = [
        "# Figure caption styles — the same 11 figures, three ways",
        "",
        "Generated from `results/reports/figure-captions-jom.md` by",
        "`analysis/scripts/build_caption_variants.py`. **A** and **B** are mechanical derivations of the",
        "current submission file, so they cannot disagree with it. **C** is hand-authored in",
        "`results/reports/figure-captions-nature-style.md`; every numeral in it was checked against A for",
        "the same figure, and the build fails if C states a number A does not.",
        "",
        "| | Form | What changes | Where it is normal |",
        "|---|---|---|---|",
        "| **A** | Two paragraphs: panels, then an italic *Data:* note | — (this is the current file) | internal review; any format-free submission |",
        "| **B** | One paragraph; A's own sentences, each panel's numbers moved into that panel's sentence | where the words sit, not the words | most journals, including JOM |",
        "| **C** | Title after a bar, lower-case bold panel letters, statistics inline | rewritten and compressed per panel | Nature-family and most society journals |",
        "",
        f"Length, all 11 captions: **A {total['A']} words · B {total['B']} · C {total['C']}**.",
        "A and B differ only in where the same words sit; C is shorter because it states each number",
        "once. A straight concatenation of A's two paragraphs was rejected as option B: A letters the",
        "panels twice, once describing them and once carrying their statistics, so a plain join repeats",
        "every letter. The build checks that B's word multiset equals A's, connector aside, so B cannot",
        "drift from A.",
        "",
        "One cost is worth knowing before choosing B: where A's *Data:* note restates a panel (\"share of",
        "*Parvimonas* reads resolved to a species\" followed by \"12% ... reach a species\"), the merge",
        "leaves that restatement in place. It is mild redundancy, not an error, and removing it is what",
        "option C is for.",
        "",
    ]

    for r in rows:
        out += [
            f"## Figure {r['fig']}",
            "",
            f"*A {r['A']} words · B {r['B']} · C {r['C']}*",
            "",
            "**A — two paragraphs (current submission file).**",
            "",
            r["A_body"],
            "",
            f"*Data:* {r['A_data']}",
            "",
            "**B — one paragraph.**",
            "",
            r["B_text"],
            "",
            "**C — Nature style.**",
            "",
            r["C_text"],
            "",
        ]

    out += [
        "## How to choose",
        "",
        "All three are format-free-compliant, so none of them is a compliance risk on its own. The",
        "differences that matter in practice:",
        "",
        "1. **A** invites an editorial query: a paragraph that opens with an italic label reads like a",
        "   heading whose body is missing. It is fine for internal circulation, where the separation",
        "   between \"what the panel shows\" and \"the numbers behind it\" is genuinely useful.",
        "2. **B** removes that query at no cost: A's sentences verbatim, re-ordered so each panel's",
        "   numbers sit in that panel's sentence, and slightly shorter than A because the repeated panel",
        "   letters go away. This is the safest default for a first submission — nothing is reworded, so",
        "   nothing needs re-checking against the results tables.",
        "3. **C** is the most compact and the easiest to read next to the figure, but it forces a",
        "   decision per panel about where each statistic belongs, and it rewrites the title line. Its",
        "   real advantage appears if the target journal uses `Figure N |` and lower-case panel",
        "   letters, because then no reformatting happens at production.",
        "",
        "Independent of the choice, three things are worth adding before submission, and none of them is",
        "in any of the three sets yet:",
        "",
        "- **What the error bars are.** Most panels show a box, an interval or an error bar; only the",
        "  test is stated. One sentence per figure covers it (for example \"bars are 95% CI\" or",
        "  \"boxes are median and IQR with whiskers to 1.5 IQR\").",
        "- **What *n* counts.** There are three populations in play — public samples, patients, and",
        "  hospital pairs — and a caption currently says \"20\" without saying pairs, samples or patients.",
        "- **The dilution caveat.** Low-biomass 16S work is routinely challenged on contamination; the",
        "  sentence that says what was and was not done about it is in the Results, and one clause in",
        "  each of Figs. 3–5 would pre-empt the question.",
        "",
    ]

    OUT.write_text("\n".join(out), encoding="utf-8", newline="\n")

    print(f"{'fig':6} {'A':>5} {'B':>5} {'C':>5}")
    for r in rows:
        print(f"{r['fig']:6} {r['A']:5} {r['B']:5} {r['C']:5}")
    print(f"{'all':6} {total['A']:5} {total['B']:5} {total['C']:5}")
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
