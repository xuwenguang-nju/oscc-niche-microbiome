#!/usr/bin/env python3
"""Figure 1 — study design: cohort x niche coverage, analysis flow, evidence map.

Why this figure exists: the paper's central constraint is that the six cohorts sample
different, overlapping sets of ecological niches, and only one of them can answer the
source-attribution question. Readers need to see that before reading any result, otherwise
"replicated in 1 of 2 cohorts" reads as weakness rather than as a design fact.

Panels
  a  cohort x niche coverage matrix (the hook: what data exist, and the gap)
  b  five-stage analysis flow, with the cohorts feeding each stage
  c  evidence map: per conclusion, how many cohorts could test it and how many confirmed it

Each panel is a separate helper so a panel can be re-rendered alone while iterating.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
STYLE_DIR = _FIGSTYLE
if STYLE_DIR.is_dir():
    sys.path.insert(0, str(STYLE_DIR))
from runtime import apply_figure_style, panel_crops  # noqa: E402
from figlib import panel_letter_fixed as panel_letter  # noqa: E402


def locate_root() -> Path:
    here = Path(__file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / "results" / "tables").is_dir() and (cand / "analysis").is_dir():
            return cand
    return here.parents[4]


ROOT = locate_root()
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures"

NICHE_GROUPS = [
    ("tumour\n(deep)", ["tumor"]),
    ("tumour\nsurface", ["cancer_surface"]),
    ("adjacent /\nparacancerous",
     ["adjacent_normal_tissue", "paracancerous", "matched_normal"]),
    ("lesion-side\nmucosa", ["lesion_side_mucosa"]),
    ("contralateral /\ndistant mucosa", ["contralateral_normal_mucosa"]),
    ("saliva", ["saliva", "preop_saliva", "postop_saliva"]),
    ("tongue\ncoat", ["tongue_coat"]),
    ("plaque", ["subgingival_plaque"]),
]
# platform came from cohort_params.yaml; country came from the ENA read_run
# `country` field. Three labels here were previously wrong (412445 said USA but its
# ENA records say China; 666746 said HiSeq but the runs are MiSeq; 597251 said MiSeq
# but the runs are HiSeq 2500). "Korea" is carried from the cohort notes for
# 1415268, whose ENA records leave country unset.
COHORT_LABEL = {
    "PRJNA866676": "PRJNA866676\nIon Torrent, China",
    "PRJNA1415359": "PRJNA1415359\nMiSeq, Korea",
    "PRJNA412445": "PRJNA412445\nMiSeq, China",
    "PRJNA666746": "PRJNA666746\nMiSeq, India",
    "PRJNA1415268": "PRJNA1415268\nMiSeq, Korea",
    "PRJNA597251": "PRJNA597251\nHiSeq, China",
    # the hospital cohort, drawn from its own manifest rather than the public one. The vendor
    # project number is deliberately not shown: it identifies the sequencing job, not the cohort.
    "HOSPITAL": "Hospital\nNovaSeq, China",
}

# (label, cohorts that could test it, cohorts that confirmed it). These are the counts the figure
# draws and the caption repeats, so they must match the source tables, not an earlier draft.
# "Parvimonas enriched" is testable in the four pooled public cohorts plus the hospital cohort
# (s5_meta_analysis.csv n_cohorts = 4; s5_score_per_cohort.csv adds the hospital arm), and every
# one of them is positive, hence 5 of 5. The caption states which five.
EVIDENCE = [
    ("Tumour differs from normal", 5, 5),
    ("Adjacent tissue is closest", 2, 1),
    ("Gradient towards saliva", 1, 1),
    ("$\\it{Parvimonas}$ enriched (genus)", 5, 5),
    ("$\\it{P.\\ micra}$ resolved (species)", 1, 1),
    ("Holds on RNA-seq", 1, 0),
    ("Hidden taxa (genus-unassigned)", 2, 0),
]

STAGES = [
    ("1  Identifiability", "Can source fractions\nbe estimated at all?", "simulation, 5 depths"),
    ("2  Source attribution", "Which niche is the\ntumour closest to?", "866676 (6 niches)"),
    ("3  Replication", "Does the ranking hold\nin another cohort?", "1415359 (n=31 pairs)"),
    ("4  Species resolution", "Is the signal\none species?", "866676 + 412445"),
    ("5  Platform check", "Does it hold on\ngenome-derived calls?", "TCGA-HNSC / TCMA"),
]


SUBSET_DIR = ROOT / "analysis" / "modules" / "02-dada2" / "input"


def analysed_manifest() -> pd.DataFrame:
    """Manifest restricted to the runs that actually entered the DADA2 analysis.

    The manifest lists every run in the archive, which is NOT the same as the runs
    that were analysed. PRJNA1415359 is the only cohort where the two differ: 241
    runs were fetched, but DADA2 processed 125 (32 tumour + 32 adjacent + 31 pre-op
    saliva + 30 post-op saliva, from 32 patients), per the subset list under
    analysis/modules/02-dada2/input/. The 116 excluded runs are 80 saliva-only
    donors who contributed no tissue and 36 saliva runs from 18 patients with no
    tumour sample; none of them enters any downstream analysis.

    Counting archive runs made panel A overstate this cohort by a factor of two
    (row total 241 vs 125, grand total 898 vs 782) and made panel F halve its
    apparent pairing depth (130 donors with a median of 1 niche, instead of 32
    patients with a median of 4). For the other five cohorts the archive runs and
    the analysed samples are the same set, so no subset file exists for them.
    """
    man = pd.read_csv(TAB / "oscc_niche_sample_manifest.csv", encoding="utf-8-sig")
    keep = pd.Series(True, index=man.index)
    for f in sorted(SUBSET_DIR.glob("subset_*.txt")):
        study = f.stem.replace("subset_", "")
        runs = set(pd.read_csv(f, sep="\t")["run"].astype(str))
        keep &= ~(man["study"].eq(study) & ~man["run"].astype(str).isin(runs))
    return man[keep].copy()


def panel_a(ax) -> None:
    man = analysed_manifest()
    ct = pd.crosstab(man.study, man.niche)

    def breadth(c: str) -> int:
        return sum(1 for _, srcs in NICHE_GROUPS
                   if any(s in ct.columns and c in ct.index and ct.loc[c, s] > 0
                          for s in srcs))

    # the in-house cohort is not in the manifest; it is appended by hand below, so exclude it here
    cohorts = sorted([c for c in COHORT_LABEL if c != "HOSPITAL"],
                     key=lambda c: (-breadth(c), c))

    M = np.zeros((len(cohorts), len(NICHE_GROUPS)), dtype=float)
    for j, (_, srcs) in enumerate(NICHE_GROUPS):
        for s in srcs:
            if s in ct.columns:
                M[:, j] += ct.reindex(cohorts)[s].fillna(0).to_numpy()

    # The in-house paired cohort is not in the public manifest, so its row is appended by hand.
    # It samples tumour and one distant normal tissue: no adjacent mucosa, no saliva.
    cohorts = cohorts + ["HOSPITAL"]
    M = np.vstack([M, np.zeros((1, len(NICHE_GROUPS)))])
    IN_HOUSE = {"tumor": 20, "adjacent /\nparacancerous": 20}
    # distant normal mucosa is closest to the "contralateral mucosa" column by definition of the
    # niche groups, but it is NOT contralateral: it is a distant mucosal site. Put it in its own
    # column so the figure does not overstate what was sampled.
    M[-1, [g for g, _ in NICHE_GROUPS].index("tumour\n(deep)")] = 20
    M[-1, [g for g, _ in NICHE_GROUPS].index("contralateral /\ndistant mucosa")] = 20

    cmap = plt.get_cmap("Blues")
    vmax = M.max()
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if v > 0:
                ax.add_patch(plt.Rectangle((j, i), 1, 1,
                                           facecolor=cmap(0.20 + 0.65 * v / vmax),
                                           edgecolor="white", linewidth=1.1))
                ax.text(j + 0.5, i + 0.5, f"{int(v)}", ha="center", va="center",
                        fontsize=6.1, color="white" if v / vmax > 0.42 else "0.15")
            else:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="#f2f2f2",
                                           edgecolor="white", linewidth=1.1))

    ax.set_xlim(0, M.shape[1])
    ax.set_ylim(M.shape[0], 0)
    ax.set_xticks(np.arange(M.shape[1]) + 0.5)
    ax.set_xticklabels([g for g, _ in NICHE_GROUPS], fontsize=6.2)
    ax.set_yticks(np.arange(M.shape[0]) + 0.5)
    # seven rows in a fixed band: the two-line labels need slightly tighter leading than
    # 1.05, but the row pitch itself is what fixes the overlap, so it is set by the grid ratio
    ax.set_yticklabels([COHORT_LABEL[c] for c in cohorts], fontsize=4.5,
                   linespacing=0.98)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    # bracket the six-niche cohort on the RIGHT (the left is occupied by row labels)
    top = cohorts.index("PRJNA866676")
    xw = M.shape[1]
    ax.annotate("", xy=(xw + 0.14, top + 0.07), xytext=(xw + 0.14, top + 0.93),
                arrowprops=dict(arrowstyle="-", color="#b2182b", lw=2.2))
    ax.text(xw + 0.30, top + 0.5, "6 niches,\nsame patients", rotation=0,
            ha="left", va="center", fontsize=5.3, color="#b2182b",
            fontweight="bold", linespacing=1.3)
    ax.set_xlim(0, xw + 1.55)

    ax.set_title("Cohort × niche coverage — only one cohort samples all six niches in "
                 "the same patients", fontsize=7.2, pad=7)
    panel_letter(ax, "A", case="upper")


BOX_PAD = 0.005        # FancyBboxPatch pad: the drawn border sits this far outside the patch


def _arrow(ax, right_edge, gap, y) -> None:
    """Compact arrow inside a gutter: the boxes carry a FancyBbox pad, so the arrow has to sit
    strictly between the two drawn borders rather than at the patch coordinates."""
    lo = right_edge + BOX_PAD + 0.0015
    hi = right_edge + gap - BOX_PAD - 0.0015
    ax.add_patch(FancyArrowPatch((lo, y), (hi, y),
                                 arrowstyle="-|>,head_width=0.55,head_length=0.9",
                                 mutation_scale=5, color="0.42",
                                 linewidth=0.8, shrinkA=0, shrinkB=0))


def panel_b(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    n = len(STAGES)
    # Five boxes plus an end-cap. The gutter has to hold an arrow, so it is set to roughly
    # 3 mm at final size; the few millimetres that costs each box is paid back by the shorter
    # labels, which no longer need the width.
    gap, x0 = 0.0180, 0.004
    cap_w = 0.138
    w = (1.000 - cap_w - (n - 1) * gap - gap - x0) / n
    # One box height and one rule for where the three text rows sit, so every box has the same
    # whitespace above, between and below its rows. The rows are placed from the box extent
    # rather than by eye, which is what left the gap under the two-line question too large.
    BOX_Y0, BOX_H = -0.28, 1.20
    BOX_TOP = BOX_Y0 + BOX_H
    ROW_PAD = 0.10                    # clear space inside the border, top and bottom
    ROW_H = 0.035                     # half-height of a single line of text
    Y_TITLE = BOX_TOP - ROW_PAD - ROW_H
    Y_COHORT = BOX_Y0 + ROW_PAD + ROW_H
    Y_QUESTION = (Y_TITLE + Y_COHORT) / 2.0
    for i, (name, question, cohorts) in enumerate(STAGES):
        x = x0 + i * (w + gap)
        cxx = x + w / 2
        # The boxes extend BELOW the axes (clip_on=False) because this axes is axis("off"):
        # without the overhang the same inter-row gap would read as dead space between B and C.
        ax.add_patch(FancyBboxPatch((x, BOX_Y0), w, BOX_H,
                                    boxstyle="round,pad=0.005,rounding_size=0.010",
                                    facecolor="#eef4fa", edgecolor="#8fb4d4", linewidth=0.8,
                                    clip_on=False))
        ax.text(cxx, Y_TITLE, name, ha="center", va="center",
                fontsize=5.6, fontweight="bold", color="#1f4e79")
        ax.text(cxx, Y_QUESTION, question, ha="center", va="center",
                fontsize=4.9, color="0.22", linespacing=1.40)
        ax.text(cxx, Y_COHORT, cohorts, ha="center", va="center",
                fontsize=4.6, color="0.42", style="italic")
        # arrow into the NEXT box; the last stage's arrow into the end-cap is drawn after the
        # loop, because the cap is a separate patch with its own x
        if i < len(STAGES) - 1:
            _arrow(ax, x + w, gap, Y_QUESTION)
    # end-cap: the in-house paired cohort closes the flow, drawn in the accent colour because it
    # is the only part of the design that is new data rather than re-analysis
    xc = x0 + n * (w + gap)
    # the arrow from stage 5 into the cap: same gutter, same rule as inside the loop
    _arrow(ax, x0 + (n - 1) * (w + gap) + w, gap, Y_QUESTION)
    ax.add_patch(FancyBboxPatch((xc, BOX_Y0), cap_w, BOX_H,
                                boxstyle="round,pad=0.005,rounding_size=0.010",
                                facecolor="#fbeaea", edgecolor="#b2182b", linewidth=0.9,
                                clip_on=False))
    # the cap is narrow, so the heading gets its own line and a little breathing room above
    # the question block instead of sitting on top of it
    ax.text(xc + cap_w / 2, Y_TITLE, "6  Our cohort", ha="center", va="center",
            fontsize=5.2, fontweight="bold", color="#b2182b")
    ax.text(xc + cap_w / 2, Y_QUESTION, "Does it repeat in\npaired tissue we\nsequenced?",
            ha="center", va="center", fontsize=4.3, color="#7a1c1c", linespacing=1.26)
    # the cap's footer is two lines where the stage boxes carry one, so it needs the extra
    # line height subtracted from the shared baseline rule
    ax.text(xc + cap_w / 2, Y_COHORT + ROW_H * 0.9, "hospital cohort\n20 paired samples",
            ha="center", va="center", fontsize=4.3, color="#7a1c1c", style="italic",
            linespacing=1.18)
    ax.set_title("Analysis flow — five questions, then the one cohort that could answer "
                 "a sixth", fontsize=7.2, pad=6)
    panel_letter(ax, "B", case="upper")


def panel_c(ax) -> None:
    labels = [e[0] for e in EVIDENCE]
    testable = np.array([e[1] for e in EVIDENCE], dtype=float)
    confirmed = np.array([e[2] for e in EVIDENCE], dtype=float)
    y = np.arange(len(labels))[::-1]

    ax.barh(y, testable, height=0.54, color="#cfe0ee", edgecolor="#8fb4d4",
            linewidth=0.7, label="cohorts that could test it", zorder=2)
    ax.barh(y, confirmed, height=0.27, color="#2c7fb8", edgecolor="white",
            linewidth=0.5, label="cohorts that confirmed it", zorder=3)

    for yi, t, c in zip(y, testable, confirmed):
        # When testable and confirmed coincide, the two counts would land on the same point.
        # Print them as "c of t" once, outside the bar, instead of stacking two digits.
        if c and c == t:
            ax.text(t + 0.12, yi, f"{int(c)}/{int(t)}", va="center", ha="left",
                    fontsize=5.4, color="#2c7fb8", fontweight="bold")
        else:
            ax.text(t + 0.12, yi, f"{int(t)}", va="center", ha="left",
                    fontsize=5.4, color="0.45")
            if c:
                # below the bar's centre line and just past its own end, never on the bar
                ax.text(c + 0.12, yi - 0.26, f"{int(c)}", va="center", ha="left",
                        fontsize=5.4, color="#2c7fb8", fontweight="bold")
        if t == 1:
            ax.text(2.30, yi, "single cohort — cannot be replicated",
                    va="center", ha="left", fontsize=4.9, color="#b2182b",
                    style="italic")

    ax.set_yticks(y)
    # two-line labels, centred so the block has one clean left edge instead of a ragged
    # one; the asterisks mark the species name for italics below
    # rows fill the band symmetrically, so the space above the first bar equals the space
    # below the last one and nothing is wasted; labels are now short enough for one line
    ax.set_ylim(-0.75, len(labels) - 0.25)
    # Species names must be italic. set_fontstyle() on a tick label is dropped when the tick
    # formatter redraws, so the one emphasised label already carries mathtext.
    ax.set_yticklabels([lab.replace("*", "") for lab in labels], fontsize=5.9)
    ax.set_xlabel("number of independent cohorts", fontsize=6.3)
    ax.set_xlim(0, 6.1)
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.legend(loc="lower right", frameon=False, fontsize=5.8,
              handlelength=1.3, labelspacing=0.3, borderaxespad=0.25)
    ax.set_title("Evidence map — replication depth differs by conclusion, and is stated "
                 "per conclusion", fontsize=7.2, pad=7)
    panel_letter(ax, "C", case="upper")


TECH = {
    # cohort: (platform, primer trimming needed?, truncLen, median non-chimera reads/sample)
    # `trunc` is carried for reference only — panel E renders cohort, platform, primer status and
    # reads, not truncLen. PRJNA412445 was 236/236 in a pre-run estimate; the run used 245/240
    # (run_status.json), corrected here so a future edit cannot inherit the stale value.
    "PRJNA866676": ("Ion Torrent", "none (SE)", "350", 56828),
    "PRJNA412445": ("MiSeq",          "none",    "245/240",   11686),
    "PRJNA1415359": ("MiSeq",         "trim 17/21",         "280/225",   83207),
    "PRJNA1415268": ("MiSeq",         "trim 17/21",         "230/225",   18184),
    "PRJNA597251": ("HiSeq 2500",     "none",    "287/273",   11197),
    "PRJNA666746": ("MiSeq",          "none",    "240/240*",  26164),
}


def panel_tech(ax) -> None:
    """E - per-cohort technical settings, all measured rather than defaulted."""
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.05, len(TECH) + 0.10)
    ax.axis("off")

    XC, XP, XR, XB, BW = 0.000, 0.115, 0.330, 0.520, 0.400
    hdr = len(TECH) - 0.34
    for x, lab in [(XC, "cohort"), (XP, "platform"), (XR, "5' primers"),
                   (XB, "median reads / sample")]:
        ax.text(x, hdr, lab, fontsize=4.8, color="0.35", va="bottom", ha="left")

    lo, hi = np.log10(8000), np.log10(100000)
    for i, (cohort, (plat, prim, trunc, reads)) in enumerate(TECH.items()):
        y = len(TECH) - 1 - i
        ax.text(XC, y, cohort.replace("PRJNA", ""), fontsize=4.9, va="center")
        ax.text(XP, y, plat, fontsize=4.4, va="center", color="0.25")
        trimmed = prim.startswith("trim")
        ax.text(XR, y, prim, fontsize=4.3, va="center",
                color="#b2182b" if trimmed else "0.45",
                fontweight="bold" if trimmed else "normal")
        frac = (np.log10(reads) - lo) / (hi - lo)
        w = BW * max(frac, 0.03)
        ax.barh(y, w, left=XB, height=0.5,
                color="#b2182b" if trimmed else "#7a9cc6",
                edgecolor="0.3", linewidth=0.4, zorder=3)
        ax.text(XB + w + 0.010, y, f"{reads/1000:.0f}k", fontsize=4.3, va="center",
                color="0.3")
    ax.text(0.0, -0.88, "all amplicons V3-V4;  red = 5' primer trimmed;  "
            "666746 downsampled", fontsize=4.2, color="0.35", style="italic")
    ax.set_title("Technical settings, measured per cohort", fontsize=7.2,
                 multialignment="center", pad=5)
    panel_letter(ax, "D", case="upper")


def panel_pairing(ax) -> None:
    """F - how many niches each patient contributed, per cohort."""
    d = analysed_manifest()
    g = d.groupby(["study", "patient_id"])["niche"].nunique().reset_index()
    tab = (g.groupby(["study", "niche"]).size().unstack(fill_value=0)
           .reindex(["PRJNA866676", "PRJNA412445", "PRJNA1415359",
                     "PRJNA1415268", "PRJNA597251", "PRJNA666746"]))
    order = [c for c in tab.columns]
    cmap = plt.get_cmap("YlGnBu")
    ys = np.arange(len(tab))[::-1]
    left = np.zeros(len(tab))
    for j, k in enumerate(order):
        v = tab[k].values
        col = cmap(0.22 + 0.13 * (j + 1))
        ax.barh(ys, v, left=left, height=0.55, color=col, edgecolor="white",
                linewidth=0.4, zorder=3, label=f"{k}")
        for y, l, vv in zip(ys, left, v):
            # tested `vv`, the segment's OWN value. It used to test `l`, the running
            # left edge, so once a row had been pushed right by earlier segments every
            # tiny segment inherited a label ("2"/"3" printed on slivers in 412445).
            if vv >= 12:
                ax.text(l + vv / 2, y, str(int(vv)), ha="center", va="center",
                        fontsize=4.2, color="white" if col[0] < 0.55 else "0.15")
        left += v
    ax.set_yticks(ys)
    ax.set_yticklabels([s.replace("PRJNA", "") for s in tab.index], fontsize=4.9)
    # the legend hangs below the axis and carries its own title, so an x-label directly
    # under the tick labels collided with it. The axis name moves into the panel title line.
    ax.set_xlabel("")
    # Legend sits in the empty upper-right of the axes. Below the axes it pushed the
    # canvas down and left a large blank band under panel E, which has no such overhang.
    ax.set_ylim(-0.6, len(tab) - 0.15)
    # The legend sat inside the upper right, which worked while the x axis reached 125
    # (it was counting archive runs). Now that the axis ends near the largest cohort,
    # the longest bar runs under it, so the key moves below the axis as one row.
    ax.legend(frameon=False, fontsize=4.1, loc="upper center",
              bbox_to_anchor=(0.5, -0.30), handlelength=0.8, labelspacing=0.12,
              borderaxespad=0.0, ncol=6, columnspacing=0.55,
              title="patients by niches per patient", title_fontsize=4.3)
    ax.set_title("Pairing depth differs by cohort", fontsize=7.2,
                 multialignment="center", pad=5)
    panel_letter(ax, "E", case="upper")


def panel_overlap(ax) -> None:
    """Which of our cohorts also appear in the concurrent 13-cohort meta-analysis.

    This panel exists because a reviewer's first question is what is new. The honest answer is a
    table: we share 2 of 6 public cohorts with Na et al. 2026, none of the four cohorts that carry
    our source-attribution or paired-replication results, and they have no in-house samples.
    """
    O = pd.read_csv(TAB / "fig1_cohort_overlap_with_na2026.csv", encoding="utf-8-sig")
    O = O.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(O))
    colours = ["#7b9ab5" if s == "public" else "#b2182b" for s in O["Source"]]
    ax.barh(y, O["Samples (n)"], color=colours, height=0.62, zorder=3)
    for i, r in O.iterrows():
        ax.text(r["Samples (n)"] + 6, i, f"{int(r['Samples (n)'])}",
                va="center", fontsize=5.0, color="0.2")
    ax.set_yticks(y)
    ax.set_yticklabels([c for c in O["Cohort"]], fontsize=5.2)
    # niches per patient inside the bar: a second label line crowded the y tick labels
    for i, (k, v) in enumerate(zip(O["Niches per patient (max)"], O["Samples (n)"])):
        if v >= 90:
            ax.text(v - 6, i, f"{int(k)}", ha="right", va="center", fontsize=4.7,
                    color="white", zorder=5)
    ax.set_xlabel("samples analysed", fontsize=6.0)
    ax.set_xlim(0, 372)

    ax.tick_params(axis="x", labelsize=5.6)
    # mark which rows are shared with Na et al.
    for i, r in O.iterrows():
        if r["Also in Na 2026"] == "yes":
            ax.text(352, i, "also in Na 2026", ha="right", va="center",
                    fontsize=4.6, color="#b2182b", style="italic", zorder=5)
    ax.set_title("Only 2 of 6 public cohorts overlap a concurrent meta-analysis,\n"
                 "and every cohort here pairs tumour with the patient's own tissue "
                 "(bar label = tissue types/patient)",
                 fontsize=6.2, multialignment="center", pad=5)
    ax.text(0.985, 0.05,
            f"we add {int(O.loc[O['Also in Na 2026']=='no','Samples (n)'].sum())} samples "
            f"not in that pool, incl. 20 paired in-house",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=4.8, color="0.25")
    panel_letter(ax, "F", case="upper")


def main() -> None:
    apply_figure_style(sizes=(8.0, 7.1, 6.2))
    fig = plt.figure(figsize=(7.09, 9.55))          # 180 mm wide
    gs = fig.add_gridspec(5, 2, width_ratios=[1.42, 0.92], left=0.134, right=0.985,
                          height_ratios=[84, 34, 46, 50, 48], hspace=0.86, wspace=0.34)

    panel_a(fig.add_subplot(gs[0, :]))
    panel_b(fig.add_subplot(gs[1, :]))
    panel_c(fig.add_subplot(gs[2, :]))
    panel_tech(fig.add_subplot(gs[3, 0]))
    panel_pairing(fig.add_subplot(gs[3, 1]))
    panel_overlap(fig.add_subplot(gs[4, :]))

    out = FIG / "manuscript-v6" / "main" / "Fig1-design-and-cohort-overlap.png"
    fig.savefig(out, dpi=330, bbox_inches="tight", facecolor="white")
    print("wrote", out)
    # panel_crops is off: its PNGs are diagnostic only and would linger in .cache


if __name__ == "__main__":
    main()
