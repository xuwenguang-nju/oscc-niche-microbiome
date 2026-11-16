#!/usr/bin/env python3
"""F4 — does the source ranking replicate, and is the signal one species?

Ten panels on a rigid 3-column x 4-row grid (see the layout block below), letters A-J
strictly left to right within every row:

A  source ranking in both paired cohorts (Bray-Curtis, mean rank to own tumour)
B  P. micra across the six OSCC niches
C  the ten replicating genera x six niches, as a heat map (replaces supplementary S2)
D  reads in ASVs that never reach genus level — the blind spot, quantified
E  median source distance per study and source, as a matrix
F  family-level test of those genus-less reads (the one hit is LOWER in tumour)
G  how much of Parvimonas can be resolved to a species at all
H  per-patient distance to the closest source, both cohorts
I  the surface proxy: lesion-side vs contralateral mucosa (the only available stand-in)
J  P. micra paired log2 fold change, tumour vs comparator, both cohorts

Reads:  results/tables/s8a_rank_tests.csv, s8b_species_summary.csv,
        s8b_species_by_niche.csv, s8a_distances_long.csv, s8a_surface_proxy.csv,
        s8a_source_distance_summary.csv, s8c_genusless_read_share.csv,
        s8c_hidden_family_tests.csv, s8b_species_paired_tests.csv
Writes: figures/main/F4-replication-and-species.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "14-figure-layout" / "scripts"))
from figlib import panel_letter_fixed, project_root  # noqa: E402

ROOT = project_root(Path(__file__).resolve().parent)
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures"
# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
exec(compile(open(_FIGSTYLE / "runtime.py", encoding="utf-8").read(),
              "runtime.py", "exec"))

apply_figure_style(sizes=(7.0, 6.2, 5.4))

rank = pd.read_csv(TAB / "s8a_rank_tests.csv")
sp_sum = pd.read_csv(TAB / "s8b_species_summary.csv")
sp_niche = pd.read_csv(TAB / "s8b_species_by_niche.csv")
dist = pd.read_csv(TAB / "s8a_distances_long.csv")
proxy = pd.read_csv(TAB / "s8a_surface_proxy.csv")
sdsum = pd.read_csv(TAB / "s8a_source_distance_summary.csv")
gl = pd.read_csv(TAB / "s8c_genusless_read_share.csv")
hid = pd.read_csv(TAB / "s8c_hidden_family_tests.csv")
sp_pair = pd.read_csv(TAB / "s8b_species_paired_tests.csv")

# The three pooled panels (K, L, M) come from the cross-cohort meta-analysis. They are drawn
# here rather than in a separate figure so the replication story and the pooled summary sit
# in one place; the remaining meta panels stay in the supplement.
meta = pd.read_csv(TAB / "s5_meta_analysis.csv")
eff = pd.read_csv(TAB / "s5_effect_by_cohort_full.csv")     # carries per-cohort se and n
nov = pd.read_csv(TAB / "s11_novelty_assessment.csv")
K_COH = 3
REP = meta[(meta.n_cohorts >= K_COH) & (meta.pooled_q_bh < 0.05) &
           meta.direction_consistent].copy()
M_STUDIES = ["PRJNA1415359", "PRJNA1415268", "PRJNA597251", "PRJNA666746"]
M_SLAB = {"PRJNA1415359": "1415359", "PRJNA1415268": "1415268",
        "PRJNA597251": "597251", "PRJNA666746": "666746"}
M_SCOL = {"PRJNA1415359": "#1f4e79", "PRJNA1415268": "#3182bd",
        "PRJNA597251": "#a6611a", "PRJNA666746": "#7b3294"}
M_LARGE_LOG = 0.92 / np.log(2)   # |d_z| = 0.8 converted, now on the log2 axis

BLUE, RED, ORANGE, GREEN = "#2166ac", "#b2182b", "#a6611a", "#1b7837"
SCO = {"PRJNA1415359": BLUE, "PRJNA412445": RED, "PRJNA866676": "#7b3294"}
# City names were hard-coded here earlier and two were wrong (412445 is not from
# Pennsylvania; 866676 is from Jinan, not Qingdao). The accession cannot be wrong
# and the city carries no scientific weight, so label by accession instead.
SLAB = {"PRJNA1415359": "1415359\n(31 pairs)", "PRJNA412445": "412445\n(8 pairs)",
        "PRJNA866676": "866676\n(313 runs)"}
SRC = {"adjacent_normal_tissue": "adjacent normal",
       "contralateral_normal_mucosa": "contralateral",
       "lesion_side_mucosa": "lesion-side mucosa",
       "preop_saliva": "pre-op saliva", "saliva": "saliva",
       "subgingival_plaque": "subgingival plaque", "postop_saliva": "post-op saliva"}

# ------------------------------------------------------------------ rigid layout
# A 3-column x 4-row grid with fixed seat spans, so panel FRAMES and panel LETTERS
# land on the same vertical and horizontal guide lines:
#
#     A | B  B          A 1 column, B 2 columns
#     C  C | D          C 2 columns, D 1 column
#     E | F | G         one column each
#     H | I | J         one column each
#
# Three vertical guide lines (column start + gutter) carry every panel's left edge
# and every letter: x0 (A C E H), x1 (B F I), x2 (D G J). Right edges fall on x1
# (A E H), x2 (C) and x3 (B D G J), and all four rows share one band height.
#
# The two things that broke this before are both gone:
#   * _indent(ax, frac) pushed each panel rightwards by a fraction OF ITS OWN cell,
#     so every panel got a different left edge and a different width;
#   * fig.colorbar(..., ax=ax) silently steals width from the parent axes, so the
#     two panels with colour bars (C, E) were narrower than their neighbours.
# Label room is now one ABSOLUTE gutter per column, sized to the widest y label that
# column must carry; colour bars get axes of their own (see hbar_below).
FW, FH = 7.09, 9.55                     # nominal inches; delivered width set later
XM_L, XM_R = 0.29, 0.06                 # XM_L + GUT (0.66) = 0.95 in, the shared
                                        # left edge of the first column, so panel
                                        # letters line up with every other figure
GUT = 0.66                              # y-label room, identical for every column
GAPS = (0.60, 0.92, 0.62, 0.74)         # between rows 0/1 .. 3/4; rows 1 and 2 own
                                        # colour bars, so those gaps carry bar + label
#   row 1 carries C's colour bar (0.36 in) and row 2 carries E's (0.50 in); each gap
#   must also hold the next row's panel letter + title (~0.28 in)
TOP, BOT = FH - 0.36, 0.52              # top of row 0, bottom of row 3
COLW = (FW - XM_L - XM_R) / 3.0
X = [XM_L + i * COLW for i in range(4)]         # column boundaries x0..x3
ONE = COLW - GUT                        # width of a 1-column panel
TWO = 2 * COLW - GUT                    # width of a 2-column panel (B, C)
ROW_H = (TOP - BOT - sum(GAPS)) / 5.0
ROW_Y, _y = [], TOP
for _g in (*GAPS, 0.0):
    ROW_Y.append(_y - ROW_H)
    _y -= ROW_H + _g

SEAT = {                                # (row, col) -> column span, end-exclusive
    (0, 0): (0, 1), (0, 1): (1, 3),
    (1, 0): (0, 2), (1, 1): (2, 3),
    (2, 0): (0, 1), (2, 1): (1, 2), (2, 2): (2, 3),
    (3, 0): (0, 1), (3, 1): (1, 2), (3, 2): (2, 3),
    (4, 0): (0, 1), (4, 1): (1, 2), (4, 2): (2, 3),
}

fig = plt.figure(figsize=(FW, FH))
axes = np.empty((5, 3), dtype=object)
for (_r, _c), (_c0, _c1) in SEAT.items():
    # every column gets the same label gutter, so a 1-column panel has the same
    # width wherever it sits and a 2-column panel is exactly ONE + COLW. With
    # per-column gutters the frames still landed on three guides, but A came out
    # 1.48 in against D's 1.78 in -- the asymmetry being fixed here.
    axes[_r, _c] = fig.add_axes([(X[_c0] + GUT) / FW, ROW_Y[_r] / FH,
                                 (TWO if _c1 - _c0 == 2 else ONE) / FW, ROW_H / FH])


def hbar_below(ax, im, label, width_frac=0.55, drop=0.20, h=0.062):
    """Horizontal colour bar in its OWN axes, hung `drop` inches under `ax`.

    fig.colorbar(im, ax=ax) shrinks the parent axes to make room, which would put
    that panel's width out of line with its neighbours -- exactly what this layout
    exists to prevent.
    """
    box = ax.get_position()
    cax = fig.add_axes([box.x0, box.y0 - (drop + h) / FH,
                        box.width * width_frac, h / FH])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.4)
    cb.set_label(label, fontsize=4.2, labelpad=1.2)
    cb.ax.tick_params(labelsize=4.0, length=1.8, pad=1.0)
    return cb


def lbl(ax, letter, text):
    ax.set_title(text, fontsize=6.2, multialignment="center", pad=5)
    panel_letter_fixed(ax, letter, case="upper")


# ------------------------------------------------------------------ A rank replication
ax = axes[0, 0]
r = rank[rank.metric == "bray"]
sources = r.groupby("source")["mean_rank"].mean().sort_values().index
ys = np.arange(len(sources))
for i, st in enumerate(["PRJNA1415359", "PRJNA412445"]):
    d = r[r.study == st].set_index("source").reindex(sources)
    ax.scatter(d["mean_rank"], ys - 0.16 + i * 0.32, s=22, color=SCO[st],
               edgecolor="0.25", linewidth=0.35, zorder=3, label=SLAB[st].split("\n")[0])
# The two cohorts have different numbers of candidate sources (2 in PRJNA1415359, 4 in
# PRJNA412445), so no single chance rank describes both. Each study gets its own reference line,
# drawn as a short vertical tick on its own row offset rather than a full-height rule.
for i, st in enumerate(["PRJNA1415359", "PRJNA412445"]):
    n_src = int(r[r.study == st]["source"].nunique())
    chance = n_src / 2 + 0.5
    ax.plot([chance, chance], [ys.min() - 0.42 + i * 0.05, ys.max() + 0.42 - (1 - i) * 0.05],
            color=SCO[st], lw=0.9, ls=(0, (4, 1.6)), zorder=2, alpha=0.9)
# One line, kept short enough to stay inside the panel width: the earlier wording ran past the
# right edge of A and, when wrapped, collided with panel C's title.
ax.text(0.0, -0.245, "dashed: chance rank 1.5 / 2.5 per cohort",
        transform=ax.transAxes, fontsize=4.2, color="0.35", va="top", ha="left",
        clip_on=False)
ax.set_yticks(ys)
ax.set_yticklabels([SRC[s] for s in sources], fontsize=4.9)
ax.invert_yaxis()
ax.set_xlabel("mean rank to own tumour", fontsize=5.8)
ax.legend(frameon=False, fontsize=4.5, loc="upper center", ncol=2,
          bbox_to_anchor=(0.5, 1.005), handletextpad=0.2, columnspacing=1.0,
          labelspacing=0.2, borderaxespad=0.0)
lbl(ax, "A", "Ranking replicates")

# ------------------------------------------------------------------ B species resolution
ax = axes[2, 2]
d = sp_sum.copy()
xs = np.arange(len(d))
ax.bar(xs, d["median_resolved_frac_of_genus"] * 100, color=[SCO[s] for s in d.study],
       edgecolor="0.3", lw=0.4, width=0.55, zorder=3)
for x, v in zip(xs, d["median_resolved_frac_of_genus"] * 100):
    ax.text(x, v + 1.4, f"{v:.0f}%", ha="center", fontsize=5.2, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels([SLAB[s] for s in d.study], fontsize=4.9)
ax.set_ylim(0, 100)
ax.set_ylabel("% of Parvimonas reads\nresolved to species", fontsize=5.8)
ax.text(0.97, 0.97, "remainder stays\nunresolved", transform=ax.transAxes,
        ha="right", va="top", fontsize=4.5, color="0.35", style="italic")
lbl(ax, "G", "Only a fifth resolves")

# ------------------------------------------------------------------ C P. micra by niche
ax = axes[0, 1]
d = sp_niche[sp_niche.study == "PRJNA866676"].copy()
ORDER = ["tumor", "paracancerous", "cancer_surface", "saliva", "tongue_coat",
         "contralateral_normal_mucosa"]
NICE = {"tumor": "tumour (deep)", "paracancerous": "paracancerous",
        "cancer_surface": "cancer surface", "saliva": "saliva",
        "tongue_coat": "tongue coat",
        "contralateral_normal_mucosa": "contralateral"}
d = d.set_index("niche").reindex([n for n in ORDER if n in set(d.niche)])
xs = np.arange(len(d))
ax.bar(xs, d["Parvimonas micra_pct_of_bacteria"], color="#1b7837",
       edgecolor="0.3", lw=0.4, width=0.62, zorder=3)
for x, v in zip(xs, d["Parvimonas micra_pct_of_bacteria"]):
    # label only the bars tall enough to carry a number; the two tiny niches stay
    # unlabelled rather than printing "-0.003"/"0.003" on top of each other
    if v >= 0.01:
        ax.text(x, v + 0.006, f"{v:.3f}", ha="center", fontsize=4.6, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels([NICE[i] for i in d.index], rotation=38, ha="right", fontsize=4.7)
ax.set_ylabel("P. micra (% of bacterial reads)", fontsize=5.8)
ax.set_ylim(0, float(d["Parvimonas micra_pct_of_bacteria"].max()) * 1.24)
lbl(ax, "B", "One species carries it")

# ------------------------------------------------------------------ D per-patient distance
ax = axes[3, 0]
adj = {"PRJNA1415359": "adjacent_normal_tissue", "PRJNA412445": "contralateral_normal_mucosa"}
groups, vals, cols = [], [], []
for st, src in adj.items():
    v = dist[(dist.study == st) & (dist.source == src)]["bray"].dropna().to_numpy()
    groups.append(v); vals.append(v); cols.append(SCO[st])
bp = ax.boxplot(vals, positions=np.arange(len(vals)), widths=0.42, patch_artist=True,
                medianprops=dict(color="0.15", lw=1.0),
                boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5),
                capprops=dict(lw=0.5), flierprops=dict(ms=1.5, mfc="0.6", mec="none"))
for patch, c in zip(bp["boxes"], cols):
    patch.set_facecolor(c); patch.set_alpha(0.45); patch.set_edgecolor("0.3")
for i, v in enumerate(vals):
    ax.scatter(np.full(len(v), i) + np.random.default_rng(0).normal(0, 0.045, len(v)), v,
               s=3.4, color=cols[i], alpha=0.75, linewidths=0, zorder=3)
LBL_H = [s.replace("PRJNA", "") + "\n"
         + ("adjacent" if s == "PRJNA1415359" else "contralateral")
         for s in adj]
ax.set_xticks(np.arange(len(vals)))
ax.set_xticklabels(LBL_H, fontsize=4.2)
ax.set_xlabel("closest available source per patient\n"
              + "  ·  ".join(f"n={len(g)}" for g in groups), fontsize=5.0)
ax.set_ylabel("Bray-Curtis to that source", fontsize=5.8)
lbl(ax, "H", "Per-patient distances")

# ------------------------------------------------------------------ E surface proxy
ax = axes[3, 1]
d = proxy.iloc[0]
ax.barh([0], [d["median_diff"]], color=ORANGE, edgecolor="0.3", lw=0.4, height=0.42, zorder=3)
ax.axvline(0, color="0.3", lw=0.8, zorder=4)
ax.set_yticks([0])
ax.set_yticklabels(["lesion-side −\ncontralateral"], fontsize=4.8)
ax.set_xlabel("median Aitchison difference\n(lesion-side − contralateral)", fontsize=5.4)
ax.text(0.5, 0.80, f"n={int(d['n_pairs'])} pairs,  p = {d['p']:.2f}",
        transform=ax.transAxes, ha="center", fontsize=4.9, color="0.25")
ax.set_ylim(-0.55, 0.62)
ax.set_xlim(min(-0.01, d["median_diff"] * 1.5), max(0.01, d["median_diff"] * 1.5))
ax.text(0.5, 0.06, "no stand-in for the surface gradient", transform=ax.transAxes,
        fontsize=4.6, ha="center", color=RED, style="italic")
lbl(ax, "I", "Surface proxy is null")

# ------------------------------------------------------------------ F genus-less read share
ax = axes[1, 1]
for i, st in enumerate(["PRJNA866676", "PRJNA412445"]):
    v = gl[gl.study == st]["genus_unassigned_read_share"].dropna().to_numpy() * 100
    v = np.sort(v)
    ax.plot(v, np.linspace(0, 100, len(v)), color=SCO.get(st, GREEN), lw=1.4,
            zorder=3, label=f"{st.replace('PRJNA','')}")
ax.set_xlabel("% of reads\nin genus-less ASVs", fontsize=5.4)
ax.set_ylabel("% of samples", fontsize=5.6)      # this is a cumulative distribution over samples
ax.legend(frameon=False, fontsize=4.7, loc="lower right",
          handlelength=1.1, labelspacing=0.2, borderaxespad=0.5)
ax.axvline(10, color="0.55", lw=0.8, ls=(0, (4, 1.6)), zorder=2)
lbl(ax, "D", "The blind spot is small")

# ------------------------------------------------------------------ G hidden family test
ax = axes[2, 1]
d = hid.copy().sort_values("q")
ys = np.arange(len(d))
cols = [RED if q < 0.05 else "#bdbdbd" for q in d["q"]]
ax.barh(ys, d["median_clr_diff"], color=cols, edgecolor="0.3", lw=0.4, height=0.6, zorder=3)
ax.axvline(0, color="0.3", lw=0.8, zorder=4)
ax.set_yticks(ys)
ax.set_yticklabels([f"{f}\n({s.replace('PRJNA','')})" for f, s in zip(d.family, d.study)],
                   fontsize=4.5)
ax.invert_yaxis()
ax.set_xlabel("median CLR diff.\n(tumour − normal)", fontsize=5.4)
for y, q, val in zip(ys, d["q"], d["median_clr_diff"]):
    # text sits ON the bar when the bar is long enough, so it must contrast with
    # the bar fill rather than default to grey-on-dark-red (invisible).
    # Long bars get white text inside, but the label must actually fit: at
    # 4.2 pt, "q=0.005" is ~0.34 data units wide, so require comfortably more.
    # Inside means "just within the bar end": for a positive bar that is right-
    # aligned at val-eps, for a negative bar left-aligned at val+eps.
    inside = abs(val) > 0.46
    dark_bar = q < 0.05          # only the crimson bar carries white text
    if inside and val >= 0:
        xpos, ha = val - 0.03, "right"
    elif inside:
        xpos, ha = val + 0.03, "left"
    elif val >= 0:
        xpos, ha = 0.02, "left"
    else:
        xpos, ha = -0.02, "right"
    col = "white" if (inside and dark_bar) else ("0.15" if inside else "0.25")
    ax.text(xpos, y, f"q={q:.3f}" if q >= 0.001 else f"q={q:.1e}", va="center",
            ha=ha, fontsize=4.2, color=col, zorder=5)
lbl(ax, "F", "Nothing hides there")

# ------------------------------------------------------------------ H source distance matrix
ax = axes[2, 0]
m = sdsum.pivot(index="source", columns="study", values="bray_median")
m = m.reindex([s for s in SRC if s in m.index])
im = ax.imshow(m.values, aspect="auto", cmap="YlOrRd", vmin=0.2, vmax=0.9)
ax.set_xticks(np.arange(m.shape[1]))
ax.set_xticklabels([s.replace("PRJNA", "") for s in m.columns], fontsize=4.9)
ax.set_yticks(np.arange(m.shape[0]))
ax.set_yticklabels([SRC[s] for s in m.index], fontsize=4.7)
for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        val = m.values[i, j]
        if np.isfinite(val):
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=4.4,
                    color="white" if val > 0.62 else "0.2")
ax.set_xlabel("cohort", fontsize=5.8)
# No colour bar here. E prints every cell value and H carries the same Bray-Curtis
# scale with a bar; a second one could only be placed by overlapping the row below at
# this row height. The unit goes in the axis label.
ax.set_xlabel("median Bray-Curtis to that source", fontsize=5.0)
lbl(ax, "E", "Same ordering")

# ------------------------------------------------------------------ I species effect size
ax = axes[3, 2]
d = sp_pair.copy()
ys = np.arange(len(d))
ax.barh(ys, d["median_log2fc"], color=[SCO[s] for s in d.study], edgecolor="0.3",
        lw=0.4, height=0.5, zorder=3)
ax.axvline(0, color="0.3", lw=0.8, zorder=4)
ax.set_yticks(ys)
ax.set_yticklabels([f"{s.replace('PRJNA', '')}\n{n} pairs"
                    for s, n in zip(d.study, d.n_pairs)], fontsize=4.6)
# The input is s8b_species_paired_tests.csv, whose taxon is *Parvimonas micra*: this panel is
# species-level, and the axis has to say so. It previously read "$\it{Parvimonas}$" (genus),
# which contradicted both the input and the caption.
ax.set_xlabel("median log$_2$ fold change, $\\it{P.\\ micra}$\n(tumour / comparator)",
              fontsize=5.4)
for y, (fc, p, f) in enumerate(zip(d["median_log2fc"], d["p_paired_wilcoxon"],
                                   d["frac_up_in_tumour"])):
    # p at two decimals printed 0.0044 as "p=0.00"; small values need three significant figures
    ptxt = f"p={p:.3f}" if p >= 0.001 else "p<0.001"
    ax.text(fc + (0.04 if fc >= 0 else -0.04), y,
            f"{fc:+.2f}  {ptxt}  {f*100:.0f}% up",
            va="center", ha="left" if fc >= 0 else "right", fontsize=4.4, color="0.25")
ax.set_xlim(0, max(d["median_log2fc"]) * 2.1)
lbl(ax, "J", "Both increase")

# ------------------------------------------------- J genus x niche heatmap (ex-S2)
# This replaces supplementary S2, which drew the same data as ten separate bar
# panels. A single heatmap carries the same 10 x 6 matrix in one cell.
ax = axes[1, 0]
gx = pd.read_csv(TAB / "s7c_s5_genera_by_niche.csv")
NICHE_COLS = [("tumor", "tumour"), ("paracancerous", "para-"),
              ("cancer_surface", "surface"), ("contralateral_normal_mucosa", "contra-"),
              ("saliva", "saliva"), ("tongue_coat", "tongue")]
M = np.log10(gx[[c for c, _ in NICHE_COLS]].values.astype(float) + 1e-4)
order_j = np.argsort(-M[:, 0])          # most abundant in tumour first
M = M[order_j]
names = gx["genus"].values[order_j]
dirs = gx["s5_direction"].values[order_j]
im = ax.imshow(M, cmap="Blues", aspect="auto", vmin=-4, vmax=np.log10(3.5))
ax.set_xticks(np.arange(len(NICHE_COLS)))
ax.set_xticklabels([lab for _, lab in NICHE_COLS], fontsize=4.2,
                   rotation=32, ha="right", rotation_mode="anchor")
ax.set_yticks(np.arange(len(names)))
ARR = {"up_in_tumour": "\u2191", "down_in_tumour": "\u2193"}
ax.set_yticklabels([ARR.get(d, "") + " " + g for g, d in zip(names, dirs)],
                   fontsize=4.2)
ax.set_xticks(np.arange(-0.5, len(NICHE_COLS), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(names), 1), minor=True)
ax.grid(which="minor", color="white", lw=0.5)
ax.tick_params(which="minor", length=0)
ax.tick_params(axis="both", length=2)
# a vertical bar here sat in the gutter between C and D and its ticks merged with
# D's y-axis; putting it under the heatmap keeps that gutter clear.
hbar_below(ax, im, "log$_{10}$ rel. abundance", width_frac=0.62, drop=0.26, h=0.050)
lbl(ax, "C", "Ten genera x six niches")


# ------------------------------------------------------------------ K selection funnel
# The volcano itself needs a wide cell: in 39 mm its leader-line labels collide. It is
# therefore panel A of Fig. S2, at full width. What stays here is the count that the
# volcano is there to make: of all testable genera, how few replicate.
ax = axes[4, 0]
n_all = int((meta.n_cohorts >= 2).sum())
n_two = int((meta.n_cohorts >= 2).sum())
n_rep = int(len(REP))
steps = [("testable\nin \u22652 cohorts", n_all),
         ("at FDR < 0.05,\nconsistent", n_rep),
         ("higher in\ntumour", int((REP.pooled_effect > 0).sum())),
         ("lower in\ntumour", int((REP.pooled_effect < 0).sum()))]
xs = np.arange(len(steps))
ax.bar(xs, [v for _, v in steps], color=[BLUE, "0.45", RED, "#8fb8dd"], width=0.62,
       edgecolor="0.25", linewidth=0.35)
for x, (lab, v) in zip(xs, steps):
    ax.text(x, v + n_all * 0.03, f"{v}", ha="center", va="bottom", fontsize=4.6)
ax.set_xticks(xs)
ax.set_xticklabels([lab for lab, _ in steps], fontsize=4.0)
ax.set_ylabel("genera", fontsize=5.2)
ax.set_ylim(0, n_all * 1.18)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
lbl(ax, "K", "Few of 267 replicate")

# ------------------------------------------------------------------ L forest, one unit
# Per-cohort mean log-ratio and the pooled log-ratio share ONE axis and one unit. An
# earlier version plotted per-cohort Cohen's d_z against a pooled log-ratio, which put a
# standardised and an unstandardised quantity on the same scale.
ax = axes[4, 1]
# Effects are converted to log2 so this panel, Table S1 and panel J of Fig. 5 share one unit.
# The meta-analysis pools natural-log ratios; a factor of ln2 sat between the figure's axis
# label and the numbers in the table before this conversion.
L2 = np.log(2)
REP = REP.assign(pooled_effect=REP.pooled_effect / L2, pooled_se=REP.pooled_se / L2)
eff = eff.assign(effect=eff.effect / L2, se=eff.se / L2)
gen = REP.genus.tolist()
ys = np.arange(len(gen))
for i, g in enumerate(gen):
    d = eff[(eff.genus == g) & (eff.study.isin(M_STUDIES))]
    for _, r in d.iterrows():
        ax.plot([r["effect"] - 1.96 * r["se"], r["effect"] + 1.96 * r["se"]], [i, i],
                color=M_SCOL[r["study"]], lw=0.7, alpha=0.9, zorder=2)
        ax.scatter(r["effect"], i, s=8, color=M_SCOL[r["study"]], edgecolor="0.25",
                   linewidth=0.2, zorder=3)
    row = REP[REP.genus == g].iloc[0]
    ax.scatter(row["pooled_effect"], i, marker="D", s=17, color="0.12", zorder=5)
    ax.plot([row["pooled_effect"] - 1.96 * row["pooled_se"],
             row["pooled_effect"] + 1.96 * row["pooled_se"]], [i, i],
            color="0.12", lw=0.9, zorder=4)
# reference band: |d_z| = 0.8 expressed on this axis, i.e. the conventional "large effect" mark
ax.axvspan(-M_LARGE_LOG, M_LARGE_LOG, color="0.93", zorder=0)
ax.axvline(0, color="0.35", lw=0.8, zorder=2)
ax.set_yticks(ys)
ax.set_yticklabels(gen, fontsize=3.9, style="italic")
ax.set_xlabel("mean log$_2$ ratio (± 95% CI; diamonds pooled)", fontsize=4.4, labelpad=1.5)
ax.set_ylim(-3.4, len(gen) - 0.3)
hd = [plt.Line2D([], [], marker="o", ls="", color=M_SCOL[k], ms=2.6, label=M_SLAB[k])
      for k in M_SCOL]
hd.append(plt.Line2D([], [], marker="D", ls="", color="0.12", ms=3.0, label="pooled"))
# Below the axis, not inside it: at this row height an in-axes legend sat on top of the
# bottom three genera, which are the ones the figure is about.
ax.legend(handles=hd, frameon=False, fontsize=3.6, loc="upper center",
          bbox_to_anchor=(0.5, -0.42), ncol=3, handletextpad=0.16,
          columnspacing=0.5, labelspacing=0.10, borderaxespad=0.0)
lbl(ax, "L", f"{len(gen)} genera, {K_COH}+ cohorts")

# ------------------------------------------------------------------ M heterogeneity
# Direction replicates even where magnitude does not.
ax = axes[4, 2]
i2 = REP.sort_values("I2_pct")
ax.scatter(i2.I2_pct, np.arange(len(i2)), s=14, color=BLUE, edgecolor="0.25",
           linewidth=0.25, zorder=3)
ax.set_yticks(np.arange(len(i2)))
ax.set_yticklabels(i2.genus, fontsize=3.9, style="italic")
ax.axvline(50, color="0.6", lw=0.8, ls=":", zorder=1)
ax.set_xlim(-4, 100)
ax.set_xlabel("heterogeneity I² (%)", fontsize=4.6)
lbl(ax, "M", "Direction, not magnitude")

out = FIG / "manuscript-v6" / "main" / "Fig4-replication-and-meta-analysis.png"
fig.savefig(out, dpi=330, facecolor="white")
print(f"wrote {out.relative_to(ROOT)}")
try:
    panel_crops(fig)
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
