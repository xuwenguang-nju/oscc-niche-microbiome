#!/usr/bin/env python3
"""F7 — the platform check: can RNA-seq-mined microbiome data stand in for 16S?

Twelve panels (4x3). The answer is nuanced and drives the paper's strongest methodological
claim: the two platforms AGREE on broad genus structure but DISAGREE on specific taxa and
on host coupling, so a platform-specific result cannot be carried across.

A  genus abundance ranked by TCMA vs by amplicon 16S, with Spearman rho per cohort
B  top-20 overlap (Jaccard) between the two platforms
C  genera where the two platforms disagree most (rank gap)
D  TCMA vs 16S abundance scatter for the worst-disagreeing genera
E  direction check: TCMA log2 FC against the 16S pooled effect
F  what the TCMA "microbiome" signal is actually made of (oral vs background)
G  host coupling: number of significant genus-programme pairs per assay
H  the same as a heat map, genus x assay (min q)
I  Parvimonas coupling: the null case, 0 of 200 tests
J  distribution of coupling |R| across all 2000 tests
K  why Parvimonas is not detectable in TCMA: abundance / prevalence
L  the genus-level agreement, summarised across all three amplicon cohorts

Reads:  results/tables/s11_genus_concordance.csv, s11_genus_concordance_disagreement.csv,
        s6a_direction_check.csv, s6b_signal_composition.csv, s6v4_genus_summary.csv,
        s6v4_coupling_long.csv, s6v4_parvimonas_top.csv, s6v2_tcma_genus.csv
Writes: figures/manuscript-v6/main/Fig6-platform-concordance.png
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
FIG = ROOT / "figures" / "manuscript-v6" / "main"
# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
exec(compile(open(_FIGSTYLE / "runtime.py", encoding="utf-8").read(),
              "runtime.py", "exec"))

apply_figure_style(sizes=(7.0, 6.2, 5.4))

conc = pd.read_csv(TAB / "s11_genus_concordance.csv")
dis = pd.read_csv(TAB / "s11_genus_concordance_disagreement.csv")
dchk = pd.read_csv(TAB / "s6a_direction_check.csv")
sig = pd.read_csv(TAB / "s6b_signal_composition.csv")
gs = pd.read_csv(TAB / "s6v4_genus_summary.csv")
coup = pd.read_csv(TAB / "s6v4_coupling_long.csv")
tcma = pd.read_csv(TAB / "s6v2_tcma_genus.csv", index_col=0)

BLUE, RED, ORANGE, GREEN = "#2166ac", "#b2182b", "#a6611a", "#1b7837"
ASSAYS = ["Methylation27K", "RNAseq", "RPPA", "miRNA"]
ALAB = {"Methylation27K": "Methyl", "RNAseq": "RNAseq", "RPPA": "RPPA", "miRNA": "miRNA"}

fig = plt.figure(figsize=(7.09, 9.05))
gsx = fig.add_gridspec(4, 3, left=0.134, right=0.985, top=0.975, bottom=0.042,
                       wspace=0.60, hspace=0.66)
axes = np.array([[fig.add_subplot(gsx[i, j]) for j in range(3)] for i in range(4)])


def lbl(ax, letter, text):
    ax.set_title(text, fontsize=6.2, multialignment="center", pad=5)
    panel_letter_fixed(ax, letter, case="upper")


# ------------------------------------------------------------------ A rank agreement
ax = axes[0, 0]
c = conc.set_index("cohort")
xs = np.arange(len(c))
ax.bar(xs, c["spearman_rho_filtered"], color=BLUE, edgecolor="0.3", lw=0.4,
       width=0.55, zorder=3)
for x, v in zip(xs, c["spearman_rho_filtered"]):
    ax.text(x, v + 0.022, f"{v:.2f}", ha="center", fontsize=5.2, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels([i.replace("PRJNA", "") for i in c.index], rotation=42, ha="right",
                   fontsize=4.9)
ax.set_ylabel("Spearman rho, genus abundance\nTCMA vs amplicon 16S", fontsize=5.6)
ax.set_ylim(0, 1)
lbl(ax, "A", "Broad structure agrees")

# ------------------------------------------------------------------ B top-20 overlap
ax = axes[0, 1]
ax.bar(xs, c["top20_jaccard"], color=GREEN, edgecolor="0.3", lw=0.4, width=0.55, zorder=3)
for x, (v, n) in enumerate(zip(c["top20_jaccard"], c["top20_overlap_n"])):
    ax.text(x, v + 0.02, f"{n}/20", ha="center", fontsize=5.0, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels([i.replace("PRJNA", "") for i in c.index], rotation=42, ha="right",
                   fontsize=4.9)
ax.set_ylabel("Jaccard of top-20 genera", fontsize=5.6)
ax.set_ylim(0, 1)
lbl(ax, "B", "but not the leading taxa")

# ------------------------------------------------------------------ C rank gap
ax = axes[0, 2]
d = dis.copy()
d = d.reindex(d["rank_gap"].abs().sort_values(ascending=False).index).head(12)
ys = np.arange(len(d))
cols = [RED if g > 0 else BLUE for g in d["rank_gap"]]
ax.barh(ys, d["rank_gap"], color=cols, edgecolor="0.3", lw=0.35, height=0.6, zorder=3)
ax.axvline(0, color="0.3", lw=0.8, zorder=4)
ax.set_yticks(ys)
# The accession was truncated to its last four digits here, which turned PRJNA866676 into "6676"
# while panels A and L print the full number. Same label everywhere now, without the PRJNA prefix.
ax.set_yticklabels([f"{g} ({c.replace('PRJNA', '')})"
                    for g, c in zip(d["genus"], d["cohort"])],
                   fontsize=4.0, style="italic")
ax.invert_yaxis()
ax.set_xlabel("rank gap (TCMA rank − 16S rank)", fontsize=5.6)
lbl(ax, "C", "Specific genera diverge")

# ------------------------------------------------------------------ D abundance scatter
ax = axes[0, 0 + 2] if False else axes[1, 0]
d = dis.dropna(subset=["tcma_mean_pct", "cohort_mean_pct"])
ax.scatter(d["cohort_mean_pct"], d["tcma_mean_pct"], s=16, color=BLUE, alpha=0.6,
           edgecolor="none", zorder=3)
lim = float(max(d["cohort_mean_pct"].max(), d["tcma_mean_pct"].max())) * 1.1
ax.plot([0, lim], [0, lim], color="0.55", lw=0.9, ls=(0, (4, 1.6)), zorder=2)
for _, r in d.reindex(d["rank_gap"].abs().sort_values(ascending=False).index).head(4).iterrows():
    ax.annotate(r["genus"], (r["cohort_mean_pct"], r["tcma_mean_pct"]),
                textcoords="offset points", xytext=(3, 3), fontsize=4.3, color="0.3",
                style="italic")
ax.set_xlabel("mean % in amplicon 16S", fontsize=5.6)
ax.set_ylabel("mean % in TCMA", fontsize=5.6)
lbl(ax, "D", "Off-diagonal genera")

# ------------------------------------------------------------------ E direction check
ax = axes[1, 1]
d = dchk.dropna(subset=["pooled_effect", "rna_cancer_minus_control"])
agree = d["direction_agrees"] == True                                   # noqa: E712
ax.scatter(d.loc[agree, "pooled_effect"], d.loc[agree, "rna_cancer_minus_control"],
           s=22, color=GREEN, edgecolor="0.25", linewidth=0.3, zorder=3, label="agrees")
ax.scatter(d.loc[~agree, "pooled_effect"], d.loc[~agree, "rna_cancer_minus_control"],
           s=22, color=RED, edgecolor="0.25", linewidth=0.3, zorder=3, label="opposite")
ax.axhline(0, color="0.5", lw=0.8, zorder=2)
ax.axvline(0, color="0.5", lw=0.8, zorder=2)
ax.legend(frameon=False, fontsize=4.4, loc="lower right", handletextpad=0.2,
          labelspacing=0.25, borderaxespad=0.4)
# Place labels with a small greedy anti-overlap pass: nine genera sit within a few
# hundredths of the origin, so fixed offsets always stack some of them.
_labels = [(r["genus"][:12], r["pooled_effect"], r["rna_cancer_minus_control"])
           for _, r in d.iterrows()]
_taken: list[tuple[float, float]] = []
_cands = [(5, 4), (5, -8), (-34, 4), (-34, -8), (5, 12), (-34, 12), (5, -16), (-34, -16)]
for _txt, _x, _y in _labels:
    for _dx, _dy in _cands:
        _hit = False
        for _px, _py in _taken:
            if abs(_dx - _px) < 26 and abs(_dy - _py) < 7:
                _hit = True
                break
        if not _hit:
            _taken.append((_dx, _dy))
            ax.annotate(_txt, (_x, _y), textcoords="offset points",
                        xytext=(_dx, _dy), fontsize=3.9, color="0.35",
                        arrowprops=dict(arrowstyle="-", lw=0.4, color="0.62",
                                        shrinkA=0, shrinkB=1.5))
            break
ax.set_xlabel("16S pooled effect (Cohen's $d_z$)", fontsize=5.6)
# The substrate must be named: this panel reads s6a_direction_check.csv, which is an
# independent total-RNA-seq oral cohort (GSE227919), NOT TCMA. Saying "TCMA" here is what
# let E and Fig. S2A be swapped for several revisions.
ax.set_ylabel("RNA-seq log$_2$ cancer − control", fontsize=5.6)
lbl(ax, "E", "Direction often flips (RNA-seq)")

# ------------------------------------------------------------------ F signal composition
ax = axes[1, 2]
bc = sig.groupby("bucket")["mean_rel_abundance"].sum().sort_values(ascending=False)
COLC = {"oral_typical": GREEN, "environmental_background": "#bdbdbd", "other": ORANGE}
ax.bar(np.arange(len(bc)), bc.values * 100, color=[COLC.get(b, BLUE) for b in bc.index],
       edgecolor="0.3", lw=0.4, width=0.55, zorder=3)
for x, v in enumerate(bc.values * 100):
    ax.text(x, v + 1.4, f"{v:.0f}%", ha="center", fontsize=5.2, color="0.25")
ax.set_xticks(np.arange(len(bc)))
ax.set_xticklabels([b.replace("_", "\n") for b in bc.index], fontsize=4.7)
ax.set_ylabel("% of TCMA microbial reads", fontsize=5.6)
ax.set_ylim(0, max(bc.values * 100) * 1.24)
ax.set_xlabel("taxon provenance", fontsize=5.6)
lbl(ax, "F", "Half is background")

# ------------------------------------------------------------------ G significant pairs
ax = axes[2, 0]
cnt = gs.groupby("assay")["q05"].sum().reindex(ASSAYS).fillna(0)
ax.bar(np.arange(len(cnt)), cnt.values, color=[BLUE, ORANGE, GREEN, "#7b3294"],
       edgecolor="0.3", lw=0.4, width=0.55, zorder=3)
for x, v in enumerate(cnt.values):
    ax.text(x, v + 0.35, f"{int(v)}", ha="center", fontsize=5.2, color="0.25")
ax.set_xticks(np.arange(len(cnt)))
ax.set_xticklabels([ALAB[a] for a in cnt.index], fontsize=5.0)
ax.set_ylabel("genus–feature pairs, q<0.05", fontsize=5.6)
ax.set_ylim(0, max(cnt.values) * 1.28)
ax.text(0.97, 0.93, f"{int(cnt.sum())}/2000 in total", transform=ax.transAxes,
        fontsize=4.6, ha="right", va="top", color="0.3")
lbl(ax, "G", "Coupling is assay-specific")

# ------------------------------------------------------------------ H coupling heat map
ax = axes[2, 1]
M = gs.pivot(index="genus", columns="assay", values="q05").reindex(columns=ASSAYS)
M = M.reindex(gs.groupby("genus")["q05"].sum().sort_values(ascending=False).index)
M = M.fillna(0)
im = ax.imshow(M.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=max(1, M.values.max()))
ax.set_xticks(np.arange(len(ASSAYS)))
ax.set_xticklabels([ALAB[a] for a in ASSAYS], fontsize=4.8)
ax.set_yticks(np.arange(len(M)))
ax.set_yticklabels(M.index, fontsize=4.3, style="italic")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = int(M.values[i, j])
        if v:
            ax.text(j, i, str(v), ha="center", va="center", fontsize=3.9,
                    color="white" if v > M.values.max() * 0.6 else "0.15")
cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.03)
cb.ax.tick_params(labelsize=4.2)
cb.set_label("pairs q<0.05", fontsize=4.5)
lbl(ax, "H", "Only Prevotella is shared")

# ------------------------------------------------------------------ I Parvimonas null
ax = axes[2, 2]
pt = pd.read_csv(TAB / "s6v4_parvimonas_top.csv")
pt = pt.reindex(pt["q"].sort_values().index).head(10)
ys = np.arange(len(pt))
ax.barh(ys, -np.log10(pt["q"].clip(lower=1e-4)), color="#bdbdbd", edgecolor="0.3",
        lw=0.35, height=0.6, zorder=3)
ax.axvline(-np.log10(0.05), color=RED, lw=1.0, ls=(0, (4, 1.6)), zorder=4)
ax.set_yticks(ys)
ax.set_yticklabels([f"{f[:20]}" for f in pt["feature"]], fontsize=4.0)
ax.invert_yaxis()
ax.set_xlabel("−log10 q (best Parvimonas couplings)", fontsize=5.6)
ax.set_xlim(0, 1.45)
ax.text(0.97, 0.28, f"0/200 pass FDR\nmin q = {pt['q'].min():.3f}",
        transform=ax.transAxes, fontsize=4.4, ha="right", va="top", color=RED)
lbl(ax, "I", "Parvimonas couples with nothing")

# ------------------------------------------------------------------ J |R| distribution
ax = axes[3, 0]
ax.hist(coup["R"].abs(), bins=30, color=BLUE, edgecolor="white", lw=0.3, zorder=3)
pv = coup[coup.genus == "Parvimonas"]["R"].abs()
ax.hist(pv, bins=30, color=RED, alpha=0.75, edgecolor="none", zorder=4)
ax.axvline(pv.max(), color=RED, lw=1.0, ls=(0, (4, 1.6)), zorder=5)
ax.set_xlabel("|R| of genus–feature coupling", fontsize=5.6)
ax.set_ylabel("tests", fontsize=5.6)
ax.text(0.97, 0.93, f"Parvimonas max |R| = {pv.max():.2f}", transform=ax.transAxes,
        fontsize=4.5, ha="right", va="top", color=RED)
lbl(ax, "J", "Its signal is not extreme")

# ------------------------------------------------------------------ K why it is invisible
ax = axes[3, 1]
d = dchk.dropna(subset=["median_rel_abundance", "pct_samples_absent"])
ax.scatter(d["median_rel_abundance"], d["pct_samples_absent"], s=22, color=BLUE,
           edgecolor="0.25", linewidth=0.3, zorder=3)
# The four genera with a zero median sit in one column, so their labels are stacked with a fixed
# vertical pitch and a leader line back to the point; without it the names overlap each other.
_zero = d[d["median_rel_abundance"] == 0].sort_values("pct_samples_absent", ascending=False)
for k, (_, r) in enumerate(_zero.iterrows()):
    ax.annotate(r["genus"][:16], (0.0, float(r["pct_samples_absent"])),
                textcoords="offset points", xytext=(16, 0), fontsize=3.9, color="0.35",
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.4, color="0.65", shrinkA=1, shrinkB=2))
for _, r in d[d["median_rel_abundance"] > 0].iterrows():
    ax.annotate(r["genus"][:16], (float(r["median_rel_abundance"]), float(r["pct_samples_absent"])),
                textcoords="offset points", xytext=(5, 4), fontsize=3.9, color="0.35", ha="left")
# A log x-axis silently dropped every genus whose median abundance is zero - which includes
# *Parvimonas*, the one genus this figure is about. Five of the nine have a zero median, so the
# axis is linear and x = 0 means "the median sample carries no reads assigned to this genus";
# those genera stack in one column at x = 0 and are labelled individually.
ax.set_xlim(-0.0006, max(d["median_rel_abundance"]) * 1.30)
# Both axes come from the same total-RNA-seq direction-check cohort (GSE227919): the builder writes
# median_rel_abundance and pct_samples_absent from that cohort's relative-abundance table, not from
# TCMA. The axis used to say "in TCMA", which would have made this a cross-substrate panel it is not.
ax.set_xlabel("median relative abundance in the RNA-seq cohort (%)\n"
              "0 = half the samples carry none", fontsize=5.4)
ax.set_ylabel("% samples where absent", fontsize=5.6)
lbl(ax, "K", "Present, just not coupled")

# ------------------------------------------------------------------ L agreement summary
ax = axes[3, 2]
v = c["spearman_rho_filtered"].values
ax.plot(np.arange(len(v)), v, color=BLUE, lw=1.4, marker="o", ms=3.6, zorder=3)
ax.fill_between(np.arange(len(v)), 0, v, color=BLUE, alpha=0.12, zorder=2)
ax.axhline(np.mean(v), color=ORANGE, lw=1.0, ls=(0, (4, 1.6)), zorder=4)
ax.text(0.03, np.mean(v) + 0.035, f"mean rho = {np.mean(v):.2f}", fontsize=4.8,
        color=ORANGE, va="bottom")
ax.set_xticks(np.arange(len(v)))
ax.set_xticklabels([i.replace("PRJNA", "") for i in c.index], rotation=42, ha="right",
                   fontsize=4.9)
ax.set_ylim(0, 1)
ax.set_ylabel("Spearman rho per cohort", fontsize=5.6)
ax.set_xlabel("amplicon cohort", fontsize=5.6)
lbl(ax, "L", "Agreement is reproducible")

out = FIG / "Fig6-platform-concordance.png"
fig.savefig(out, dpi=330, facecolor="white")
print(f"wrote {out.relative_to(ROOT)}")
try:
    panel_crops(fig)
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
