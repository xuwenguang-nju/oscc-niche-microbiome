#!/usr/bin/env python3
"""F3 — which ecological niche is the intratumoral community closest to? (PRJNA866676)

Nine panels, moving from "the niches are distinguishable at all" to "which one is the
source", then to the two falsification checks.

A  source ranking: mean rank of distance to that patient's own tumour
B  surface-to-depth gradient: deep tumour vs cancer surface, distance to saliva
C  niche separation (PERMANOVA) and the ordering of pairwise distances
D  the six mean profiles (top discriminative genera), as a heat map
E  genera that discriminate the niches (effect size)
F  reagent-contaminant burden per niche — the confound that had to be checked
G  the ranking before vs after decontamination (sensitivity)
H  paired distances per patient for the closest contrast (tumour vs paracancerous)
I  how many patients place each niche first

Reads:  results/tables/s7_source_ranking.csv, s7_source_ranking_tests.csv,
        s7_backflow_gradient.csv, s7_pairwise_niche_distances.csv, s7_summary.json,
        s7_niche_composition.csv, s7_discriminative_genera.csv,
        s7b_contaminant_burden.csv, s7b_before_after.csv, s7_paired_distances.csv
Writes: figures/main/F3-source-attribution.png
"""

from __future__ import annotations

import json
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

rank = pd.read_csv(TAB / "s7_source_ranking.csv")
rank_t = pd.read_csv(TAB / "s7_source_ranking_tests.csv")
back = pd.read_csv(TAB / "s7_backflow_gradient.csv")
pair = pd.read_csv(TAB / "s7_pairwise_niche_distances.csv")
summ = json.load(open(TAB / "s7_summary.json"))
comp = pd.read_csv(TAB / "s7_niche_composition.csv", index_col=0)
disc = pd.read_csv(TAB / "s7_discriminative_genera.csv")
cont = pd.read_csv(TAB / "s7b_contaminant_burden.csv")
ba = pd.read_csv(TAB / "s7b_before_after.csv")
pdist = pd.read_csv(TAB / "s7_paired_distances.csv")

BLUE, RED, ORANGE = "#2166ac", "#b2182b", "#a6611a"
NICE = {"paracancerous": "paracancerous", "cancer_surface": "cancer surface",
        "saliva": "saliva", "contralateral_normal_mucosa": "contralateral\nnormal",
        "tongue_coat": "tongue coat", "tumor": "tumour (deep)"}
NCOL = {"paracancerous": "#a6611a", "cancer_surface": "#e66101", "saliva": "#2166ac",
        "contralateral_normal_mucosa": "#1b7837", "tongue_coat": "#92c5de",
        "tumor": "#7b3294"}

fig = plt.figure(figsize=(7.09, 7.15))
gs = fig.add_gridspec(3, 3, left=0.134, right=0.985, top=0.962, bottom=0.072,
                      wspace=0.95, hspace=0.66)
axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(3)] for i in range(3)])


def lbl(ax, letter, text):
    ax.set_title(text, fontsize=6.2, multialignment="center", pad=5)
    panel_letter_fixed(ax, letter, case="upper")


# ------------------------------------------------------------------ A source ranking
ax = axes[0, 0]
r = rank[rank.source_niche != "tumor"]
order = r.groupby("source_niche")["rank"].mean().sort_values().index
mu = r.groupby("source_niche")["rank"].mean().reindex(order)
ptab = rank_t.set_index("source_niche")["p_vs_others"].to_dict()
xs = np.arange(len(order))
ax.barh(xs, mu.values, color=[NCOL[x] for x in order], edgecolor="0.3", lw=0.4,
        height=0.62, zorder=3)
ax.axvline(3, color="0.4", lw=0.9, ls=(0, (4, 1.6)), zorder=4)
for i, (k, v) in enumerate(mu.items()):
    p = ptab.get(k, np.nan)
    mark = "ns" if not np.isfinite(p) else ("***" if p < 0.001 else
                                            "**" if p < 0.01 else
                                            "*" if p < 0.05 else "ns")
    ax.text(v + 0.05, i, f"{v:.2f} {mark}", va="center", fontsize=4.9, color="0.25")
ax.set_yticks(xs)
ax.set_yticklabels([NICE.get(x, x) for x in order], fontsize=5.4)
ax.invert_yaxis()
ax.set_xlim(0, max(mu.values) * 1.32)
ax.set_xlabel("mean rank of distance to own tumour", fontsize=5.8)
ax.text(3.10, 0.42, "random (rank 3)", fontsize=4.6, color="0.45",
        va="bottom", ha="left")
lbl(ax, "A", "Adjacent tissue is closest")

# ------------------------------------------------------------------ B backflow gradient
ax = axes[0, 1]
d = back.dropna(subset=["d_tumour_saliva", "d_surface_saliva"])
for _, r in d.iterrows():
    ax.plot([0, 1], [r["d_tumour_saliva"], r["d_surface_saliva"]],
            color="0.72", lw=0.7, zorder=2)
md = d["d_tumour_saliva"].mean()
ms = d["d_surface_saliva"].mean()
ax.plot([0, 1], [md, ms], color=RED, lw=2.2, marker="o", ms=4.2, zorder=4)
ax.set_xticks([0, 1])
ax.set_xticklabels(["deep\ntumour", "cancer\nsurface"], fontsize=5.4)
ax.set_ylabel("Bray-Curtis distance to saliva", fontsize=5.8)
ax.set_xlim(-0.28, 1.28)
frac = float((d["d_surface_saliva"] < d["d_tumour_saliva"]).mean())
ax.text(0.5, d[["d_tumour_saliva", "d_surface_saliva"]].to_numpy().max() * 0.99,
        f"{frac*100:.0f}% closer (n={len(d)})",
        ha="center", fontsize=5.2, color=RED)
lbl(ax, "B", "Surface is closer to saliva")

# ------------------------------------------------------------------ C pairwise distances
ax = axes[0, 2]
p = pair.copy()
p["label"] = p.niche_a.map(NICE).str.replace("\n", " ") + " / " + \
             p.niche_b.map(NICE).str.replace("\n", " ")
p = p.sort_values("mean_bray_curtis", ascending=False).head(7)
ys = np.arange(len(p))
ax.barh(ys, p["mean_bray_curtis"], color="#9ecae1", edgecolor="0.3", lw=0.4,
        height=0.62, zorder=3)
ax.set_yticks(ys)
ax.set_yticklabels(p["label"], fontsize=4.7)
ax.invert_yaxis()
ax.set_xlabel("mean Bray-Curtis between niches", fontsize=5.8)
pf = summ["permanova"]["blocked"]
ax.set_xlabel(f"mean Bray-Curtis between niches\n"
              f"PERMANOVA F={pf['pseudo_F']:.1f}, R²={pf['r2']:.2f}, p={pf['p']:.4f}",
              fontsize=5.6)
lbl(ax, "C", "Niches are distinguishable")

# ------------------------------------------------------------------ D profile heat map
ax = axes[1, 0]
top = disc.groupby("genus")["delta_rel_abundance"].max().sort_values(ascending=False).head(12).index
M = comp.loc[:, [g for g in top if g in comp.columns]]
order_n = ["tumor", "paracancerous", "cancer_surface",
           "contralateral_normal_mucosa", "saliva", "tongue_coat"]
M = M.reindex(order_n)
Z = M.div(M.max(axis=0), axis=1)
im = ax.imshow(Z.values, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
ax.set_xticks(np.arange(Z.shape[1]))
ax.set_xticklabels([g[:13] for g in Z.columns], rotation=55, ha="right", fontsize=4.4)
ax.set_yticks(np.arange(len(order_n)))
ax.set_yticklabels([NICE[n].replace("\n", " ") for n in order_n], fontsize=4.9)

cb = fig.colorbar(im, ax=ax, fraction=0.024, pad=0.030)
cb.ax.tick_params(labelsize=4.4)
cb.set_label("relative to niche max", fontsize=4.6)
lbl(ax, "D", "Six mean profiles")

# ------------------------------------------------------------------ E discriminative genera
ax = axes[1, 1]
e = disc.groupby("genus")["delta_rel_abundance"].max().sort_values(ascending=False).head(10)
ys = np.arange(len(e))[::-1]
ax.barh(ys, e.values, color="#7b3294", edgecolor="0.3", lw=0.4, height=0.62, zorder=3)
ax.set_yticks(ys)
ax.set_yticklabels(e.index, fontsize=4.9)
ax.set_xlabel("max enrichment over other niches\n(top hit is a reagent contaminant)",
              fontsize=5.4)
lbl(ax, "E", "What distinguishes niches")

# ------------------------------------------------------------------ F contaminant burden
ax = axes[1, 2]
c = cont.copy().sort_values("contaminant_pct", ascending=False)
ys = np.arange(len(c))
cols = [RED if v > 20 else "#9ecae1" for v in c["contaminant_pct"]]
ax.barh(ys, c["contaminant_pct"], color=cols, edgecolor="0.3", lw=0.4, height=0.62, zorder=3)
ax.set_yticks(ys)
ax.set_yticklabels([NICE.get(x, x).replace("\n", " ") for x in c.niche], fontsize=5.0)
ax.invert_yaxis()
for y, v in zip(ys, c["contaminant_pct"]):
    ax.text(v + 0.9, y, f"{v:.1f}%", va="center", fontsize=4.8, color="0.25")
ax.set_xlim(0, max(c["contaminant_pct"]) * 1.26)
ax.set_xlabel("reagent-contaminant share (%)", fontsize=5.8)
lbl(ax, "F", "Contamination is uneven")

# ------------------------------------------------------------------ G before vs after
ax = axes[2, 0]
nn = ["paracancerous", "cancer_surface", "saliva",
      "contralateral_normal_mucosa", "tongue_coat"]
b = ba.set_index("analysis")
before = [b.loc["before", f"rank_{n}"] for n in nn]
after = [b.loc["after", f"rank_{n}"] for n in nn]
xs = np.arange(len(nn))
ax.plot(xs, before, color="0.55", lw=1.4, marker="o", ms=3.0, label="before", zorder=3)
ax.plot(xs, after, color=BLUE, lw=1.6, marker="s", ms=3.0, label="after", zorder=4)
ax.axhline(3, color="0.5", lw=0.8, ls=(0, (4, 1.6)), zorder=2)
ax.set_xticks(xs)
ax.set_xticklabels([NICE[n].replace("\n", "\n") for n in nn], fontsize=4.7, rotation=30,
                   ha="right")
ax.set_ylabel("mean rank to own tumour", fontsize=5.8)
ax.legend(frameon=False, fontsize=4.9, loc="upper left", handlelength=1.1,
          labelspacing=0.2, borderaxespad=0.2)
lbl(ax, "G", "Holds after decontamination")

# ------------------------------------------------------------------ H paired contrast
ax = axes[2, 1]
d = pdist[pdist.source_niche == "paracancerous"].dropna(subset=["bray_curtis"])
ax.hist(d["bray_curtis"], bins=16, color="#a6611a", edgecolor="white", lw=0.5, zorder=3)
ax.axvline(d["bray_curtis"].median(), color=RED, lw=1.2, ls=(0, (4, 1.6)), zorder=4)
ax.set_xlabel("Bray-Curtis, tumour vs its own paracancerous", fontsize=5.8)
ax.set_ylabel("patients", fontsize=5.8)
# reserve a clear band above the tallest bar so the two notes never touch a bar
top = ax.get_ylim()[1]
ax.set_ylim(0, top * 1.34)
ax.text(d["bray_curtis"].median() + 0.014, top * 1.26,
        f"median {d['bray_curtis'].median():.2f}", color=RED, fontsize=4.8,
        ha="left", va="top")
ax.text(0.03, 0.93, f"n={len(d)} pairs", transform=ax.transAxes,
        fontsize=4.9, va="top", color="0.3")
lbl(ax, "H", "Distance to adjacent tissue")

# ------------------------------------------------------------------ I first-place count
ax = axes[2, 2]
first = rank[rank["rank"] == 1].groupby("source_niche").size()
first = first.reindex([n for n in NICE if n != "tumor"]).fillna(0)
ys = np.arange(len(first))
ax.barh(ys, first.values, color=[NCOL[x] for x in first.index], edgecolor="0.3",
        lw=0.4, height=0.62, zorder=3)
for y, v in zip(ys, first.values):
    ax.text(v + 0.25, y, f"{int(v)}", va="center", fontsize=5.0, color="0.25")
ax.set_yticks(ys)
ax.set_yticklabels([NICE[x].replace("\n", " ") for x in first.index], fontsize=5.0)
ax.invert_yaxis()
ax.set_xlim(0, max(first.values) * 1.28)
ax.set_xlabel("patients ranking it first", fontsize=5.8)
lbl(ax, "I", "How often each wins")

out = FIG / "manuscript-v6" / "main" / "Fig3-six-niche-attribution.png"
fig.savefig(out, dpi=330, facecolor="white")
print(f"wrote {out.relative_to(ROOT)}")
try:
    panel_crops(fig)
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
