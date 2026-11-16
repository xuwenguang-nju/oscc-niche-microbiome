"""Figure: niche attribution of the OSCC intratumoral microbiome (PRJNA1415359).

Two panels, matching the two readouts that survived S3:

  a  distance from each tumour to each niche, paired per patient -> the ranking
     (H1). Shows that tumour sits closest to adjacent normal mucosa and furthest
     from saliva, with every patient contributing a paired comparison.
  b  the reason H2 fails: effect size against adjacent normal is uniformly small.
     A volcano of mean log-ratio vs Cohen's dz, with the large-effect band marked,
     so the null is visibly "small effect" rather than "underpowered".

Reads:
  results/tables/s4b_distance_own_patient_PRJNA1415359.csv
  results/tables/s4a_residual_taxa_PRJNA1415359.csv
Writes:
  figures/s4-niche-attribution-PRJNA1415359.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures"
# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
exec(compile(open(_FIGSTYLE / "runtime.py", encoding="utf-8").read(),
              "runtime.py", "exec"))

apply_figure_style(sizes=(9, 7.5, 6.5))

d = pd.read_csv(TAB / "s4b_distance_own_patient_PRJNA1415359.csv", index_col=0)
res = pd.read_csv(TAB / "s4a_residual_taxa_PRJNA1415359.csv")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3))

# ---------------- panel a: paired distance by niche ----------------
ax = axes[0]
NICHE_LABEL = {"adjacent_normal_tissue": "adjacent\nnormal mucosa",
               "preop_saliva": "pre-operative\nsaliva"}
cols = {"adjacent_normal_tissue": "#1f4e79", "preop_saliva": "#a6611a"}
order = ["adjacent_normal_tissue", "preop_saliva"]
xx = np.arange(len(order))

for i, patient in enumerate(d.index):
    vals = [d.loc[patient, n] for n in order]
    if any(pd.isna(v) for v in vals):
        continue
    ax.plot(xx, vals, color="0.72", lw=0.6, alpha=0.55, zorder=1)
for j, n in enumerate(order):
    v = d[n].dropna().values
    ax.scatter(np.full(v.size, j) + np.linspace(-0.06, 0.06, v.size), v,
               s=13, color=cols[n], zorder=3, edgecolors="white", linewidths=0.4)
    ax.plot([j - 0.18, j + 0.18], [np.median(v)] * 2, color=cols[n], lw=2.2, zorder=4)

ax.set_xticks(xx)
ax.set_xticklabels([NICHE_LABEL[n] for n in order])
ax.set_ylabel("distance (Aitchison)", labelpad=8)
ax.set_title("Tumour sits closer to adjacent mucosa\nthan to saliva, patient by patient", multialignment="center", loc="center")
panel_letter(ax, "A", case="upper")
ax.margins(0.08)
ax.set_ylim(bottom=0)
ax.text(0.5, -0.20, "paired Wilcoxon p = 0.010, n = 31", transform=ax.transAxes,
        fontsize=6.6, va="top", ha="center", color="0.3")

# ---------------- panel b: effect size vs adjacent normal ----------------
ax = axes[1]
x = res["cohens_dz_vs_normal"].values
y = res["mean_logratio_vs_normal"].values
# grey: no sizeable effect; blue: large effect band would be |dz| >= 0.8
big = np.abs(x) >= 0.8
ax.axvspan(-0.8, 0.8, color="0.90", zorder=0)
ax.axvline(0, color="0.65", lw=0.8, zorder=1)
ax.axhline(0, color="0.65", lw=0.8, zorder=1)
ax.scatter(x[~big], y[~big], s=14, color="#6b8fb0", alpha=0.85,
           edgecolors="white", linewidths=0.3, zorder=3)
if big.any():
    ax.scatter(x[big], y[big], s=22, color="#b2182b", zorder=4)

# label the three largest effects, both directions
# The three largest effects sit within 0.2 log units of each other, so naming each
# one individually would need colliding leaders. The panel's claim is about the
# whole cloud, so the cluster is annotated once with the actual largest |dz|.
idx = np.argsort(-np.abs(y))[:3]
ax.scatter(x[idx], y[idx], s=46, facecolors="none", edgecolors="#b2182b",
           linewidths=0.9, zorder=6)
i_max = idx[np.argmax(np.abs(x[idx]))]
ax.annotate(f"3 largest effects,\nall |$d_z$| < {abs(x[i_max]):.2f}",
            xy=(float(np.median(x[idx])), float(np.median(y[idx]))),
            xytext=(0.04, 0.86), textcoords="axes fraction",
            fontsize=6.3, ha="left", va="top", color="#b2182b", zorder=7,
            arrowprops=dict(arrowstyle="-", color="#b2182b", lw=0.5,
                            shrinkA=2, shrinkB=4,
                            connectionstyle="arc3,rad=-0.25"))

ax.set_xlabel("Cohen's $d_z$ (tumour vs adjacent normal)")
ax.set_ylabel("mean log-ratio\n(tumour vs adjacent normal)")
ax.set_title("No genus reaches a large effect, so the\nnull is small effect, not low power", multialignment="center", loc="center")
panel_letter(ax, "B", case="upper")
ax.margins(0.10)

fig.tight_layout(w_pad=2.4)
out = FIG / "s4-niche-attribution-PRJNA1415359.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")

# ---- prove the render ----
rend = fig.canvas.get_renderer()
tick_texts = set()
for ax_ in fig.axes:
    for t in (ax_.get_xticklabels(which="both") + ax_.get_yticklabels(which="both")):
        tick_texts.add(t.get_text())
labels = [(t, t.get_window_extent(rend)) for t in fig.findobj(matplotlib.text.Text)
          if t.get_text().strip() and t.get_visible() and t.get_text() not in tick_texts]
frames = [(s, s.get_window_extent(rend)) for ax_ in fig.axes
          for s in ax_.spines.values() if s.get_visible()]
hits = [(a, b) for i, (a, ba) in enumerate(labels)
        for b, bb in labels[i + 1:] if ba.overlaps(bb)]
hits += [(t, s) for t, bt in labels for s, bs in frames if bt.overlaps(bs)]
print("file:", out.name, out.stat().st_size, "bytes")
print("non-tick text collisions:", len(hits))
for h in hits[:8]:
    print("   ", getattr(h[0], "get_text", lambda: type(h[0]).__name__)()[:36], "|",
          getattr(h[1], "get_text", lambda: type(h[1]).__name__)()[:36])

save_panel_crops(out, panel_crops(fig))
print("crops: .cache/figure-style/")
