#!/usr/bin/env python3
"""F2 — what can and cannot be recovered from a niche-resolved source-attribution design.

Nine panels, each answering a different identifiability question. The figure exists to
justify one methodological decision: the study reports DISTANCES AND RANKS, never mixing
proportions, because the proportion estimator is not identifiable on this data.

A  mixing fraction is not identifiable (the null sits at the signal level)
B  the recoverable quantity is detectability, and it depends on depth
C  every objective is biased; L2 is blind, chi2 saturates, L1 carries the smallest bias
D  patient-level vs cohort-level spread: cohorts, not patients, are the unit
E  reference leakage (25% of the tumour-specific taxa also present in the sources) is tolerated
F  the minimum detectable tumour-specific FRACTION falls with depth
G  the excess over the null depends on how distinct the two source profiles are
H  residual-taxa call: sensitivity and specificity are flat and high across the range
I  what actually buys detection is sequencing depth, not source separation

Reads:  results/tables/s3_simulation_calibration.csv, s3_objective_comparison.csv,
        s3_perpatient_variability.csv, s3_simulation_mdc.csv
Writes: figures/main/F2-identifiability.png
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

cal = pd.read_csv(TAB / "s3_simulation_calibration.csv")
obj = pd.read_csv(TAB / "s3_objective_comparison.csv")
ppv = pd.read_csv(TAB / "s3_perpatient_variability.csv")
mdc = pd.read_csv(TAB / "s3_simulation_mdc.csv")

DEPTHS = ["PRJNA412445_23k", "PRJNA597251_43k", "PRJNA866676_105k",
          "PRJNA1415359_130k", "PRJNA666746_968k"]
DCOL = dict(zip(DEPTHS, ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08306b"]))
CORE = "PRJNA1415359_130k"
BLUE, RED, ORANGE = "#2166ac", "#b2182b", "#a6611a"

fig = plt.figure(figsize=(7.09, 7.45))
gs = fig.add_gridspec(3, 3, left=0.134, right=0.985, top=0.968, bottom=0.085,
                      wspace=0.52, hspace=0.62)
axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(3)] for i in range(3)])


def lbl(ax, letter, text):
    ax.set_title(text, fontsize=6.2, multialignment="center", pad=5)
    panel_letter_fixed(ax, letter, case="upper")


def core(scenario):
    d = cal[(cal.scenario == scenario) & (cal.depth_label == CORE)]
    return d.sort_values("p_novel_true")


# ------------------------------------------------------------------ A
ax = axes[0, 0]
d = core("moderate_c80")
ax.plot(d["p_novel_true"], d["u_null_median"], color=ORANGE, lw=1.4, marker="s", ms=3.0,
        ls=(0, (4, 1.6)), zorder=3)
ax.plot(d["p_novel_true"], d["u_t_median"], color=BLUE, lw=1.9, marker="o", ms=3.2, zorder=4)
ax.text(0.055, 0.425, "source sample\n(the null)", color=ORANGE, fontsize=5.2,
        ha="left", va="bottom", linespacing=1.15)
ax.text(0.130, 0.470, "tumour", color=BLUE, fontsize=5.4, ha="center", va="bottom")
ax.set_xlabel("true tumour-specific fraction", fontsize=5.8)
ax.set_ylabel("estimated unexplained fraction", fontsize=5.8)
ax.set_xlim(-0.012, 0.315)
lbl(ax, "A", "Mixing fraction is not identifiable")

# ------------------------------------------------------------------ B
ax = axes[0, 1]
for lab in DEPTHS:
    dd = cal[(cal.scenario == "moderate_c80") & (cal.depth_label == lab)]
    dd = dd.sort_values("p_novel_true")
    ax.plot(dd["p_novel_true"], dd["n_scorable_novel_median"], color=DCOL[lab], lw=1.6,
            marker="o", ms=2.4, zorder=3, label=f"{lab.split('_')[-1]}")
ax.legend(title="reads/run", frameon=False, fontsize=4.8, title_fontsize=5.0,
          loc="lower right", handlelength=1.0, labelspacing=0.18, borderaxespad=0.2)
ax.set_xlabel("true tumour-specific fraction", fontsize=5.8)
ax.set_ylabel("tumour-specific genera\ndetected", fontsize=5.8)
lbl(ax, "B", "The limit is detectability")

# ------------------------------------------------------------------ C
ax = axes[0, 2]
d = obj[obj.depth_label == CORE].sort_values("p_novel_true")
OC = {"l1": "#7b3294", "l2": "#008837", "chi2": "#e66101"}
ON = {"l1": "L1", "l2": "L2", "chi2": r"$\chi^2$"}
for o in ["l1", "l2", "chi2"]:
    dd = d[d.objective == o]
    ax.plot(dd["p_novel_true"], dd["bias_vs_truth"], color=OC[o], lw=1.6, marker="o",
            ms=2.6, label=ON[o], zorder=3)
ax.axhline(0, color="0.75", lw=0.6, ls=":", zorder=1)
ax.legend(frameon=False, fontsize=5.0, loc="center left", handlelength=1.1,
          labelspacing=0.25, borderaxespad=0.4)
ax.set_xlabel("true tumour-specific fraction", fontsize=5.8)
ax.set_ylabel("bias vs truth", fontsize=5.8)
lbl(ax, "C", "Every objective is biased")

# ------------------------------------------------------------------ D
ax = axes[1, 0]
d = ppv[(ppv.scenario == "moderate_c80") & (ppv.depth_label == CORE)]
d = d.sort_values("p_novel_true")
ax.plot(d["p_novel_true"], d["patient_sd_mean"], color=RED, lw=1.6, marker="o", ms=2.8,
        label="1 patient", zorder=3)
ax.plot(d["p_novel_true"], d["cohort_sd"], color=BLUE, lw=1.6, marker="s", ms=2.4,
        label="cohort (n=20)", zorder=4)
ax.legend(frameon=False, fontsize=5.0, loc="upper right", handlelength=1.1,
          labelspacing=0.2, borderaxespad=0.25)
ax.set_xlabel("true tumour-specific fraction", fontsize=5.8)
ax.set_ylabel("SD of the estimate", fontsize=5.8)
lbl(ax, "D", "Cohort-level, not per patient")

# ------------------------------------------------------------------ E
ax = axes[1, 1]
for scen, col, mk, name in [("moderate_c80", BLUE, "o", "clean sources"),
                            ("moderate_c80_leak25", RED, "s", "25% leakage")]:
    d = core(scen)
    ax.plot(d["p_novel_true"], d["sens_median"], color=col, lw=1.5, marker=mk, ms=2.4,
            label=name, zorder=3)
ax.set_ylim(0.55, 1.05)
ax.legend(frameon=False, fontsize=5.0, loc="lower right", handlelength=1.1,
          labelspacing=0.2, borderaxespad=0.25)
ax.set_xlabel("true tumour-specific fraction", fontsize=5.8)
ax.set_ylabel("sensitivity of residual call", fontsize=5.8)
lbl(ax, "E", "Robust to reference leakage")

# ------------------------------------------------------------------ F
ax = axes[1, 2]
d = mdc[mdc.scenario == "moderate_c80"].set_index("depth_label").reindex(DEPTHS)
xs = np.arange(len(DEPTHS))
vals = d["mdc_taxa_sens80_spec90"].to_numpy() * 100.0
ax.bar(xs, vals, color=[DCOL[x] for x in DEPTHS], edgecolor="0.3", linewidth=0.4,
       width=0.66, zorder=3)
for x, v in zip(xs, vals):
    ax.text(x, v + 0.2, f"{v:.0f}%", ha="center", fontsize=5.0, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels([x.split("_")[-1] for x in DEPTHS], rotation=42, ha="right",
                   rotation_mode="anchor", fontsize=5.0)
ax.set_xlabel("reads per run", fontsize=5.8)
ax.set_ylabel("min. detectable tumour-\nspecific fraction (%)", fontsize=5.8)
ax.set_ylim(0, vals.max() * 1.32)
lbl(ax, "F", "Depth sets the detection floor")

# ------------------------------------------------------------------ G
ax = axes[2, 0]
SC = [("distinct_c20", "distinct", "#1b7837"),
      ("moderate_c80", "moderate", BLUE),
      ("moderate_c80_leak25", "leaky", ORANGE),
      ("near_identical_c400", "near-identical", RED)]
for scen, name, col in SC:
    d = core(scen)
    ax.plot(d["p_novel_true"], d["excess_median"], color=col, lw=1.5, marker="o", ms=2.2,
            label=name, zorder=3)
ax.axhline(0, color="0.75", lw=0.6, ls=":", zorder=1)
ax.legend(frameon=False, fontsize=4.8, loc="upper left", handlelength=1.1,
          labelspacing=0.16, borderaxespad=0.2)
ax.set_xlabel("true tumour-specific fraction", fontsize=5.8)
ax.set_ylabel("excess over the null", fontsize=5.8)
lbl(ax, "G", "Depends on source distinctness")

# ------------------------------------------------------------------ H
ax = axes[2, 1]
d = core("moderate_c80")
ax.plot(d["p_novel_true"], d["sens_median"], color=RED, lw=1.6, marker="o", ms=2.6,
        label="sensitivity", zorder=3)
ax.plot(d["p_novel_true"], d["spec_median"], color=BLUE, lw=1.6, marker="s", ms=2.4,
        label="specificity", zorder=4)
ax.set_ylim(0.55, 1.05)
ax.legend(frameon=False, fontsize=5.0, loc="lower right", handlelength=1.1,
          labelspacing=0.2, borderaxespad=0.25)
ax.set_xlabel("true tumour-specific fraction", fontsize=5.8)
ax.set_ylabel("call accuracy", fontsize=5.8)
lbl(ax, "H", "Which taxa, not how much")

# ------------------------------------------------------------------ I
ax = axes[2, 2]
d = mdc[mdc.scenario == "moderate_c80"].set_index("depth_label").reindex(DEPTHS)
vals = d["scorable_novel_at_p_novel_0.10"].to_numpy()
ax.plot(xs, vals, color=BLUE, lw=1.5, marker="o", ms=3.2, zorder=3)
for x, v in zip(xs, vals):
    ax.text(x, v + 0.28, f"{v:.0f}", ha="center", fontsize=5.0, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels([x.split("_")[-1] for x in DEPTHS], rotation=42, ha="right",
                   rotation_mode="anchor", fontsize=5.0)
ax.set_ylim(vals.min() - 1.1, vals.max() + 1.5)
ax.set_xlabel("reads per run", fontsize=5.8)
ax.set_ylabel("genera detected at 10%\ntumour-specific fraction", fontsize=5.8)
lbl(ax, "I", "Depth, not separation")

out = FIG / "manuscript-v6" / "main" / "Fig2-identifiability.png"
fig.savefig(out, dpi=330, facecolor="white")
print(f"wrote {out.relative_to(ROOT)}")
try:
    panel_crops(fig)
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
