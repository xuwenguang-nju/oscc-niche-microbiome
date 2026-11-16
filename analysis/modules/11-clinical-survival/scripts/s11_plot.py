#!/usr/bin/env python3
"""F6 — clinical and prognostic context in TCGA-HNSC (TCMA microbiome).

Nine panels. Read them as a NEGATIVE result stated properly: the microbiome metrics do not
associate with the standard clinical variables after FDR, they attenuate when stage enters
the model, and none adds predictive value over clinical variables alone.

A  the clinical association screen: -log10 q for 3 metrics x 6 variables (all null)
B  alpha diversity does not track pathological stage
C  Shannon diversity KM (median split), the one nominal log-rank
D  genus richness KM, the other nominal log-rank
E  Parvimonas KM: flat
F  univariable HRs per SD, all near 1
G  the same effect attenuates after stage adjustment
H  likelihood-ratio test: no metric improves the clinical model
I  Harrell's C-index: clinical only vs clinical + metric

Reads:  results/tables/s11_clinical_assoc.csv, s11_survival_univariable.csv,
        s11_survival_multivariable.csv, s11_survival_lrt.csv, s11_survival_km.csv,
        s11_survival_stage_sensitivity.csv, master.csv
Writes: figures/main/F6-clinical-prognosis.png
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

ca = pd.read_csv(TAB / "s11_clinical_assoc.csv")
uni = pd.read_csv(TAB / "s11_survival_univariable.csv")
multi = pd.read_csv(TAB / "s11_survival_multivariable.csv")
lrt = pd.read_csv(TAB / "s11_survival_lrt.csv")
km = pd.read_csv(TAB / "s11_survival_km.csv").set_index("metric")
sens = pd.read_csv(TAB / "s11_survival_stage_sensitivity.csv")
master = pd.read_csv(TAB / "master.csv")   # one row per SAMPLE; panel C-E need patients

BLUE, RED, ORANGE, GREEN = "#2166ac", "#b2182b", "#a6611a", "#1b7837"
NAMED = {"shannon": "Shannon", "richness": "Richness", "parvimonas_rel": "Parvimonas"}
VARLAB = {"stage": "stage", "grade_num": "grade", "site_group": "site",
          "smoker": "smoker", "sex": "sex", "hpv": "HPV"}

fig = plt.figure(figsize=(7.09, 9.30))
# 4 x 3 grid: nine clinical/prognostic panels plus two site-stratification panels
# (previously in supplementary S5). Panel K spans the last two columns of row 3,
# so the grid has no empty cell; axes[3,1] and axes[3,2] are therefore not created.
gs = fig.add_gridspec(4, 3, left=0.134, right=0.985, top=0.975, bottom=0.050,
                      wspace=0.66, hspace=0.72)
axes = np.array([[fig.add_subplot(gs[i, j]) if not (i == 3 and j >= 1) else None
                     for j in range(3)] for i in range(4)])


def lbl(ax, letter, text):
    ax.set_title(text, fontsize=6.2, multialignment="center", pad=5)
    panel_letter_fixed(ax, letter, case="upper")


def to_patient_level(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse samples to one row per patient, mirroring s11_survival.R exactly.

    master.csv is one row per SAMPLE: 177 rows for 155 patients, and 156 rows for the
    134 patients that have Shannon plus survival. Drawing the KM curves from the sample
    rows therefore counts 22 patients twice and treats their two samples as independent
    observations -- pseudo-replication, which makes the log-rank test anti-conservative.
    It mattered here: richness came out at p = 0.021 on 156 samples but p = 0.042 on the
    134 patients, i.e. the sample-level curve doubled the apparent significance of the
    only nominal result in the figure.

    134 is the right unit. Every other statistic in this figure -- HRs, 95% CIs, the
    likelihood-ratio tests, the C-indices -- is already computed by s11_survival.R on
    one row per patient, so the curves have to use the same 134 or the panel and its
    printed p value describe different datasets.

    The aggregation reproduces the R rule verbatim (mean of the metrics, first non-NA
    value of the clinical fields); validated against s11_survival_km.csv, which it
    reproduces to three decimals (67/67 per group, chi-square 0.976 / 4.132 / 0.082).
    """
    metrics = [c for c in ["shannon", "richness", "parvimonas_rel"] if c in df.columns]
    clin = ["os_time", "os_event", "age", "sex", "stage", "grade_num",
            "site_group", "smoker", "hpv"]
    clin = [c for c in clin if c in df.columns]
    agg = df.groupby("case_id")[metrics].mean()
    first = df.groupby("case_id")[clin].agg(
        lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan)
    out = agg.join(first)
    return out[out["os_time"].notna() & (out["os_time"] > 0)]


def km_curve(df, col, cut):
    """Kaplan-Meier step data for a median split (no lifelines dependency)."""
    out = {}
    for name, mask in [("high", df[col] > cut), ("low", df[col] <= cut)]:
        d = df[mask]
        t = np.sort(d["os_time"].dropna().values)
        e = d.loc[d["os_time"].notna()].sort_values("os_time")["os_event"].values
        at_risk, surv, curve = len(t), 1.0, [(0.0, 1.0)]
        for ti, ei in zip(t, e):
            if at_risk > 0:
                surv *= (1 - ei / at_risk)
            curve.append((ti, surv))
            at_risk -= 1
        out[name] = (np.array([c[0] for c in curve]), np.array([c[1] for c in curve]))
    return out


patients = to_patient_level(master)      # 134 patients for the KM panels (see helper)

# ------------------------------------------------------------------ A association screen
ax = axes[0, 0]
M = ca.pivot_table(index="metric", columns="variable", values="q_bh", aggfunc="min")
M = M.reindex(index=["shannon", "richness", "parvimonas_rel"],
              columns=["stage", "grade_num", "site_group", "smoker", "sex", "hpv"])
V = -np.log10(M.values)
im = ax.imshow(V, aspect="auto", cmap="YlOrRd", vmin=0, vmax=max(2.0, np.nanmax(V)))
ax.set_xticks(np.arange(M.shape[1]))
ax.set_xticklabels([VARLAB[c] for c in M.columns], rotation=60, ha="right", fontsize=4.4)
ax.set_yticks(np.arange(M.shape[0]))
ax.set_yticklabels([NAMED[i] for i in M.index], fontsize=4.9)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M.values[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.2,
                    color="0.15")
cb = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.03)
cb.ax.tick_params(labelsize=4.2)
cb.set_label("−log10 q", fontsize=4.5)
ax.set_xlabel("clinical variable", fontsize=5.6)
lbl(ax, "A", "No clinical association")

# ------------------------------------------------------------------ B diversity vs stage
ax = axes[0, 1]
d = master.dropna(subset=["stage", "shannon"])
groups = [d.loc[d.stage == s, "shannon"].dropna().values for s in [1, 2, 3, 4]]
bp = ax.boxplot(groups, positions=np.arange(4), widths=0.5, patch_artist=True,
                medianprops=dict(color="0.15", lw=1.0), boxprops=dict(lw=0.5),
                whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5),
                flierprops=dict(ms=1.4, mfc="0.6", mec="none"))
for patch in bp["boxes"]:
    patch.set_facecolor("#cfe0ee"); patch.set_edgecolor("0.35")
rng = np.random.default_rng(0)
for i, g in enumerate(groups):
    ax.scatter(np.full(len(g), i) + rng.normal(0, 0.05, len(g)), g, s=2.6,
               color=BLUE, alpha=0.55, linewidths=0, zorder=3)
ax.set_xticks(np.arange(4))
ax.set_xticklabels(["I", "II", "III", "IV"], fontsize=5.0)
ax.set_xlabel("pathological stage", fontsize=5.6)
ax.set_ylabel("Shannon diversity", fontsize=5.6)
pv = ca[(ca.metric == "shannon") & (ca.variable == "stage")]["p"].iloc[0]
ax.text(0.03, 0.95, f"Kruskal–Wallis p = {pv:.2f}", transform=ax.transAxes,
        fontsize=4.8, va="top", color="0.3")
lbl(ax, "B", "Diversity is stage-flat")

# ------------------------------------------------------------------ C/D/E KM panels
for ax, col, letter, title in [(axes[0, 2], "shannon", "C", "Shannon"),
                               (axes[1, 0], "richness", "D", "Richness"),
                               (axes[1, 1], "parvimonas_rel", "E", "Parvimonas")]:
    d = patients.dropna(subset=[col, "os_time", "os_event"])
    cut = float(d[col].median())
    cur = km_curve(d, col, cut)
    for (name, (t, s)), c, ls in zip(cur.items(), [RED, BLUE], ["-", "--"]):
        ax.step(t / 30.4375, s, where="post", color=c, lw=1.4, ls=ls,
                label=f"{name} (n={int((d[col] > cut).sum() if name=='high' else (d[col] <= cut).sum())})")
    ax.set_ylim(0, 1.03)
    # follow-up runs to 3,981 days (130.8 months); the old hard limit of 92 months cut the curves
    # while four patients were still at risk, hiding the tail for no stated reason
    ax.set_xlim(0, float((d["os_time"] / 30.4375).max()) * 1.02)
    ax.set_xlabel("months since surgery", fontsize=5.6)
    ax.set_ylabel("overall survival", fontsize=5.6)
    ax.legend(frameon=False, fontsize=4.3, loc="lower left", handlelength=1.2,
              labelspacing=0.18, borderaxespad=0.25)
    ax.text(0.96, 0.93, f"log-rank p = {km.loc[col, 'p_logrank']:.2f}",
            transform=ax.transAxes, fontsize=4.7, ha="right", va="top", color="0.3")
    lbl(ax, letter, f"{title} KM")

# ------------------------------------------------------------------ F univariable HR
ax = axes[1, 2]
ys = np.arange(len(uni))
ax.errorbar(uni["HR_per_SD"], ys, xerr=[uni["HR_per_SD"] - uni["lower95"],
                                        uni["upper95"] - uni["HR_per_SD"]],
            fmt="o", ms=3.6, color=BLUE, ecolor="0.55", elinewidth=0.9, capsize=1.6,
            zorder=3)
ax.axvline(1, color="0.35", lw=0.9, ls=(0, (4, 1.6)), zorder=2)
ax.set_yticks(ys)
ax.set_yticklabels([NAMED[m] for m in uni.metric], fontsize=5.0)
ax.invert_yaxis()
ax.set_xlabel("HR per SD (95% CI)", fontsize=5.6)
for y, (hr, c) in enumerate(zip(uni["HR_per_SD"], uni["concordance"])):
    ax.text(1.32, y, f"C={c:.2f}", va="center", fontsize=4.3, color="0.35")
ax.set_xlim(0.55, 1.45)
lbl(ax, "F", "Univariable effects small")

# ------------------------------------------------------------------ G attenuation
ax = axes[2, 0]
mm = multi[multi.term.str.endswith("_z")]
ys = np.arange(len(mm))
ax.scatter(uni["HR_per_SD"], np.arange(len(uni)), s=26, color=ORANGE, zorder=3,
           label="univariable", edgecolor="0.25", linewidth=0.3)
# the multivariable model holds stage, age and sex, so the label cannot say "stage-adjusted" alone
ax.scatter(mm["HR"], ys, s=26, marker="s", color=BLUE, zorder=4,
           label="adjusted", edgecolor="0.25", linewidth=0.3)
for y in ys:
    ax.plot([uni["HR_per_SD"].iloc[y], mm["HR"].iloc[y]], [y, y], color="0.6", lw=0.8,
            zorder=2)
ax.axvline(1, color="0.35", lw=0.9, ls=(0, (4, 1.6)), zorder=2)
ax.set_yticks(ys)
ax.set_yticklabels([NAMED[m] for m in mm.metric], fontsize=5.0)
ax.set_xlabel("HR per SD", fontsize=5.6)
ax.legend(frameon=False, fontsize=4.3, loc="upper left", handletextpad=0.2,
          labelspacing=0.2, borderaxespad=0.25, ncol=2)
ax.set_xlim(0.55, 1.15)
ax.set_ylim(len(ys) - 0.4, -2.0)
lbl(ax, "G", "Adjustment moves one of three toward 1")

# ------------------------------------------------------------------ H LRT
ax = axes[2, 1]
ys = np.arange(len(lrt))
ax.barh(ys, lrt["LRT_chisq"], color="#bdbdbd", edgecolor="0.3", lw=0.4, height=0.55,
        zorder=3)
ax.axvline(3.841, color=RED, lw=1.0, ls=(0, (4, 1.6)), zorder=4)
ax.set_yticks(ys)
ax.set_yticklabels([NAMED[m] for m in lrt.metric], fontsize=5.0)
ax.invert_yaxis()
ax.set_xlabel("likelihood-ratio χ²", fontsize=5.6)
for y, (c, p) in enumerate(zip(lrt["LRT_chisq"], lrt["LRT_p"])):
    ax.text(c + 0.08, y, f"p={p:.2f}", va="center", fontsize=4.4, color="0.3")
ax.set_xlim(0, lrt["LRT_chisq"].max() * 1.6)
ax.text(0.97, 0.06, "dashed = p 0.05", transform=ax.transAxes, fontsize=4.3,
        ha="right", color="0.4", style="italic")
lbl(ax, "H", "No added information")

# ------------------------------------------------------------------ I C-index
ax = axes[2, 2]
base = float(lrt["loglik_clin"].iloc[0])
xs = np.arange(len(uni))
ax.bar(xs, uni["concordance"], color=BLUE, edgecolor="0.3", lw=0.4, width=0.55, zorder=3)
ax.axhline(0.5, color=RED, lw=1.0, ls=(0, (4, 1.6)), zorder=4)
for x, v in zip(xs, uni["concordance"]):
    ax.text(x, v + 0.008, f"{v:.2f}", ha="center", fontsize=4.8, color="0.25")
ax.set_xticks(xs)
ax.set_xticklabels([NAMED[m] for m in uni.metric], rotation=42, ha="right", fontsize=4.8)
ax.set_ylim(0.45, 0.62)
ax.set_ylabel("Harrell's C-index", fontsize=5.6)
ax.text(0.5, 0.06, "dashed = chance", transform=ax.transAxes, fontsize=4.3,
        ha="center", color="0.4", style="italic")
lbl(ax, "I", "Barely above chance")

# ------------------------------------------------------------------ J/K site stratification
# These two panels answer the obvious reviewer question - "did you check anatomic
# site?" - and used to sit in supplementary S5.
strat = pd.read_csv(TAB / "s6v3_stratified_tests.csv")
pooled_dir = pd.read_csv(TAB / "s6v2_direction_check.csv")

ax = axes[3, 0]
pv = strat[strat["genus"] == "Parvimonas"].dropna(subset=["log2fc"]).copy()
# "oral" pools the six explicit oral subsites; "oral*" is the GDC "other and ill-defined sites
# in lip, oral cavity and pharynx" bucket, whose subsite is not resolvable. The dagger marks the
# second one, and both captions use the same wording so the panel and the text cannot drift.
pv["label"] = pv["site_group"].str.replace("*", "\u2020", regex=False)
pv = pv.sort_values("log2fc")
yy = np.arange(len(pv))
SITE_COL = {"oral": "#1f78b4", "oral*": "#8da0cb", "larynx": "#d95f02", "pharynx": "#9e9e9e"}
cols_j = [SITE_COL.get(g, "#9e9e9e") for g in pv["site_group"]]
ax.barh(yy, pv["log2fc"], color=cols_j, edgecolor="0.3", lw=0.4, height=0.6, zorder=3)
ax.axvline(0, color="0.4", lw=0.9, zorder=4)
ax.set_yticks(yy)
ax.set_yticklabels([f"{l}\n(n={int(t)}T/{int(n)}N)"
                    for l, t, n in zip(pv["label"], pv["n_tumour"], pv["n_normal"])],
                   fontsize=4.4)
for i, (v, pv_) in enumerate(zip(pv["log2fc"], pv["p_paired"])):
    txt = f"{v:+.2f}" + (f" p={pv_:.2f}" if pd.notna(pv_) else "")
    ax.text(v + (0.08 if v >= 0 else -0.08), i, txt, fontsize=4.2, va="center",
            ha="left" if v >= 0 else "right", color="0.25")
ax.set_xlabel("Parvimonas log$_2$ (tumour / normal)", fontsize=5.4)
ax.margins(x=0.42, y=0.10)
ax.tick_params(axis="x", labelsize=4.6, pad=1.5)
lbl(ax, "J", "Site does not rescue it")

ax = fig.add_subplot(gs[3, 1:])   # spans the last two columns; axes[3,2] is not created
tn = strat.drop_duplicates("site_group").set_index("site_group")[
    ["n_tumour", "n_normal", "n_pairs"]].reindex(["oral", "oral*", "larynx", "pharynx"])
xx = np.arange(len(tn))
ax.bar(xx - 0.21, tn["n_tumour"], 0.40, color="#7b3294", zorder=3, label="tumour")
ax.bar(xx + 0.21, tn["n_normal"], 0.40, color="#a6dba0", zorder=3, label="normal")
for i, pr in enumerate(tn["n_pairs"]):
    ax.text(i, max(tn["n_tumour"].iloc[i], tn["n_normal"].iloc[i]) + 2.5,
            f"{int(pr)} pr", ha="center", fontsize=4.4, color="0.3")
ax.set_xticks(xx)
ax.set_xticklabels([s.replace("*", "\u2020") for s in tn.index], fontsize=4.8,
                   rotation=28, ha="right")
ax.set_ylabel("samples", fontsize=5.6)
ax.margins(y=0.36)
ax.legend(frameon=False, fontsize=4.6, loc="upper right", ncol=1,
          handlelength=0.9, labelspacing=0.2, borderaxespad=0.3)
lbl(ax, "K", "Testability by site")

out = FIG / "manuscript-v6" / "supplementary" / "FigS3-clinical-negative.png"
fig.savefig(out, dpi=330, facecolor="white")
print(f"wrote {out.relative_to(ROOT)}")
try:
    panel_crops(fig)
except Exception as e:                                        # noqa: BLE001
    print("panel_crops skipped:", e)
