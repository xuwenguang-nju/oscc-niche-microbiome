#!/usr/bin/env python3
"""Draw Fig. 9: an independent primary OSCC cohort validates the framework.

Supersedes F9, F10 and F11 (all ten panels kept) in one figure, following the
decision of 2026-10-28 to present the primary cohort as a single, visible
validation rather than as three small satellite figures.

Row 1 (A-C) The shift is real: paired ordination, per-patient ratio, F. nucleatum
Row 2 (D-G) The shift survives a 50-fold depth reduction
Row 3 (H-J) The shift survives the choice of taxonomic rank and replicates the
            public cross-cohort meta-analysis at whole-profile level

Cohort : FW2024-420, 20 paired OSCC tumour vs DISTANT normal mucosa (author-confirmed).
         No adjacent-mucosa sample exists, so this figure cannot address localisation.
Input  : results/tables/fw2024_420/*.csv, results/derived/fw2024_420/genus_counts.csv,
         data/raw/FW2024-420/results/05.beta_diversity.zip (PCoA)
Output : figures/main/F9-primary-cohort-validation.png   (176 mm wide)
         figures/submission_tiff/F9-primary-cohort-validation.tif (600 dpi, RGB, LZW)

Run: python analysis/scripts/fw2024_fig9_validation.py
"""
from __future__ import annotations

import importlib.util
import re
import zipfile
from io import StringIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]

_spec = importlib.util.spec_from_file_location(
    "down", ROOT / "analysis/scripts/fw2024_downsampling.py")
down = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(down)

_spec_sd = importlib.util.spec_from_file_location(
    "sd", ROOT / "analysis/scripts/score_definitions.py")
sd = importlib.util.module_from_spec(_spec_sd); _spec_sd.loader.exec_module(sd)
TAB = ROOT / "results/tables/fw2024_420"
DER = ROOT / "results/derived/fw2024_420"
RZ = ROOT / "data/raw/FW2024-420/results"
OUT_PNG = ROOT / "figures/manuscript-v6/main/Fig5-hospital-paired-cohort.png"
OUT_TIF = ROOT / "figures/manuscript-v6/main/Fig5-hospital-paired-cohort.tif"

TARGET_MM = 176.0
FOC_T, FOC_N, MUT = "C3", "C0", "#b9b9b9"
RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

pid = lambda s: re.sub(r"^[NTY]", "", s)
grp = lambda s: "T" if s[0] in "TY" else "N"


def rd_zip(zpath, member):
    with zipfile.ZipFile(RZ / zpath) as zf:
        b = zf.read(member)
    for enc in ("utf-8-sig", "utf-8", "gbk"):
        try:
            return pd.read_csv(StringIO(b.decode(enc)), sep="\t")
        except Exception:
            continue
    raise RuntimeError(member)


def load():
    pcoa = rd_zip("05.beta_diversity.zip",
                  "05.beta_diversity/beta_PCoA/T-vs-N/T-vs-N.braycurtis.PCoA.xls")
    pcoa["pid"] = pcoa["sample"].map(pid)
    pcoa["grp"] = pcoa["sample"].map(grp)
    rat = pd.read_csv(TAB / "consortium_ratio_per_patient.csv",
                      encoding="utf-8-sig", index_col=0)
    fn = pd.read_csv(TAB / "fusobacterium_species_per_sample.csv", encoding="utf-8-sig")
    T = pd.read_csv(TAB / "downsampling_for_manuscript.csv",
                    encoding="utf-8-sig").sort_values("depth")
    counts = pd.read_csv(DER / "genus_counts.csv", index_col=0, encoding="utf-8-sig")
    rel = counts.div(counts.sum(axis=0), axis=1)
    PG = pd.read_csv(TAB / "downsampling_per_genus.csv", encoding="utf-8-sig")
    med = PG.groupby(["genus", "depth"])["p"].median().unstack()
    sig = med < 0.05
    info = pd.DataFrame({"prev": (counts > 0).mean(axis=1) * 100,
                         "ab": rel.mean(axis=1) * 100})
    info["retained"] = sig[20000]
    info["retained_500"] = sig[500]
    # composite score from the published taxonomy sets (29 vs 16), computed independently of any
    # differential test in this cohort. The arms are imported: this file used to carry its own
    # 26 + 16 copy, and that copy is how the figure came to use a set Table S3 does not publish.
    _c = counts.div(counts.sum(axis=0), axis=1)
    _ana = [g for g in sd.MAIN_ANA if g in _c.index]
    _aer = [g for g in sd.MAIN_AER if g in _c.index]
    SCORE = pd.DataFrame({"ana": _c.reindex(_ana).sum(),
                           "aer": _c.reindex(_aer).sum()})
    SCORE["lr"] = np.log10((SCORE.ana + 1e-9) / (SCORE.aer + 1e-9))
    SCORE["grp"] = [grp(s_) for s_ in SCORE.index]
    SCORE["pid"] = [pid(s_) for s_ in SCORE.index]
    rt = pd.read_csv(TAB / "rank_level_tests.csv", encoding="utf-8-sig")
    wp = pd.read_csv(TAB / "whole_profile_vs_public.csv", encoding="utf-8-sig")
    return pcoa, rat, fn, T, info, rt, wp, SCORE


RED_F = "#b2182b"
BLUE_F = "#2166ac"
GREY_F = "#9a9a9a"
ORANGE_F = "#e08214"

# canonical bacterial pathway names, short enough to print on a small panel
PATH_NAMES = {
    "ko03010": "Ribosome",
    "ko03430": "Mismatch repair",
    "ko03440": "Homologous recomb.",
    "ko00970": "Aminoacyl-tRNA",
}
# the four the concurrent meta-analysis (Na et al. 2026) also reports
PATH_SHARED = {"ko03010", "ko03430", "ko03440", "ko00970"}


def draw():
    pcoa, rat, fn, T, info, rt, wp, SCORE = load()

    fig = plt.figure(figsize=(7.09, 9.95))
    # Each row gets its own wspace: rows with long y-labels need more room.
    outer = fig.add_gridspec(5, 1, hspace=0.86, height_ratios=[1, 1, 1, 1, 0.84],
                             left=0.134, right=0.985, top=0.950, bottom=0.050)
    r1 = outer[0].subgridspec(1, 3, wspace=0.55)
    r2 = outer[1].subgridspec(1, 4, wspace=1.05)
    r3 = outer[2].subgridspec(1, 12, wspace=0.62)
    r4 = outer[3].subgridspec(1, 12, wspace=0.30)
    r5 = outer[4].subgridspec(1, 2, wspace=0.36)
    axes = [
        fig.add_subplot(r1[0, 0]), fig.add_subplot(r1[0, 1]), fig.add_subplot(r1[0, 2]),
        fig.add_subplot(r2[0, 0]), fig.add_subplot(r2[0, 1]),
        fig.add_subplot(r2[0, 2]), fig.add_subplot(r2[0, 3]),
        fig.add_subplot(r3[0, 0:3]), fig.add_subplot(r3[0, 4:9]),
        fig.add_subplot(r3[0, 9:12]),
        fig.add_subplot(r4[0, 0:7]), fig.add_subplot(r4[0, 7:12]),
        fig.add_subplot(r5[0, 0]), fig.add_subplot(r5[0, 1]),
    ]
    x = T.depth.values

    # ---------- A: paired PCoA, arrows normal -> tumour
    ax = axes[0]
    px = pcoa.pivot_table(index="pid", columns="grp", values="PCo1")
    py = pcoa.pivot_table(index="pid", columns="grp", values="PCo2")
    ax.scatter(px["N"], py["N"], s=15, color=FOC_N, edgecolor="white", lw=0.4,
               label="distant normal", zorder=3)
    for p_ in px.index:
        ax.annotate("", xy=(px.loc[p_, "T"], py.loc[p_, "T"]),
                    xytext=(px.loc[p_, "N"], py.loc[p_, "N"]),
                    arrowprops=dict(arrowstyle="-|>", color=FOC_T, lw=0.55,
                                    shrinkA=1.5, shrinkB=1.5, alpha=.9), zorder=2)
    ax.scatter(px["T"], py["T"], s=15, color=FOC_T, edgecolor="white", lw=0.4,
               label="tumour", zorder=4)
    ax.axvline(0, color="#dddddd", lw=0.7, zorder=0)
    ax.set_xlabel("PCoA1 (Bray–Curtis)")
    ax.set_ylabel("PCoA2")
    ax.set_title("Normal $\\rightarrow$ tumour", loc="left")
    ax.legend(frameon=False, fontsize=5.6, loc="lower right")

    # ---------- B: per-patient anaerobe:aerobe ratio
    # Panel B uses the taxonomy-defined score (strict), not the original set, which was
    # selected on this cohort's own differential result and would make the AUC circular.
    ax = axes[1]
    pr = SCORE.pivot_table(index="pid", columns="grp", values="lr")
    for p_ in pr.index:
        ax.plot([0, 1], [pr.loc[p_, "N"], pr.loc[p_, "T"]], color="#c9c9c9", lw=0.7, zorder=1)
    ax.scatter(np.zeros(len(pr)), pr["N"], s=15, color=FOC_N, edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(np.ones(len(pr)), pr["T"], s=15, color=FOC_T, edgecolor="white", lw=0.4, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["distant\nnormal", "tumour"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("log$_{10}$ ratio")
    from sklearn.metrics import roc_auc_score as _auc
    _y = np.r_[np.zeros(len(pr)), np.ones(len(pr))]
    _s = np.r_[pr["N"].values, pr["T"].values]
    # two lines: at 15/20 the single-line form ran into panel C's title
    ax.set_title("%d/20 up, AUC %.2f" % (int((pr["T"] > pr["N"]).sum()),
                                         _auc(_y, _s)), loc="left", fontsize=6.4)

    # ---------- C: F. nucleatum, species level
    ax = axes[2]
    fv = fn.pivot_table(index="pid", columns="grp", values="fn")
    for p_ in fv.index:
        ax.plot([0, 1], [np.log10(fv.loc[p_, "N"] + 0.01), np.log10(fv.loc[p_, "T"] + 0.01)],
                color="#c9c9c9", lw=0.7, zorder=1)
    ax.scatter(np.zeros(len(fv)), np.log10(fv["N"] + 0.01), s=15, color=FOC_N,
               edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(np.ones(len(fv)), np.log10(fv["T"] + 0.01), s=15, color=FOC_T,
               edgecolor="white", lw=0.4, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["distant\nnormal", "tumour"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("log$_{10}$ % of reads")
    # d is computed, not hard-coded: the tumour/normal label correction for P16 moved it
    _dv = fv["T"] - fv["N"]
    _d = float(_dv.mean() / _dv.std(ddof=1)) if _dv.std(ddof=1) > 0 else float("nan")
    ax.set_title("$\\it{F.\\ nucleatum}$ %d/20, $d$ = %.2f"
                 % (int((fv["T"] > fv["N"]).sum()), _d), loc="left", fontsize=6.4)

    # ---------- D: genera detected with depth
    ax = axes[3]
    ax.plot(x, T.genera_detected, "-o", color=MUT, ms=3.0, lw=1.2, label="genera detected")
    ax.plot(x, T.nominal_p05, "-s", color=FOC_T, ms=3.0, lw=1.2, label="nominal p < 0.05")
    ax.plot(x, T.fdr05, "-^", color="#7a0177", ms=3.0, lw=1.2, label="FDR < 0.05")
    ax.set_xscale("log")
    ax.set_xlabel("reads per sample")
    ax.set_ylabel("number of genera")
    ax.set_title("Discovery vs depth", loc="left")
    ax.legend(frameon=False, fontsize=5.4, loc="center left")

    # ---------- E: key signals are depth-robust
    ax = axes[4]
    ax.plot(x, T.parv_p, "-o", color=FOC_N, ms=3.0, lw=1.2, label="$\\it{Parvimonas}$")
    ax.plot(x, T.score_p, "-s", color=FOC_T, ms=3.0, lw=1.2,
            label="score")
    ax.axhline(0.05, color="#444", lw=0.8, ls=":")
    ax.annotate("p = 0.05", xy=(x.max(), 0.05), xytext=(-2, 5), textcoords="offset points",
                ha="right", va="bottom", fontsize=5.0, color="#444")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(6e-4, 6e-1)
    ax.set_xlabel("reads per sample")
    ax.set_ylabel("paired Wilcoxon $p$")
    ax.set_title("Signals robust", loc="left")
    ax.legend(frameon=False, fontsize=5.0, loc="upper right", bbox_to_anchor=(0.99, 0.99),
              handlelength=1.2, labelspacing=0.3, borderaxespad=0.2)

    # ---------- F: AUC of the same score panel B reports, and of Parvimonas alone
    ax = axes[5]
    # The two series sit within 0.01 of each other, so an identical solid line with filled markers
    # let the second hide the first. They are drawn as dashed/filled to keep both legible without
    # offsetting either one: the flatness relative to 0.5 is the finding and must not be disturbed.
    ax.plot(x, T.score_auc, ls="--", marker="o", mfc="none", color=FOC_T, ms=3.4, lw=1.2,
            zorder=3)
    # the a-priori up set is the single genus Parvimonas, so its AUC is that of the genus alone
    ax.plot(x, T.prior_up_auc, ls="-", marker="s", color=FOC_N, ms=2.8, lw=1.0, zorder=2)
    ax.axhline(0.5, color="#444", lw=0.8, ls=":")
    # Values carried in the legend rather than as annotations: the two series sit within 0.01 of
    # each other, so free-floating labels overlapped the title and each other. The numbers are read
    # from the plotted series; the previous version hard-coded "0.90 / 0.81" beside lines drawn at
    # 0.86 / 0.78.
    ax.legend(handles=[
        Line2D([], [], color=FOC_T, marker="o", mfc="none", ls="--", ms=3.4, lw=1.2,
               label="score  %.2f" % T.score_auc.iloc[-1]),
        Line2D([], [], color=FOC_N, marker="s", ms=3.0, lw=1.2,
               label="$\\it{Parvimonas}$  %.2f" % T.prior_up_auc.iloc[-1])],
        frameon=False, fontsize=4.8, loc="upper center", handlelength=1.6, labelspacing=0.35,
        borderaxespad=0.1, ncol=1)
    ax.set_xscale("log"); ax.set_ylim(0.45, 1.32)
    ax.set_xlabel("reads per sample")
    ax.set_ylabel("AUC")
    ax.set_title("Separation holds", loc="left")

    # ---------- G: which genera are lost, and why
    ax = axes[6]
    lost = info[info.retained & ~info.retained_500]
    kept = info[info.retained & info.retained_500]
    ax.scatter(kept.prev, kept.ab, s=15, color=FOC_N, alpha=.85,
               label="retained (n=%d)" % len(kept))
    ax.scatter(lost.prev, lost.ab, s=20, color=FOC_T, marker="^",
               label="lost (n=%d)" % len(lost))
    ax.set_yscale("log")
    ax.set_xlabel("prevalence (% of 40)")
    ax.set_ylabel("abundance (%)")
    ax.set_title("Lost = low abundance", loc="left")
    ax.legend(frameon=False, fontsize=5.4, loc="lower right")

    # ---------- H: signal survives at every taxonomic rank
    ax = axes[7]
    cnt = rt.groupby("rank").agg(tested=("taxon", "size"),
                                 fdr=("q", lambda s: (s < 0.05).sum())).reindex(RANKS)
    y = np.arange(len(RANKS))
    ax.hlines(y, cnt.tested, cnt.fdr, color="#d5d5d5", lw=1.0, zorder=1)
    ax.scatter(cnt.tested, y, s=20, facecolor="white", edgecolor="#8a8a8a", lw=0.9,
               zorder=3, label="tested")
    ax.scatter(cnt.fdr, y, s=24, color=FOC_T, edgecolor="white", lw=0.4, zorder=4,
               label="FDR < 0.05")
    for i, f_ in enumerate(cnt.fdr):
        ax.annotate("%d" % int(f_), (f_, i), xytext=(-6, 0), textcoords="offset points",
                    ha="right", va="center", fontsize=5.4, color=FOC_T)
    ax.set_yticks(y); ax.set_yticklabels(RANKS, fontsize=6)
    # The left limit has to sit below 1: phylum passes exactly one test, and with the old limit of
    # 1.5 that dot and its number were clipped outside the axes, so the top row looked like zero
    # FDR hits while the caption (correctly) said 1 of 20.
    ax.set_xscale("log"); ax.set_xlim(0.62, 1200)
    ax.set_xlabel("taxa (log scale)")
    ax.set_title("Every rank", loc="left")
    ax.legend(frameon=False, fontsize=5.2, loc="lower right",
              bbox_to_anchor=(1.02, 0.02))

    # ---------- I: phylum-level change
    ax = axes[8]
    ph = rt[rt["rank"] == "Phylum"].copy()
    ph["delta_pct"] = (ph.mean_T - ph.mean_N) * 100
    PH = ph[(ph.mean_N > 0.005) | (ph.mean_T > 0.005)].sort_values("delta_pct")
    yy = np.arange(len(PH))
    ax.barh(yy, PH.delta_pct, color=[FOC_T if d > 0 else FOC_N for d in PH.delta_pct],
            height=0.52, zorder=2)
    ax.axvline(0, color="#555", lw=0.8, zorder=3)
    ax.set_yticks(yy)
    SHORT = {"Actinobacteriota": "Actinobact.", "Bacteroidota": "Bacteroidota",
             "Campylobacterota": "Campylobact.", "Firmicutes": "Firmicutes",
             "Fusobacteriota": "Fusobact.", "Patescibacteria": "Patescibact.",
             "Proteobacteria": "Proteobact.", "Spirochaetota": "Spirochaet."}
    # names sit on the y axis for every row, so the eye does not jump from one side
    # of the zero line to the other depending on the sign of the bar
    ax.set_yticklabels([SHORT.get(nm, nm.replace("_", " ")) for nm in PH.taxon],
                       fontsize=4.5, color="#333")
    ax.tick_params(axis="y", length=0, pad=1.2)
    for i, (d, q) in enumerate(zip(PH.delta_pct, PH.q)):
        star = "**" if q < 0.01 else ("*" if q < 0.05 else "")
        ax.annotate("%+.1f%s" % (d, star), (d, i), xytext=(3 if d > 0 else -3, 0),
                    textcoords="offset points", ha="left" if d > 0 else "right",
                    va="center", fontsize=5.0, color="#333")
    ax.set_xlabel("$\\Delta$ abundance (points)")
    ax.set_title("Opposite phyla", loc="left")
    # x-limits follow the data: the old (-24, 18) left ~13 points of empty axis
    # to the left of the longest bar
    # right limit leaves room for the value label at the bar tip without pushing
    # it into panel J's y-axis label
    ax.set_xlim(-13.2, 11.0)

    # ---------- J: whole-profile replication against the public meta-analysis
    ax = axes[9]
    hi = (wp.prev_N >= 0.5) & (wp.prev_T >= 0.5)
    # pooled_effect is a natural-log ratio; log2FC is base 2. Plotting them against each other
    # put a factor of ln2 between the axes, so the identity line was drawn between two different
    # units. Convert the pooled effect to log2 so both axes and the identity line agree.
    wp = wp.assign(pooled_log2=wp.pooled_effect / np.log(2))
    ax.scatter(wp.loc[~hi, "pooled_log2"], wp.loc[~hi, "log2FC"], s=7,
               color="#c4c4c4", edgecolor="none", label="low prevalence")
    ax.scatter(wp.loc[hi, "pooled_log2"], wp.loc[hi, "log2FC"], s=9,
               color=FOC_T, edgecolor="white", lw=0.25, label="prevalence $\\geq$ 50%")
    # Two problems met here. The axes were pinned to -2.2..2.2, which clipped 26 of the 120 genera
    # and hid the focus genus itself (Parvimonas, log2FC = 2.9). But simply following the data lets
    # ~24 genera with 0-5% prevalence set the scale: their log2 ratios reach +/-17 because they are
    # driven by one or two reads, not by biology. So clip the view to a range that shows the
    # informative majority and say so in the caption, rather than silently truncating.
    CLIP = 4.0
    xlim = [wp.pooled_log2.min() - 0.15, wp.pooled_log2.max() + 0.15]
    ylim = [-CLIP, CLIP]
    ax.plot([-6, 6], [-6, 6], ls=":", color="#666", lw=0.8)
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    n_out = int((wp.log2FC.abs() > CLIP).sum())
    ax.text(0.03, 0.005, "y clipped at \u00b1%d: %d low-prevalence genera outside"
            % (CLIP, n_out), transform=ax.transAxes, fontsize=4.0, color="0.45", va="bottom")
    rho_all = stats.spearmanr(wp.pooled_log2, wp.log2FC)
    rho_hi = stats.spearmanr(wp.loc[hi, "pooled_log2"], wp.loc[hi, "log2FC"])
    ax.annotate("all $\\rho$ = %.2f\nprevalent $\\rho$ = %.2f" % (rho_all.statistic, rho_hi.statistic),
                xy=(0.03, 0.97), xycoords="axes fraction", ha="left", va="top", fontsize=5.6)
    ax.set_xlabel("public pooled effect (log$_2$ ratio)", labelpad=13)
    ax.set_ylabel("this cohort, log$_2$ ratio")
    ax.set_title("Whole profile replicates", loc="left")
    # Place the legend just under the tick labels (axes y ~ -0.17) and drop the x label below it,
    # so the two never collide in this short panel.
    ax.legend(frameon=False, fontsize=4.8, loc="upper center", bbox_to_anchor=(0.5, -0.155),
              ncol=2, handletextpad=0.3, columnspacing=1.0)

    # ---------- K: the same axis, cohort by cohort (all six cohorts)
    REC = pd.read_csv(ROOT / "results/tables/s5_score_per_cohort.csv", encoding="utf-8-sig")
    REC = REC[REC.definition == "main"].drop(columns="definition")
    REC["normal_type"] = ["peri-tumoural","peri-tumoural","peri-tumoural",
                          "peri-tumoural","contralateral (distant)","distant"]
    REC["lr_mean"] = REC.effect; REC["lr_se"] = REC.se
    REC = REC.iloc[::-1].reset_index(drop=True)
    POOL = pd.read_csv(ROOT / "results/tables/s5_score_pooled.csv",
                       encoding="utf-8-sig")
    POOL = POOL[POOL.definition == "main"].iloc[0]
    ax = axes[10]
    yy = np.arange(len(REC))
    ax.axvline(0, color="#555", lw=0.8, zorder=1)
    ax.errorbar(REC.lr_mean, yy, xerr=1.96 * REC.lr_se, fmt="none", ecolor="#666",
                elinewidth=0.9, capsize=2.0, zorder=3)
    for i, nt in enumerate(REC.normal_type):
        ax.plot(REC.lr_mean.iloc[i], i, "o", ms=4.2,
                color=FOC_T if not str(nt).startswith("peri") else FOC_N, zorder=4)
    ax.set_yticks(yy); ax.set_yticklabels(REC.study, fontsize=5.8)
    hi6 = float((REC.lr_mean + 1.96 * REC.lr_se).max())
    lo6 = float((REC.lr_mean - 1.96 * REC.lr_se).min())
    ax.set_xlim(lo6 - 0.55, hi6 + 0.45); ax.set_ylim(-0.7, len(REC) - 0.3)
    ax.set_xlabel("$\\Delta\\log_{10}$ (anaerobe : aerobe)", fontsize=6.0)
    ax.set_title("Same direction in all six cohorts", loc="left")
    # the note sat on top of the FW2024-420 point; it now sits under the axis, where
    # the x-label leaves clear space
    ax.annotate("score as in panel B (29 vs 16 genera)", xy=(0.02, -0.30),
                xycoords="axes fraction", fontsize=4.4, va="top")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4.0, color=FOC_N,
                              label="peri-tumoural"),
                       Line2D([], [], marker="o", ls="", ms=4.0, color=FOC_T,
                              label="distant")],
              frameon=False, fontsize=4.6, loc="lower left", handletextpad=0.3,
              borderaxespad=0.2)

    # ---------- L: pooled, under all three genus-set definitions
    ax = axes[11]
    PD = pd.read_csv(ROOT / "results/tables/s5_score_pooled.csv",
                     encoding="utf-8-sig").iloc[::-1].reset_index(drop=True)
    yb = np.arange(len(PD))
    ax.axvline(0, color="#555", lw=0.8, ls=":", zorder=1)
    # `hi` ALREADY is the 95% upper bound (pooled + 1.96*SE), so the half-width is
    # (hi - pooled) as it stands. An earlier version wrote xerr = 1.96*(hi-pooled),
    # which multiplied the half-width by 1.96 a second time and drew bars ~2x too
    # wide -- they crossed zero for all four definitions while the tabulated CIs do
    # not. Kept asymmetric so the bounds are exactly the tabulated ones.
    ax.errorbar(PD.pooled, yb, xerr=[PD.pooled - PD.lo, PD.hi - PD.pooled], fmt="none",
                ecolor="#666", elinewidth=0.9, capsize=2.2, zorder=3)
    for i, tag in enumerate(PD.definition):
        ax.plot(PD.pooled.iloc[i], i, "o", ms=4.4,
                color="#888" if tag == "original" else "#222", zorder=4,
                mfc="white" if tag == "original" else "#222")
    ax.set_yticks(yb)
    ax.set_yticklabels([{"original": "original$^{a}$", "full": "full",
                         "main": "main", "eHOMD": "eHOMD"}[t] for t in PD.definition],
                       fontsize=5.6)
    # keep zero visible inside the frame; the widest bar (eHOMD, up to 0.816) must fit
    ax.set_xlim(-0.10, float((PD.hi).max()) + 0.35)
    ax.set_ylim(-0.55, len(PD) - 0.45)
    ax.set_xlabel("pooled $\\Delta\\log_{10}$ (anaerobe : aerobe)", fontsize=6.0)
    ax.set_title("Same direction under all four", loc="left")
    ax.annotate("bars = 95% CI (DerSimonian-Laird, z). $^{a}$sets chosen from this\ncohort's own result (circular; shown for transparency, not used for\ninference). eHOMD = only genera with an eHOMD annotation. All six\ncohorts positive under every definition; eHOMD is the widest interval\nand does not exclude zero under a t(k-2) interval (Table S3).",
                xy=(0.02, -0.34), xycoords="axes fraction", fontsize=4.1,
                va="top", ha="left")

    # ---------- M: predicted function, volcano view
    # PICRUSt2 pathways in the paired design. This row is a bound on the taxonomic result, not an
    # extra line of evidence: panel N shows the whole signal collapses once composition is
    # adjusted for.
    kg = pd.read_csv(ROOT / "results/tables/fw2024_420/picrust2_kegg_pathway_tests.csv",
                     encoding="utf-8-sig")
    kg["neglog10q"] = -np.log10(kg["q"].clip(lower=1e-12))
    qcut = -np.log10(0.05)

    ax = axes[12]
    sig = kg["q"] < 0.05
    ax.axhline(qcut, color="0.45", lw=0.6, ls=(0, (4, 2)), zorder=1)
    ax.axvline(0, color="0.55", lw=0.6, zorder=1)
    ax.scatter(kg.loc[~sig, "delta"], kg.loc[~sig, "neglog10q"], s=4.5, color="#cfcfcf",
               lw=0, zorder=2, label="q ≥ 0.05")
    ax.scatter(kg.loc[sig, "delta"], kg.loc[sig, "neglog10q"], s=6, color=ORANGE_F,
               lw=0, zorder=3, label="q < 0.05")
    lab = kg[kg["pathway"].isin(PATH_SHARED) & sig].copy()
    ax.scatter(lab["delta"], lab["neglog10q"], s=22, facecolor="none", edgecolor=RED_F,
               lw=0.9, zorder=5, label="also in Na et al. 2026")
    PIN = {"ko03440": (0.70, 0.95), "ko03430": (0.70, 0.85),
           "ko00970": (0.70, 0.75), "ko03010": (0.70, 0.65)}
    for _, r in lab.iterrows():
        if r["pathway"] not in PIN:
            continue
        lx, ly = PIN[r["pathway"]]
        ax.annotate(PATH_NAMES[r["pathway"]], xy=(r["delta"], r["neglog10q"]),
                    xycoords="data", xytext=(lx, ly), textcoords="axes fraction",
                    ha="left", va="center", fontsize=4.2, color="0.15", zorder=7,
                    arrowprops=dict(arrowstyle="-", lw=0.35, color="0.55",
                                    shrinkA=1, shrinkB=2.5))
    ax.set_xlabel("pathway abundance,\ntumour − normal", fontsize=5.4)
    ax.set_ylabel(r"$-\log_{10} q$", labelpad=2.0, fontsize=5.6)
    ax.set_title("%d of %d pathways; more replication, less biosynthesis" % (int(sig.sum()), len(kg)),
                 loc="left", fontsize=6.0)
    ax.legend(frameon=False, fontsize=4.0, loc="lower right", handletextpad=0.3,
              borderaxespad=0.3, labelspacing=0.22)
    ax.margins(x=0.52, y=0.10)
    # letter drawn by the loop below

    # ---------- N: the same pathways after adjusting for composition
    ax = axes[13]
    bb = kg[kg["q"] < 0.05].copy()
    ax.axhline(qcut, color="0.45", lw=0.6, ls=(0, (4, 2)), zorder=1)
    ax.axhline(1.0, color=GREY_F, lw=0.5, ls=":", zorder=1)
    for _, r in bb.iterrows():
        ax.plot([0, 1], [r["q"], r["q_after_consortium"]], color="0.82", lw=0.4,
                zorder=2, solid_capstyle="round")
    ax.scatter(np.zeros(len(bb)), bb["q"], s=10, facecolor="none", edgecolor=RED_F,
               lw=0.7, zorder=3, label="before")
    ax.scatter(np.ones(len(bb)), bb["q_after_consortium"], s=11, color=RED_F, lw=0,
               zorder=4, label="after")
    ax.set_yscale("log")
    ax.set_ylim(1.5e-3, 2.6)
    ax.set_xlim(-0.4, 1.4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["q", "q |\nconsortium"], fontsize=4.6)
    ax.set_ylabel("q (log scale)", fontsize=5.6)
    ax.set_title("None survives (0 of %d)" % len(bb), loc="left")
    ax.annotate("q = %.2f–%.2f" % (bb["q_after_consortium"].min(),
                                   bb["q_after_consortium"].max()),
                xy=(0.30, 0.62), xycoords="axes fraction", ha="center", va="center",
                fontsize=4.4, color="0.25", zorder=6)
    ax.legend(frameon=False, fontsize=4.0, loc="lower right", handletextpad=0.3,
              borderaxespad=0.3, labelspacing=0.22)


    # Letters sit a fixed number of inches left of each axes' left edge (figure
    # coordinates). The old per-panel transAxes offsets put the first column of every
    # figure at a different x, so panel letters did not line up down a page.
    from figlib import panel_letter_fixed  # noqa: E402
    for a, L in zip(axes, "ABCDEFGHIJKLMN"):
        panel_letter_fixed(a, L, pad_in=0.30, dy=1.0, case="upper")
        a.margins(0.10)
    for a_ in axes:
        a_.tick_params(labelsize=6)
    return fig


def main():
    # Render ONCE at 600 dpi and derive everything from those pixels. An earlier version saved the
    # PNG at 330 dpi and then re-saved it with a retagged dpi, which cannot change the pixel
    # dimensions: the "deliverable" ended up 188 mm wide instead of 176 mm. Both outputs are now
    # written from the same array, so their physical sizes agree by construction.
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    tmp = ROOT / ".cache/fig5_600.png"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    fig = draw()
    fig.savefig(tmp, dpi=600, bbox_inches="tight")
    with Image.open(tmp) as im:
        arr = im.convert("RGB").copy()
    tmp.unlink(missing_ok=True)
    dpi = round(arr.width / (TARGET_MM / 25.4))
    arr.save(OUT_PNG, format="PNG", dpi=(dpi, dpi))
    arr.save(OUT_TIF, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
    print("PNG %s  %s  %d dpi = %.1f x %.1f mm"
          % (OUT_PNG.relative_to(ROOT), arr.size, dpi,
             arr.width / dpi * 25.4, arr.height / dpi * 25.4))
    print("TIF %s  %.2f MB" % (OUT_TIF.relative_to(ROOT), OUT_TIF.stat().st_size / 1e6))


if __name__ == "__main__":
    main()
