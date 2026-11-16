#!/usr/bin/env python3
"""S2 - host coupling, anatomic site, tongue subsite and the platform signature.

One figure replacing what used to be three separate supplements (old S4 host-coupling
panels, old S5 site stratification, old S6 tongue subsite). They share a single
question - "does the TCGA-HNSC signal depend on where the sample came from, and on
which assay was used?" - so they belong together rather than in three places.

Layout 2 x 12 (four panels of three columns, then three of four): A-B platform
direction and Parvimonas coupling, C-D oral-only restriction and the Parvimonas sign
flip across groupings, E-F tongue subsite and per-assay coupling, G the rank-agreement
scatter.

FIVE PANELS WERE REMOVED from the original twelve, because each re-drew a panel that is
already in a main figure from the same table and the same statistic: the composition bar
(now Fig. 6L and Fig. 7F - it appeared three times), Parvimonas-by-site (Fig. 6J),
testability-by-site (Fig. 6K), and rank agreement by cohort (Fig. 7A). A sixth, the top-8
rank-gap bar chart, was a nested subset of Fig. 7C's top 12 and added no data. The panels
were re-lettered A-G so the sequence has no gaps.

Note on A: it is NOT a duplicate of Fig. 7E. A uses s6v2_direction_check.csv, whose log2FC
comes from TCMA tumour versus normal abundances in TCGA WXS/WGS (9 genera, 6/9 agreeing
with 16S). Fig. 7E uses s6a_direction_check.csv, an independent total-RNA-seq oral cohort
(GSE227919; 9 testable genera, 5/9 agreeing). Different substrate, different genus set.
The column name "rna_direction" inside s6v2 is a legacy misnomer - see the warning in
s6v2_host_coupling.py. It has previously caused A and Fig. 7E to be swapped; the tables
must be read from this file, not from that column name.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT_SRC = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT_SRC / "14-figure-layout" / "scripts"))
from figlib import panel_letter_fixed, project_root  # noqa: E402

ROOT = project_root(Path(__file__).resolve().parent)
TAB = ROOT / "results" / "tables"

# figure-style helper, vendored beside the analysis code
_FIGSTYLE = next(p for p in Path(__file__).resolve().parents
              if (p / "figure_style" / "runtime.py").is_file()) / "figure_style"
exec(compile(open(_FIGSTYLE / "runtime.py", encoding="utf-8").read(),
              "runtime.py", "exec"))

# --------------------------------------------------------------------- data
# Tables needed by the SEVEN retained panels only. Three tables that used to be read
# here are gone with their panels: s6b_signal_composition.csv (was C), and the two
# s11_genus_concordance files were kept only for L, which still reads the disagreement
# table below.
dir_tab = pd.read_csv(TAB / "s6v2_direction_check.csv")
corr = pd.read_csv(TAB / "s6v4_parvimonas_top.csv")
strat = pd.read_csv(TAB / "s6v3_stratified_tests.csv")
pooled = pd.read_csv(TAB / "s6v2_direction_check.csv")
tongue = pd.read_csv(TAB / "s6v5_subsite_tests.csv")

BLUE, RED, ORANGE, GREEN = "#2166ac", "#b2182b", "#a6611a", "#1b7837"
GREY = "#9e9e9e"


def lbl(ax, letter, text):
    ax.set_title(text, fontsize=6.2, multialignment="center", loc="center", pad=5)
    panel_letter_fixed(ax, letter, case="upper")


def main() -> None:
    apply_figure_style(sizes=(7.0, 6.2, 5.4))
    fig = plt.figure(figsize=(7.09, 5.05))
    # Three rows of 2 + 3 + 2 on a 12-column lattice; letters run A-G strictly left to
    # right, top to bottom. Two rows packed all seven panels into half the height and
    # neighbouring panels collided; every panel now gets at least four columns.
    # Row 1 is deliberately uneven: B carries y tick labels of the form
    # "Methylation27K: ZNF354A", which need ~1.3 in of gutter to B's left, so A and B
    # split the width evenly instead of squeezing A.
    # bottom raised from 0.075: the two-line tick labels in E need the room, or their n is clipped
    gs = fig.add_gridspec(3, 12, left=0.134, right=0.985, top=0.945, bottom=0.115,
                          wspace=1.55, hspace=1.10)
    SPAN = {(0, 0): (0, 6), (0, 1): (6, 12),
            (1, 0): (0, 4), (1, 1): (4, 8), (1, 2): (8, 12),
            (2, 0): (0, 6), (2, 1): (6, 12)}
    axes = {k: fig.add_subplot(gs[k[0], a:b]) for k, (a, b) in SPAN.items()}

    # ---------------------------------------------------------- A direction check
    ax = axes[(0, 0)]
    d = dir_tab[dir_tab["rna_direction"].isin(["up_in_tumour", "down_in_tumour"])].copy()
    d = d.dropna(subset=["log2fc"]).sort_values("log2fc")
    y = np.arange(len(d))
    ax.barh(y, d["log2fc"], height=0.62, zorder=3,
            color=[RED if bool(a) else BLUE for a in d["direction_agrees"]])
    ax.axvline(0, color="0.4", lw=0.9, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(d["genus"], fontsize=4.6)
    ax.set_xlabel("log$_2$ (tumour / normal)", fontsize=5.4)
    # the legend used to sit inside the plotting area at the far right, where B's long
    # "Methylation27K: …" tick labels reached in; it now sits under the axis
    ax.text(0.0, -0.19, "red = same as 16S\nblue = opposite", transform=ax.transAxes,
            fontsize=4.4, color="0.35", va="top", ha="left")
    ax.margins(y=0.05)
    n_agree = int(d["direction_agrees"].sum())
    # Name the substrate in the title. A and B-G ALL use TCMA, whose microbial reads come from
    # TCGA WXS/WGS (PathSeq); A asks the direction question and B-G the coupling question. The
    # only RNA-seq comparison in the paper is Fig. 7E (GSE227919, 5/9). Leaving the substrate
    # implicit let A and Fig. 7E be swapped, so the title now states it.
    lbl(ax, "A", f"WXS/WGS vs 16S: {n_agree}/{len(d)}")

    # ---------------------------------------------------------- B Parvimonas coupling
    ax = axes[(0, 1)]
    top = corr.sort_values("q").head(12).iloc[::-1]
    yy = np.arange(len(top))
    ax.barh(yy, top["R"], height=0.62, zorder=3,
            color=[RED if q < 0.05 else GREY for q in top["q"]])
    ax.axvline(0, color="0.4", lw=0.9, zorder=4)
    ax.set_yticks(yy)
    # labels are shortened to "assay · feature" so they fit the gutter B owns; the full
    # names are in the data table. Previously they reached into A.
    SHORT_A = {"Methylation27K": "Met27K", "RNAseq": "RNA", "RPPA": "RPPA", "miRNA": "miR"}
    ax.set_yticklabels([f"{SHORT_A.get(a, a)} · {f}" for a, f in
                        zip(top["assay"].astype(str), top["feature"].astype(str))],
                       fontsize=3.8)
    ax.tick_params(axis="y", length=0, pad=1.2)
    ax.set_xlabel("Spearman $\\rho$ with Parvimonas", fontsize=5.4)
    ax.margins(x=0.22, y=0.06)
    lbl(ax, "B", "Nothing couples")

    # ---------------------------------------------------------- C three genera, oral vs pooled
    ax = axes[(1, 0)]
    rep = ["Actinomyces", "Lancefieldella", "Veillonella"]
    rows = []
    for g in rep:
        qp = pooled.loc[pooled["genus"] == g, "q_paired"]
        qo = strat.loc[(strat["genus"] == g) & (strat["site_group"] == "oral"), "q_paired"]
        rows.append({"genus": g,
                     "pooled": float(qp.iloc[0]) if len(qp) and pd.notna(qp.iloc[0]) else np.nan,
                     "oral": float(qo.iloc[0]) if len(qo) and pd.notna(qo.iloc[0]) else np.nan})
    B = pd.DataFrame(rows)
    xx = np.arange(len(B))
    w = 0.36
    ax.bar(xx - w / 2, -np.log10(B["pooled"]), w, color=GREY, zorder=3, label="pooled")
    ax.bar(xx + w / 2, -np.log10(B["oral"]), w, color="#1f78b4", zorder=3, label="oral only")
    ax.axhline(-np.log10(0.05), color=RED, lw=1.0, ls=(0, (4, 2)), zorder=4)
    ax.set_xticks(xx)
    ax.set_xticklabels([g[:9] for g in B["genus"]], fontsize=4.2, rotation=28, ha="right")
    ax.set_ylabel("$-\\log_{10}$ q", fontsize=5.4)
    # one horizontal row above the axes: inside the frame the key covered both the tall
    # oral-only bar and the y tick labels
    ax.legend(frameon=False, fontsize=4.2, loc="upper center", ncol=2,
              bbox_to_anchor=(0.5, -0.30), handlelength=0.8, columnspacing=0.9,
              labelspacing=0.15, borderaxespad=0.0)
    ax.margins(y=0.30)
    lbl(ax, "C", "Oral-only sharpens")

    # ---------------------------------------------------------- D Parvimonas across groupings
    ax = axes[(1, 1)]
    # Site groups only. This panel used to also carry the two tongue subsets, which are a
    # different stratification of a different table; mixing them put five bars on one axis and
    # made "by site group" untrue. The tongue subsets belong to panel E, and the pooled
    # comparison (all tumour vs all normal) is already shown in panel A, so neither is repeated.
    rows = []
    for sg in ["oral", "oral*", "larynx"]:
        r = strat[(strat.genus == "Parvimonas") & (strat.site_group == sg)]
        if len(r) and pd.notna(r["log2fc"].iloc[0]):
            rows.append((sg.replace("*", "\u2020"), float(r["log2fc"].iloc[0]),
                         float(r["p_paired"].iloc[0]) if pd.notna(r["p_paired"].iloc[0]) else np.nan,
                         int(r["n_tumour"].iloc[0]), int(r["n_normal"].iloc[0])))
    G = pd.DataFrame(rows, columns=["label", "lfc", "p", "nT", "nN"])
    yy = np.arange(len(G))
    ax.barh(yy, G["lfc"], height=0.6, zorder=3,
            color=[GREEN if v > 0 else ORANGE for v in G["lfc"]])
    ax.axvline(0, color="0.4", lw=0.9, zorder=4)
    ax.set_yticks(yy)
    # n per group, as in Fig. S3J, so the reader can see at once which bars rest on two normals
    ax.set_yticklabels([f"{l}\n({t}T/{n}N)" for l, t, n in zip(G["label"], G["nT"], G["nN"])],
                       fontsize=4.2)
    ax.tick_params(axis="y", length=0, pad=1.2)
    # value labels go at the bar tip on the side away from the zero line, so they never
    # sit on top of a bar
    for i, (v, p_) in enumerate(zip(G["lfc"], G["p"])):
        txt = f"{v:+.2f}" + (f" p={p_:.2f}" if pd.notna(p_) else " n/a")
        ax.text(v + (0.06 if v >= 0 else -0.06), i, txt, fontsize=4.0, va="center",
                ha="left" if v >= 0 else "right", color="0.25")
    ax.set_xlabel("Parvimonas log$_2$ (tumour / normal)", fontsize=5.4)
    ax.margins(x=0.55, y=0.22)
    lbl(ax, "D", "Sign flips across site groups")

    # ---------------------------------------------------------- E tongue subsites
    ax = axes[(1, 2)]
    sub = tongue[tongue.genus == "Parvimonas"].dropna(subset=["log2fc"]).copy()
    xs = np.arange(len(sub))
    bars = ax.bar(xs, sub["log2fc"], width=0.55, color=GREEN, zorder=3)
    ax.axhline(0, color="0.4", lw=0.9, zorder=4)
    ax.set_xticks(xs)
    # the n goes in the tick label: the last bar rests on a single pair, which the reader has to
    # see, and the value label of that bar was being clipped at the axis edge before
    ax.set_xticklabels([f"{s.replace('_', ' ')}\n({int(t)}T/{int(n)}N)"
                        for s, t, n in zip(sub["subset"], sub["n_tumour"], sub["n_normal"])],
                       fontsize=4.0)
    # "Parvimonas log2" alone, because the two-line form was wide enough to run into panel D's
    # tick labels; the axis meaning is the same as D's and is given once in the caption.
    ax.set_ylabel("Parvimonas log$_2$", fontsize=5.4)
    for x, v, p_ in zip(xs, sub["log2fc"], sub["p_paired"]):
        # labels sit on the far side of the baseline from the bar, so a long negative bar cannot
        # push its own text past the axis limit
        txt = f"{v:+.2f}" + (f"\np={p_:.2f}" if pd.notna(p_) else "\nn/a")
        ax.text(x, v + (0.10 if v >= 0 else -0.10), txt, ha="center",
                va="bottom" if v >= 0 else "top", fontsize=3.9, color="0.25")
    ax.margins(y=0.34)
    # "oral others" is not a tongue subsite, so the title cannot say "Tongue subsites" for a panel
    # that also carries it; the subset names themselves carry what each bar is.
    lbl(ax, "E", "Tongue subsets and oral remainder")

    # ---------------------------------------------------------- F min q per assay
    # No Parvimonas pair reaches q < 0.05 in any assay, so counting significant pairs
    # gives an empty panel. Plotting the smallest q per assay shows the same fact and
    # keeps the panel informative.
    ax = axes[(2, 0)]
    # The assay is labelled "Methylation27K" in the table. This list used to read
    # "methylation", which never matched, so the methylation assay was silently
    # dropped -- the very assay that carries 32 of the 34 significant couplings in
    # Fig. 7G. Omitting it weakened the panel's own claim (that nothing couples in ANY
    # assay), so the real label is used and all four assays are plotted.
    _ORDER = ["Methylation27K", "RNAseq", "RPPA", "miRNA"]
    minq = corr.groupby("assay")["q"].min().reindex(
        [a for a in _ORDER if a in set(corr["assay"])]).dropna()
    xs = np.arange(len(minq))
    ax.bar(xs, minq.values, width=0.55, color=ORANGE, zorder=3)
    ax.axhline(0.05, color=RED, lw=1.0, ls=(0, (4, 2)), zorder=4)
    ax.set_xticks(xs)
    ax.set_xticklabels(list(minq.index), fontsize=4.4, rotation=28, ha="right")
    ax.set_ylabel("min q per assay", fontsize=5.2, labelpad=2.0)
    for x, v in zip(xs, minq.values):
        ax.text(x, v + 0.006, f"{v:.2f}", ha="center", fontsize=4.2, color="0.25")
    ax.margins(y=0.24)
    lbl(ax, "F", "Assay-specific coupling")

    # ---------------------------------------------------------- G concordance scatter
    ax = axes[(2, 1)]
    dis = pd.read_csv(TAB / "s11_genus_concordance_disagreement.csv")
    ax.scatter(dis["tcma_rank"], dis["cohort_rank"], s=9, color=BLUE,
               edgecolor="0.3", lw=0.3, zorder=3)
    lim = [0, max(dis["tcma_rank"].max(), dis["cohort_rank"].max()) + 8]
    ax.plot(lim, lim, color="0.55", lw=0.9, ls=(0, (4, 1.6)), zorder=2)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("TCMA genus rank", fontsize=5.4)
    ax.set_ylabel("amplicon genus rank", fontsize=5.4)
    lbl(ax, "G", "Concordance is reproducible")

    # In the v6 set this figure is Fig. S4. The path used to say FigS3, which left a stale
    # orphan PNG beside the real file on every export run.
    out = ROOT / "figures" / "manuscript-v6" / "supplementary" / "FigS4-tcma-host-coupling.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=330, facecolor="white")
    print(f"wrote {out.name} ({out.stat().st_size:,} bytes)")
    try:
        panel_crops(fig)
    except Exception as e:                                       # noqa: BLE001
        print("panel_crops skipped:", e)


if __name__ == "__main__":
    main()
