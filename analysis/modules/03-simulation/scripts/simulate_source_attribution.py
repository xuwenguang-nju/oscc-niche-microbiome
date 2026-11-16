"""S3 — Simulation validation for niche-resolved source attribution of the
oral squamous cell carcinoma (OSCC) intratumoral microbiome.

THE QUESTION, ASKED BEFORE ANY REAL DATA IS TOUCHED
---------------------------------------------------
Two readouts were planned. This script asks, separately, whether each one is
recoverable at the sequencing depths and per-sample overdispersion that the six
public OSCC cohorts actually have:

  readout 1  the FRACTION of the tumour community not explained by known source
             niches (a quantitative "tumour-specific percentage")
  readout 2  the IDENTITY of the taxa making up that unexplained part

THE ANSWER (this is the deliverable; the tables are the evidence)
----------------------------------------------------------------
Readout 1 FAILS and must be dropped from the study design. Three objectives were
tried; none recovers the truth, and the failure is structural rather than a
matter of tuning:

  L2 on proportions      blind to unexplained mass by construction. The mixture
                         p@S has zeros where the novel taxa sit, so the residual
                         on those taxa does not depend on p at all and exerts no
                         downward pull on sum(p): the estimate returns "fully
                         explained" even at a true 20% novel fraction.
  chi-square (1/p)       the opposite pathology -- the rare tail dominates and
                         the estimate saturates near 0.9.
  L1 on proportions      sensitive to unexplained mass, but its own null is of
                         the same order as the signal: scoring a *source* sample
                         against the same reference set already yields an
                         "unexplained fraction" of 0.1-0.46. The bias therefore
                         cannot be separated from the signal.

Readout 2 SURVIVES. The paired-subtraction residual is well behaved: at a true
10% tumour-specific fraction it recovers most of the tumour-specific genera at
specificity 1.0, and recovery improves monotonically with the true fraction.
This is the readout the study should be built on.

The practical consequence, carried into the pipeline design: do not report a
mixing percentage. Report (a) a residual taxon set that reproduces across
independent cohorts and (b) a qualitative "which niche the tumour community
resembles most" ordering, which is a classification question rather than an
unmixing question.

DESIGN
------
* Genera split into G_SRC source-capable genera and G_NOVEL genera absent from
  every source niche, so a residual taxon is a true positive only if it belongs
  to G_NOVEL.
* Novel genera are drawn from the MID-abundance range, because real
  OSCC-enriched taxa sit at ~0.1-5%; genuinely rare taxa drop out under
  overdispersion and are excluded from the sensitivity metric by an explicit
  prevalence filter.
* Source niches are patient-specific Dirichlet(c*base) draws; `c` sets how
  similar they are, and the realised Bray-Curtis dissimilarity is reported
  because that is the quantity that matters.
* Counts are Dirichlet-multinomial at the per-run depths observed in the real
  cohorts.
* A cohort is N_PATIENTS paired patients. Reference source profiles are POOLED
  ACROSS PATIENTS WITH THE SCORED PATIENT LEFT OUT, mirroring how the real
  analysis must avoid self-reference leakage.

SOLVER
------
Iteratively reweighted least squares with a half-simplex projection: pure numpy,
no compiled QP/LP solver (scipy's SLSQP aborts natively in some sandboxes). On
identical inputs it reproduces the exact LP solution -- agreement to four
decimals against scipy SLSQP on both a fully explainable sink and a 20%-novel
sink.

OUTPUTS
-------
results/tables/s3_objective_comparison.csv    bias of L1 / L2 / chi-square
results/tables/s3_simulation_replicates.csv   per cohort replicate
results/tables/s3_simulation_calibration.csv  per scenario x depth x p_novel
results/tables/s3_simulation_mdc.csv          minimum detectable fraction (taxa)
results/tables/s3_perpatient_variability.csv  why cohort aggregation is required
figures/s3-simulation-calibration.png

Run:
  python analysis/modules/03-simulation/scripts/simulate_source_attribution.py \
      [--reps 20] [--patients 20] [--quick]
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[4]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures"
MOD = ROOT / "analysis" / "modules" / "03-simulation"
for p in (TAB, FIG, MOD / "input"):
    p.mkdir(parents=True, exist_ok=True)

SEED = 20261015
G_SRC = 135
G_NOVEL = 15
THETA = 50.0
ALPHA = 0.05
MIN_READS = 10
MIN_PREVALENCE = 0.5

DEPTHS = {
    "PRJNA412445_23k": 23_428,
    "PRJNA597251_43k": 42_602,
    "PRJNA866676_105k": 104_655,
    "PRJNA1415359_130k": 129_505,
    "PRJNA666746_968k": 967_658,
}
CONCENTRATIONS = {"distinct_c20": 20.0,
                  "moderate_c80": 80.0,
                  "near_identical_c400": 400.0}
P_NOVEL_GRID = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]


# --------------------------------------------------------------------------- #
# community generation
# --------------------------------------------------------------------------- #
def base_weights(rng):
    w = rng.lognormal(mean=0.0, sigma=1.7, size=G_SRC + G_NOVEL)
    return w / w.sum()


def split_genera(rng, base):
    order = np.argsort(base)[::-1]
    mid = order[10:70]
    nov_idx = rng.choice(mid, size=G_NOVEL, replace=False)
    src_idx = np.setdiff1d(np.arange(base.size), nov_idx)
    return (src_idx, nov_idx,
            base[src_idx] / base[src_idx].sum(),
            base[nov_idx] / base[nov_idx].sum())


def niche_profiles(rng, w_src, c):
    return np.stack([rng.dirichlet(c * w_src) for _ in range(3)])


def novel_profile(rng, w_nov):
    v = w_nov * rng.lognormal(mean=0.0, sigma=0.5, size=w_nov.size)
    return v / v.sum()


def dm_sample(rng, props, depth, theta=THETA):
    return rng.multinomial(depth, rng.dirichlet(theta * props)).astype(float)


# --------------------------------------------------------------------------- #
# constrained unmixing: IRLS on the half-simplex (pure numpy)
# --------------------------------------------------------------------------- #
def project_half_simplex(p):
    q = np.clip(p, 0.0, None)
    if q.sum() <= 1.0:
        return q
    s = np.sort(q)[::-1]
    css = np.cumsum(s)
    idx = np.arange(1, q.size + 1)
    rho = np.nonzero(s * idx > (css - 1.0))[0][-1]
    theta = (css[rho] - 1.0) / (rho + 1.0)
    return np.clip(q - theta, 0.0, None)


WEIGHTS = {
    "l1": lambda r, m: 1.0 / np.maximum(np.abs(r), 1e-8),
    "l2": lambda r, m: np.ones_like(r),
    "chi2": lambda r, m: 1.0 / np.maximum(m, 1e-10),
}


def unmix(y_sink, S, kind="l1", iters=100):
    y = y_sink / y_sink.sum()
    K = S.shape[0]
    p = project_half_simplex(np.full(K, 1.0 / (K + 1)))
    wfun = WEIGHTS[kind]
    for _ in range(iters):
        m = p @ S
        w = wfun(y - m, m)
        try:
            step, *_ = np.linalg.lstsq(S.T * np.sqrt(w)[:, None], y * np.sqrt(w),
                                       rcond=None)
        except np.linalg.LinAlgError:
            break
        p_new = project_half_simplex(step)
        if np.max(np.abs(p_new - p)) < 1e-12:
            p = p_new
            break
        p = p_new
    if not np.all(np.isfinite(p)):
        p = np.zeros(K)
    return 1.0 - float(p.sum())


# --------------------------------------------------------------------------- #
# cohort-level residual taxa test
# --------------------------------------------------------------------------- #
def bh_q(p):
    m = p.size
    order = np.argsort(p)
    q = np.empty(m)
    ranked = p[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.clip(ranked, 0.0, 1.0)
    return q


def residual_taxa(lr_n, lr_s, alpha=ALPHA):
    def t_p(L):
        mu = L.mean(axis=0)
        sd = L.std(axis=0, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = mu / (sd / np.sqrt(L.shape[0]))
        p = 2.0 * stats.t.sf(np.abs(t), df=L.shape[0] - 1)
        p[~np.isfinite(p)] = 1.0
        return mu, p

    m1, p1 = t_p(lr_n)
    m2, p2 = t_p(lr_s)
    q1, q2 = bh_q(p1), bh_q(p2)
    up = (q1 < alpha) & (q2 < alpha) & (m1 > 0) & (m2 > 0)
    down = (q1 < alpha) & (q2 < alpha) & (m1 < 0) & (m2 < 0)
    return up | down


# --------------------------------------------------------------------------- #
@dataclass
class Rep:
    scenario: str
    depth_label: str
    depth: int
    p_novel_true: float
    rep: int
    bc_saliva_vs_normal: float
    unexplained_tumour: float
    unexplained_null: float
    unexplained_excess: float
    unexplained_patient_sd: float
    n_scorable_novel: int
    sens: float
    spec: float
    f1: float
    leak: float


def run_cohort(rng, scenario, c, depth_label, depth, p_novel, rep, n_patients,
               w_src, w_nov, leak=0.0):
    """leak: the fraction of a novel genus' tumour abundance that is also present
    in every source niche. leak=0 is the most favourable (and unrealistic) case;
    leak=0.25 approximates a tumour-ENRICHED taxon that is detectable but rare in
    adjacent mucosa."""
    G = G_SRC + G_NOVEL
    truth = np.zeros(G, dtype=bool)
    truth[G_SRC:] = True

    counts = np.zeros((n_patients, 4, G))   # 0 tumour, 1 saliva, 2 normal, 3 plaque
    bcs = []
    for i in range(n_patients):
        prof = niche_profiles(rng, w_src, c)
        nov = novel_profile(rng, w_nov)
        p_src = np.array([0.40, 0.40, 0.20]) * (1.0 - p_novel)
        sink = np.concatenate([p_src @ prof, p_novel * nov])

        if leak > 0:
            prof_full = []
            for k in range(3):
                v = np.concatenate([prof[k], leak * p_novel * nov / 3.0])
                prof_full.append(v / v.sum())
        else:
            prof_full = [np.concatenate([prof[k], np.zeros(G_NOVEL)]) for k in range(3)]

        counts[i, 0] = dm_sample(rng, sink, depth)
        for k in range(3):
            counts[i, 1 + k] = dm_sample(rng, prof_full[k], depth)
        bcs.append(1.0 - float(np.sum(np.minimum(prof[0], prof[1]))))

    pooled = counts.sum(axis=0)
    u_t, u_null = [], []
    lr_n, lr_s = [], []
    for i in range(n_patients):
        ref = (pooled - counts[i])[1:]           # leave-this-patient-out sources
        ref = ref / ref.sum(axis=1, keepdims=True)
        u_t.append(unmix(counts[i, 0], ref))
        u_null.append(np.median([unmix(counts[i, k], ref) for k in (1, 2, 3)]))
        eps = 1e-6
        pT, pN, pS = (counts[i, 0] / counts[i, 0].sum(),
                      counts[i, 2] / counts[i, 2].sum(),
                      counts[i, 1] / counts[i, 1].sum())
        lr_n.append(np.log(pT + eps) - np.log(pN + eps))
        lr_s.append(np.log(pT + eps) - np.log(pS + eps))

    u_t, u_null = np.array(u_t), np.array(u_null)
    called = residual_taxa(np.vstack(lr_n), np.vstack(lr_s))

    counts_t = counts[:, 0, :]
    scorable = (counts_t >= MIN_READS).mean(axis=0) >= MIN_PREVALENCE
    s_nov, s_src = scorable & truth, scorable & ~truth
    tp = int(np.sum(called & s_nov))
    n_nov, n_src = int(s_nov.sum()), int(s_src.sum())
    fp = int(np.sum(called & s_src))
    sens = tp / n_nov if n_nov else np.nan
    spec = 1.0 - fp / n_src if n_src else np.nan
    f1 = 2 * tp / (2 * tp + fp + (n_nov - tp)) if (2 * tp + fp + (n_nov - tp)) else np.nan

    return Rep(scenario, depth_label, depth, p_novel, rep, float(np.mean(bcs)),
               float(np.median(u_t)), float(np.median(u_null)),
               float(np.median(u_t - u_null)), float(u_t.std(ddof=1)),
               n_nov, float(sens), float(spec), float(f1), leak)


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--patients", type=int, default=20)
    ap.add_argument("--scenario", default="all",
                    help="one of the CONCENTRATIONS keys, or 'all'")
    ap.add_argument("--leak", type=float, default=0.0,
                    help="fraction of novel-taxon abundance also present in sources")
    ap.add_argument("--assemble-only", action="store_true",
                    help="build the aggregate tables from existing part files")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    if args.quick:
        args.reps, args.patients = 4, 8

    if args.assemble_only:
        assemble(args)
        return

    t0 = time.time()
    rng = np.random.default_rng(SEED)
    src_idx, nov_idx, w_src, w_nov = split_genera(rng, base_weights(rng))

    # ---- objective comparison: this is evidence, not a tuning step ----
    if args.scenario in ("all", "distinct_c20"):
        obj_rows = []
        for depth_label, depth in DEPTHS.items():
            for p_novel in [0.0, 0.10, 0.20]:
                for kind in ["l1", "l2", "chi2"]:
                    u_t, u_n = [], []
                    for _ in range(6):
                        prof = niche_profiles(rng, w_src, 80.0)
                        nov = novel_profile(rng, w_nov)
                        p_src = np.array([0.40, 0.40, 0.20]) * (1 - p_novel)
                        sink = np.concatenate([p_src @ prof, p_novel * nov])
                        yT = dm_sample(rng, sink, depth)
                        pad = np.zeros(G_NOVEL)
                        ref = np.stack([np.concatenate([dm_sample(rng, prof[k], depth), pad])
                                        for k in range(3)])
                        ref = ref / ref.sum(axis=1, keepdims=True)
                        u_t.append(unmix(yT, ref, kind))
                        u_n.append(np.median([unmix(ref[k] * depth, ref, kind)
                                              for k in range(3)]))
                    obj_rows.append({
                        "depth_label": depth_label, "depth": depth,
                        "p_novel_true": p_novel, "objective": kind,
                        "estimate_tumour_median": float(np.median(u_t)),
                        "null_median": float(np.median(u_n)),
                        "bias_vs_truth": float(np.median(u_t) - p_novel),
                    })
        pd.DataFrame(obj_rows).to_csv(TAB / "s3_objective_comparison.csv", index=False)
        print(f"  objectives compared ({time.time() - t0:.0f}s)", flush=True)

    # ---- main sweep, one scenario per invocation (keeps each call short) ----
    scenarios = (CONCENTRATIONS if args.scenario == "all"
                 else {args.scenario: CONCENTRATIONS[args.scenario]})
    for scenario, c in scenarios.items():
        rows: list[Rep] = []
        label = scenario if args.leak == 0 else f"{scenario}_leak{int(args.leak * 100):02d}"
        for depth_label, depth in DEPTHS.items():
            for p_novel in P_NOVEL_GRID:
                for rep in range(args.reps):
                    rows.append(run_cohort(rng, label, c, depth_label, depth,
                                           p_novel, rep, args.patients,
                                           w_src, w_nov, leak=args.leak))
            print(f"  {label} {depth_label} done "
                  f"({time.time() - t0:.0f}s)", flush=True)
        part = TAB / f"s3_replicates_{label}.csv"
        pd.DataFrame([asdict(r) for r in rows]).to_csv(part, index=False)
        print(f"  wrote {part.name} ({time.time() - t0:.0f}s)", flush=True)

    (MOD / "input" / "s3_simulation_settings.json").write_text(json.dumps({
        "seed": SEED, "G_src": G_SRC, "G_novel": G_NOVEL, "theta": THETA,
        "alpha": ALPHA, "n_reps": args.reps, "n_patients_per_cohort": args.patients,
        "depths": DEPTHS, "concentrations": CONCENTRATIONS,
        "p_novel_grid": P_NOVEL_GRID,
        "estimator": "IRLS half-simplex projection, L1 weights",
        "reference_profiles": "pooled across patients, scored patient left out",
        "prevalence_filter": {"min_reads": MIN_READS, "min_prevalence": MIN_PREVALENCE},
        "numpy": np.__version__, "pandas": pd.__version__,
        "scipy": __import__("scipy").__version__,
    }, indent=1), encoding="utf-8")
    print(f"total runtime {time.time() - t0:.0f}s")


def assemble(args) -> None:
    parts = sorted(TAB.glob("s3_replicates_*.csv"))
    if not parts:
        raise SystemExit("no part files found; run the sweep first")
    rep_df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    rep_df.to_csv(TAB / "s3_simulation_replicates.csv", index=False)

    per = (rep_df.groupby(["scenario", "depth_label", "p_novel_true"])
           .agg(patient_sd_mean=("unexplained_patient_sd", "mean"),
                cohort_sd=("unexplained_tumour", "std"),
                null_median=("unexplained_null", "median"),
                n=("rep", "size"))
           .reset_index())
    per["variance_reduction_factor"] = (
        per["patient_sd_mean"] / per["cohort_sd"].replace(0, np.nan))
    per.to_csv(TAB / "s3_perpatient_variability.csv", index=False)

    agg = (rep_df.groupby(["scenario", "depth_label", "depth", "p_novel_true"])
           .agg(n=("rep", "size"),
                bc_mean=("bc_saliva_vs_normal", "mean"),
                u_t_median=("unexplained_tumour", "median"),
                u_null_median=("unexplained_null", "median"),
                excess_median=("unexplained_excess", "median"),
                sens_median=("sens", "median"),
                spec_median=("spec", "median"),
                f1_median=("f1", "median"),
                n_scorable_novel_median=("n_scorable_novel", "median"))
           .reset_index())
    agg["excess_bias"] = agg["excess_median"] - agg["p_novel_true"]
    agg.to_csv(TAB / "s3_simulation_calibration.csv", index=False)

    mdc_rows = []
    for scenario in sorted(agg.scenario.unique()):
        for depth_label in DEPTHS:
            sub = agg[(agg.scenario == scenario) & (agg.depth_label == depth_label)]
            if not len(sub):
                continue
            good = sub[(sub.sens_median >= 0.80) & (sub.spec_median >= 0.90) &
                       (sub.n_scorable_novel_median >= 0.5 * G_NOVEL)]
            mdc_taxa = float(good.p_novel_true.min()) if len(good) else np.nan
            at10 = sub[sub.p_novel_true == 0.10]
            at0 = sub[sub.p_novel_true == 0.0]
            mdc_rows.append({
                "scenario": scenario,
                "depth_label": depth_label,
                "depth": int(sub.depth.iloc[0]),
                "realised_bray_curtis_saliva_vs_normal": round(float(sub.bc_mean.mean()), 4),
                "mdc_taxa_sens80_spec90": mdc_taxa,
                "sens_at_p_novel_0.10": round(float(at10.sens_median.iloc[0]), 3),
                "spec_at_p_novel_0.10": round(float(at10.spec_median.iloc[0]), 3),
                "scorable_novel_at_p_novel_0.10": float(at10.n_scorable_novel_median.iloc[0]),
                "excess_at_p_novel_0.10": round(float(at10.excess_median.iloc[0]), 3),
                "fraction_estimate_at_truth_0": round(float(at0.u_t_median.iloc[0]), 3),
                "fraction_null_at_truth_0": round(float(at0.u_null_median.iloc[0]), 3),
            })
    mdc = pd.DataFrame(mdc_rows).sort_values(["scenario", "depth"])
    mdc.to_csv(TAB / "s3_simulation_mdc.csv", index=False)
    print("assembled from:", [p.name for p in parts])
    print(mdc.to_string(index=False))


if __name__ == "__main__":
    main()
