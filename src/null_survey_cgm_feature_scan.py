"""
Scan the survey features that were null against BOTH glycemic targets for any
association with a wider set of CGM-derived features.

Motivation: src/survey_glycemic_spearman_bootstrap.py tested every survey
feature against lab HbA1c and mean CGM. Both are *central tendency* measures --
they say where glucose sits on average and nothing about how it moves. Sixteen
survey features came back null against both. A feature that tracks glycemic
variability, hypo/hyperglycaemic exposure, or excursion frequency rather than
average level would have been invisible to that design.

This script asks whether any of those 16 survivors correlate with a CGM feature
that captures shape rather than level.

Design notes:

  - Candidate CGM features are wear-normalised. The cgm_<range>_minutes columns
    are excluded because they are frac_time multiplied by wear duration, and
    days_tracked ranges from 2 to 15 here, so they conflate glycemic exposure
    with how long the sensor was on. total_spikes is excluded for the same
    reason in favour of avg_spikes_per_day.

  - Conditional within-range means are mostly excluded. cgm_sev_hypo_mean_mg_dl
    and cgm_sev_hyper_mean_mg_dl are defined only for the 31% / 38% of
    participants who ever entered those ranges, so a correlation there is
    contaminated by selection into the subsample, and sev_hyper is additionally
    censored at the 401 mg/dL sensor ceiling. Only cgm_normal_mean_mg_dl (99.3%
    coverage) and cgm_mod_hyper_mean_mg_dl (89.9%) are retained.

  - Time-below-range and time-above-range are derived here because the consensus
    ADA/ATTD targets are stated on the pooled ranges, not the four separate bins.

  - Participants flagged cgm_quality_suspect are dropped. Per
    src/add_cgm_ranges_to_full_person_aggregates.py these are failed sensors
    (pinned at the floor, or under one day of wear), not physiology.

  - This is an explicitly exploratory scan, so Benjamini-Hochberg runs across the
    ENTIRE survey x CGM grid rather than within a column. Testing 16 features
    against 22 outcomes and correcting per-outcome would let roughly one false
    positive per column through by construction.

  - Two-stage: Spearman + global FDR over the whole grid first, then the same
    10,000-resample BCa bootstrap used by the original script, run only on the
    survivors. The bootstrap is expensive and uninformative for the ~97% of
    cells that are obviously null.

  - Any survivor is re-checked within meta_study_group. Diabetes severity drives
    both glycemia and survey responses, and it already manufactured two spurious
    "discordant" results in the original analysis (dm_selfcare_healthcare,
    dm_selfcare_footcare), so an unstratified hit here is not believable on its own.

Run from repo root:
    .venv/bin/python src/null_survey_cgm_feature_scan.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap, rankdata, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"
HBA1C_RESULTS = ROOT / "outputs" / "survey_glycemia" / "survey_hba1c_spearman_bootstrap.csv"
CGM_RESULTS = ROOT / "outputs" / "survey_glycemia" / "survey_mean_cgm_spearman_bootstrap.csv"
OUT_PATH = ROOT / "outputs" / "survey_glycemia" / "null_survey_cgm_feature_scan.csv"

N_BOOT = 10_000
SEED = 42
CI_LEVEL = 0.95
ALPHA = 0.05

SEVERITY_ORDER = {
    "healthy": 0,
    "pre_diabetes_lifestyle_controlled": 1,
    "oral_medication_and_or_non_insulin_injectable_medication_controlled": 2,
    "insulin_dependent": 3,
}

# label -> (column, family). Derived columns are built in load_data().
CGM_FEATURES: list[tuple[str, str]] = [
    # Consensus time-in-range family: where glucose sits, as a share of wear.
    ("cgm_sev_hypo_frac_time", "time in range"),
    ("cgm_mod_hypo_frac_time", "time in range"),
    ("cgm_normal_frac_time", "time in range"),
    ("cgm_mod_hyper_frac_time", "time in range"),
    ("cgm_sev_hyper_frac_time", "time in range"),
    ("cgm_time_below_range_frac", "time in range"),
    ("cgm_time_above_range_frac", "time in range"),
    # Excursion frequency: how often glucose crosses into a range. This is the
    # axis a mean cannot see -- two people with identical means differ here.
    ("cgm_sev_hypo_windows_per_day", "excursion frequency"),
    ("cgm_mod_hypo_windows_per_day", "excursion frequency"),
    ("cgm_normal_windows_per_day", "excursion frequency"),
    ("cgm_mod_hyper_windows_per_day", "excursion frequency"),
    ("cgm_sev_hyper_windows_per_day", "excursion frequency"),
    # Same, counting only episodes that persist >= 15 min (clinical event).
    ("cgm_sev_hypo_sustained_windows_per_day", "sustained excursions"),
    ("cgm_mod_hypo_sustained_windows_per_day", "sustained excursions"),
    ("cgm_normal_sustained_windows_per_day", "sustained excursions"),
    ("cgm_mod_hyper_sustained_windows_per_day", "sustained excursions"),
    ("cgm_sev_hyper_sustained_windows_per_day", "sustained excursions"),
    # Spike pipeline: postprandial-style rises, independent of the range bins.
    ("avg_spikes_per_day", "spike pipeline"),
    ("avg_spike_duration_minutes", "spike pipeline"),
    ("avg_cgm_per_spike_mg_dl", "spike pipeline"),
    # Level within a range, for the two ranges with near-complete coverage.
    ("cgm_normal_mean_mg_dl", "within-range level"),
    ("cgm_mod_hyper_mean_mg_dl", "within-range level"),
]

# Reported next to each result for context but NOT part of the FDR grid: these
# are the two targets already tested, and every survey feature here is null
# against them by construction.
REFERENCE_TARGETS = ["hba1c_percent", "mean_cgm_overall_mg_dl"]


def spearman_axis(x: np.ndarray, y: np.ndarray, axis: int = -1) -> np.ndarray:
    """Spearman rho along `axis`, vectorised over bootstrap resamples.

    Ranks are recomputed inside each resample; resampling with replacement
    creates ties that change them. Lifted from
    src/survey_glycemic_spearman_bootstrap.py so the CIs are comparable.
    """
    rx = rankdata(x, axis=axis).astype(float)
    ry = rankdata(y, axis=axis).astype(float)
    rx -= rx.mean(axis=axis, keepdims=True)
    ry -= ry.mean(axis=axis, keepdims=True)
    num = (rx * ry).sum(axis=axis)
    den = np.sqrt((rx**2).sum(axis=axis) * (ry**2).sum(axis=axis))
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def bootstrap_ci(x: np.ndarray, y: np.ndarray) -> tuple[float, float, str]:
    """BCa 95% CI for Spearman rho, falling back to percentile if BCa fails."""
    for method in ("BCa", "percentile"):
        try:
            with np.errstate(invalid="ignore", divide="ignore"):
                res = bootstrap(
                    (x, y),
                    spearman_axis,
                    n_resamples=N_BOOT,
                    confidence_level=CI_LEVEL,
                    method=method,
                    paired=True,
                    random_state=SEED,
                    vectorized=True,
                )
            low = float(res.confidence_interval.low)
            high = float(res.confidence_interval.high)
            if np.isfinite(low) and np.isfinite(high):
                return low, high, method
        except Exception:
            continue
    return float("nan"), float("nan"), "failed"


def null_survey_features() -> pd.DataFrame:
    """Survey features that failed FDR against BOTH HbA1c and mean CGM."""
    hba1c = pd.read_csv(HBA1C_RESULTS)
    cgm = pd.read_csv(CGM_RESULTS)
    merged = hba1c.merge(
        cgm, on=["feature", "survey", "direction"], suffixes=("_hba1c", "_cgm")
    )
    null = merged[
        ~merged["significant_fdr_0_05_hba1c"] & ~merged["significant_fdr_0_05_cgm"]
    ]
    return null[
        ["feature", "survey", "direction", "spearman_rho_hba1c", "spearman_rho_cgm"]
    ].reset_index(drop=True)


def load_data() -> pd.DataFrame:
    df = pd.read_csv(DATA_PATH)
    n0 = len(df)

    # Pooled consensus ranges. frac_time across the five bins sums to 1, so
    # these are exact rather than approximations.
    df["cgm_time_below_range_frac"] = (
        df["cgm_sev_hypo_frac_time"] + df["cgm_mod_hypo_frac_time"]
    )
    df["cgm_time_above_range_frac"] = (
        df["cgm_mod_hyper_frac_time"] + df["cgm_sev_hyper_frac_time"]
    )

    suspect = df["cgm_quality_suspect"].fillna(False).astype(bool)
    df = df[~suspect].copy()
    print(f"Dropped {int(suspect.sum())} cgm_quality_suspect participants "
          f"({n0} -> {len(df)})")

    df["severity"] = df["meta_study_group"].map(SEVERITY_ORDER)
    return df


def scan(df: pd.DataFrame, survey: pd.DataFrame) -> pd.DataFrame:
    """Stage 1: Spearman over the full survey x CGM grid."""
    rows = []
    for _, s in survey.iterrows():
        for cgm_col, family in CGM_FEATURES:
            pair = df[[s["feature"], cgm_col]].apply(
                pd.to_numeric, errors="coerce"
            ).dropna()
            n = len(pair)
            if n < 3 or pair[s["feature"]].nunique() < 2 or pair[cgm_col].nunique() < 2:
                continue
            res = spearmanr(pair[s["feature"]], pair[cgm_col])
            rows.append(
                {
                    "survey_feature": s["feature"],
                    "survey": s["survey"],
                    "direction": s["direction"],
                    "cgm_feature": cgm_col,
                    "cgm_family": family,
                    "n": n,
                    "spearman_rho": float(res.correlation),
                    "p_value": float(res.pvalue),
                    "rho_vs_hba1c": s["spearman_rho_hba1c"],
                    "rho_vs_mean_cgm": s["spearman_rho_cgm"],
                }
            )

    results = pd.DataFrame(rows)
    reject, p_fdr, _, _ = multipletests(
        results["p_value"].to_numpy(), alpha=ALPHA, method="fdr_bh"
    )
    results["p_value_fdr_bh"] = p_fdr
    results["significant_fdr_0_05"] = reject
    return results.sort_values("p_value").reset_index(drop=True)


def add_bootstrap(df: pd.DataFrame, results: pd.DataFrame) -> pd.DataFrame:
    """Stage 2: BCa CI for the FDR survivors only."""
    results = results.copy()
    for col in ("ci_low", "ci_high"):
        results[col] = np.nan
    results["ci_method"] = ""

    hits = results.index[results["significant_fdr_0_05"]]
    print(f"\nBootstrapping {len(hits)} survivor(s) at n_boot={N_BOOT:,}...")
    for i in hits:
        r = results.loc[i]
        pair = df[[r["survey_feature"], r["cgm_feature"]]].apply(
            pd.to_numeric, errors="coerce"
        ).dropna()
        low, high, method = bootstrap_ci(
            pair[r["survey_feature"]].to_numpy(), pair[r["cgm_feature"]].to_numpy()
        )
        results.loc[i, ["ci_low", "ci_high", "ci_method"]] = [low, high, method]
        print(f"  {r['survey_feature']:<28} x {r['cgm_feature']:<42} "
              f"rho={r['spearman_rho']:+.3f} CI[{low:+.3f}, {high:+.3f}]")

    results["ci_excludes_zero"] = np.where(
        results["ci_low"].notna() & (results["ci_low"] * results["ci_high"] > 0),
        "\u2713",
        "\u2717",
    )
    return results


def stratified_check(df: pd.DataFrame, results: pd.DataFrame) -> None:
    """Re-run each survivor within diabetes severity stratum."""
    hits = results[results["significant_fdr_0_05"]]
    if hits.empty:
        return
    print("\n" + "=" * 100)
    print("SEVERITY-STRATIFIED RE-CHECK (survivors only)")
    print("=" * 100)
    for _, r in hits.iterrows():
        sf, cf = r["survey_feature"], r["cgm_feature"]
        print(f"\n{sf} x {cf}   pooled rho={r['spearman_rho']:+.3f}")
        conf = df[["severity", cf]].dropna()
        rc = spearmanr(conf["severity"], conf[cf])
        print(f"  {cf} vs severity: rho={rc.correlation:+.3f} p={rc.pvalue:.2e}")
        conf = df[["severity", sf]].dropna()
        rs = spearmanr(conf["severity"], conf[sf])
        print(f"  {sf} vs severity: rho={rs.correlation:+.3f} p={rs.pvalue:.2e}")
        for group, sub in df.groupby("meta_study_group"):
            pair = sub[[sf, cf]].dropna()
            if len(pair) < 30:
                continue
            rr = spearmanr(pair[sf], pair[cf])
            print(f"    {group[:45]:<47} n={len(pair):<5} "
                  f"rho={rr.correlation:+.3f} p={rr.pvalue:.3f}")


def main() -> None:
    survey = null_survey_features()
    df = load_data()

    print(f"\nSurvey features null against BOTH targets: {len(survey)}")
    print(survey[["feature", "survey", "direction"]].to_string(index=False))
    print(f"\nCandidate CGM features: {len(CGM_FEATURES)}")
    for fam in dict.fromkeys(f for _, f in CGM_FEATURES):
        cols = [c for c, f in CGM_FEATURES if f == fam]
        print(f"  {fam} ({len(cols)}): {', '.join(cols)}")
    print(f"\nGrid: {len(survey)} x {len(CGM_FEATURES)} = "
          f"{len(survey) * len(CGM_FEATURES)} tests, BH-FDR across ALL of them")

    results = scan(df, survey)
    n_sig = int(results["significant_fdr_0_05"].sum())
    print(f"\nTests run: {len(results)}")
    print(f"Significant after global BH-FDR at alpha={ALPHA}: {n_sig}/{len(results)}")
    print(f"Expected false positives at alpha if all null: "
          f"{ALPHA * len(results):.1f} uncorrected")

    results = add_bootstrap(df, results)

    print("\n" + "=" * 100)
    print("TOP 20 BY RAW p (for context; FDR verdict in the significant column)")
    print("=" * 100)
    show = results.head(20).copy()
    show["rho"] = show["spearman_rho"].map(lambda v: f"{v:+.3f}")
    show["p"] = show["p_value"].map(lambda v: f"{v:.2e}")
    show["p_FDR"] = show["p_value_fdr_bh"].map(lambda v: f"{v:.3f}")
    show["sig"] = np.where(show["significant_fdr_0_05"], "\u2713", "\u2717")
    print(
        show[["survey_feature", "cgm_feature", "cgm_family", "n", "rho", "p", "p_FDR", "sig"]]
        .to_string(index=False)
    )

    print("\n--- best hit per survey feature (raw p) ---")
    best = results.loc[results.groupby("survey_feature")["p_value"].idxmin()]
    best = best.sort_values("p_value")
    b = best.copy()
    b["rho"] = b["spearman_rho"].map(lambda v: f"{v:+.3f}")
    b["p_FDR"] = b["p_value_fdr_bh"].map(lambda v: f"{v:.3f}")
    b["sig"] = np.where(b["significant_fdr_0_05"], "\u2713", "\u2717")
    print(b[["survey_feature", "cgm_feature", "n", "rho", "p_FDR", "sig"]].to_string(index=False))

    stratified_check(df, results)

    results.to_csv(OUT_PATH, index=False)
    print(f"\nWrote: {OUT_PATH}  ({len(results)} rows)")


if __name__ == "__main__":
    main()
