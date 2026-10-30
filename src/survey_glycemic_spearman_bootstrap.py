"""
Spearman correlation between every survey feature and each glycemic outcome
(lab HbA1c and mean CGM), each with a bootstrap 95% confidence interval and an
FDR-corrected p-value.

Features are the survey columns added to full_person_aggregates.csv by
src/add_surveys_to_full_person_aggregates.py and documented in
docs/survey_features.md. Of the 46 survey columns, 44 are numeric and tested
here; sdoh_education_label and sdoh_housing_status_label are text renderings of
sdoh_education_code and sdoh_housing_status_code, which are already included.

For each feature we report:
  - n (complete pairs after dropping NaNs)
  - Spearman rho and its asymptotic p-value
  - bootstrap 95% CI for rho (BCa, resampling person rows with replacement,
    falling back to the percentile method if BCa cannot be computed)
  - FDR (Benjamini-Hochberg) adjusted p-value across all tested features,
    corrected within each target rather than across both
  - direction metadata, because the diabetes self-care scale is scored so that
    higher = BETTER while the other three scales are higher = worse, and the
    sign of rho cannot be interpreted without knowing which

The bootstrap statistic is a vectorised Spearman (rank each resample, then
Pearson along the resample axis) so that scipy can evaluate all 10,000
resamples in one batched call. It is exact to floating point against
scipy.stats.spearmanr, verified on continuous, binary, and tie-heavy features,
and reproduces the same BCa interval roughly 3x faster.

Note the two targets are not interchangeable. HbA1c is a single lab draw
covering 2182 of 2245 participants; mean CGM is a multi-day wearable average
covering all 2245, so CGM tests run at slightly larger n.

Run from repo root:
    .venv/bin/python src/survey_glycemic_spearman_bootstrap.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap, rankdata, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"
OUTPUTS = ROOT / "outputs"

# target column -> output filename
TARGETS = {
    "hba1c_percent": "survey_hba1c_spearman_bootstrap.csv",
    "mean_cgm_overall_mg_dl": "survey_mean_cgm_spearman_bootstrap.csv",
}

N_BOOT = 10_000
SEED = 42
CI_LEVEL = 0.95
ALPHA = 0.05

WORSE = "higher = worse"
BETTER = "higher = better"
NOMINAL = "nominal, sign not interpretable"
QC = "completeness check, not a construct"

# (feature, survey, direction) in survey order.
FEATURES: list[tuple[str, str, str]] = [
    ("cesd10_total", "CES-D-10 depression", WORSE),
    ("cesd10_depressed", "CES-D-10 depression", WORSE),
    ("cesd10_depressed_affect", "CES-D-10 depression", WORSE),
    ("cesd10_somatic", "CES-D-10 depression", WORSE),
    ("cesd10_anhedonia", "CES-D-10 depression", WORSE),
    ("cesd10_items_answered", "CES-D-10 depression", QC),

    ("paid5_total", "PAID-5 diabetes distress", WORSE),
    ("paid5_distress", "PAID-5 diabetes distress", WORSE),
    ("paid5_items_answered", "PAID-5 diabetes distress", QC),

    ("diet_risk_score", "Dietary questionnaire", WORSE),
    ("diet_sugary_drinks", "Dietary questionnaire", WORSE),
    ("diet_fruit_veg_deficit", "Dietary questionnaire", WORSE),
    ("diet_fat_processed", "Dietary questionnaire", WORSE),
    ("diet_items_answered", "Dietary questionnaire", QC),

    ("dm_selfcare_total", "Diabetes self-care", BETTER),
    ("dm_selfcare_exercise", "Diabetes self-care", BETTER),
    ("dm_selfcare_diet", "Diabetes self-care", BETTER),
    ("dm_selfcare_healthcare", "Diabetes self-care", BETTER),
    ("dm_selfcare_footcare", "Diabetes self-care", BETTER),
    ("dm_selfcare_items_answered", "Diabetes self-care", QC),

    ("sdoh_education_code", "PhenX SDOH", BETTER),
    ("sdoh_education_years", "PhenX SDOH", BETTER),
    ("sdoh_college_graduate", "PhenX SDOH", BETTER),
    ("sdoh_english_first_language", "PhenX SDOH", NOMINAL),
    ("sdoh_housing_status_code", "PhenX SDOH", NOMINAL),
    ("sdoh_housing_instability", "PhenX SDOH", WORSE),
    ("sdoh_housing_problem_count", "PhenX SDOH", WORSE),
    ("sdoh_any_housing_problem", "PhenX SDOH", WORSE),
    ("sdoh_insurance_type_count", "PhenX SDOH", NOMINAL),
    ("sdoh_uninsured", "PhenX SDOH", WORSE),
    ("sdoh_rx_cost_barrier_count", "PhenX SDOH", WORSE),
    ("sdoh_rx_cost_barrier_self", "PhenX SDOH", WORSE),
    ("sdoh_care_access_barrier_count", "PhenX SDOH", WORSE),
    ("sdoh_no_usual_care_place", "PhenX SDOH", WORSE),
    ("sdoh_hospitalized_overnight", "PhenX SDOH", WORSE),
    ("sdoh_urgent_care_visits_code", "PhenX SDOH", WORSE),
    ("sdoh_med_discrimination_mean", "PhenX SDOH", WORSE),
    ("sdoh_med_discrimination_any", "PhenX SDOH", WORSE),
    ("sdoh_job_loss_likely", "PhenX SDOH", WORSE),
    ("sdoh_adversity_index", "PhenX SDOH", WORSE),

    ("food_insecurity_score", "PhenX food insecurity", WORSE),
    ("food_insecure", "PhenX food insecurity", WORSE),
    ("food_insecurity_hvs_positive", "PhenX food insecurity", WORSE),
    ("food_insecurity_items_answered", "PhenX food insecurity", QC),
]

# Text renderings of *_code columns already tested above.
EXCLUDED_TEXT = ["sdoh_education_label", "sdoh_housing_status_label"]


def spearman_axis(x: np.ndarray, y: np.ndarray, axis: int = -1) -> np.ndarray:
    """Spearman rho along `axis`, vectorised over bootstrap resamples.

    Ranks are recomputed within each resample, which is what makes this
    equivalent to scipy.stats.spearmanr rather than an approximation: resampling
    with replacement creates ties that change the ranks.
    """
    rx = rankdata(x, axis=axis).astype(float)
    ry = rankdata(y, axis=axis).astype(float)
    rx -= rx.mean(axis=axis, keepdims=True)
    ry -= ry.mean(axis=axis, keepdims=True)
    num = (rx * ry).sum(axis=axis)
    den = np.sqrt((rx**2).sum(axis=axis) * (ry**2).sum(axis=axis))
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def bootstrap_ci(target: np.ndarray, feat: np.ndarray) -> tuple[float, float, str]:
    """BCa 95% CI for Spearman rho, falling back to percentile if BCa fails.

    BCa needs a finite acceleration estimate, which it cannot produce for very
    sparse binary features where many jackknife replicates are identical.
    """
    for method in ("BCa", "percentile"):
        try:
            with np.errstate(invalid="ignore", divide="ignore"):
                res = bootstrap(
                    (target, feat),
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


def analyze(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Run the full bootstrap + FDR pipeline for one target column."""
    rows = []
    for i, (col, survey, direction) in enumerate(FEATURES, start=1):
        pair = df[[target, col]].apply(pd.to_numeric, errors="coerce").dropna()
        n = len(pair)
        print(f"  [{i:2d}/{len(FEATURES)}] {col:<32} n={n:<5}", end="", flush=True)

        if n < 3 or pair[col].nunique() < 2:
            print(" skipped (degenerate)")
            rows.append(
                {
                    "feature": col, "survey": survey, "direction": direction,
                    "n": n, "spearman_rho": np.nan, "ci_low": np.nan,
                    "ci_high": np.nan, "p_value": np.nan, "ci_method": "n/a",
                }
            )
            continue

        t = pair[target].to_numpy()
        f = pair[col].to_numpy()
        res = spearmanr(t, f)
        low, high, method = bootstrap_ci(t, f)
        print(f" rho={res.correlation:+.3f}  CI[{low:+.3f}, {high:+.3f}]")

        rows.append(
            {
                "feature": col, "survey": survey, "direction": direction,
                "n": n, "spearman_rho": float(res.correlation),
                "ci_low": low, "ci_high": high,
                "p_value": float(res.pvalue), "ci_method": method,
            }
        )

    results = pd.DataFrame(rows)

    valid = results["p_value"].notna()
    results["p_value_fdr_bh"] = np.nan
    results["significant_fdr_0_05"] = False
    if valid.any():
        reject, p_fdr, _, _ = multipletests(
            results.loc[valid, "p_value"].to_numpy(), alpha=ALPHA, method="fdr_bh"
        )
        results.loc[valid, "p_value_fdr_bh"] = p_fdr
        results.loc[valid, "significant_fdr_0_05"] = reject

    results["significant"] = np.where(
        results["significant_fdr_0_05"], "\u2713", "\u2717"
    )
    # A CI excluding zero is the bootstrap's own verdict, reported alongside the
    # FDR verdict because the two can disagree for sparse or skewed features.
    results["ci_excludes_zero"] = np.where(
        results["ci_low"].notna() & (results["ci_low"] * results["ci_high"] > 0),
        "\u2713",
        "\u2717",
    )

    results = results.sort_values("p_value", na_position="last").reset_index(drop=True)
    results = results[
        [
            "feature", "survey", "direction", "n", "spearman_rho",
            "ci_low", "ci_high", "ci_excludes_zero", "p_value",
            "p_value_fdr_bh", "significant_fdr_0_05", "significant", "ci_method",
        ]
    ]
    n_sig = int(results["significant_fdr_0_05"].sum())
    print(f"Significant after Benjamini-Hochberg FDR at alpha={ALPHA}: "
          f"{n_sig}/{int(valid.sum())}")
    if (results["ci_method"] == "percentile").any():
        fell_back = results.loc[results["ci_method"] == "percentile", "feature"]
        print(f"Percentile CI fallback used for: {list(fell_back)}")

    show = results.copy()
    show["95% CI"] = [
        f"[{lo:+.3f}, {hi:+.3f}]" if np.isfinite(lo) else "n/a"
        for lo, hi in zip(show["ci_low"], show["ci_high"])
    ]
    show["rho"] = show["spearman_rho"].map(lambda v: f"{v:+.3f}")
    show["p"] = show["p_value"].map(lambda v: f"{v:.2e}")
    show["p_FDR"] = show["p_value_fdr_bh"].map(lambda v: f"{v:.2e}")
    print()
    print(
        show[["feature", "n", "rho", "95% CI", "p", "p_FDR", "significant", "direction"]]
        .to_string(index=False)
    )
    return results


def main() -> None:
    df = pd.read_csv(DATA_PATH)
    missing = [f for f, _, _ in FEATURES if f not in df.columns]
    assert not missing, f"features absent from {DATA_PATH.name}: {missing}"

    print(f"Testing {len(FEATURES)} survey features against {len(TARGETS)} targets")
    print(f"Excluded text duplicates: {EXCLUDED_TEXT}")
    print(f"(n_boot={N_BOOT}, method=BCa, ci_level={CI_LEVEL}, seed={SEED})")

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    for target, filename in TARGETS.items():
        assert target in df.columns, f"target absent: {target}"
        print("\n" + "=" * 100)
        print(f"TARGET: {target}  (non-null: {df[target].notna().sum()}/{len(df)})")
        print("=" * 100)
        results = analyze(df, target)
        out = OUTPUTS / filename
        results.to_csv(out, index=False)
        print(f"\nWrote: {out}")


if __name__ == "__main__":
    main()
