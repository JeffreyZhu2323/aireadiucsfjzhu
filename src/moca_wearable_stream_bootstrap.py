"""
Spearman correlation between MoCA score and the person-level MEAN of the five
non-physical-activity wearable streams, each with a bootstrap 95% CI (BCa),
plus Benjamini-Hochberg FDR-adjusted p-values.

Features (person-level means from outputs/wearable_stream_means.csv):
  - mean_heartrate_bpm
  - mean_oxygen_saturation_pct
  - mean_stress_level
  - mean_respiratory_rate_bpm
  - mean_sleep_hours

Analyses produced:
  1. OVERALL                         -> moca_wearable_stream_bootstrap.csv
  2. by diabetes status              -> moca_wearable_stream_bootstrap_by_diabetes.csv
  3. by age band                     -> moca_wearable_stream_bootstrap_by_age.csv
  4. by diabetes status x age band   -> moca_wearable_stream_bootstrap_by_diabetes_age.csv

For every (analysis, stratum) we report per feature:
  n, spearman_rho, p_value (asymptotic), bootstrap 95% CI (BCa, paired resample
  of person rows), FDR-adjusted p, and the FDR significance flag. FDR is applied
  across the five features WITHIN each stratum, matching the repo's spike-analysis
  convention.

Run from repo root (after wearable_stream_aggregates.py):
    .venv/bin/python src/moca_wearable_stream_bootstrap.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path.home() / "aireadi"
MEANS_PATH = ROOT / "outputs" / "wearable_stream_means.csv"
FEATURES_PATH = ROOT / "outputs" / "eda_advanced" / "tables" / "master_patient_features.csv"

OUT_OVERALL = ROOT / "outputs" / "moca_wearable_stream_bootstrap.csv"
OUT_BY_DIABETES = ROOT / "outputs" / "moca_wearable_stream_bootstrap_by_diabetes.csv"
OUT_BY_AGE = ROOT / "outputs" / "moca_wearable_stream_bootstrap_by_age.csv"
OUT_BY_DIABETES_AGE = ROOT / "outputs" / "moca_wearable_stream_bootstrap_by_diabetes_age.csv"

TARGET = "moca_score"
FEATURES = [
    "mean_heartrate_bpm",
    "mean_oxygen_saturation_pct",
    "mean_stress_level",
    "mean_respiratory_rate_bpm",
    "mean_sleep_hours",
]

N_BOOT = 10_000
SEED = 42
CI_LEVEL = 0.95
MIN_N = 10  # minimum complete pairs to attempt a correlation in a stratum

DIABETES_COL = "meta_study_group"
DIABETES_GROUPS = [
    "healthy",
    "pre_diabetes_lifestyle_controlled",
    "oral_medication_and_or_non_insulin_injectable_medication_controlled",
    "insulin_dependent",
]

AGE_COL = "meta_age"
AGE_BINS = [40, 55, 70, 200]
AGE_LABELS = ["40-54", "55-69", "70+"]


def _spearman_stat(x: np.ndarray, y: np.ndarray) -> float:
    return spearmanr(x, y).correlation


def _boot_ci(target: np.ndarray, feat: np.ndarray) -> tuple[float, float]:
    """BCa bootstrap 95% CI for Spearman rho, falling back to percentile."""
    for method in ("BCa", "percentile"):
        try:
            boot = bootstrap(
                (target, feat),
                _spearman_stat,
                n_resamples=N_BOOT,
                confidence_level=CI_LEVEL,
                method=method,
                paired=True,
                random_state=SEED,
            )
            return float(boot.confidence_interval.low), float(boot.confidence_interval.high)
        except Exception:
            continue
    return float("nan"), float("nan")


def _analyze_subset(sub: pd.DataFrame, strata_keys: dict) -> pd.DataFrame:
    """Spearman + BCa CI per feature within one stratum; BH-FDR across features."""
    rows: list[dict] = []
    for col in FEATURES:
        tmp = sub[[TARGET, col]].copy()
        tmp[TARGET] = pd.to_numeric(tmp[TARGET], errors="coerce")
        tmp[col] = pd.to_numeric(tmp[col], errors="coerce")
        tmp = tmp.dropna()
        n = len(tmp)

        row = {**strata_keys, "feature": col, "n": n}

        if n < MIN_N or tmp[TARGET].nunique() < 2 or tmp[col].nunique() < 2:
            row.update(
                spearman_rho=float("nan"),
                p_value=float("nan"),
                ci_low=float("nan"),
                ci_high=float("nan"),
            )
            rows.append(row)
            continue

        target = tmp[TARGET].to_numpy()
        feat = tmp[col].to_numpy()
        res = spearmanr(target, feat)
        ci_low, ci_high = _boot_ci(target, feat)
        row.update(
            spearman_rho=float(res.correlation),
            p_value=float(res.pvalue),
            ci_low=ci_low,
            ci_high=ci_high,
        )
        rows.append(row)

    res_df = pd.DataFrame(rows)
    valid = res_df["p_value"].notna()
    res_df["p_value_fdr_bh"] = float("nan")
    res_df["significant_fdr_0_05"] = False
    pv = res_df.loc[valid, "p_value"].to_numpy()
    if len(pv) > 0:
        reject, pv_fdr, _, _ = multipletests(pv, alpha=0.05, method="fdr_bh")
        res_df.loc[valid, "p_value_fdr_bh"] = pv_fdr
        res_df.loc[valid, "significant_fdr_0_05"] = reject
    return res_df


def _prep() -> pd.DataFrame:
    means = pd.read_csv(MEANS_PATH)
    means["participant_id"] = pd.to_numeric(means["participant_id"], errors="coerce")

    feats = pd.read_csv(FEATURES_PATH)
    feats["participant_id"] = pd.to_numeric(feats["participant_id"], errors="coerce")

    keep = ["participant_id", TARGET, AGE_COL, DIABETES_COL]
    keep = [c for c in keep if c in feats.columns]
    df = means.merge(feats[keep], on="participant_id", how="inner")

    df[AGE_COL] = pd.to_numeric(df[AGE_COL], errors="coerce")
    df["age_stratum"] = pd.cut(df[AGE_COL], bins=AGE_BINS, right=False, labels=AGE_LABELS)
    return df


def run_overall(df: pd.DataFrame) -> pd.DataFrame:
    return _analyze_subset(df, {"stratum": "ALL"})


def run_by_diabetes(df: pd.DataFrame) -> pd.DataFrame:
    parts = [_analyze_subset(df[df[DIABETES_COL] == g], {"diabetes_group": g}) for g in DIABETES_GROUPS]
    return pd.concat(parts, ignore_index=True)


def run_by_age(df: pd.DataFrame) -> pd.DataFrame:
    parts = [_analyze_subset(df[df["age_stratum"] == b], {"age_stratum": b}) for b in AGE_LABELS]
    return pd.concat(parts, ignore_index=True)


def run_by_diabetes_age(df: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for g in DIABETES_GROUPS:
        for b in AGE_LABELS:
            sub = df[(df[DIABETES_COL] == g) & (df["age_stratum"] == b)]
            parts.append(_analyze_subset(sub, {"diabetes_group": g, "age_stratum": b}))
    return pd.concat(parts, ignore_index=True)


def main() -> None:
    df = _prep()
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)

    print(f"(n_boot={N_BOOT}, method=BCa, ci_level={CI_LEVEL}, seed={SEED}, min_n={MIN_N})")
    print(f"merged persons with MoCA + wearable means: {df[TARGET].notna().sum()}\n")

    for label, runner, out_path in [
        ("OVERALL", run_overall, OUT_OVERALL),
        ("BY DIABETES STATUS", run_by_diabetes, OUT_BY_DIABETES),
        ("BY AGE", run_by_age, OUT_BY_AGE),
        ("BY DIABETES x AGE", run_by_diabetes_age, OUT_BY_DIABETES_AGE),
    ]:
        results = runner(df)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(out_path, index=False)
        print(f"===== {label} -> {out_path.name} =====")
        print(results.to_string(index=False))
        print()


if __name__ == "__main__":
    main()
