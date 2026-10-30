"""
Spearman correlation between MoCA score and person-level MEAN GLUCOSE
(overall mean of all CGM readings across all time points), with a bootstrap
95% CI (BCa), across the same stratifications used for the wearable-stream
analysis:

  1. OVERALL                         -> moca_glucose_bootstrap.csv
  2. by diabetes status              -> moca_glucose_bootstrap_by_diabetes.csv
  3. by age band                     -> moca_glucose_bootstrap_by_age.csv
  4. by diabetes status x age band   -> moca_glucose_bootstrap_by_diabetes_age.csv

Feature: mean_glucose (from master_patient_features.csv; verified to equal the
per-person mean of all raw CGM readings in outputs/glycemic_phenotype/glucose_raw.csv).

Multiple-comparison note: the wearable analysis applied BH-FDR across the FIVE
features within each stratum. Here there is only ONE feature, so within-stratum
FDR is a no-op. Instead we apply BH-FDR ACROSS THE STRATA within each
stratification family (4 diabetes groups; 3 age bands; 12 diabetes x age cells).
The raw p-value is also reported. OVERALL has a single test (FDR == raw p).

Run from repo root:
    .venv/bin/python src/moca_glucose_bootstrap.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
FEATURES_PATH = ROOT / "outputs" / "eda_advanced" / "tables" / "master_patient_features.csv"

OUT_OVERALL = ROOT / "outputs" / "cognition" / "moca_glucose_bootstrap.csv"
OUT_BY_DIABETES = ROOT / "outputs" / "cognition" / "moca_glucose_bootstrap_by_diabetes.csv"
OUT_BY_AGE = ROOT / "outputs" / "cognition" / "moca_glucose_bootstrap_by_age.csv"
OUT_BY_DIABETES_AGE = ROOT / "outputs" / "cognition" / "moca_glucose_bootstrap_by_diabetes_age.csv"

TARGET = "moca_score"
FEATURE = "mean_glucose"

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


def _corr_one(sub: pd.DataFrame, strata_keys: dict) -> dict:
    """Spearman + BCa CI for mean_glucose vs MoCA within one stratum."""
    tmp = sub[[TARGET, FEATURE]].copy()
    tmp[TARGET] = pd.to_numeric(tmp[TARGET], errors="coerce")
    tmp[FEATURE] = pd.to_numeric(tmp[FEATURE], errors="coerce")
    tmp = tmp.dropna()
    n = len(tmp)

    row = {**strata_keys, "feature": FEATURE, "n": n}
    if n < MIN_N or tmp[TARGET].nunique() < 2 or tmp[FEATURE].nunique() < 2:
        row.update(spearman_rho=float("nan"), p_value=float("nan"), ci_low=float("nan"), ci_high=float("nan"))
        return row

    target = tmp[TARGET].to_numpy()
    feat = tmp[FEATURE].to_numpy()
    res = spearmanr(target, feat)
    ci_low, ci_high = _boot_ci(target, feat)
    row.update(
        spearman_rho=float(res.correlation),
        p_value=float(res.pvalue),
        ci_low=ci_low,
        ci_high=ci_high,
    )
    return row


def _add_fdr(df: pd.DataFrame) -> pd.DataFrame:
    """BH-FDR across the strata (rows) of one stratification family."""
    df = df.copy()
    df["p_value_fdr_bh"] = float("nan")
    df["significant_fdr_0_05"] = False
    valid = df["p_value"].notna()
    pv = df.loc[valid, "p_value"].to_numpy()
    if len(pv) > 0:
        reject, pv_fdr, _, _ = multipletests(pv, alpha=0.05, method="fdr_bh")
        df.loc[valid, "p_value_fdr_bh"] = pv_fdr
        df.loc[valid, "significant_fdr_0_05"] = reject
    return df


def _prep() -> pd.DataFrame:
    keep = ["participant_id", TARGET, FEATURE, AGE_COL, DIABETES_COL]
    df = pd.read_csv(FEATURES_PATH, usecols=lambda c: c in keep)
    df[AGE_COL] = pd.to_numeric(df[AGE_COL], errors="coerce")
    df["age_stratum"] = pd.cut(df[AGE_COL], bins=AGE_BINS, right=False, labels=AGE_LABELS)
    return df


def run_overall(df: pd.DataFrame) -> pd.DataFrame:
    return _add_fdr(pd.DataFrame([_corr_one(df, {"stratum": "ALL"})]))


def run_by_diabetes(df: pd.DataFrame) -> pd.DataFrame:
    rows = [_corr_one(df[df[DIABETES_COL] == g], {"diabetes_group": g}) for g in DIABETES_GROUPS]
    return _add_fdr(pd.DataFrame(rows))


def run_by_age(df: pd.DataFrame) -> pd.DataFrame:
    rows = [_corr_one(df[df["age_stratum"] == b], {"age_stratum": b}) for b in AGE_LABELS]
    return _add_fdr(pd.DataFrame(rows))


def run_by_diabetes_age(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g in DIABETES_GROUPS:
        for b in AGE_LABELS:
            sub = df[(df[DIABETES_COL] == g) & (df["age_stratum"] == b)]
            rows.append(_corr_one(sub, {"diabetes_group": g, "age_stratum": b}))
    return _add_fdr(pd.DataFrame(rows))


def main() -> None:
    df = _prep()
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)

    print(f"(feature={FEATURE}, n_boot={N_BOOT}, method=BCa, ci_level={CI_LEVEL}, seed={SEED}, min_n={MIN_N})")
    print(f"persons with MoCA + mean_glucose: {df[[TARGET, FEATURE]].dropna().shape[0]}")
    print("FDR: BH across strata within each family (single feature)\n")

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
