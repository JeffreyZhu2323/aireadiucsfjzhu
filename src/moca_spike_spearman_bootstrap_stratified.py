"""
STRATIFIED Spearman correlation between MoCA score and the three person-level
glucose-spike aggregate features, each with a bootstrap 95% CI.

Same statistic/engine as ``moca_spike_spearman_bootstrap.py`` (Spearman rho,
BCa bootstrap CI, BH-FDR across the three features), but computed WITHIN strata:

  1. by diabetes status   (meta_study_group; 4 groups)  -> ..._by_diabetes.csv
  2. by age               (meta_age; 40-54 / 55-69 / 70+) -> ..._by_age.csv
  3. by diabetes x age    (all group x age-band combos)   -> ..._by_diabetes_age.csv

FDR correction is applied across the three features *within* each stratum.

Run from repo root:
    .venv/bin/python src/moca_spike_spearman_bootstrap_stratified.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"

OUT_BY_DIABETES = ROOT / "outputs" / "cognition" / "moca_spike_spearman_bootstrap_by_diabetes.csv"
OUT_BY_AGE = ROOT / "outputs" / "cognition" / "moca_spike_spearman_bootstrap_by_age.csv"
OUT_BY_DIABETES_AGE = ROOT / "outputs" / "cognition" / "moca_spike_spearman_bootstrap_by_diabetes_age.csv"

TARGET = "moca_score"
FEATURES = [
    "avg_spike_duration_minutes",
    "avg_spikes_per_day",
    "avg_cgm_per_spike_mg_dl",
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


def _analyze_subset(sub: pd.DataFrame, strata_keys: dict) -> list[dict]:
    """Spearman + BCa bootstrap CI for each feature within one stratum.

    FDR (BH) is applied across the features that yielded a valid p-value.
    """
    rows: list[dict] = []
    raw_pvals: list[float] = []

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
            raw_pvals.append(float("nan"))
            continue

        target = tmp[TARGET].to_numpy()
        feat = tmp[col].to_numpy()

        res = spearmanr(target, feat)
        rho = float(res.correlation)
        pval = float(res.pvalue)

        try:
            boot = bootstrap(
                (target, feat),
                _spearman_stat,
                n_resamples=N_BOOT,
                confidence_level=CI_LEVEL,
                method="BCa",
                paired=True,
                random_state=SEED,
            )
            ci_low = float(boot.confidence_interval.low)
            ci_high = float(boot.confidence_interval.high)
        except Exception:
            # BCa can fail on degenerate resamples; fall back to percentile.
            boot = bootstrap(
                (target, feat),
                _spearman_stat,
                n_resamples=N_BOOT,
                confidence_level=CI_LEVEL,
                method="percentile",
                paired=True,
                random_state=SEED,
            )
            ci_low = float(boot.confidence_interval.low)
            ci_high = float(boot.confidence_interval.high)

        row.update(
            spearman_rho=rho,
            p_value=pval,
            ci_low=ci_low,
            ci_high=ci_high,
        )
        rows.append(row)
        raw_pvals.append(pval)

    res_df = pd.DataFrame(rows)
    valid = res_df["p_value"].notna()
    pv = res_df.loc[valid, "p_value"].to_numpy()
    res_df["p_value_fdr_bh"] = float("nan")
    res_df["significant_fdr_0_05"] = False
    if len(pv) > 0:
        reject, pv_fdr, _, _ = multipletests(pv, alpha=0.05, method="fdr_bh")
        res_df.loc[valid, "p_value_fdr_bh"] = pv_fdr
        res_df.loc[valid, "significant_fdr_0_05"] = reject

    return res_df.to_dict("records")


def _prep(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df[AGE_COL] = pd.to_numeric(df[AGE_COL], errors="coerce")
    df["age_stratum"] = pd.cut(df[AGE_COL], bins=AGE_BINS, right=False, labels=AGE_LABELS)
    return df


def run_by_diabetes(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for grp in DIABETES_GROUPS:
        sub = df[df[DIABETES_COL] == grp]
        rows.extend(_analyze_subset(sub, {"diabetes_group": grp}))
    return pd.DataFrame(rows)


def run_by_age(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for band in AGE_LABELS:
        sub = df[df["age_stratum"] == band]
        rows.extend(_analyze_subset(sub, {"age_stratum": band}))
    return pd.DataFrame(rows)


def run_by_diabetes_age(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for grp in DIABETES_GROUPS:
        for band in AGE_LABELS:
            sub = df[(df[DIABETES_COL] == grp) & (df["age_stratum"] == band)]
            rows.extend(
                _analyze_subset(sub, {"diabetes_group": grp, "age_stratum": band})
            )
    return pd.DataFrame(rows)


def main() -> None:
    df = _prep(pd.read_csv(DATA_PATH))

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)

    print(f"(n_boot={N_BOOT}, method=BCa, ci_level={CI_LEVEL}, seed={SEED}, min_n={MIN_N})\n")

    for label, runner, out_path in [
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
