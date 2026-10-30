"""
Age + diabetes-status JOINTLY-STRATIFIED bootstrap of impaired-vs-unimpaired
physical-activity means.

For each (age group x diabetes group) cell -- age groups 40-54/55-69/70+ from
meta_age, diabetes groups from meta_study_group -- and each MoCA cutoff
(26/24/23; impaired = moca_score < cutoff), we report ONLY the bootstrapped
group means:

  - mean physical-activity value for impaired and unimpaired, each with a 95%
    percentile CI from a stratified bootstrap (resample each group with
    replacement, sizes fixed), and
  - the mean difference (impaired - unimpaired) with its 95% CI and a flag for
    whether that CI excludes 0.

Uses the shared, vectorized ``bootstrap_means`` engine (same N_BOOT and seed as
the non-stratified analysis). Joint cells are small; any group with n<3 is
reported without a bootstrap.

Run from repo root:
    .venv/bin/python src/moca_impairment_activity_bootstrap_by_age_diabetes.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from moca_impairment_activity_tests import (
    CUTOFFS,
    N_BOOT,
    PA_METRICS,
    SEED,
    TARGET,
    _empty_mean_boot,
    bootstrap_means,
    load_data,
)

AGE = "meta_age"
AGE_BINS = [40, 55, 70, 200]
AGE_LABELS = ["40-54", "55-69", "70+"]

STRATUM_COL = "meta_study_group"
DIAB_STRATA = [
    "healthy",
    "pre_diabetes_lifestyle_controlled",
    "oral_medication_and_or_non_insulin_injectable_medication_controlled",
    "insulin_dependent",
]
SHORT = {
    "healthy": "healthy",
    "pre_diabetes_lifestyle_controlled": "pre_diabetes",
    "oral_medication_and_or_non_insulin_injectable_medication_controlled": "oral_med",
    "insulin_dependent": "insulin_dep",
}

OUT_PATH = (
    Path(__file__).resolve().parents[1]
    / "outputs"
    / "moca_impairment_activity_bootstrap_by_age_diabetes.csv"
)


def run() -> pd.DataFrame:
    df = load_data()
    df[AGE] = pd.to_numeric(df[AGE], errors="coerce")
    df["age_stratum"] = pd.cut(df[AGE], bins=AGE_BINS, right=False, labels=AGE_LABELS)

    rng = np.random.default_rng(SEED)
    rows: list[dict] = []

    for age_stratum in AGE_LABELS:
        for diab in DIAB_STRATA:
            cell = df[
                (df["age_stratum"] == age_stratum)
                & (df[STRATUM_COL] == diab)
                & df[TARGET].notna()
            ].copy()
            for cutoff in CUTOFFS:
                cell["impaired"] = cell[TARGET] < cutoff
                for col, desc in PA_METRICS.items():
                    sub = cell[[col, "impaired"]].dropna()
                    imp = sub.loc[sub["impaired"], col]
                    unimp = sub.loc[~sub["impaired"], col]
                    n_imp, n_unimp = len(imp), len(unimp)

                    row = {
                        "age_stratum": age_stratum,
                        "diabetes_group": diab,
                        "cutoff": cutoff,
                        "metric": col,
                        "description": desc,
                        "n_impaired": n_imp,
                        "n_unimpaired": n_unimp,
                    }
                    if n_imp < 3 or n_unimp < 3:
                        row.update(_empty_mean_boot())
                    else:
                        row.update(
                            bootstrap_means(
                                imp.to_numpy(dtype=float),
                                unimp.to_numpy(dtype=float),
                                rng,
                            )
                        )
                    rows.append(row)

    return pd.DataFrame(rows)


def main() -> None:
    results = run()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(OUT_PATH, index=False)

    pd.set_option("display.width", 300)
    pd.set_option("display.max_columns", 60)
    print(f"Wrote: {OUT_PATH}")
    print(f"Bootstrap: {N_BOOT} stratified resamples, seed={SEED}, 95% percentile CIs (means only)")
    print("Strata: age (40-54/55-69/70+) x diabetes (healthy/pre_diabetes/oral_med/insulin_dep)\n")

    def _ci(r, lo, hi):
        if pd.isna(r[lo]):
            return "[n/a]"
        return f"[{r[lo]:,.1f}, {r[hi]:,.1f}]"

    for age_stratum in AGE_LABELS:
        for diab in DIAB_STRATA:
            for cutoff in CUTOFFS:
                block = results[
                    (results["age_stratum"] == age_stratum)
                    & (results["diabetes_group"] == diab)
                    & (results["cutoff"] == cutoff)
                ].copy()
                n_imp = int(block["n_impaired"].iloc[0])
                n_unimp = int(block["n_unimpaired"].iloc[0])
                print(
                    f"===== Age {age_stratum} | {SHORT[diab]} | cutoff {cutoff}: "
                    f"impaired n={n_imp} vs unimpaired n={n_unimp} ====="
                )
                if n_imp < 3 or n_unimp < 3:
                    print("  (skipped: a group has n<3)\n")
                    continue
                block["mean_impaired_95ci"] = block.apply(
                    lambda r: _ci(r, "mean_impaired_ci_low", "mean_impaired_ci_high"), axis=1
                )
                block["mean_unimpaired_95ci"] = block.apply(
                    lambda r: _ci(r, "mean_unimpaired_ci_low", "mean_unimpaired_ci_high"), axis=1
                )
                block["mean_diff_95ci"] = block.apply(
                    lambda r: _ci(r, "mean_diff_ci_low", "mean_diff_ci_high"), axis=1
                )
                show = block[
                    [
                        "metric",
                        "mean_impaired",
                        "mean_impaired_95ci",
                        "mean_unimpaired",
                        "mean_unimpaired_95ci",
                        "mean_diff",
                        "mean_diff_95ci",
                        "mean_diff_ci_excludes_0",
                    ]
                ]
                print(show.to_string(index=False))
                print()


if __name__ == "__main__":
    main()
