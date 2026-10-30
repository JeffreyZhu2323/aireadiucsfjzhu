"""
Diabetes-status-STRATIFIED repeat of the two MoCA / physical-activity analyses.

Identical design to src/moca_pa_age_stratified.py, but strata are the four
diabetes study groups (meta_study_group):
  - healthy
  - pre_diabetes_lifestyle_controlled
  - oral_medication_and_or_non_insulin_injectable_medication_controlled
  - insulin_dependent

Within each diabetes group we run, independently:

  PART A -- numerical correlation between MoCA score and physical-activity metrics
    * Spearman rho + p (primary), Pearson r + p (secondary)
    * Benjamini-Hochberg FDR correction across metrics within the group
    Output: outputs/cognition/moca_pa_correlation_by_diabetes.csv

  PART B -- impaired vs unimpaired group differences (Mann-Whitney U)
    cutoffs 26 / 24 / 23 ; impaired = moca_score < cutoff
    metrics: total_movement, physical_activity_num_days,
             average_active_calories_kcal, average_daily_activity
    * group medians, Mann-Whitney U, two-sided p, rank-biserial effect size
    * FDR-BH correction across the 4 metrics within each (group, cutoff)
    Output: outputs/cognition/moca_impairment_activity_by_diabetes.csv

Run from repo root: .venv/bin/python src/moca_pa_diabetes_stratified.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
FEATURES_PATH = ROOT / "outputs" / "eda_advanced" / "tables" / "master_patient_features.csv"
ACTIVITY_MANIFEST = (
    ROOT
    / "data"
    / "manifests"
    / "jeffreyzhuaireadi"
    / "7bc8c89b-8457-41fc-917c-b5571ff90ea0"
    / "dataset"
    / "wearable_activity_monitor"
    / "manifest.tsv"
)
CORR_OUT = ROOT / "outputs" / "cognition" / "moca_pa_correlation_by_diabetes.csv"
DIFF_OUT = ROOT / "outputs" / "cognition" / "moca_impairment_activity_by_diabetes.csv"

TARGET = "moca_score"
STRATUM_COL = "meta_study_group"

# order strata from healthy -> most severe; short labels for printing
STRATA = [
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

CUTOFFS = [26, 24, 23]

CORR_METRICS = {
    "total_movement": "Total movement over study period",
    "day_movement": "Daytime movement (total - night)",
    "night_movement": "Nighttime movement",
    "average_daily_activity": "Mean daily steps (manifest)",
    "average_active_calories_kcal": "Mean active calories/day (manifest)",
    "physical_activity_num_days": "Days of activity coverage (manifest)",
}

DIFF_METRICS = {
    "total_movement": "Total movement over study period",
    "physical_activity_num_days": "Days of activity coverage (manifest)",
    "average_active_calories_kcal": "Mean active calories/day (manifest)",
    "average_daily_activity": "Mean daily steps (manifest)",
}


def load_data() -> pd.DataFrame:
    feats = pd.read_csv(FEATURES_PATH)
    feats["participant_id"] = pd.to_numeric(feats["participant_id"], errors="coerce")

    manifest = pd.read_csv(ACTIVITY_MANIFEST, sep="\t").rename(
        columns={"person_id": "participant_id"}
    )
    keep = [
        "participant_id",
        "average_daily_activity",
        "average_active_calories_kcal",
        "physical_activity_num_days",
    ]
    df = feats.merge(manifest[keep], on="participant_id", how="left")

    df["day_movement"] = pd.to_numeric(df["total_movement"], errors="coerce") - pd.to_numeric(
        df["night_movement"], errors="coerce"
    )
    df[TARGET] = pd.to_numeric(df[TARGET], errors="coerce")
    for col in set(CORR_METRICS) | set(DIFF_METRICS):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def correlation_by_group(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for stratum in STRATA:
        sdf = df[df[STRATUM_COL] == stratum]
        for col, desc in CORR_METRICS.items():
            pair = sdf[[TARGET, col]].dropna()
            n = len(pair)
            if n < 5:
                rows.append(
                    {
                        "diabetes_group": stratum,
                        "metric": col,
                        "description": desc,
                        "n": n,
                        "spearman_rho": np.nan,
                        "spearman_p": np.nan,
                        "pearson_r": np.nan,
                        "pearson_p": np.nan,
                    }
                )
                continue
            rho, p_s = spearmanr(pair[TARGET], pair[col])
            r, p_p = pearsonr(pair[TARGET], pair[col])
            rows.append(
                {
                    "diabetes_group": stratum,
                    "metric": col,
                    "description": desc,
                    "n": n,
                    "spearman_rho": float(rho),
                    "spearman_p": float(p_s),
                    "pearson_r": float(r),
                    "pearson_p": float(p_p),
                }
            )

    res = pd.DataFrame(rows)
    res["spearman_p_fdr_bh"] = np.nan
    res["significant_fdr_0_05"] = False
    for stratum in STRATA:
        m = (res["diabetes_group"] == stratum) & res["spearman_p"].notna()
        if m.any():
            reject, p_fdr, _, _ = multipletests(
                res.loc[m, "spearman_p"].to_numpy(), alpha=0.05, method="fdr_bh"
            )
            res.loc[m, "spearman_p_fdr_bh"] = p_fdr
            res.loc[m, "significant_fdr_0_05"] = reject
    return res


def differences_by_group(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for stratum in STRATA:
        sdf = df[(df[STRATUM_COL] == stratum) & df[TARGET].notna()].copy()
        for cutoff in CUTOFFS:
            sdf["impaired"] = sdf[TARGET] < cutoff
            for col, desc in DIFF_METRICS.items():
                sub = sdf[[col, "impaired"]].dropna()
                imp = sub.loc[sub["impaired"], col]
                unimp = sub.loc[~sub["impaired"], col]
                n_imp, n_unimp = len(imp), len(unimp)
                row = {
                    "diabetes_group": stratum,
                    "cutoff": cutoff,
                    "metric": col,
                    "description": desc,
                    "n_impaired": n_imp,
                    "n_unimpaired": n_unimp,
                    "median_impaired": float(imp.median()) if n_imp else np.nan,
                    "median_unimpaired": float(unimp.median()) if n_unimp else np.nan,
                    "mann_whitney_u": np.nan,
                    "p_value": np.nan,
                    "rank_biserial": np.nan,
                }
                if n_imp >= 3 and n_unimp >= 3:
                    u, p = mannwhitneyu(imp, unimp, alternative="two-sided")
                    row["mann_whitney_u"] = float(u)
                    row["p_value"] = float(p)
                    row["rank_biserial"] = float(2.0 * u / (n_imp * n_unimp) - 1.0)
                rows.append(row)

    res = pd.DataFrame(rows)
    res["p_value_fdr_bh"] = np.nan
    res["significant_fdr_0_05"] = False
    for stratum in STRATA:
        for cutoff in CUTOFFS:
            m = (
                (res["diabetes_group"] == stratum)
                & (res["cutoff"] == cutoff)
                & res["p_value"].notna()
            )
            if m.any():
                reject, p_fdr, _, _ = multipletests(
                    res.loc[m, "p_value"].to_numpy(), alpha=0.05, method="fdr_bh"
                )
                res.loc[m, "p_value_fdr_bh"] = p_fdr
                res.loc[m, "significant_fdr_0_05"] = reject
    return res


def main() -> None:
    df = load_data()
    pd.set_option("display.width", 240)
    pd.set_option("display.max_columns", 40)

    print("Diabetes group sizes:")
    print(df[STRATUM_COL].value_counts().reindex(STRATA).to_string())
    print()

    corr = correlation_by_group(df)
    corr.to_csv(CORR_OUT, index=False)
    print(f"Wrote: {CORR_OUT}\n")
    print("===== PART A: Spearman correlation (MoCA vs PA metric), by diabetes group =====")
    for stratum in STRATA:
        print(f"\n--- {SHORT[stratum]} ---")
        block = corr[corr["diabetes_group"] == stratum].sort_values("spearman_p")
        print(
            block[
                ["metric", "n", "spearman_rho", "spearman_p", "spearman_p_fdr_bh", "significant_fdr_0_05"]
            ].to_string(index=False)
        )

    diff = differences_by_group(df)
    diff.to_csv(DIFF_OUT, index=False)
    print(f"\n\nWrote: {DIFF_OUT}\n")
    print("===== PART B: impaired vs unimpaired (Mann-Whitney U), by diabetes group & cutoff =====")
    for stratum in STRATA:
        for cutoff in CUTOFFS:
            block = diff[(diff["diabetes_group"] == stratum) & (diff["cutoff"] == cutoff)]
            n_imp = int(block["n_impaired"].iloc[0])
            n_unimp = int(block["n_unimpaired"].iloc[0])
            print(
                f"\n--- {SHORT[stratum]} | cutoff {cutoff} "
                f"(impaired n={n_imp} vs unimpaired n={n_unimp}) ---"
            )
            print(
                block[
                    [
                        "metric",
                        "median_impaired",
                        "median_unimpaired",
                        "p_value",
                        "p_value_fdr_bh",
                        "significant_fdr_0_05",
                        "rank_biserial",
                    ]
                ].to_string(index=False)
            )


if __name__ == "__main__":
    main()
