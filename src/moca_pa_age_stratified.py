"""
Age-STRATIFIED repeat of the two MoCA / physical-activity analyses.

Within each age stratum (40-54, 55-69, 70+) we run, independently:

  PART A -- numerical correlation between MoCA score and physical-activity metrics
    * Spearman rho + p (primary, nonparametric)
    * Pearson r + p (secondary)
    * Benjamini-Hochberg FDR correction across metrics within the stratum
    Output: outputs/moca_pa_correlation_by_age.csv

  PART B -- impaired vs unimpaired group differences (Mann-Whitney U)
    cutoffs 26 / 24 / 23 ; impaired = moca_score < cutoff
    metrics: total_movement, physical_activity_num_days,
             average_active_calories_kcal, average_daily_activity
    * group medians, Mann-Whitney U, two-sided p, rank-biserial effect size
    * FDR-BH correction across the 4 metrics within each (stratum, cutoff)
    Output: outputs/moca_impairment_activity_by_age.csv

Stratifying by age removes age as a confounder by analysing within
age-homogeneous subgroups (rather than adjusting for it in a model).

Run from repo root: .venv/bin/python src/moca_pa_age_stratified.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path.home() / "aireadi"
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
CORR_OUT = ROOT / "outputs" / "moca_pa_correlation_by_age.csv"
DIFF_OUT = ROOT / "outputs" / "moca_impairment_activity_by_age.csv"

TARGET = "moca_score"
AGE = "meta_age"

AGE_BINS = [40, 55, 70, 200]
AGE_LABELS = ["40-54", "55-69", "70+"]

CUTOFFS = [26, 24, 23]

# full PA metric set for the correlation analysis
CORR_METRICS = {
    "total_movement": "Total movement over study period",
    "day_movement": "Daytime movement (total - night)",
    "night_movement": "Nighttime movement",
    "average_daily_activity": "Mean daily steps (manifest)",
    "average_active_calories_kcal": "Mean active calories/day (manifest)",
    "physical_activity_num_days": "Days of activity coverage (manifest)",
}

# the four metrics requested for the impaired-vs-unimpaired difference test
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
    df[AGE] = pd.to_numeric(df[AGE], errors="coerce")
    for col in set(CORR_METRICS) | set(DIFF_METRICS):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["age_stratum"] = pd.cut(df[AGE], bins=AGE_BINS, right=False, labels=AGE_LABELS)
    return df


def correlation_by_age(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for stratum in AGE_LABELS:
        sdf = df[df["age_stratum"] == stratum]
        for col, desc in CORR_METRICS.items():
            pair = sdf[[TARGET, col]].dropna()
            n = len(pair)
            if n < 5:
                rows.append(
                    {
                        "age_stratum": stratum,
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
                    "age_stratum": stratum,
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
    for stratum in AGE_LABELS:
        m = (res["age_stratum"] == stratum) & res["spearman_p"].notna()
        if m.any():
            reject, p_fdr, _, _ = multipletests(
                res.loc[m, "spearman_p"].to_numpy(), alpha=0.05, method="fdr_bh"
            )
            res.loc[m, "spearman_p_fdr_bh"] = p_fdr
            res.loc[m, "significant_fdr_0_05"] = reject
    return res


def differences_by_age(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for stratum in AGE_LABELS:
        sdf = df[(df["age_stratum"] == stratum) & df[TARGET].notna()].copy()
        for cutoff in CUTOFFS:
            sdf["impaired"] = sdf[TARGET] < cutoff
            for col, desc in DIFF_METRICS.items():
                sub = sdf[[col, "impaired"]].dropna()
                imp = sub.loc[sub["impaired"], col]
                unimp = sub.loc[~sub["impaired"], col]
                n_imp, n_unimp = len(imp), len(unimp)
                row = {
                    "age_stratum": stratum,
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
    for stratum in AGE_LABELS:
        for cutoff in CUTOFFS:
            m = (
                (res["age_stratum"] == stratum)
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

    print("Age stratum sizes:")
    print(df["age_stratum"].value_counts().reindex(AGE_LABELS).to_string())
    print()

    corr = correlation_by_age(df)
    corr.to_csv(CORR_OUT, index=False)
    print(f"Wrote: {CORR_OUT}\n")
    print("===== PART A: Spearman correlation (MoCA vs PA metric), by age =====")
    for stratum in AGE_LABELS:
        print(f"\n--- Age {stratum} ---")
        block = corr[corr["age_stratum"] == stratum].sort_values("spearman_p")
        print(
            block[
                ["metric", "n", "spearman_rho", "spearman_p", "spearman_p_fdr_bh", "significant_fdr_0_05"]
            ].to_string(index=False)
        )

    diff = differences_by_age(df)
    diff.to_csv(DIFF_OUT, index=False)
    print(f"\n\nWrote: {DIFF_OUT}\n")
    print("===== PART B: impaired vs unimpaired (Mann-Whitney U), by age & cutoff =====")
    for stratum in AGE_LABELS:
        for cutoff in CUTOFFS:
            block = diff[(diff["age_stratum"] == stratum) & (diff["cutoff"] == cutoff)]
            n_imp = int(block["n_impaired"].iloc[0])
            n_unimp = int(block["n_unimpaired"].iloc[0])
            print(f"\n--- Age {stratum} | cutoff {cutoff} (impaired n={n_imp} vs unimpaired n={n_unimp}) ---")
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
