"""
Age x Diabetes-status CROSS-STRATIFIED repeat of the two MoCA / physical-activity
analyses.

Identical design to src/moca_pa_age_stratified.py and
src/moca_pa_diabetes_stratified.py, but strata are the 12 combinations of:
  age band : 40-54 / 55-69 / 70+
  diabetes : healthy / pre_diabetes / oral_med / insulin_dep

Within each (age band x diabetes group) cell we run, independently:

  PART A -- Spearman (primary) + Pearson correlation of MoCA vs PA metrics,
            FDR-BH corrected across metrics within the cell.
            Output: outputs/moca_pa_correlation_by_age_diabetes.csv

  PART B -- impaired (moca < cutoff) vs unimpaired Mann-Whitney U for
            total_movement, physical_activity_num_days,
            average_active_calories_kcal, average_daily_activity
            at cutoffs 26 / 24 / 23 ; FDR-BH across the 4 metrics within
            each (cell, cutoff).
            Output: outputs/moca_impairment_activity_by_age_diabetes.csv

NOTE: cross-stratifying fragments the sample (12 cells); the impaired-vs-unimpaired
tests in particular are underpowered (small impaired n at strict cutoffs). Cells
with <3 in either group are reported with NaN test stats.

Run from repo root: .venv/bin/python src/moca_pa_age_diabetes_stratified.py
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
CORR_OUT = ROOT / "outputs" / "moca_pa_correlation_by_age_diabetes.csv"
DIFF_OUT = ROOT / "outputs" / "moca_impairment_activity_by_age_diabetes.csv"

TARGET = "moca_score"
AGE = "meta_age"
AGE_BINS = [40, 55, 70, 200]
AGE_LABELS = ["40-54", "55-69", "70+"]

DG_SHORT = {
    "healthy": "healthy",
    "pre_diabetes_lifestyle_controlled": "pre_diab",
    "oral_medication_and_or_non_insulin_injectable_medication_controlled": "oral_med",
    "insulin_dependent": "insulin_dep",
}
DG_ORDER = ["healthy", "pre_diab", "oral_med", "insulin_dep"]

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
    df[AGE] = pd.to_numeric(df[AGE], errors="coerce")
    for col in set(CORR_METRICS) | set(DIFF_METRICS):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["age_stratum"] = pd.cut(df[AGE], bins=AGE_BINS, right=False, labels=AGE_LABELS)
    df["diabetes_group"] = df["meta_study_group"].map(DG_SHORT)
    return df


def iter_cells():
    for age in AGE_LABELS:
        for dg in DG_ORDER:
            yield age, dg


def correlation_by_cell(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for age, dg in iter_cells():
        cell = df[(df["age_stratum"] == age) & (df["diabetes_group"] == dg)]
        for col, desc in CORR_METRICS.items():
            pair = cell[[TARGET, col]].dropna()
            n = len(pair)
            rec = {
                "age_stratum": age,
                "diabetes_group": dg,
                "metric": col,
                "description": desc,
                "n": n,
                "spearman_rho": np.nan,
                "spearman_p": np.nan,
                "pearson_r": np.nan,
                "pearson_p": np.nan,
            }
            if n >= 5:
                rho, p_s = spearmanr(pair[TARGET], pair[col])
                r, p_p = pearsonr(pair[TARGET], pair[col])
                rec.update(
                    spearman_rho=float(rho),
                    spearman_p=float(p_s),
                    pearson_r=float(r),
                    pearson_p=float(p_p),
                )
            rows.append(rec)

    res = pd.DataFrame(rows)
    res["spearman_p_fdr_bh"] = np.nan
    res["significant_fdr_0_05"] = False
    for age, dg in iter_cells():
        m = (
            (res["age_stratum"] == age)
            & (res["diabetes_group"] == dg)
            & res["spearman_p"].notna()
        )
        if m.any():
            reject, p_fdr, _, _ = multipletests(
                res.loc[m, "spearman_p"].to_numpy(), alpha=0.05, method="fdr_bh"
            )
            res.loc[m, "spearman_p_fdr_bh"] = p_fdr
            res.loc[m, "significant_fdr_0_05"] = reject
    return res


def differences_by_cell(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for age, dg in iter_cells():
        cell = df[
            (df["age_stratum"] == age)
            & (df["diabetes_group"] == dg)
            & df[TARGET].notna()
        ].copy()
        for cutoff in CUTOFFS:
            cell["impaired"] = cell[TARGET] < cutoff
            for col, desc in DIFF_METRICS.items():
                sub = cell[[col, "impaired"]].dropna()
                imp = sub.loc[sub["impaired"], col]
                unimp = sub.loc[~sub["impaired"], col]
                n_imp, n_unimp = len(imp), len(unimp)
                rec = {
                    "age_stratum": age,
                    "diabetes_group": dg,
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
                    rec["mann_whitney_u"] = float(u)
                    rec["p_value"] = float(p)
                    rec["rank_biserial"] = float(2.0 * u / (n_imp * n_unimp) - 1.0)
                rows.append(rec)

    res = pd.DataFrame(rows)
    res["p_value_fdr_bh"] = np.nan
    res["significant_fdr_0_05"] = False
    for age, dg in iter_cells():
        for cutoff in CUTOFFS:
            m = (
                (res["age_stratum"] == age)
                & (res["diabetes_group"] == dg)
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
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)

    print("Cell sizes (age x diabetes):")
    print(
        pd.crosstab(df["age_stratum"], df["diabetes_group"]).reindex(columns=DG_ORDER).to_string()
    )
    print()

    corr = correlation_by_cell(df)
    corr.to_csv(CORR_OUT, index=False)
    print(f"Wrote: {CORR_OUT}")

    diff = differences_by_cell(df)
    diff.to_csv(DIFF_OUT, index=False)
    print(f"Wrote: {DIFF_OUT}\n")

    # ---- PART A summary: significant correlations after FDR ----
    print("===== PART A: Spearman correlations significant after FDR (per cell) =====")
    sigA = corr[corr["significant_fdr_0_05"]]
    if len(sigA):
        print(
            sigA[
                ["age_stratum", "diabetes_group", "metric", "n", "spearman_rho", "spearman_p_fdr_bh"]
            ].to_string(index=False)
        )
    else:
        print("  (none survive FDR in any cell)")

    # strongest |rho| per cell regardless of significance
    print("\n-- Strongest |rho| per cell (for context) --")
    corr_valid = corr[corr["spearman_rho"].notna()].copy()
    corr_valid["abs_rho"] = corr_valid["spearman_rho"].abs()
    idx = corr_valid.groupby(["age_stratum", "diabetes_group"])["abs_rho"].idxmax()
    top = corr_valid.loc[idx].sort_values(["age_stratum", "diabetes_group"])
    print(
        top[
            ["age_stratum", "diabetes_group", "metric", "n", "spearman_rho", "spearman_p", "spearman_p_fdr_bh"]
        ].to_string(index=False)
    )

    # ---- PART B summary: significant differences after FDR + raw-significant ----
    print("\n\n===== PART B: impaired vs unimpaired (Mann-Whitney) =====")
    sigB = diff[diff["significant_fdr_0_05"]]
    print("Significant after FDR (per cell x cutoff):")
    if len(sigB):
        print(
            sigB[
                [
                    "age_stratum",
                    "diabetes_group",
                    "cutoff",
                    "metric",
                    "n_impaired",
                    "n_unimpaired",
                    "p_value",
                    "p_value_fdr_bh",
                    "rank_biserial",
                ]
            ].to_string(index=False)
        )
    else:
        print("  (none survive FDR in any cell)")

    rawB = diff[(diff["p_value"] < 0.05) & (~diff["significant_fdr_0_05"])]
    print("\nRaw p < 0.05 but NOT significant after FDR (suggestive only):")
    if len(rawB):
        print(
            rawB[
                [
                    "age_stratum",
                    "diabetes_group",
                    "cutoff",
                    "metric",
                    "n_impaired",
                    "n_unimpaired",
                    "p_value",
                    "p_value_fdr_bh",
                    "rank_biserial",
                ]
            ].to_string(index=False)
        )
    else:
        print("  (none)")


if __name__ == "__main__":
    main()
