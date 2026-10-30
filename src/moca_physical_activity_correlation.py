"""
Test for correlation between MoCA cognitive score and physical-activity metrics.

Physical-activity metrics tested:
  - total_movement              (feature table; total movement over study period)
  - day_movement                (derived: total_movement - night_movement)
  - night_movement              (feature table; nighttime movement)
  - average_daily_activity      (activity manifest; mean daily steps)
  - average_active_calories_kcal(activity manifest; mean active calories)
  - physical_activity_num_days  (activity manifest; days of activity coverage)

For each metric we report:
  - Spearman rho + p           (primary; nonparametric, robust to MoCA's bounded/skewed scale)
  - Pearson r + p              (secondary; linear)
  - age-adjusted partial Spearman (controls for meta_age, which strongly drives MoCA)
  - Benjamini-Hochberg FDR-corrected Spearman p-values across all metrics

Run from repo root: .venv/bin/python src/moca_physical_activity_correlation.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
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
OUT_PATH = ROOT / "outputs" / "cognition" / "moca_physical_activity_correlation.csv"

TARGET = "moca_score"
AGE = "meta_age"

# label -> human description
PA_METRICS = {
    "total_movement": "Total movement over study period",
    "day_movement": "Daytime movement (total - night)",
    "night_movement": "Nighttime movement",
    "average_daily_activity": "Mean daily steps (manifest)",
    "average_active_calories_kcal": "Mean active calories/day (manifest)",
    "physical_activity_num_days": "Days of activity coverage (manifest)",
}


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> tuple[float, float, int]:
    """Spearman partial correlation of x and y controlling for z.

    Computed as the Pearson partial correlation on rank-transformed data.
    Returns (partial_rho, p_value, n).
    """
    d = pd.concat([x, y, z], axis=1).dropna()
    n = len(d)
    if n < 5:
        return float("nan"), float("nan"), n
    rx = d.iloc[:, 0].rank()
    ry = d.iloc[:, 1].rank()
    rz = d.iloc[:, 2].rank()

    rxy = np.corrcoef(rx, ry)[0, 1]
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]

    denom = np.sqrt((1 - rxz**2) * (1 - ryz**2))
    if denom == 0:
        return float("nan"), float("nan"), n
    pr = (rxy - rxz * ryz) / denom

    # t-test on partial correlation, dof = n - 3 (one covariate)
    from scipy.stats import t as t_dist

    dof = n - 3
    if dof <= 0 or abs(pr) >= 1:
        return float(pr), float("nan"), n
    t_stat = pr * np.sqrt(dof / (1 - pr**2))
    p = 2 * t_dist.sf(abs(t_stat), dof)
    return float(pr), float(p), n


def main() -> None:
    feats = pd.read_csv(FEATURES_PATH)
    feats["participant_id"] = pd.to_numeric(feats["participant_id"], errors="coerce")

    manifest = pd.read_csv(ACTIVITY_MANIFEST, sep="\t")
    manifest = manifest.rename(columns={"person_id": "participant_id"})
    manifest_cols = [
        "participant_id",
        "average_daily_activity",
        "average_active_calories_kcal",
        "physical_activity_num_days",
    ]
    df = feats.merge(manifest[manifest_cols], on="participant_id", how="left")

    # derived metric
    df["day_movement"] = pd.to_numeric(df["total_movement"], errors="coerce") - pd.to_numeric(
        df["night_movement"], errors="coerce"
    )

    df[TARGET] = pd.to_numeric(df[TARGET], errors="coerce")
    df[AGE] = pd.to_numeric(df[AGE], errors="coerce")

    rows = []
    for col, desc in PA_METRICS.items():
        s = pd.to_numeric(df[col], errors="coerce")
        pair = pd.concat([df[TARGET], s], axis=1).dropna()
        n = len(pair)
        if n < 5:
            rows.append(
                {
                    "metric": col,
                    "description": desc,
                    "n": n,
                    "spearman_rho": np.nan,
                    "spearman_p": np.nan,
                    "pearson_r": np.nan,
                    "pearson_p": np.nan,
                    "age_adj_partial_spearman_rho": np.nan,
                    "age_adj_partial_p": np.nan,
                    "n_age_adj": 0,
                }
            )
            continue

        rho, p_s = spearmanr(pair[TARGET], pair[col])
        r, p_p = pearsonr(pair[TARGET], pair[col])
        pr, p_pr, n_adj = partial_spearman(df[TARGET], s, df[AGE])

        rows.append(
            {
                "metric": col,
                "description": desc,
                "n": n,
                "spearman_rho": float(rho),
                "spearman_p": float(p_s),
                "pearson_r": float(r),
                "pearson_p": float(p_p),
                "age_adj_partial_spearman_rho": pr,
                "age_adj_partial_p": p_pr,
                "n_age_adj": n_adj,
            }
        )

    results = pd.DataFrame(rows)

    # Benjamini-Hochberg FDR across the Spearman p-values
    mask = results["spearman_p"].notna()
    if mask.any():
        reject, p_fdr, _, _ = multipletests(
            results.loc[mask, "spearman_p"].to_numpy(), alpha=0.05, method="fdr_bh"
        )
        results.loc[mask, "spearman_p_fdr_bh"] = p_fdr
        results.loc[mask, "significant_fdr_0_05"] = reject
    else:
        results["spearman_p_fdr_bh"] = np.nan
        results["significant_fdr_0_05"] = False

    results = results.sort_values("spearman_p", na_position="last").reset_index(drop=True)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(OUT_PATH, index=False)

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 30)
    print(f"Wrote: {OUT_PATH}\n")
    show = results[
        [
            "metric",
            "n",
            "spearman_rho",
            "spearman_p",
            "spearman_p_fdr_bh",
            "significant_fdr_0_05",
            "pearson_r",
            "age_adj_partial_spearman_rho",
            "age_adj_partial_p",
        ]
    ]
    print(show.to_string(index=False))


if __name__ == "__main__":
    main()
