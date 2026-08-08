"""
Same design as diabetes_group_differences.py, for the 7 wearable mean features
plus one new engineered feature: mean_cgm_overall_mg_dl -- each person's mean
glucose across ALL CGM readings (not just readings inside a detected spike
episode). Cleaning matches src/glucose_spike_detection.py and
src/glucose_summary_stratified.py: Dexcom "Low"/"High" sentinel strings are
mapped to 39/401 mg/dL before averaging.

Wearable features: average_daily_activity, average_active_calories_kcal,
mean_heartrate_bpm, mean_oxygen_saturation_pct, mean_stress_level,
mean_respiratory_rate_bpm, mean_sleep_hours.

See diabetes_group_differences.py for full methodology notes.

Outputs:
  outputs/diabetes_group_differences_wearable_plus_mean_cgm.csv
  outputs/diabetes_group_differences_wearable_plus_mean_cgm_dunn.csv

Run from repo root:
    .venv/bin/python src/diabetes_group_differences_wearable_plus_mean_cgm.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import diabetes_group_differences as base

ROOT = Path.home() / "aireadi"
GLUCOSE_RAW = ROOT / "outputs" / "glucose_raw.csv"
FULL_PATH = ROOT / "outputs" / "full_person_aggregates.csv"
OUT_OMNIBUS = ROOT / "outputs" / "diabetes_group_differences_wearable_plus_mean_cgm.csv"
OUT_DUNN = ROOT / "outputs" / "diabetes_group_differences_wearable_plus_mean_cgm_dunn.csv"

NEW_COL = "mean_cgm_overall_mg_dl"

WEARABLE_MEAN_METRICS = {
    "average_daily_activity": base.WEARABLE_METRICS["average_daily_activity"],
    "average_active_calories_kcal": base.WEARABLE_METRICS["average_active_calories_kcal"],
    "mean_heartrate_bpm": base.WEARABLE_METRICS["mean_heartrate_bpm"],
    "mean_oxygen_saturation_pct": base.WEARABLE_METRICS["mean_oxygen_saturation_pct"],
    "mean_stress_level": base.WEARABLE_METRICS["mean_stress_level"],
    "mean_respiratory_rate_bpm": base.WEARABLE_METRICS["mean_respiratory_rate_bpm"],
    "mean_sleep_hours": base.WEARABLE_METRICS["mean_sleep_hours"],
}


def compute_mean_cgm_overall() -> pd.DataFrame:
    readings = pd.read_csv(GLUCOSE_RAW, usecols=["person_id", "glucose_mg_dl"], dtype=str)
    readings["glucose_mg_dl"] = readings["glucose_mg_dl"].replace(
        {"Low": 39, "LOW": 39, "High": 401, "HIGH": 401}
    )
    readings["glucose_mg_dl"] = pd.to_numeric(readings["glucose_mg_dl"], errors="coerce")
    readings["person_id"] = pd.to_numeric(readings["person_id"], errors="coerce").astype("Int64")
    readings = readings.dropna(subset=["person_id", "glucose_mg_dl"])
    person_mean = (
        readings.groupby("person_id")["glucose_mg_dl"].mean().rename(NEW_COL).reset_index()
    )
    return person_mean


def main() -> None:
    print("Computing overall mean CGM per person from raw readings...")
    mean_cgm = compute_mean_cgm_overall()

    full = pd.read_csv(FULL_PATH)
    full["person_id"] = pd.to_numeric(full["person_id"], errors="coerce").astype("Int64")
    df = full.merge(mean_cgm, on="person_id", how="left")

    base.METRICS = {
        NEW_COL: ("glucose", "Mean CGM overall, all readings (mg/dL)"),
        **{m: ("wearable", d) for m, d in WEARABLE_MEAN_METRICS.items()},
    }

    # Reuse base.run()'s body by temporarily swapping the module-level frame it reads from.
    orig_read_csv = pd.read_csv
    tmp_path = ROOT / "outputs" / "_tmp_wearable_plus_mean_cgm_input.csv"
    df.to_csv(tmp_path, index=False)
    base.DATA_PATH = tmp_path
    try:
        omnibus, dunn = base.run()
    finally:
        tmp_path.unlink(missing_ok=True)

    OUT_OMNIBUS.parent.mkdir(parents=True, exist_ok=True)
    omnibus.to_csv(OUT_OMNIBUS, index=False)
    dunn.to_csv(OUT_DUNN, index=False)

    pd.set_option("display.width", 240)
    pd.set_option("display.max_columns", 60)

    print(f"Wrote: {OUT_OMNIBUS}")
    print(f"Wrote: {OUT_DUNN}\n")

    print("===== Omnibus (Kruskal-Wallis): wearable means + overall mean CGM =====")
    show_cols = [
        "family",
        "metric",
        "median_healthy",
        "median_pre_diabetes",
        "median_oral_med",
        "median_insulin",
        "kruskal_H",
        "kruskal_p",
        "kruskal_p_fdr_bh",
        "epsilon_squared",
        "significant_fdr_0_05",
    ]
    with pd.option_context("display.float_format", lambda v: f"{v:,.4g}"):
        print(omnibus[show_cols].to_string(index=False))

    if not dunn.empty:
        print("\n===== Dunn's post-hoc (BH-adjusted) for FDR-significant metrics =====")
        dshow = dunn[
            ["metric", "group_a", "group_b", "z", "p_value", "p_value_bh", "significant_bh_0_05"]
        ]
        with pd.option_context("display.float_format", lambda v: f"{v:,.4g}"):
            print(dshow.to_string(index=False))


if __name__ == "__main__":
    main()
