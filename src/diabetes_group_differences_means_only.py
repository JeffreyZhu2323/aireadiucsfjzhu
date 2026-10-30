"""
Same design as diabetes_group_differences.py, restricted to the *mean/average*
glucose and wearable features only (drops total_spikes, total_movement,
night_movement, day_movement -- these are totals/exposure-driven, not means).

See diabetes_group_differences.py for full methodology notes.

Outputs:
  outputs/group_differences/diabetes_group_differences_means_only.csv
  outputs/group_differences/diabetes_group_differences_means_only_dunn.csv

Run from repo root:
    .venv/bin/python src/diabetes_group_differences_means_only.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import diabetes_group_differences as base

ROOT = Path(__file__).resolve().parents[1]
OUT_OMNIBUS = ROOT / "outputs" / "group_differences" / "diabetes_group_differences_means_only.csv"
OUT_DUNN = ROOT / "outputs" / "group_differences" / "diabetes_group_differences_means_only_dunn.csv"

GLUCOSE_MEAN_METRICS = {
    "avg_cgm_per_spike_mg_dl": base.GLUCOSE_METRICS["avg_cgm_per_spike_mg_dl"],
    "avg_spike_duration_minutes": base.GLUCOSE_METRICS["avg_spike_duration_minutes"],
    "avg_spikes_per_day": base.GLUCOSE_METRICS["avg_spikes_per_day"],
}
WEARABLE_MEAN_METRICS = {
    "average_daily_activity": base.WEARABLE_METRICS["average_daily_activity"],
    "average_active_calories_kcal": base.WEARABLE_METRICS["average_active_calories_kcal"],
    "mean_heartrate_bpm": base.WEARABLE_METRICS["mean_heartrate_bpm"],
    "mean_oxygen_saturation_pct": base.WEARABLE_METRICS["mean_oxygen_saturation_pct"],
    "mean_stress_level": base.WEARABLE_METRICS["mean_stress_level"],
    "mean_respiratory_rate_bpm": base.WEARABLE_METRICS["mean_respiratory_rate_bpm"],
    "mean_sleep_hours": base.WEARABLE_METRICS["mean_sleep_hours"],
}

# Patch the module-level METRICS dict that `run()` reads, then reuse it as-is.
base.METRICS = {
    **{m: ("glucose", d) for m, d in GLUCOSE_MEAN_METRICS.items()},
    **{m: ("wearable", d) for m, d in WEARABLE_MEAN_METRICS.items()},
}


def main() -> None:
    omnibus, dunn = base.run()
    OUT_OMNIBUS.parent.mkdir(parents=True, exist_ok=True)
    omnibus.to_csv(OUT_OMNIBUS, index=False)
    dunn.to_csv(OUT_DUNN, index=False)

    pd.set_option("display.width", 240)
    pd.set_option("display.max_columns", 60)

    print(f"Wrote: {OUT_OMNIBUS}")
    print(f"Wrote: {OUT_DUNN}\n")

    print("===== Omnibus (Kruskal-Wallis) by metric, means only =====")
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
