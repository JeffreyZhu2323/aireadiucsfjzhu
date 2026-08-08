"""
Merge wearable activity metrics into the person-level glucose aggregate table.

Starts from the existing per-person glucose-spike aggregate table
(outputs/full_person_aggregates.csv: one row per person_id, with glucose-spike
aggregates + moca_score + study-group/age meta) and left-merges the wearable
activity metrics onto person_id, producing a single combined per-person table.

Wearable metrics merged in (all on person_id):

  Physical activity (from the wearable_activity_monitor manifest):
    - average_daily_activity          mean daily steps
    - average_active_calories_kcal     mean active calories / day

  Movement (from master_patient_features.csv):
    - total_movement                   total movement over study period
    - night_movement                   nighttime movement
    - day_movement                     total_movement - night_movement (derived)

  Cleaned wearable vitals (from wearable_stream_means.csv; sentinel-cleaned means):
    - mean_heartrate_bpm
    - mean_oxygen_saturation_pct
    - mean_stress_level
    - mean_respiratory_rate_bpm
    - mean_sleep_hours

Note: physical_activity_num_days and the *_n_valid / *_n_nights coverage-count
columns are intentionally excluded.

Overwrites outputs/full_person_aggregates.csv in place.

Run from repo root:
    .venv/bin/python src/merge_wearable_into_full_person_aggregates.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path.home() / "aireadi"
FULL_PATH = ROOT / "outputs" / "full_person_aggregates.csv"
STREAM_MEANS_PATH = ROOT / "outputs" / "wearable_stream_means.csv"
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

# Existing (glucose-block) columns kept up front, in order.
GLUCOSE_COLS = [
    "person_id",
    "avg_cgm_per_spike_mg_dl",
    "avg_spike_duration_minutes",
    "avg_spikes_per_day",
    "total_spikes",
    "days_tracked",
]
# Cognitive / meta columns kept at the end, in order.
TAIL_COLS = ["moca_score", "meta_study_group", "meta_age"]

# Wearable columns inserted between the glucose block and the tail.
PA_COLS = ["average_daily_activity", "average_active_calories_kcal"]
MOVEMENT_COLS = ["total_movement", "night_movement", "day_movement"]
STREAM_COLS = [
    "mean_heartrate_bpm",
    "mean_oxygen_saturation_pct",
    "mean_stress_level",
    "mean_respiratory_rate_bpm",
    "mean_sleep_hours",
]


def _to_person_id(df: pd.DataFrame, src: str) -> pd.DataFrame:
    """Return df with an integer ``person_id`` merge key derived from ``src``."""
    out = df.copy()
    out["person_id"] = pd.to_numeric(out[src], errors="coerce").astype("Int64")
    return out


def main() -> None:
    full = pd.read_csv(FULL_PATH)
    full["person_id"] = pd.to_numeric(full["person_id"], errors="coerce").astype("Int64")
    n0, cols0 = len(full), list(full.columns)

    # ---- physical activity from the manifest ----
    manifest = pd.read_csv(ACTIVITY_MANIFEST, sep="\t")
    manifest = _to_person_id(manifest, "person_id")
    pa = manifest[["person_id", *PA_COLS]].copy()

    # ---- movement from master_patient_features ----
    feats = pd.read_csv(FEATURES_PATH)
    feats = _to_person_id(feats, "participant_id")
    feats["total_movement"] = pd.to_numeric(feats["total_movement"], errors="coerce")
    feats["night_movement"] = pd.to_numeric(feats["night_movement"], errors="coerce")
    feats["day_movement"] = feats["total_movement"] - feats["night_movement"]
    movement = feats[["person_id", *MOVEMENT_COLS]].copy()

    # ---- cleaned wearable vitals ----
    streams = pd.read_csv(STREAM_MEANS_PATH)
    streams = _to_person_id(streams, "participant_id")
    stream_cols_present = [c for c in STREAM_COLS if c in streams.columns]
    vitals = streams[["person_id", *stream_cols_present]].copy()

    merged = (
        full.merge(pa, on="person_id", how="left")
        .merge(movement, on="person_id", how="left")
        .merge(vitals, on="person_id", how="left")
    )

    ordered = [
        *GLUCOSE_COLS,
        *PA_COLS,
        *MOVEMENT_COLS,
        *stream_cols_present,
        *TAIL_COLS,
    ]
    ordered = [c for c in ordered if c in merged.columns]
    extra = [c for c in merged.columns if c not in ordered]
    merged = merged[ordered + extra]

    assert len(merged) == n0, f"row count changed: {n0} -> {len(merged)}"

    merged.to_csv(FULL_PATH, index=False)

    added = [c for c in merged.columns if c not in cols0]
    print(f"Wrote: {FULL_PATH}  (rows={len(merged)}, cols={len(merged.columns)})")
    print(f"Added {len(added)} wearable columns: {added}")
    cov = merged[added].notna().mean().mul(100).round(1)
    print("\nNon-null coverage of added columns (%):")
    print(cov.to_string())


if __name__ == "__main__":
    main()
