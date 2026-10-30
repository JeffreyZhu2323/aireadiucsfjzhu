"""
Add per-person CGM range features to outputs/cohort/full_person_aggregates.csv.

Readings are classified into five mutually exclusive glycaemic ranges:

    sev_hypo   Severe hypoglycaemia      < 54 mg/dL
    mod_hypo   Moderate hypoglycaemia    54 - 69 mg/dL
    normal     Normal                    70 - 180 mg/dL
    mod_hyper  Moderate hyperglycaemia   181 - 250 mg/dL
    sev_hyper  Severe hyperglycaemia     > 250 mg/dL

Five metrics are emitted per range (25 columns):

    cgm_<range>_minutes          total time in range, minutes
    cgm_<range>_frac_time        share of the person's valid readings in range
    cgm_<range>_mean_mg_dl       mean glucose across readings in range
    cgm_<range>_windows_per_day  contiguous episodes in range, per tracked day
    cgm_<range>_sustained_windows_per_day
                                 as above, counting only episodes lasting at
                                 least MIN_EPISODE_MINUTES

Five data-quality columns are emitted alongside them, so unreliable
participants can be filtered by rule rather than by hardcoded person_id:

    cgm_n_readings           usable readings after cleaning and dedupe
    cgm_low_sentinel_frac    share of readings reported as Dexcom "Low"
    cgm_high_sentinel_frac   share of readings reported as Dexcom "High"
    cgm_short_wear           fewer than one full day of readings
    cgm_quality_suspect      short wear, or a sensor implausibly pinned at the floor

The two sentinel fractions mean different things and are deliberately left as
raw numbers rather than folded into one flag. A high "Low" share is a
malfunction signal: three participants here are 77-100% "Low" across only
28-75 readings, and nobody is genuinely below 54 mg/dL for two days. A high
"High" share is instead a CENSORING signal: five participants exceed 50%
"High" but have full 10-11 day wear, so their severe hyperglycaemia is real
while their cgm_sev_hyper_mean_mg_dl is capped at the 401 ceiling and
understates true exposure. Filtering those people out would discard real
signal, so cgm_quality_suspect does not key on the "High" share -- threshold
cgm_high_sentinel_frac yourself if censoring matters to your analysis.

Time accounting. The Dexcom G6 cadence here is a near-perfect 5 minutes
(99.94% of consecutive gaps are exactly 5.00 min), so each valid reading is
credited with EPOCH_MINUTES = 5 minutes of time. Counting readings and
integrating timestamps therefore agree to within rounding, and this avoids
attributing multi-hour sensor dropouts to whichever range preceded them.
frac_time is the standard consensus time-in-range quantity: readings in range
divided by all valid readings.

Windows. A window is a maximal run of consecutive readings in the same range.
A run is broken by a range change or by a sampling gap exceeding
MAX_GAP_MINUTES, so a sensor dropout never bridges two separate episodes.
Runs are counted strictly: a single reading crossing a range boundary does
split a window. This differs on purpose from the dip tolerance in
src/glucose_spike_detection.py, where a lone sub-threshold reading is folded
back into a spike -- here the five ranges are exhaustive, so every boundary
crossing is a genuine range transition rather than noise around one cutoff.

Because strict counting lets glucose hovering on a boundary register as many
separate episodes, the sustained_windows_per_day variant counts only episodes
of at least MIN_EPISODE_MINUTES (15 min = 3 readings), the usual clinical
convention for calling an excursion an event. It filters which episodes are
COUNTED and deliberately does not alter minutes, frac_time, or mean_mg_dl --
time spent in a brief excursion was still time spent in that range. The two
columns answer different questions: how often glucose touches a range, versus
how often it settles there.

Cleaning matches src/glucose_spike_detection.py and
src/add_mean_cgm_overall_to_full_person_aggregates.py: the Dexcom "Low"/"High"
out-of-sensor-range sentinels are mapped to 39/401 mg/dL, which places them in
severe hypoglycaemia and severe hyperglycaemia respectively.

Inserted after days_tracked so the CGM-derived block stays contiguous.

Overwrites outputs/cohort/full_person_aggregates.csv in place.

Run from repo root:
    .venv/bin/python src/add_cgm_ranges_to_full_person_aggregates.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GLUCOSE_RAW = ROOT / "outputs" / "glycemic_phenotype" / "glucose_raw.csv"
FULL_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"

AFTER_COL = "days_tracked"

# Each reading represents one 5-minute sampling epoch.
EPOCH_MINUTES = 5.0
# A gap longer than this breaks a window; matches src/glucose_spike_detection.py.
MAX_GAP_MINUTES = 10.0

RANGE_KEYS = ["sev_hypo", "mod_hypo", "normal", "mod_hyper", "sev_hyper"]
# Half-open [lo, hi) edges. All glucose values in this export are integers, so
# these reproduce the inclusive bounds in the docstring exactly; the half-open
# form additionally leaves no gap for a hypothetical fractional value.
RANGE_EDGES = [-np.inf, 54, 70, 181, 251, np.inf]

METRIC_SUFFIXES = [
    "minutes",
    "frac_time",
    "mean_mg_dl",
    "windows_per_day",
    "sustained_windows_per_day",
]
RANGE_COLS = [f"cgm_{key}_{suffix}" for key in RANGE_KEYS for suffix in METRIC_SUFFIXES]

# Clinical convention counts an excursion as an event only once it persists.
MIN_EPISODE_MINUTES = 15.0
MIN_EPISODE_READINGS = int(round(MIN_EPISODE_MINUTES / EPOCH_MINUTES))

# One full day of readings at the 5-minute cadence.
MIN_READINGS_FULL_DAY = 288
# Sustained floor-pinning above this share is a failed sensor, not physiology.
MAX_PLAUSIBLE_LOW_SENTINEL_FRAC = 0.5

QUALITY_COLS = [
    "cgm_n_readings",
    "cgm_low_sentinel_frac",
    "cgm_high_sentinel_frac",
    "cgm_short_wear",
    "cgm_quality_suspect",
]

# Quality context is written before the features it qualifies.
NEW_COLS = QUALITY_COLS + RANGE_COLS


def load_readings() -> pd.DataFrame:
    """Load, clean, and range-classify every CGM reading."""
    readings = pd.read_csv(
        GLUCOSE_RAW,
        usecols=["person_id", "timestamp_utc", "glucose_mg_dl"],
        dtype=str,
    )
    n_raw = len(readings)

    # Record sentinel provenance before the strings are replaced by 39/401.
    readings["is_low_sentinel"] = readings["glucose_mg_dl"].isin(["Low", "LOW"])
    readings["is_high_sentinel"] = readings["glucose_mg_dl"].isin(["High", "HIGH"])

    readings["glucose_mg_dl"] = pd.to_numeric(
        readings["glucose_mg_dl"].replace({"Low": 39, "LOW": 39, "High": 401, "HIGH": 401}),
        errors="coerce",
    )
    readings["person_id"] = pd.to_numeric(readings["person_id"], errors="coerce")
    readings["timestamp_utc"] = pd.to_datetime(
        readings["timestamp_utc"], utc=True, errors="coerce", format="ISO8601"
    )
    readings = readings.dropna(subset=["person_id", "timestamp_utc", "glucose_mg_dl"])
    readings["person_id"] = readings["person_id"].astype(int)

    # Dedupe in file order, before sorting, so which row survives does not
    # depend on sort stability. Mirrored in
    # src/add_mean_cgm_overall_to_full_person_aggregates.py.
    n_dedup = readings.duplicated(["person_id", "timestamp_utc"]).sum()
    readings = readings.drop_duplicates(["person_id", "timestamp_utc"], keep="first")
    readings = readings.sort_values(["person_id", "timestamp_utc"]).reset_index(drop=True)

    readings["range"] = pd.cut(
        readings["glucose_mg_dl"], bins=RANGE_EDGES, right=False, labels=RANGE_KEYS
    )
    assert readings["range"].notna().all(), "range bins failed to cover some readings"

    print(f"  rows read: {n_raw:,}")
    print(f"  dropped as unparseable: {n_raw - len(readings) - n_dedup:,}")
    print(f"  dropped as duplicate (person_id, timestamp_utc): {n_dedup:,}")
    print(f"  usable readings: {len(readings):,} across {readings['person_id'].nunique():,} persons")
    return readings


def build_windows(readings: pd.DataFrame) -> pd.DataFrame:
    """One row per contiguous in-range episode, with its length in readings.

    A window starts whenever the person changes, the range changes, or the
    sampling gap since the previous reading exceeds MAX_GAP_MINUTES. Because
    ``readings`` is sorted and index-reset, each window is a contiguous block
    of rows, so lengths are just the gaps between successive start positions.
    """
    gap_minutes = readings.groupby("person_id")["timestamp_utc"].diff().dt.total_seconds() / 60.0
    starts_window = (
        (readings["person_id"] != readings["person_id"].shift())
        | (readings["range"] != readings["range"].shift())
        | (gap_minutes > MAX_GAP_MINUTES)
    )

    start_pos = np.flatnonzero(starts_window.to_numpy())
    lengths = np.diff(np.append(start_pos, len(readings)))
    return pd.DataFrame(
        {
            "person_id": readings["person_id"].to_numpy()[start_pos],
            "range": pd.Categorical(
                readings["range"].to_numpy()[start_pos], categories=RANGE_KEYS
            ),
            "n_readings": lengths,
        }
    )


def compute_quality_flags(readings: pd.DataFrame) -> pd.DataFrame:
    """Per-person CGM data-quality context, computed post-dedupe."""
    per_person = readings.groupby("person_id").agg(
        cgm_n_readings=("glucose_mg_dl", "size"),
        cgm_low_sentinel_frac=("is_low_sentinel", "mean"),
        cgm_high_sentinel_frac=("is_high_sentinel", "mean"),
    )
    per_person["cgm_short_wear"] = per_person["cgm_n_readings"] < MIN_READINGS_FULL_DAY
    per_person["cgm_quality_suspect"] = per_person["cgm_short_wear"] | (
        per_person["cgm_low_sentinel_frac"] > MAX_PLAUSIBLE_LOW_SENTINEL_FRAC
    )
    return per_person


def compute_range_features(readings: pd.DataFrame) -> pd.DataFrame:
    """Build the 20 per-person range columns."""
    grouped = readings.groupby(["person_id", "range"], observed=False)

    counts = grouped.size().unstack("range")
    means = grouped["glucose_mg_dl"].mean().unstack("range")

    windows = build_windows(readings)
    by_range = windows.groupby(["person_id", "range"], observed=False)
    all_windows = by_range.size().unstack("range").reindex(counts.index).fillna(0)
    sustained_windows = (
        windows[windows["n_readings"] >= MIN_EPISODE_READINGS]
        .groupby(["person_id", "range"], observed=False)
        .size()
        .unstack("range")
        .reindex(counts.index)
        .fillna(0)
    )

    days_tracked = readings.groupby("person_id")["timestamp_utc"].apply(
        lambda ts: ts.dt.tz_convert("UTC").dt.date.nunique()
    )

    minutes = counts * EPOCH_MINUTES
    frac_time = counts.div(counts.sum(axis=1), axis=0)
    windows_per_day = all_windows.div(days_tracked, axis=0)
    sustained_per_day = sustained_windows.div(days_tracked, axis=0)

    out = pd.DataFrame(index=counts.index)
    for key in RANGE_KEYS:
        out[f"cgm_{key}_minutes"] = minutes[key]
        out[f"cgm_{key}_frac_time"] = frac_time[key]
        # Undefined, not zero, when the person never entered the range.
        out[f"cgm_{key}_mean_mg_dl"] = means[key]
        out[f"cgm_{key}_windows_per_day"] = windows_per_day[key]
        out[f"cgm_{key}_sustained_windows_per_day"] = sustained_per_day[key]

    return out.reset_index()


def validate(features: pd.DataFrame) -> None:
    frac_cols = [f"cgm_{key}_frac_time" for key in RANGE_KEYS]
    frac_sum = features[frac_cols].sum(axis=1)
    assert np.allclose(frac_sum, 1.0), f"fractions do not sum to 1 (min {frac_sum.min()})"

    for key in RANGE_KEYS:
        never = features[f"cgm_{key}_minutes"] == 0
        assert features.loc[never, f"cgm_{key}_mean_mg_dl"].isna().all(), (
            f"cgm_{key}_mean_mg_dl should be NaN when the range was never entered"
        )
        assert (features.loc[~never, f"cgm_{key}_windows_per_day"] > 0).all(), (
            f"cgm_{key}_windows_per_day should be positive when time was spent in range"
        )
        # Sustained is a subset of all windows, so it may be 0 where the strict
        # count is not -- every episode can be shorter than the minimum.
        assert (
            features[f"cgm_{key}_sustained_windows_per_day"]
            <= features[f"cgm_{key}_windows_per_day"] + 1e-9
        ).all(), f"cgm_{key}_sustained_windows_per_day exceeds the unfiltered count"
        assert features.loc[never, f"cgm_{key}_sustained_windows_per_day"].eq(0).all(), (
            f"cgm_{key}_sustained_windows_per_day should be 0 when the range was never entered"
        )
    print("  validation passed: fractions sum to 1, empty ranges are NaN-mean, "
          "sustained <= all windows")


def main() -> None:
    # round_trip parsing keeps columns this script does not own bit-identical.
    # pandas' default float reader is up to 1 ULP off from strtod, which would
    # otherwise let untouched values drift as scripts rewrite this file in place.
    full = pd.read_csv(FULL_PATH, float_precision="round_trip")
    full["person_id"] = pd.to_numeric(full["person_id"], errors="coerce").astype("Int64")
    n0 = len(full)

    cols0 = [c for c in full.columns if c not in NEW_COLS]
    full = full[cols0]

    print("Loading CGM readings from glucose_raw.csv...")
    readings = load_readings()

    print("Computing range features...")
    features = compute_range_features(readings)
    validate(features)

    print("Computing data-quality flags...")
    quality = compute_quality_flags(readings).reset_index()
    features = quality.merge(features, on="person_id", how="outer")
    features["person_id"] = features["person_id"].astype("Int64")

    merged = full.merge(features, on="person_id", how="left")
    assert len(merged) == n0, f"row count changed: {n0} -> {len(merged)}"

    ordered: list[str] = []
    for c in cols0:
        ordered.append(c)
        if c == AFTER_COL:
            ordered.extend(NEW_COLS)
    if AFTER_COL not in cols0:
        ordered.extend(NEW_COLS)
    merged = merged[ordered]

    merged.to_csv(FULL_PATH, index=False)

    print(f"\nWrote: {FULL_PATH}  (rows={len(merged)}, cols={len(merged.columns)})")
    covered = merged["cgm_normal_frac_time"].notna()
    print(f"Persons with range features: {covered.sum()}/{len(merged)} ({covered.mean() * 100:.1f}%)")

    print("\n=== Data quality ===")
    for col in ["cgm_short_wear", "cgm_quality_suspect"]:
        flagged = merged[col].fillna(False).astype(bool)
        print(f"{col}: {flagged.sum()}/{len(merged)} flagged")
    censored = merged["cgm_high_sentinel_frac"] > 0.10
    print(f'cgm_high_sentinel_frac > 0.10 (severe-hyper mean censored): {censored.sum()} persons')
    suspect = merged.loc[merged["cgm_quality_suspect"].fillna(False).astype(bool)]
    print(suspect[["person_id", *QUALITY_COLS[:3], "days_tracked"]].round(3).to_string(index=False))

    print("\n=== Range feature summary (across persons) ===")
    summary = merged[RANGE_COLS].describe().T[["count", "mean", "50%", "min", "max"]]
    print(summary.round(3).to_string())


if __name__ == "__main__":
    main()
