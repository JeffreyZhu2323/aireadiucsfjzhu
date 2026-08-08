import pandas as pd

from config import OUTDIR

CACHE_PATH = OUTDIR / "glucose_raw.csv"
SPIKES_PATH = OUTDIR / "glucose_spikes.csv"
PERSON_AGG_PATH = OUTDIR / "spike_aggregates_by_person.csv"
DAILY_PATH = OUTDIR / "spike_daily_by_cgm_day.csv"
SPIKE_THRESHOLD_MG_DL = 140.0
# A single sub-threshold reading sandwiched between two above-threshold
# readings is folded into the spike (tolerance = 1 dip reading), but only
# when the neighbouring samples are within normal CGM cadence.
MAX_DIP_READINGS = 1
MAX_GAP_MINUTES = 10.0


def load_glucose_raw_cache() -> pd.DataFrame:
    """Read the existing official CGM cache from outputs/glucose_raw.csv."""
    if not CACHE_PATH.exists():
        raise FileNotFoundError(
            f"Missing {CACHE_PATH}. This repo now treats partner outputs as official."
        )
    df = pd.read_csv(CACHE_PATH, parse_dates=["timestamp_utc"])
    print(f"Loaded cached glucose dataframe from: {CACHE_PATH}")
    return df


def compute_spike_mask(g: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Return (is_spike, gap_fill) boolean Series aligned to ``g``.

    ``g`` must be sorted by ``person_id`` then ``timestamp_utc``. A reading is a
    spike if it is above threshold, or if it is a single sub-threshold reading
    directly between two above-threshold, time-adjacent samples of the same
    person (tolerance = MAX_DIP_READINGS = 1, both gaps <= MAX_GAP_MINUTES).
    """
    above = g["glucose_mg_dl"] > SPIKE_THRESHOLD_MG_DL
    by_person = above.groupby(g["person_id"])
    prev_above = by_person.shift(1, fill_value=False)
    next_above = by_person.shift(-1, fill_value=False)
    dt_prev_min = g.groupby("person_id")["timestamp_utc"].diff().dt.total_seconds() / 60.0
    dt_next_min = dt_prev_min.shift(-1)
    close_prev = dt_prev_min <= MAX_GAP_MINUTES
    close_next = dt_next_min <= MAX_GAP_MINUTES
    gap_fill = (~above) & prev_above & next_above & close_prev & close_next
    is_spike = above | gap_fill
    return is_spike, gap_fill


def main() -> None:
    glucose_df = load_glucose_raw_cache()

    glucose_df["glucose_mg_dl"] = glucose_df["glucose_mg_dl"].replace(
        {"Low": 39, "LOW": 39, "High": 401, "HIGH": 401}
    )
    glucose_df["glucose_mg_dl"] = pd.to_numeric(glucose_df["glucose_mg_dl"], errors="coerce")
    glucose_df["timestamp_utc"] = pd.to_datetime(glucose_df["timestamp_utc"], utc=True, errors="coerce")
    glucose_df["person_id"] = pd.to_numeric(glucose_df["person_id"], errors="coerce")

    g = glucose_df.dropna(subset=["person_id", "timestamp_utc", "glucose_mg_dl"]).copy()
    g["person_id"] = g["person_id"].astype(int)
    g = g.sort_values(["person_id", "timestamp_utc"]).reset_index(drop=True)

    # Fold in a lone dip (tolerance = MAX_DIP_READINGS = 1) only when both the
    # preceding and following samples are above threshold and time-adjacent.
    is_spike, gap_fill = compute_spike_mask(g)
    g["is_spike"] = is_spike
    g["segment_id"] = g.groupby("person_id")["is_spike"].transform(lambda s: s.ne(s.shift()).cumsum())

    spike_rows = g[g["is_spike"]].copy()
    spikes_df = spike_rows.groupby(["person_id", "segment_id"], as_index=False).agg(
        spike_start_utc=("timestamp_utc", "min"),
        spike_end_utc=("timestamp_utc", "max"),
        mean_blood_sugar_mg_dl=("glucose_mg_dl", "mean"),
        max_blood_sugar_mg_dl=("glucose_mg_dl", "max"),
    )

    spikes_df["duration_minutes"] = (
        (spikes_df["spike_end_utc"] - spikes_df["spike_start_utc"]).dt.total_seconds() / 60.0
    )
    spikes_df["spike_start_date"] = spikes_df["spike_start_utc"].dt.date
    spikes_df["spike_start_time_of_day"] = spikes_df["spike_start_utc"].dt.strftime("%H:%M:%S")

    spikes_df = spikes_df[
        [
            "person_id",
            "spike_start_utc",
            "spike_start_date",
            "spike_start_time_of_day",
            "duration_minutes",
            "mean_blood_sugar_mg_dl",
            "max_blood_sugar_mg_dl",
        ]
    ].sort_values(["person_id", "spike_start_utc"])

    spikes_df.to_csv(SPIKES_PATH, index=False)

    g["cgm_date"] = g["timestamp_utc"].dt.date
    days_tracked = g.groupby("person_id")["cgm_date"].nunique()

    person_agg = spikes_df.groupby("person_id").agg(
        avg_cgm_per_spike_mg_dl=("mean_blood_sugar_mg_dl", "mean"),
        avg_spike_duration_minutes=("duration_minutes", "mean"),
        total_spikes=("spike_start_utc", "size"),
    )
    person_agg = person_agg.reindex(days_tracked.index)
    person_agg["total_spikes"] = person_agg["total_spikes"].fillna(0).astype(int)
    person_agg["days_tracked"] = days_tracked
    person_agg["avg_spikes_per_day"] = person_agg["total_spikes"] / person_agg["days_tracked"]
    person_agg = person_agg.reset_index()[
        [
            "person_id",
            "avg_cgm_per_spike_mg_dl",
            "avg_spike_duration_minutes",
            "avg_spikes_per_day",
            "total_spikes",
            "days_tracked",
        ]
    ].sort_values("person_id")
    person_agg.to_csv(PERSON_AGG_PATH, index=False)

    day_grid = g[["person_id", "cgm_date"]].drop_duplicates().rename(columns={"cgm_date": "date"})
    daily_spikes = spikes_df.groupby(["person_id", "spike_start_date"]).agg(
        mean_duration_minutes=("duration_minutes", "mean"),
        n_spikes=("spike_start_utc", "size"),
    ).reset_index().rename(columns={"spike_start_date": "date"})
    daily = day_grid.merge(daily_spikes, on=["person_id", "date"], how="left")
    daily["n_spikes"] = daily["n_spikes"].fillna(0).astype(int)
    daily = daily[["person_id", "date", "mean_duration_minutes", "n_spikes"]].sort_values(
        ["person_id", "date"]
    )
    daily.to_csv(DAILY_PATH, index=False)

    print(f"Loaded glucose rows: {len(g):,}")
    print(f"Unique people in glucose data: {g['person_id'].nunique(dropna=True)}")
    print(f"Lone dips folded into spikes: {int(gap_fill.sum()):,}")
    print(f"Detected spike episodes: {len(spikes_df):,}")
    print(f"Saved spike dataframe to: {SPIKES_PATH}")
    print(f"Saved per-person aggregates to: {PERSON_AGG_PATH}")
    print(f"Saved daily aggregates to: {DAILY_PATH}")


if __name__ == "__main__":
    main()
