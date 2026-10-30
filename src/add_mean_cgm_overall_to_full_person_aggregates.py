"""
Add mean_cgm_overall_mg_dl to outputs/cohort/full_person_aggregates.csv: each person's
mean glucose across ALL CGM readings (not just readings inside a detected
spike episode). Cleaning matches src/glucose_spike_detection.py and
src/glucose_summary_stratified.py: Dexcom "Low"/"High" sentinel strings are
mapped to 39/401 mg/dL before averaging.

Duplicate (person_id, timestamp_utc) rows are dropped before averaging, keeping
the first occurrence in file order. These are a daylight-saving artifact in the
upstream timestamp conversion, not repeated measurements: all 960 of them fall
on 2024-03-10 and 2025-03-09 -- the two US spring-forward dates in the study
window -- and each of the 80 affected participants contributes exactly 24 rows,
i.e. 12 timestamps x 2, exactly one hour. The paired values are genuinely
different readings (921 of 960 pairs disagree), so an hour of real trace has
been collapsed onto one set of timestamp labels. Keeping both would count that
hour twice; this drops 12 readings per affected person, 0.4% of their record.

The same dedupe is applied in src/add_cgm_ranges_to_full_person_aggregates.py,
so the range means reconstruct this column exactly.

Inserted into the glucose block, right after avg_cgm_per_spike_mg_dl.

Overwrites outputs/cohort/full_person_aggregates.csv in place.

Run from repo root:
    .venv/bin/python src/add_mean_cgm_overall_to_full_person_aggregates.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GLUCOSE_RAW = ROOT / "outputs" / "glycemic_phenotype" / "glucose_raw.csv"
FULL_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"

NEW_COL = "mean_cgm_overall_mg_dl"


def compute_mean_cgm_overall() -> pd.DataFrame:
    readings = pd.read_csv(
        GLUCOSE_RAW, usecols=["person_id", "timestamp_utc", "glucose_mg_dl"], dtype=str
    )
    readings["glucose_mg_dl"] = readings["glucose_mg_dl"].replace(
        {"Low": 39, "LOW": 39, "High": 401, "HIGH": 401}
    )
    readings["glucose_mg_dl"] = pd.to_numeric(readings["glucose_mg_dl"], errors="coerce")
    readings["person_id"] = pd.to_numeric(readings["person_id"], errors="coerce").astype("Int64")
    readings = readings.dropna(subset=["person_id", "timestamp_utc", "glucose_mg_dl"])

    # Deduped on the raw timestamp string in file order, which is what the
    # range script does, so both keep the identical row without depending on
    # sort stability.
    n_before = len(readings)
    readings = readings.drop_duplicates(["person_id", "timestamp_utc"], keep="first")
    print(f"Dropped {n_before - len(readings):,} duplicate (person_id, timestamp_utc) rows")

    return readings.groupby("person_id")["glucose_mg_dl"].mean().rename(NEW_COL).reset_index()


def main() -> None:
    # round_trip parsing keeps columns this script does not own bit-identical;
    # pandas' default float reader is up to 1 ULP off from strtod.
    full = pd.read_csv(FULL_PATH, float_precision="round_trip")
    full["person_id"] = pd.to_numeric(full["person_id"], errors="coerce").astype("Int64")
    n0, cols0 = len(full), list(full.columns)

    if NEW_COL in cols0:
        full = full.drop(columns=[NEW_COL])
        cols0 = [c for c in cols0 if c != NEW_COL]

    print("Computing overall mean CGM per person from raw readings...")
    mean_cgm = compute_mean_cgm_overall()

    merged = full.merge(mean_cgm, on="person_id", how="left")
    assert len(merged) == n0, f"row count changed: {n0} -> {len(merged)}"

    # Insert right after avg_cgm_per_spike_mg_dl.
    ordered = []
    for c in cols0:
        ordered.append(c)
        if c == "avg_cgm_per_spike_mg_dl":
            ordered.append(NEW_COL)
    if NEW_COL not in ordered:
        ordered.append(NEW_COL)
    merged = merged[ordered]

    merged.to_csv(FULL_PATH, index=False)

    print(f"Wrote: {FULL_PATH}  (rows={len(merged)}, cols={len(merged.columns)})")
    cov = merged[NEW_COL].notna().mean() * 100
    print(f"Non-null coverage of {NEW_COL}: {cov:.1f}%")


if __name__ == "__main__":
    main()
