"""
Add mean_cgm_overall_mg_dl to outputs/full_person_aggregates.csv: each person's
mean glucose across ALL CGM readings (not just readings inside a detected
spike episode). Cleaning matches src/glucose_spike_detection.py and
src/glucose_summary_stratified.py: Dexcom "Low"/"High" sentinel strings are
mapped to 39/401 mg/dL before averaging.

Inserted into the glucose block, right after avg_cgm_per_spike_mg_dl.

Overwrites outputs/full_person_aggregates.csv in place.

Run from repo root:
    .venv/bin/python src/add_mean_cgm_overall_to_full_person_aggregates.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path.home() / "aireadi"
GLUCOSE_RAW = ROOT / "outputs" / "glucose_raw.csv"
FULL_PATH = ROOT / "outputs" / "full_person_aggregates.csv"

NEW_COL = "mean_cgm_overall_mg_dl"


def compute_mean_cgm_overall() -> pd.DataFrame:
    readings = pd.read_csv(GLUCOSE_RAW, usecols=["person_id", "glucose_mg_dl"], dtype=str)
    readings["glucose_mg_dl"] = readings["glucose_mg_dl"].replace(
        {"Low": 39, "LOW": 39, "High": 401, "HIGH": 401}
    )
    readings["glucose_mg_dl"] = pd.to_numeric(readings["glucose_mg_dl"], errors="coerce")
    readings["person_id"] = pd.to_numeric(readings["person_id"], errors="coerce").astype("Int64")
    readings = readings.dropna(subset=["person_id", "glucose_mg_dl"])
    return readings.groupby("person_id")["glucose_mg_dl"].mean().rename(NEW_COL).reset_index()


def main() -> None:
    full = pd.read_csv(FULL_PATH)
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
