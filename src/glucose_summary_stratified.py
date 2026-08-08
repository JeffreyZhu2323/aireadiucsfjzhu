"""
Glucose summary statistics, overall and stratified by age band, diabetes
status, and age x diabetes status jointly.

Two levels of summary are reported for each stratum:

  PERSON-LEVEL (primary) -- each participant first gets one mean glucose
      value (mean of all their CGM readings), then summary stats (mean,
      sd, median, min, max, quartiles) are computed across participants.
      This is the analysis-ready number: every participant counts once,
      regardless of how many CGM readings they happened to have.
      Output: outputs/glucose_summary_person_level.csv

  READING-LEVEL (pooled) -- all individual CGM readings in the stratum are
      pooled and summarized directly. Participants with more readings
      contribute more to this number.
      Output: outputs/glucose_summary_reading_level.csv

Age bands and diabetes group labels follow the convention used elsewhere in
this repo (src/moca_pa_age_diabetes_stratified.py):
  age band : 40-54 / 55-69 / 70+
  diabetes : healthy / pre_diab / oral_med / insulin_dep

Run from repo root: .venv/bin/python src/glucose_summary_stratified.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path.home() / "aireadi"
GLUCOSE_RAW = ROOT / "outputs" / "glucose_raw.csv"
PERSON_META = ROOT / "outputs" / "full_person_aggregates.csv"

PERSON_OUT = ROOT / "outputs" / "glucose_summary_person_level.csv"
READING_OUT = ROOT / "outputs" / "glucose_summary_reading_level.csv"

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


def summarize(values: pd.Series) -> dict:
    v = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if v.size == 0:
        return {
            "n": 0, "mean": np.nan, "sd": np.nan, "median": np.nan,
            "min": np.nan, "max": np.nan, "q25": np.nan, "q75": np.nan,
            "cv_pct": np.nan,
        }
    mean = v.mean()
    sd = v.std(ddof=1) if v.size > 1 else np.nan
    return {
        "n": v.size,
        "mean": mean,
        "sd": sd,
        "median": float(np.median(v)),
        "min": v.min(),
        "max": v.max(),
        "q25": float(np.percentile(v, 25)),
        "q75": float(np.percentile(v, 75)),
        "cv_pct": (sd / mean * 100) if not np.isnan(sd) and mean != 0 else np.nan,
    }


def load_person_meta() -> pd.DataFrame:
    meta = pd.read_csv(PERSON_META, usecols=["person_id", AGE, "meta_study_group"])
    meta["age_stratum"] = pd.cut(meta[AGE], bins=AGE_BINS, right=False, labels=AGE_LABELS)
    meta["diabetes_group"] = meta["meta_study_group"].map(DG_SHORT)
    return meta[["person_id", "age_stratum", "diabetes_group"]]


def build_strata(df: pd.DataFrame) -> list[tuple[str, str, pd.Series]]:
    """Yield (stratum_type, stratum_label, boolean_mask) for every stratum."""
    strata = [("overall", "all", pd.Series(True, index=df.index))]
    for age in AGE_LABELS:
        strata.append(("age", age, df["age_stratum"] == age))
    for dg in DG_ORDER:
        strata.append(("diabetes", dg, df["diabetes_group"] == dg))
    for age in AGE_LABELS:
        for dg in DG_ORDER:
            strata.append(
                ("age_x_diabetes", f"{age} | {dg}",
                 (df["age_stratum"] == age) & (df["diabetes_group"] == dg))
            )
    return strata


def main() -> None:
    print("Loading person metadata (age, diabetes group)...")
    meta = load_person_meta()

    print("Loading glucose_raw.csv (person_id, glucose_mg_dl only; "
          "~6.2M rows)...")
    readings = pd.read_csv(
        GLUCOSE_RAW, usecols=["person_id", "glucose_mg_dl"], dtype=str
    )
    # Dexcom G6 reports "High"/"Low" for out-of-sensor-range values; map to
    # the sensor's clinical bounds, matching src/glucose_spike_detection.py.
    readings["glucose_mg_dl"] = readings["glucose_mg_dl"].replace(
        {"Low": 39, "LOW": 39, "High": 401, "HIGH": 401}
    )
    readings["glucose_mg_dl"] = pd.to_numeric(readings["glucose_mg_dl"], errors="coerce")
    readings["person_id"] = pd.to_numeric(readings["person_id"], errors="coerce")
    readings = readings.merge(meta, on="person_id", how="left")

    person_mean = (
        readings.groupby("person_id")["glucose_mg_dl"]
        .agg(n_readings="count", person_mean_glucose_mg_dl="mean")
        .reset_index()
        .merge(meta, on="person_id", how="left")
    )
    print(f"  {len(person_mean)} persons with CGM data "
          f"({person_mean['n_readings'].sum():,} total readings)")

    # ---- PERSON-LEVEL summary (mean of per-person means, across persons) ----
    strata_person = build_strata(person_mean)
    person_rows = []
    for stype, label, mask in strata_person:
        stats = summarize(person_mean.loc[mask, "person_mean_glucose_mg_dl"])
        person_rows.append({"stratum_type": stype, "stratum": label, **stats})
    person_summary = pd.DataFrame(person_rows)
    person_summary.to_csv(PERSON_OUT, index=False)
    print(f"Wrote {PERSON_OUT}")

    # ---- READING-LEVEL summary (all raw readings pooled within stratum) ----
    strata_reading = build_strata(readings)
    reading_rows = []
    for stype, label, mask in strata_reading:
        stats = summarize(readings.loc[mask, "glucose_mg_dl"])
        reading_rows.append({"stratum_type": stype, "stratum": label, **stats})
    reading_summary = pd.DataFrame(reading_rows)
    reading_summary.to_csv(READING_OUT, index=False)
    print(f"Wrote {READING_OUT}")

    print("\n=== PERSON-LEVEL SUMMARY (n = number of persons) ===")
    print(person_summary.round(2).to_string(index=False))
    print("\n=== READING-LEVEL / POOLED SUMMARY (n = number of CGM readings) ===")
    print(reading_summary.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
