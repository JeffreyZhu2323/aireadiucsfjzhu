"""
Add hba1c_percent to outputs/full_person_aggregates.csv: each person's lab
HbA1c from the OMOP measurement table, in percent.

HbA1c is drawn once per participant, so this is a single value per person, not
a summary over repeated measures. Rows are selected by measurement_concept_id
3004410 (LOINC Hemoglobin A1c/Hemoglobin.total in Blood), which maps 1:1 to the
import_hba1c source value in this export. Note that lbscat_a1c is NOT HbA1c
despite its name -- its label is "Hemoglobin - g/dL" and it is the CBC
hemoglobin.

Inserted after days_tracked, so the CGM-derived block stays contiguous and the
lab-based glycemic anchor sits next to it.

Overwrites outputs/full_person_aggregates.csv in place.

Run from repo root:
    .venv/bin/python src/add_hba1c_to_full_person_aggregates.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MEASUREMENT = ROOT / "data" / "raw" / "omop" / "measurement.csv"
FULL_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"

NEW_COL = "hba1c_percent"
AFTER_COL = "days_tracked"
HBA1C_CONCEPT_ID = 3004410


def load_hba1c() -> pd.DataFrame:
    meas = pd.read_csv(
        MEASUREMENT,
        usecols=[
            "person_id",
            "measurement_concept_id",
            "measurement_source_value",
            "value_as_number",
        ],
        dtype=str,
    )
    concept = pd.to_numeric(meas["measurement_concept_id"], errors="coerce")
    is_hba1c = (concept == HBA1C_CONCEPT_ID) & meas["measurement_source_value"].str.startswith(
        "import_hba1c"
    )
    hba1c = meas.loc[is_hba1c].copy()

    hba1c["person_id"] = pd.to_numeric(hba1c["person_id"], errors="coerce").astype("Int64")
    hba1c[NEW_COL] = pd.to_numeric(hba1c["value_as_number"], errors="coerce")
    hba1c = hba1c.dropna(subset=["person_id", NEW_COL])

    dupes = hba1c["person_id"].duplicated().sum()
    assert dupes == 0, f"expected one HbA1c per person, found {dupes} duplicated person_id"

    return hba1c[["person_id", NEW_COL]]


def main() -> None:
    full = pd.read_csv(FULL_PATH)
    full["person_id"] = pd.to_numeric(full["person_id"], errors="coerce").astype("Int64")
    n0 = len(full)

    cols0 = [c for c in full.columns if c != NEW_COL]
    full = full[cols0]

    print("Loading HbA1c from measurement table...")
    hba1c = load_hba1c()
    print(f"HbA1c values found: {len(hba1c)}")

    merged = full.merge(hba1c, on="person_id", how="left")
    assert len(merged) == n0, f"row count changed: {n0} -> {len(merged)}"

    ordered = []
    for c in cols0:
        ordered.append(c)
        if c == AFTER_COL:
            ordered.append(NEW_COL)
    if NEW_COL not in ordered:
        ordered.append(NEW_COL)
    merged = merged[ordered]

    merged.to_csv(FULL_PATH, index=False)

    print(f"Wrote: {FULL_PATH}  (rows={len(merged)}, cols={len(merged.columns)})")
    present = merged[NEW_COL].notna()
    print(f"Non-null coverage of {NEW_COL}: {present.sum()}/{len(merged)} ({present.mean() * 100:.1f}%)")
    print(merged.loc[present, NEW_COL].describe())


if __name__ == "__main__":
    main()
