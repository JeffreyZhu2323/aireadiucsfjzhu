"""
Build person-level MEAN aggregates for the five non-physical-activity wearable
streams from the Garmin Vivosmart 5 (AIREADI wearable_activity_monitor):

  - mean_heartrate_bpm          (heart_rate)
  - mean_oxygen_saturation_pct  (oxygen_saturation)
  - mean_stress_level           (stress)
  - mean_respiratory_rate_bpm   (respiratory_rate)
  - mean_sleep_hours            (sleep; mean nightly asleep hours)

Why we recompute instead of using the manifest averages:
  The manifest's precomputed averages are contaminated by device sentinel codes
  (HR value==0, stress/respiratory value in {-1,-2} meaning "no reading"/"too much
  motion"). This drags e.g. average_respiratory_rate_bpm to ~6.9 (impossible) and
  makes average_sleep_hours (~0.39) meaningless. Here we clean at the sample level
  first, then take the mean.

Cleaning (sample-level, applied before averaging):
  heart_rate         keep value in [30, 220]
  oxygen_saturation  keep value in [70, 100]
  stress             keep value in [0, 100]      (drop negative sentinels)
  respiratory_rate   keep value in [4, 40]        (drop negative sentinels)

Sleep metric:
  Episodes are (start, end, stage) intervals. We dedupe exact-duplicate records,
  keep asleep stages (light/deep/rem), assign each episode to a "night" via a
  noon-to-noon day (start - 12h), sum asleep hours per night, then take the mean
  across nights -> mean_sleep_hours.

Point-series parsing is done with a regex over the raw bytes (each record has
exactly one `"value": <number>`), which is validated to be identical to a full
json.load but far faster over ~31 GB of files.

Output: outputs/wearable_stream_means.csv  (one row per participant_id)

Run from repo root:
    .venv/bin/python src/wearable_stream_aggregates.py
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path.home() / "aireadi"
WAM = (
    ROOT
    / "data"
    / "raw"
    / "jeffreyzhuaireadi"
    / "jeffreyzhuaireadi"
    / "7bc8c89b-8457-41fc-917c-b5571ff90ea0"
    / "dataset"
    / "wearable_activity_monitor"
)
OUT_PATH = ROOT / "outputs" / "wearable_stream_means.csv"

VALUE_RE = re.compile(rb'"value":\s*(-?\d+(?:\.\d+)?)')

# stream_dir -> (json filename suffix, output mean column, (lo, hi) valid range)
POINT_STREAMS = {
    "heart_rate": ("heartrate", "mean_heartrate_bpm", (30.0, 220.0)),
    "oxygen_saturation": ("oxygensaturation", "mean_oxygen_saturation_pct", (70.0, 100.0)),
    "stress": ("stress", "mean_stress_level", (0.0, 100.0)),
    "respiratory_rate": ("respiratoryrate", "mean_respiratory_rate_bpm", (4.0, 40.0)),
}

ASLEEP_STAGES = {"light", "deep", "rem"}


def point_mean(path: Path, lo: float, hi: float) -> tuple[float, int]:
    """Mean of in-range measurement values in a point-series file.

    Returns (mean, n_valid). mean is NaN if no valid samples.
    """
    data = path.read_bytes()
    vals = np.array([float(m) for m in VALUE_RE.findall(data)], dtype=float)
    if vals.size == 0:
        return float("nan"), 0
    valid = vals[(vals >= lo) & (vals <= hi)]
    if valid.size == 0:
        return float("nan"), 0
    return float(valid.mean()), int(valid.size)


def sleep_mean(path: Path) -> tuple[float, int]:
    """Mean nightly asleep hours for one sleep file.

    Returns (mean_sleep_hours, n_nights). NaN if no asleep episodes.
    """
    with path.open() as fh:
        recs = json.load(fh)["body"]["sleep"]

    seen: set[tuple[str, str, str]] = set()
    starts: list[str] = []
    ends: list[str] = []
    for r in recs:
        state = r["sleep_stage_state"]
        if state not in ASLEEP_STAGES:
            continue
        ti = r["effective_time_frame"]["time_interval"]
        key = (ti["start_date_time"], ti["end_date_time"], state)
        if key in seen:
            continue
        seen.add(key)
        starts.append(ti["start_date_time"])
        ends.append(ti["end_date_time"])

    if not starts:
        return float("nan"), 0

    s = pd.to_datetime(pd.Series(starts))
    e = pd.to_datetime(pd.Series(ends))
    hrs = (e - s).dt.total_seconds() / 3600.0
    night = (s - pd.Timedelta(hours=12)).dt.date
    nightly = hrs.groupby(night).sum()
    if nightly.empty:
        return float("nan"), 0
    return float(nightly.mean()), int(nightly.size)


def _person_ids(stream: str) -> list[str]:
    d = WAM / stream / "garmin_vivosmart5"
    if not d.exists():
        return []
    return sorted(p.name for p in d.iterdir() if p.is_dir())


def main() -> None:
    records: dict[int, dict] = {}

    def slot(pid_str: str) -> dict:
        pid = int(pid_str)
        return records.setdefault(pid, {"participant_id": pid})

    # ---- point streams ----
    for stream, (suffix, col, (lo, hi)) in POINT_STREAMS.items():
        ids = _person_ids(stream)
        t0 = time.time()
        print(f"[{stream}] {len(ids)} persons ...", flush=True)
        for i, pid_str in enumerate(ids, 1):
            f = WAM / stream / "garmin_vivosmart5" / pid_str / f"{pid_str}_{suffix}.json"
            if not f.exists():
                continue
            m, n = point_mean(f, lo, hi)
            row = slot(pid_str)
            row[col] = m
            row[f"{col}_n_valid"] = n
            if i % 250 == 0:
                print(f"  {stream}: {i}/{len(ids)}  ({time.time()-t0:.0f}s)", flush=True)
        print(f"[{stream}] done in {time.time()-t0:.0f}s", flush=True)

    # ---- sleep ----
    ids = _person_ids("sleep")
    t0 = time.time()
    print(f"[sleep] {len(ids)} persons ...", flush=True)
    for i, pid_str in enumerate(ids, 1):
        f = WAM / "sleep" / "garmin_vivosmart5" / pid_str / f"{pid_str}_sleep.json"
        if not f.exists():
            continue
        m, n = sleep_mean(f)
        row = slot(pid_str)
        row["mean_sleep_hours"] = m
        row["mean_sleep_hours_n_nights"] = n
        if i % 250 == 0:
            print(f"  sleep: {i}/{len(ids)}  ({time.time()-t0:.0f}s)", flush=True)
    print(f"[sleep] done in {time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(sorted(records.values(), key=lambda r: r["participant_id"]))
    cols = ["participant_id"]
    for _, col, _r in [(s, c, r) for s, (suf, c, r) in POINT_STREAMS.items()]:
        cols += [col, f"{col}_n_valid"]
    cols += ["mean_sleep_hours", "mean_sleep_hours_n_nights"]
    df = df[[c for c in cols if c in df.columns]]

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_PATH, index=False)
    print(f"\nWrote: {OUT_PATH}  (rows={len(df)})")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(df.describe(include="all").to_string())


if __name__ == "__main__":
    main()
