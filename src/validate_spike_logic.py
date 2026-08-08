"""Validation harness for the single-dip spike-merge logic.

Runs three kinds of checks:

1. Explicit edge cases with hand-computed expected results.
2. The vectorised point-method mask (``compute_spike_mask``) vs an independent
   brute-force reference, over many randomised streams and a real-data sample.
3. Equivalence between the point method and the interpolated state machine
   (``detect_interpolated_episodes``): they must partition each person's stream
   into the same episodes.

Exit code is non-zero if any check fails.
"""

import random
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

from glucose_spike_detection import (
    MAX_GAP_MINUTES,
    SPIKE_THRESHOLD_MG_DL,
    compute_spike_mask,
)
from spike_method_comparison import detect_interpolated_episodes

RAW_PATH = Path.home() / "aireadi" / "outputs" / "glucose_raw.csv"
THRESH = SPIKE_THRESHOLD_MG_DL
GAP = MAX_GAP_MINUTES
T0 = datetime(2023, 1, 1, tzinfo=timezone.utc)


def brute_is_spike(times: list[datetime], glucose: list[float]) -> list[bool]:
    """Reference implementation of the single-dip fold for one person."""
    n = len(glucose)
    above = [g > THRESH for g in glucose]
    out = list(above)
    for i in range(n):
        if above[i]:
            continue
        prev_ok = i - 1 >= 0 and above[i - 1]
        next_ok = i + 1 < n and above[i + 1]
        if not (prev_ok and next_ok):
            continue
        gap_prev = (times[i] - times[i - 1]).total_seconds() / 60.0
        gap_next = (times[i + 1] - times[i]).total_seconds() / 60.0
        if gap_prev <= GAP and gap_next <= GAP:
            out[i] = True
    return out


def runs_from_mask(mask: list[bool]) -> list[tuple[int, int]]:
    """Contiguous True runs as (start_idx, end_idx) inclusive."""
    runs = []
    start = None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        elif (not v) and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask) - 1))
    return runs


def point_mask_for_stream(times: list[datetime], glucose: list[float]) -> list[bool]:
    """Run the actual vectorised mask on a single-person stream."""
    df = pd.DataFrame(
        {
            "person_id": [1] * len(times),
            "timestamp_utc": pd.to_datetime(pd.Series(times), utc=True),
            "glucose_mg_dl": glucose,
        }
    )
    is_spike, _ = compute_spike_mask(df)
    return list(is_spike)


def check_edge_cases() -> list[str]:
    failures = []

    def ts(minutes: list[int]) -> list[datetime]:
        return [T0 + timedelta(minutes=m) for m in minutes]

    cases = [
        # name, minute offsets, glucose, expected is_spike mask
        ("user example single dips", [0, 5, 10, 15, 20], [142, 138, 145, 136, 142],
         [True, True, True, True, True]),
        ("two consecutive dips split", [0, 5, 10, 15], [142, 138, 136, 145],
         [True, False, False, True]),
        ("no spike at all", [0, 5, 10], [90, 100, 110],
         [False, False, False]),
        ("all spike", [0, 5, 10], [200, 201, 202],
         [True, True, True]),
        ("dip with big gap before", [0, 40, 45], [142, 138, 145],
         [True, False, True]),
        ("dip with big gap after", [0, 5, 45], [142, 138, 145],
         [True, False, True]),
        # 140 is NOT > 140, so it is a sub-threshold dip and folds in.
        ("boundary value 140 folds as dip", [0, 5, 10], [141, 140, 142],
         [True, True, True]),
        ("leading dip cannot fold", [0, 5], [138, 145],
         [False, True]),
        ("trailing dip cannot fold", [0, 5], [145, 138],
         [True, False]),
        ("gap exactly at limit folds", [0, 10, 20], [142, 138, 145],
         [True, True, True]),
        ("gap just over limit splits", [0, 11, 21], [142, 138, 145],
         [True, False, True]),
    ]

    for name, mins, gluc, expected in cases:
        got = point_mask_for_stream(ts(mins), [float(x) for x in gluc])
        if got != expected:
            failures.append(f"[edge/point] {name}: expected {expected}, got {got}")
        ref = brute_is_spike(ts(mins), [float(x) for x in gluc])
        if ref != expected:
            failures.append(f"[edge/brute] {name}: reference {ref} != expected {expected}")

        # interpolated method must produce the same run partition
        episodes = detect_interpolated_episodes(
            1, list(zip(ts(mins), [float(x) for x in gluc])), THRESH
        )
        exp_runs = runs_from_mask(expected)
        if len(episodes) != len(exp_runs):
            failures.append(
                f"[edge/interp] {name}: expected {len(exp_runs)} episodes, got {len(episodes)}"
            )
    return failures


def random_stream(rng: random.Random) -> tuple[list[datetime], list[float]]:
    n = rng.randint(1, 40)
    times = []
    t = T0
    for _ in range(n):
        times.append(t)
        # mostly 5-min cadence, occasionally a gap or slight jitter
        step = rng.choice([5, 5, 5, 5, 6, 10, 11, 30, 60])
        t = t + timedelta(minutes=step)
    glucose = [float(rng.randint(60, 260)) for _ in range(n)]
    return times, glucose


def check_random_equivalence(n_streams: int = 5000) -> list[str]:
    failures = []
    rng = random.Random(1234)
    for k in range(n_streams):
        times, glucose = random_stream(rng)
        mask = point_mask_for_stream(times, glucose)
        ref = brute_is_spike(times, glucose)
        if mask != ref:
            failures.append(f"[random/mask] stream {k}: vectorised != brute reference")
            if len(failures) > 5:
                break
            continue
        point_runs = runs_from_mask(mask)
        episodes = detect_interpolated_episodes(1, list(zip(times, glucose)), THRESH)
        if len(point_runs) != len(episodes):
            failures.append(
                f"[random/interp] stream {k}: point {len(point_runs)} != interp {len(episodes)} episodes"
            )
            if len(failures) > 10:
                break
    return failures


def check_real_data_sample(n_people: int = 40) -> list[str]:
    failures = []
    if not RAW_PATH.exists():
        return [f"[real] missing {RAW_PATH}"]
    df = pd.read_csv(RAW_PATH)
    df["glucose_mg_dl"] = pd.to_numeric(
        df["glucose_mg_dl"].replace({"Low": 39, "LOW": 39, "High": 401, "HIGH": 401}),
        errors="coerce",
    )
    df["timestamp_utc"] = pd.to_datetime(df["timestamp_utc"], utc=True, errors="coerce")
    df["person_id"] = pd.to_numeric(df["person_id"], errors="coerce")
    df = df.dropna(subset=["person_id", "timestamp_utc", "glucose_mg_dl"]).copy()
    df["person_id"] = df["person_id"].astype(int)
    people = sorted(df["person_id"].unique())[:n_people]
    df = df[df["person_id"].isin(people)].sort_values(
        ["person_id", "timestamp_utc"]
    ).reset_index(drop=True)

    is_spike, _ = compute_spike_mask(df)
    df = df.assign(is_spike=list(is_spike))

    for pid, grp in df.groupby("person_id"):
        times = list(grp["timestamp_utc"].dt.to_pydatetime())
        glucose = [float(x) for x in grp["glucose_mg_dl"]]
        ref = brute_is_spike(times, glucose)
        vec = list(grp["is_spike"])
        if vec != ref:
            failures.append(f"[real/mask] person {pid}: vectorised != brute reference")
            continue
        point_runs = runs_from_mask(vec)
        episodes = detect_interpolated_episodes(int(pid), list(zip(times, glucose)), THRESH)
        if len(point_runs) != len(episodes):
            failures.append(
                f"[real/interp] person {pid}: point {len(point_runs)} != interp {len(episodes)}"
            )
    return failures


def main() -> None:
    all_failures = []
    print("1) edge cases...")
    f1 = check_edge_cases()
    print(f"   {'PASS' if not f1 else 'FAIL'} ({len(f1)} failures)")
    all_failures += f1

    print("2) randomised brute-force + cross-method equivalence (5000 streams)...")
    f2 = check_random_equivalence()
    print(f"   {'PASS' if not f2 else 'FAIL'} ({len(f2)} failures)")
    all_failures += f2

    print("3) real-data sample (40 people)...")
    f3 = check_real_data_sample()
    print(f"   {'PASS' if not f3 else 'FAIL'} ({len(f3)} failures)")
    all_failures += f3

    if all_failures:
        print("\nFAILURES:")
        for msg in all_failures[:50]:
            print("  -", msg)
        raise SystemExit(1)
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
