import csv
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

OUTDIR = Path(__file__).resolve().parents[1] / "outputs"
RAW_PATH = OUTDIR / "glycemic_phenotype" / "glucose_raw.csv"
POINT_PERSON_PATH = OUTDIR / "glycemic_phenotype" / "spike_aggregates_by_person.csv"
INTERP_SPIKES_PATH = OUTDIR / "glycemic_phenotype" / "glucose_spikes_interpolated.csv"
INTERP_SPIKES_MIN10_PATH = OUTDIR / "glycemic_phenotype" / "glucose_spikes_interpolated_min10.csv"
INTERP_PERSON_PATH = OUTDIR / "glycemic_phenotype" / "spike_aggregates_by_person_interpolated.csv"
INTERP_PERSON_MIN10_PATH = OUTDIR / "glycemic_phenotype" / "spike_aggregates_by_person_interpolated_min10.csv"
METHOD_COMPARE_PATH = OUTDIR / "glycemic_phenotype" / "spike_aggregate_method_comparison.csv"

SPIKE_THRESHOLD_MG_DL = 140.0
MIN_DURATION_MINUTES = 10.0
# A single sub-threshold reading between two above-threshold readings is folded
# into the spike (tolerance = 1 dip), only when both neighbours are within
# normal CGM cadence.
MAX_DIP_READINGS = 1
MAX_GAP_MINUTES = 10.0


@dataclass
class Episode:
    person_id: int
    spike_start_utc: datetime
    spike_end_utc: datetime
    duration_minutes: float
    mean_blood_sugar_mg_dl: float
    max_blood_sugar_mg_dl: float


def parse_timestamp_utc(ts: str) -> datetime:
    if ts.endswith("Z"):
        ts = ts[:-1] + "+00:00"
    return datetime.fromisoformat(ts)


def parse_glucose(value: str) -> float | None:
    if value in {"Low", "LOW"}:
        return 39.0
    if value in {"High", "HIGH"}:
        return 401.0
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def lerp_crossing_time(
    t0: datetime, g0: float, t1: datetime, g1: float, threshold: float
) -> datetime:
    if g1 == g0:
        return t1
    dt_seconds = (t1 - t0).total_seconds()
    if dt_seconds <= 0:
        return t1
    frac = (threshold - g0) / (g1 - g0)
    if frac < 0.0:
        frac = 0.0
    elif frac > 1.0:
        frac = 1.0
    return t0 + timedelta(seconds=dt_seconds * frac)


def finalize_run(
    person_id: int,
    rows: list[tuple[datetime, float]],
    start_idx: int,
    end_idx: int,
    threshold: float,
) -> Episode:
    start_ts, start_glucose = rows[start_idx]
    end_ts, end_glucose = rows[end_idx]

    if start_idx > 0:
        prev_ts, prev_glucose = rows[start_idx - 1]
        if prev_glucose <= threshold < start_glucose:
            start_ts = lerp_crossing_time(prev_ts, prev_glucose, start_ts, start_glucose, threshold)

    if end_idx < len(rows) - 1:
        next_ts, next_glucose = rows[end_idx + 1]
        if end_glucose > threshold >= next_glucose:
            end_ts = lerp_crossing_time(end_ts, end_glucose, next_ts, next_glucose, threshold)

    run_glucose = [g for _, g in rows[start_idx : end_idx + 1]]
    duration_minutes = max(0.0, (end_ts - start_ts).total_seconds() / 60.0)
    return Episode(
        person_id=person_id,
        spike_start_utc=start_ts,
        spike_end_utc=end_ts,
        duration_minutes=duration_minutes,
        mean_blood_sugar_mg_dl=sum(run_glucose) / len(run_glucose),
        max_blood_sugar_mg_dl=max(run_glucose),
    )


def detect_interpolated_episodes(
    person_id: int, rows: list[tuple[datetime, float]], threshold: float
) -> list[Episode]:
    rows = sorted(rows, key=lambda x: x[0])
    episodes: list[Episode] = []
    run_start_idx: int | None = None
    last_spike_idx: int | None = None
    dip_count = 0

    for idx, (_, glucose) in enumerate(rows):
        is_spike = glucose > threshold
        if is_spike:
            if run_start_idx is None:
                run_start_idx = idx
                last_spike_idx = idx
                dip_count = 0
                continue
            # Resuming after one or more tolerated dips: only keep them merged
            # if the recovery sample is time-adjacent, otherwise close the prior
            # run at the last real spike and start a fresh run here.
            if dip_count > 0:
                gap_min = (rows[idx][0] - rows[idx - 1][0]).total_seconds() / 60.0
                if gap_min > MAX_GAP_MINUTES:
                    episodes.append(
                        finalize_run(person_id, rows, run_start_idx, last_spike_idx, threshold)
                    )
                    run_start_idx = idx
            last_spike_idx = idx
            dip_count = 0
        else:
            if run_start_idx is None:
                continue
            gap_min = (rows[idx][0] - rows[idx - 1][0]).total_seconds() / 60.0
            if dip_count < MAX_DIP_READINGS and gap_min <= MAX_GAP_MINUTES:
                # Tolerate this lone sub-threshold reading and keep the run open.
                dip_count += 1
            else:
                episodes.append(
                    finalize_run(person_id, rows, run_start_idx, last_spike_idx, threshold)
                )
                run_start_idx = None
                last_spike_idx = None
                dip_count = 0

    if run_start_idx is not None:
        episodes.append(finalize_run(person_id, rows, run_start_idx, last_spike_idx, threshold))
    return episodes


def write_episode_csv(path: Path, episodes: list[Episode]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "person_id",
                "spike_start_utc",
                "spike_start_date",
                "spike_start_time_of_day",
                "spike_end_utc",
                "duration_minutes",
                "mean_blood_sugar_mg_dl",
                "max_blood_sugar_mg_dl",
            ],
        )
        writer.writeheader()
        for ep in episodes:
            writer.writerow(
                {
                    "person_id": ep.person_id,
                    "spike_start_utc": ep.spike_start_utc.isoformat(),
                    "spike_start_date": ep.spike_start_utc.date().isoformat(),
                    "spike_start_time_of_day": ep.spike_start_utc.strftime("%H:%M:%S"),
                    "spike_end_utc": ep.spike_end_utc.isoformat(),
                    "duration_minutes": ep.duration_minutes,
                    "mean_blood_sugar_mg_dl": ep.mean_blood_sugar_mg_dl,
                    "max_blood_sugar_mg_dl": ep.max_blood_sugar_mg_dl,
                }
            )


def person_aggregates(
    episodes: list[Episode], days_tracked: dict[int, int]
) -> list[dict[str, float | int | None]]:
    by_person: dict[int, list[Episode]] = {}
    for ep in episodes:
        by_person.setdefault(ep.person_id, []).append(ep)

    rows: list[dict[str, float | int | None]] = []
    for person_id in sorted(days_tracked):
        person_eps = by_person.get(person_id, [])
        total_spikes = len(person_eps)
        dtracked = days_tracked[person_id]
        if total_spikes == 0:
            avg_cgm = None
            avg_duration = None
        else:
            avg_cgm = sum(ep.mean_blood_sugar_mg_dl for ep in person_eps) / total_spikes
            avg_duration = sum(ep.duration_minutes for ep in person_eps) / total_spikes
        rows.append(
            {
                "person_id": person_id,
                "avg_cgm_per_spike_mg_dl": avg_cgm,
                "avg_spike_duration_minutes": avg_duration,
                "avg_spikes_per_day": total_spikes / dtracked,
                "total_spikes": total_spikes,
                "days_tracked": dtracked,
            }
        )
    return rows


def write_person_csv(path: Path, rows: list[dict[str, float | int | None]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "person_id",
                "avg_cgm_per_spike_mg_dl",
                "avg_spike_duration_minutes",
                "avg_spikes_per_day",
                "total_spikes",
                "days_tracked",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_point_method_aggregates(path: Path) -> dict[int, dict[str, float | int | None]]:
    out: dict[int, dict[str, float | int | None]] = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        days_col = "days_tracked" if "days_tracked" in reader.fieldnames else "n_cgm_days"
        for row in reader:
            person_id = int(row["person_id"])
            out[person_id] = {
                "avg_cgm_per_spike_mg_dl": (
                    None if row["avg_cgm_per_spike_mg_dl"] == "" else float(row["avg_cgm_per_spike_mg_dl"])
                ),
                "avg_spike_duration_minutes": (
                    None if row["avg_spike_duration_minutes"] == "" else float(row["avg_spike_duration_minutes"])
                ),
                "avg_spikes_per_day": float(row["avg_spikes_per_day"]),
                "total_spikes": int(float(row["total_spikes"])),
                "days_tracked": int(float(row[days_col])),
            }
    return out


def write_method_comparison_csv(
    point_rows: dict[int, dict[str, float | int | None]],
    interp_min10_rows: list[dict[str, float | int | None]],
) -> None:
    with METHOD_COMPARE_PATH.open("w", newline="") as f:
        fieldnames = [
            "person_id",
            "total_spikes_point",
            "total_spikes_interp_min10",
            "delta_total_spikes",
            "avg_spikes_per_day_point",
            "avg_spikes_per_day_interp_min10",
            "delta_avg_spikes_per_day",
            "avg_spike_duration_minutes_point",
            "avg_spike_duration_minutes_interp_min10",
            "delta_avg_spike_duration_minutes",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in interp_min10_rows:
            person_id = int(row["person_id"])
            base = point_rows[person_id]
            d_spikes = int(row["total_spikes"]) - int(base["total_spikes"])
            d_rate = float(row["avg_spikes_per_day"]) - float(base["avg_spikes_per_day"])

            base_dur = base["avg_spike_duration_minutes"]
            new_dur = row["avg_spike_duration_minutes"]
            if base_dur is None or new_dur is None:
                d_dur = None
            else:
                d_dur = float(new_dur) - float(base_dur)

            writer.writerow(
                {
                    "person_id": person_id,
                    "total_spikes_point": int(base["total_spikes"]),
                    "total_spikes_interp_min10": int(row["total_spikes"]),
                    "delta_total_spikes": d_spikes,
                    "avg_spikes_per_day_point": float(base["avg_spikes_per_day"]),
                    "avg_spikes_per_day_interp_min10": float(row["avg_spikes_per_day"]),
                    "delta_avg_spikes_per_day": d_rate,
                    "avg_spike_duration_minutes_point": base_dur,
                    "avg_spike_duration_minutes_interp_min10": new_dur,
                    "delta_avg_spike_duration_minutes": d_dur,
                }
            )


def main() -> None:
    if not RAW_PATH.exists():
        raise FileNotFoundError(f"Missing {RAW_PATH}")
    if not POINT_PERSON_PATH.exists():
        raise FileNotFoundError(f"Missing {POINT_PERSON_PATH}")

    print("Reading glucose stream and detecting interpolated spikes...")
    days_by_person: dict[int, set] = {}
    episodes_all: list[Episode] = []
    current_person: int | None = None
    person_rows: list[tuple[datetime, float]] = []

    with RAW_PATH.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            glucose = parse_glucose(row["glucose_mg_dl"])
            if glucose is None:
                continue
            try:
                person_id = int(float(row["person_id"]))
            except (TypeError, ValueError):
                continue
            ts = parse_timestamp_utc(row["timestamp_utc"])

            days_by_person.setdefault(person_id, set()).add(ts.date())
            if current_person is None:
                current_person = person_id
            if person_id != current_person:
                episodes_all.extend(
                    detect_interpolated_episodes(current_person, person_rows, SPIKE_THRESHOLD_MG_DL)
                )
                person_rows = []
                current_person = person_id
            person_rows.append((ts, glucose))

    if current_person is not None and person_rows:
        episodes_all.extend(detect_interpolated_episodes(current_person, person_rows, SPIKE_THRESHOLD_MG_DL))

    episodes_all.sort(key=lambda ep: (ep.person_id, ep.spike_start_utc))
    episodes_min10 = [ep for ep in episodes_all if ep.duration_minutes >= MIN_DURATION_MINUTES]

    write_episode_csv(INTERP_SPIKES_PATH, episodes_all)
    write_episode_csv(INTERP_SPIKES_MIN10_PATH, episodes_min10)

    days_tracked = {pid: len(days) for pid, days in days_by_person.items()}
    person_interp = person_aggregates(episodes_all, days_tracked)
    person_interp_min10 = person_aggregates(episodes_min10, days_tracked)
    write_person_csv(INTERP_PERSON_PATH, person_interp)
    write_person_csv(INTERP_PERSON_MIN10_PATH, person_interp_min10)

    point_person = read_point_method_aggregates(POINT_PERSON_PATH)
    write_method_comparison_csv(point_person, person_interp_min10)

    baseline_total = sum(int(v["total_spikes"]) for v in point_person.values())
    interp_total = sum(int(v["total_spikes"]) for v in person_interp)
    interp_min10_total = sum(int(v["total_spikes"]) for v in person_interp_min10)
    changed_people = sum(
        1
        for row in person_interp_min10
        if int(row["total_spikes"]) != int(point_person[int(row["person_id"])]["total_spikes"])
    )
    print(f"Baseline spikes (point timestamps): {baseline_total:,}")
    print(f"Interpolated spikes (no min duration): {interp_total:,}")
    print(f"Interpolated spikes (>= {MIN_DURATION_MINUTES:.0f} min): {interp_min10_total:,}")
    print(
        f"People with changed spike counts under interpolated + min-duration: "
        f"{changed_people:,} / {len(person_interp_min10):,}"
    )
    print(f"Wrote {INTERP_SPIKES_PATH}")
    print(f"Wrote {INTERP_SPIKES_MIN10_PATH}")
    print(f"Wrote {INTERP_PERSON_PATH}")
    print(f"Wrote {INTERP_PERSON_MIN10_PATH}")
    print(f"Wrote {METHOD_COMPARE_PATH}")


if __name__ == "__main__":
    main()
