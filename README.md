# UCSF AIreadi Workspace

Current workspace layout:

- `data/`: immutable source data and manifests
- `src/`: active scripts
- `outputs/`: official active outputs
- `outputs/eda_advanced/`: partner official raw/MOCA outputs
- `outputs/_archive_exploratory_2026-05-01/`: archived exploratory and legacy outputs
- `logs/`: runtime/debug logs

## Official output ownership

- Spike pipeline (your official outputs, in `outputs/`):
  - `glucose_raw.csv`
  - `glucose_spikes.csv`
  - `spike_daily_by_cgm_day.csv`
  - `spike_aggregates_by_person.csv`
  - interpolation/method variants:
    - `glucose_spikes_interpolated.csv`
    - `glucose_spikes_interpolated_min10.csv`
    - `spike_aggregates_by_person_interpolated.csv`
    - `spike_aggregates_by_person_interpolated_min10.csv`
    - `spike_aggregate_method_comparison.csv`
- Partner official overall raw/MOCA outputs:
  - `outputs/eda_advanced/tables/`
  - `outputs/eda_advanced/reports/`

## Active runs

- Glucose spike detection from official glucose cache:
  - `python3 src/glucose_spike_detection.py`
- Spike method comparison:
  - `python3 src/spike_method_comparison.py`
- Sleep interval chart:
  - `python3 src/plot_sleep_interval_summary.py`
- Sleep overlap chart:
  - `python3 src/plot_sleep_overlap_pairs.py`
