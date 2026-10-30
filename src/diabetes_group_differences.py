"""
Compare glucose-spike aggregates and wearable activity metrics across the four
diabetes-status groups (one row per person in outputs/full_person_aggregates.csv):

    healthy
    pre_diabetes_lifestyle_controlled
    oral_medication_and_or_non_insulin_injectable_medication_controlled
    insulin_dependent

Statistical design (deliberately simple, nonparametric, matches the repo's
existing Mann-Whitney/BH philosophy for the 2-group case):

  1. Omnibus per metric -- Kruskal-Wallis H test across all 4 groups.
     Rank-based one-way "ANOVA on ranks"; robust to the right-skew of these
     metrics. Answers "does this metric differ across the diabetes groups?".
  2. Effect size -- epsilon-squared = H / (n - 1). Magnitude of the group
     effect (0..1), reported alongside p because n ~ 2200 makes trivial
     differences "significant".
  3. Multiplicity across metrics -- Benjamini-Hochberg FDR applied to the
     omnibus p-values, corrected *within* each metric family (glucose, wearable).
  4. Post-hoc pairwise -- Dunn's test (tie-corrected, BH-adjusted across the 6
     pairs within a metric), run only for metrics whose omnibus is significant.
  5. Descriptives -- per-group n, median and IQR (q25, q75) for every metric.

Exposure note: days_tracked (wear time) is excluded as a comparison target, and
total_spikes is exposure-driven (more tracked days -> more total spikes); prefer
avg_spikes_per_day. total_spikes is still reported but read it with that caveat.

Outputs:
  outputs/diabetes_group_differences.csv        (omnibus + per-group descriptives)
  outputs/diabetes_group_differences_dunn.csv   (pairwise Dunn's, sig. metrics)

Run from repo root:
    .venv/bin/python src/diabetes_group_differences.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, norm, rankdata
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"
OUT_OMNIBUS = ROOT / "outputs" / "group_differences" / "diabetes_group_differences.csv"
OUT_DUNN = ROOT / "outputs" / "group_differences" / "diabetes_group_differences_dunn.csv"

GROUP_COL = "meta_study_group"

# Diabetes-severity order (used for stable, interpretable output ordering).
GROUP_ORDER = [
    "healthy",
    "pre_diabetes_lifestyle_controlled",
    "oral_medication_and_or_non_insulin_injectable_medication_controlled",
    "insulin_dependent",
]
# Short labels for compact column / row naming.
GROUP_LABEL = {
    "healthy": "healthy",
    "pre_diabetes_lifestyle_controlled": "pre_diabetes",
    "oral_medication_and_or_non_insulin_injectable_medication_controlled": "oral_med",
    "insulin_dependent": "insulin",
}

GLUCOSE_METRICS = {
    "avg_cgm_per_spike_mg_dl": "Mean CGM within spikes (mg/dL)",
    "avg_spike_duration_minutes": "Mean spike duration (min)",
    "avg_spikes_per_day": "Spikes per day",
    "total_spikes": "Total spikes (exposure-driven; see note)",
}
WEARABLE_METRICS = {
    "average_daily_activity": "Mean daily steps",
    "average_active_calories_kcal": "Mean active calories/day (kcal)",
    "total_movement": "Total movement",
    "night_movement": "Nighttime movement",
    "day_movement": "Daytime movement",
    "mean_heartrate_bpm": "Mean heart rate (bpm)",
    "mean_oxygen_saturation_pct": "Mean SpO2 (%)",
    "mean_stress_level": "Mean stress level",
    "mean_respiratory_rate_bpm": "Mean respiratory rate (bpm)",
    "mean_sleep_hours": "Mean sleep (hours/night)",
}
# metric -> (family, description)
METRICS = {
    **{m: ("glucose", d) for m, d in GLUCOSE_METRICS.items()},
    **{m: ("wearable", d) for m, d in WEARABLE_METRICS.items()},
}

ALPHA = 0.05


def epsilon_squared(h: float, n: int) -> float:
    """Kruskal-Wallis effect size epsilon^2 = H / (n - 1), range ~[0, 1]."""
    if n <= 1:
        return float("nan")
    return float(h / (n - 1))


def dunn_posthoc(groups: dict[str, np.ndarray]) -> list[dict]:
    """Tie-corrected Dunn's test for all group pairs.

    Ranks are computed once over the pooled sample (average ranks for ties).
    For groups i, j:
        z = (Rbar_i - Rbar_j) / sqrt(sigma2 * (1/n_i + 1/n_j))
        sigma2 = N(N+1)/12 - sum_t(t^3 - t) / (12 (N - 1))
    Two-sided p from the standard normal; BH-adjusted across the returned pairs.
    """
    labels = [g for g in GROUP_ORDER if g in groups and len(groups[g]) > 0]
    pooled = np.concatenate([groups[g] for g in labels])
    n_total = pooled.size
    ranks = rankdata(pooled)  # average method handles ties

    # tie correction sum_t (t^3 - t)
    _, counts = np.unique(pooled, return_counts=True)
    tie_sum = float(np.sum(counts**3 - counts))
    sigma2 = n_total * (n_total + 1) / 12.0 - tie_sum / (12.0 * (n_total - 1))

    mean_rank: dict[str, float] = {}
    sizes: dict[str, int] = {}
    offset = 0
    for g in labels:
        n_g = len(groups[g])
        mean_rank[g] = float(ranks[offset : offset + n_g].mean())
        sizes[g] = n_g
        offset += n_g

    rows: list[dict] = []
    pvals: list[float] = []
    for a_idx in range(len(labels)):
        for b_idx in range(a_idx + 1, len(labels)):
            a, b = labels[a_idx], labels[b_idx]
            se = np.sqrt(sigma2 * (1.0 / sizes[a] + 1.0 / sizes[b]))
            z = (mean_rank[a] - mean_rank[b]) / se
            p = 2.0 * norm.sf(abs(z))
            pvals.append(p)
            rows.append(
                {
                    "group_a": GROUP_LABEL[a],
                    "group_b": GROUP_LABEL[b],
                    "mean_rank_a": mean_rank[a],
                    "mean_rank_b": mean_rank[b],
                    "mean_rank_diff": mean_rank[a] - mean_rank[b],
                    "z": float(z),
                    "p_value": float(p),
                }
            )

    if pvals:
        reject, p_bh, _, _ = multipletests(pvals, alpha=ALPHA, method="fdr_bh")
        for r, pb, rej in zip(rows, p_bh, reject):
            r["p_value_bh"] = float(pb)
            r["significant_bh_0_05"] = bool(rej)
    return rows


def run() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(DATA_PATH)
    df = df[df[GROUP_COL].isin(GROUP_ORDER)].copy()

    omni_rows: list[dict] = []
    dunn_records: list[dict] = []

    for metric, (family, desc) in METRICS.items():
        if metric not in df.columns:
            continue
        sub = df[[GROUP_COL, metric]].copy()
        sub[metric] = pd.to_numeric(sub[metric], errors="coerce")
        sub = sub.dropna(subset=[metric])

        groups = {
            g: sub.loc[sub[GROUP_COL] == g, metric].to_numpy(dtype=float)
            for g in GROUP_ORDER
            if (sub[GROUP_COL] == g).any()
        }

        row: dict = {"metric": metric, "family": family, "description": desc}
        n_total = 0
        for g in GROUP_ORDER:
            lab = GROUP_LABEL[g]
            arr = groups.get(g, np.array([]))
            n_total += arr.size
            row[f"n_{lab}"] = int(arr.size)
            row[f"median_{lab}"] = float(np.median(arr)) if arr.size else np.nan
            row[f"q25_{lab}"] = float(np.percentile(arr, 25)) if arr.size else np.nan
            row[f"q75_{lab}"] = float(np.percentile(arr, 75)) if arr.size else np.nan
        row["n_total"] = int(n_total)

        usable = [a for a in groups.values() if a.size >= 3]
        if len(usable) >= 2:
            h, p = kruskal(*usable)
            row["kruskal_H"] = float(h)
            row["kruskal_p"] = float(p)
            row["epsilon_squared"] = epsilon_squared(float(h), n_total)
        else:
            row["kruskal_H"] = np.nan
            row["kruskal_p"] = np.nan
            row["epsilon_squared"] = np.nan
        omni_rows.append(row)

    omnibus = pd.DataFrame(omni_rows)

    # BH-FDR across metrics, within each family separately.
    omnibus["kruskal_p_fdr_bh"] = np.nan
    omnibus["significant_fdr_0_05"] = False
    for fam in omnibus["family"].unique():
        mask = (omnibus["family"] == fam) & omnibus["kruskal_p"].notna()
        if mask.any():
            reject, p_fdr, _, _ = multipletests(
                omnibus.loc[mask, "kruskal_p"].to_numpy(), alpha=ALPHA, method="fdr_bh"
            )
            omnibus.loc[mask, "kruskal_p_fdr_bh"] = p_fdr
            omnibus.loc[mask, "significant_fdr_0_05"] = reject

    # Dunn's post-hoc only for metrics significant after FDR.
    for _, r in omnibus[omnibus["significant_fdr_0_05"]].iterrows():
        metric = r["metric"]
        sub = df[[GROUP_COL, metric]].copy()
        sub[metric] = pd.to_numeric(sub[metric], errors="coerce")
        sub = sub.dropna(subset=[metric])
        groups = {
            g: sub.loc[sub[GROUP_COL] == g, metric].to_numpy(dtype=float)
            for g in GROUP_ORDER
            if (sub[GROUP_COL] == g).sum() >= 3
        }
        for rec in dunn_posthoc(groups):
            rec = {"metric": metric, "family": r["family"], **rec}
            dunn_records.append(rec)

    dunn = pd.DataFrame(dunn_records)
    return omnibus, dunn


def main() -> None:
    omnibus, dunn = run()
    OUT_OMNIBUS.parent.mkdir(parents=True, exist_ok=True)
    omnibus.to_csv(OUT_OMNIBUS, index=False)
    dunn.to_csv(OUT_DUNN, index=False)

    pd.set_option("display.width", 240)
    pd.set_option("display.max_columns", 60)

    print(f"Wrote: {OUT_OMNIBUS}")
    print(f"Wrote: {OUT_DUNN}\n")

    print("===== Omnibus (Kruskal-Wallis) by metric =====")
    show_cols = [
        "family",
        "metric",
        "median_healthy",
        "median_pre_diabetes",
        "median_oral_med",
        "median_insulin",
        "kruskal_H",
        "kruskal_p",
        "kruskal_p_fdr_bh",
        "epsilon_squared",
        "significant_fdr_0_05",
    ]
    with pd.option_context("display.float_format", lambda v: f"{v:,.4g}"):
        print(omnibus[show_cols].to_string(index=False))

    if not dunn.empty:
        print("\n===== Dunn's post-hoc (BH-adjusted) for FDR-significant metrics =====")
        dshow = dunn[
            ["metric", "group_a", "group_b", "z", "p_value", "p_value_bh", "significant_bh_0_05"]
        ]
        with pd.option_context("display.float_format", lambda v: f"{v:,.4g}"):
            print(dshow.to_string(index=False))


if __name__ == "__main__":
    main()
