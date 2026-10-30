"""
Compare physical-activity metrics between cognitively impaired vs unimpaired
participants, using three MoCA impairment cutoffs: 26, 24, 23.

Convention: a participant is "impaired" if moca_score < cutoff (standard MoCA
direction), "unimpaired" if moca_score >= cutoff.

Optimal test choice
-------------------
We compare a continuous physical-activity metric between two independent groups.
The activity metrics (movement, steps, calories) are right-skewed and not normal,
so the most appropriate/robust test is the Mann-Whitney U test (a.k.a. Wilcoxon
rank-sum) -- it is nonparametric and tests whether one group tends to have larger
values than the other. We report:
  - group sizes and medians (+ means)
  - Mann-Whitney U statistic and two-sided p-value
  - rank-biserial correlation as an effect size (range -1..1)
  - Benjamini-Hochberg FDR-corrected p-values across the metrics within each cutoff

Bootstrap
---------
On top of the analytic test we run a *stratified* nonparametric bootstrap: we
resample the impaired and unimpaired groups independently, with replacement,
preserving each group's size, and recompute the statistic on every replicate.
This yields percentile confidence intervals (and a bootstrap standard error) for:
  - the rank-biserial correlation (the Mann-Whitney effect size), and
  - the difference in medians (impaired - unimpaired).
A result is flagged "bootstrap significant" when the 95% CI for the rank-biserial
correlation excludes 0. The bootstrap is reproducible via a fixed seed.

Run from repo root: .venv/bin/python src/moca_impairment_activity_tests.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
FEATURES_PATH = ROOT / "outputs" / "eda_advanced" / "tables" / "master_patient_features.csv"
ACTIVITY_MANIFEST = (
    ROOT
    / "data"
    / "manifests"
    / "jeffreyzhuaireadi"
    / "7bc8c89b-8457-41fc-917c-b5571ff90ea0"
    / "dataset"
    / "wearable_activity_monitor"
    / "manifest.tsv"
)
OUT_PATH = ROOT / "outputs" / "cognition" / "moca_impairment_activity_tests.csv"

TARGET = "moca_score"
CUTOFFS = [26, 24, 23]

# Bootstrap configuration
N_BOOT = 10000
SEED = 12345
CI_LOW_PCT = 2.5
CI_HIGH_PCT = 97.5

PA_METRICS = {
    "total_movement": "Total movement over study period",
    "day_movement": "Daytime movement (total - night)",
    "night_movement": "Nighttime movement",
    "average_daily_activity": "Mean daily steps (manifest)",
    "average_active_calories_kcal": "Mean active calories/day (manifest)",
    "physical_activity_num_days": "Days of activity coverage (manifest)",
}


def load_data() -> pd.DataFrame:
    feats = pd.read_csv(FEATURES_PATH)
    feats["participant_id"] = pd.to_numeric(feats["participant_id"], errors="coerce")

    manifest = pd.read_csv(ACTIVITY_MANIFEST, sep="\t").rename(
        columns={"person_id": "participant_id"}
    )
    keep = [
        "participant_id",
        "average_daily_activity",
        "average_active_calories_kcal",
        "physical_activity_num_days",
    ]
    df = feats.merge(manifest[keep], on="participant_id", how="left")

    df["day_movement"] = pd.to_numeric(df["total_movement"], errors="coerce") - pd.to_numeric(
        df["night_movement"], errors="coerce"
    )
    df[TARGET] = pd.to_numeric(df[TARGET], errors="coerce")
    for col in PA_METRICS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _u_impaired(imp: np.ndarray, unimp_sorted: np.ndarray) -> float:
    """Mann-Whitney U for the impaired group (with tie handling), vectorized.

    U = #(imp > unimp) + 0.5 * #(imp == unimp), computed against a pre-sorted
    ``unimp_sorted`` via binary search so each bootstrap replicate is O(n log n).
    """
    left = np.searchsorted(unimp_sorted, imp, side="left")
    right = np.searchsorted(unimp_sorted, imp, side="right")
    strictly_less = left.sum()
    ties = (right - left).sum()
    return float(strictly_less) + 0.5 * float(ties)


def bootstrap_effect(
    imp: np.ndarray,
    unimp: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = N_BOOT,
) -> dict:
    """Stratified bootstrap CIs for rank-biserial correlation and median diff.

    Resamples ``imp`` and ``unimp`` independently with replacement (sizes fixed),
    recomputing the rank-biserial correlation and the median difference
    (impaired - unimpaired) on each replicate. Returns percentile CIs, bootstrap
    means, and bootstrap standard errors.
    """
    n_imp, n_unimp = len(imp), len(unimp)
    denom = n_imp * n_unimp

    rb = np.empty(n_boot, dtype=float)
    md = np.empty(n_boot, dtype=float)
    mean_imp = np.empty(n_boot, dtype=float)
    mean_unimp = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        a = rng.choice(imp, size=n_imp, replace=True)
        c = rng.choice(unimp, size=n_unimp, replace=True)
        c_sorted = np.sort(c)
        u = _u_impaired(a, c_sorted)
        rb[b] = 2.0 * u / denom - 1.0
        md[b] = np.median(a) - np.median(c)
        mean_imp[b] = a.mean()
        mean_unimp[b] = c.mean()

    mean_diff = mean_imp - mean_unimp

    rb_lo, rb_hi = np.percentile(rb, [CI_LOW_PCT, CI_HIGH_PCT])
    md_lo, md_hi = np.percentile(md, [CI_LOW_PCT, CI_HIGH_PCT])
    mi_lo, mi_hi = np.percentile(mean_imp, [CI_LOW_PCT, CI_HIGH_PCT])
    mu_lo, mu_hi = np.percentile(mean_unimp, [CI_LOW_PCT, CI_HIGH_PCT])
    mdf_lo, mdf_hi = np.percentile(mean_diff, [CI_LOW_PCT, CI_HIGH_PCT])
    return {
        "rank_biserial_boot_mean": float(rb.mean()),
        "rank_biserial_boot_se": float(rb.std(ddof=1)),
        "rank_biserial_ci_low": float(rb_lo),
        "rank_biserial_ci_high": float(rb_hi),
        "rank_biserial_ci_excludes_0": bool(rb_lo > 0 or rb_hi < 0),
        "median_diff_boot_mean": float(md.mean()),
        "median_diff_boot_se": float(md.std(ddof=1)),
        "median_diff_ci_low": float(md_lo),
        "median_diff_ci_high": float(md_hi),
        "median_diff_ci_excludes_0": bool(md_lo > 0 or md_hi < 0),
        "mean_impaired_ci_low": float(mi_lo),
        "mean_impaired_ci_high": float(mi_hi),
        "mean_unimpaired_ci_low": float(mu_lo),
        "mean_unimpaired_ci_high": float(mu_hi),
        "mean_diff": float(imp.mean() - unimp.mean()),
        "mean_diff_ci_low": float(mdf_lo),
        "mean_diff_ci_high": float(mdf_hi),
        "mean_diff_ci_excludes_0": bool(mdf_lo > 0 or mdf_hi < 0),
    }


def bootstrap_means(
    imp: np.ndarray,
    unimp: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = N_BOOT,
) -> dict:
    """Fast, fully-vectorized stratified bootstrap of GROUP MEANS.

    Resamples ``imp`` and ``unimp`` independently with replacement (sizes fixed)
    and returns each group's mean with a 95% percentile CI, plus the mean
    difference (impaired - unimpaired) with its CI. Because we only need means,
    all ``n_boot`` resamples are drawn at once via integer indexing (no Python
    loop), which is dramatically faster than the rank-based bootstrap.
    """
    n_imp, n_unimp = len(imp), len(unimp)

    idx_i = rng.integers(0, n_imp, size=(n_boot, n_imp))
    boot_imp = imp[idx_i].mean(axis=1)
    idx_u = rng.integers(0, n_unimp, size=(n_boot, n_unimp))
    boot_unimp = unimp[idx_u].mean(axis=1)
    boot_diff = boot_imp - boot_unimp

    mi_lo, mi_hi = np.percentile(boot_imp, [CI_LOW_PCT, CI_HIGH_PCT])
    mu_lo, mu_hi = np.percentile(boot_unimp, [CI_LOW_PCT, CI_HIGH_PCT])
    md_lo, md_hi = np.percentile(boot_diff, [CI_LOW_PCT, CI_HIGH_PCT])
    return {
        "mean_impaired": float(imp.mean()),
        "mean_impaired_ci_low": float(mi_lo),
        "mean_impaired_ci_high": float(mi_hi),
        "mean_unimpaired": float(unimp.mean()),
        "mean_unimpaired_ci_low": float(mu_lo),
        "mean_unimpaired_ci_high": float(mu_hi),
        "mean_diff": float(imp.mean() - unimp.mean()),
        "mean_diff_ci_low": float(md_lo),
        "mean_diff_ci_high": float(md_hi),
        "mean_diff_ci_excludes_0": bool(md_lo > 0 or md_hi < 0),
    }


_MEAN_BOOT_KEYS = [
    "mean_impaired",
    "mean_impaired_ci_low",
    "mean_impaired_ci_high",
    "mean_unimpaired",
    "mean_unimpaired_ci_low",
    "mean_unimpaired_ci_high",
    "mean_diff",
    "mean_diff_ci_low",
    "mean_diff_ci_high",
    "mean_diff_ci_excludes_0",
]


def _empty_mean_boot() -> dict:
    return {
        k: (False if k.endswith("excludes_0") else np.nan) for k in _MEAN_BOOT_KEYS
    }


_BOOT_KEYS = [
    "rank_biserial_boot_mean",
    "rank_biserial_boot_se",
    "rank_biserial_ci_low",
    "rank_biserial_ci_high",
    "rank_biserial_ci_excludes_0",
    "median_diff_boot_mean",
    "median_diff_boot_se",
    "median_diff_ci_low",
    "median_diff_ci_high",
    "median_diff_ci_excludes_0",
    "mean_impaired_ci_low",
    "mean_impaired_ci_high",
    "mean_unimpaired_ci_low",
    "mean_unimpaired_ci_high",
    "mean_diff",
    "mean_diff_ci_low",
    "mean_diff_ci_high",
    "mean_diff_ci_excludes_0",
]


def _empty_boot() -> dict:
    return {k: (False if k.endswith("excludes_0") else np.nan) for k in _BOOT_KEYS}


def run() -> pd.DataFrame:
    df = load_data()
    rng = np.random.default_rng(SEED)
    rows: list[dict] = []

    for cutoff in CUTOFFS:
        valid = df[df[TARGET].notna()].copy()
        valid["impaired"] = valid[TARGET] < cutoff

        for col, desc in PA_METRICS.items():
            sub = valid[[col, "impaired"]].dropna()
            imp = sub.loc[sub["impaired"], col]
            unimp = sub.loc[~sub["impaired"], col]
            n_imp, n_unimp = len(imp), len(unimp)

            if n_imp < 3 or n_unimp < 3:
                row = {
                    "cutoff": cutoff,
                    "metric": col,
                    "description": desc,
                    "n_impaired": n_imp,
                    "n_unimpaired": n_unimp,
                    "median_impaired": imp.median() if n_imp else np.nan,
                    "median_unimpaired": unimp.median() if n_unimp else np.nan,
                    "mean_impaired": imp.mean() if n_imp else np.nan,
                    "mean_unimpaired": unimp.mean() if n_unimp else np.nan,
                    "mann_whitney_u": np.nan,
                    "p_value": np.nan,
                    "rank_biserial": np.nan,
                    "median_diff": np.nan,
                }
                row.update(_empty_boot())
                rows.append(row)
                continue

            imp_arr = imp.to_numpy(dtype=float)
            unimp_arr = unimp.to_numpy(dtype=float)

            u, p = mannwhitneyu(imp, unimp, alternative="two-sided")
            # rank-biserial correlation; positive => impaired tend to have HIGHER values
            rank_biserial = 2.0 * u / (n_imp * n_unimp) - 1.0

            row = {
                "cutoff": cutoff,
                "metric": col,
                "description": desc,
                "n_impaired": n_imp,
                "n_unimpaired": n_unimp,
                "median_impaired": float(imp.median()),
                "median_unimpaired": float(unimp.median()),
                "mean_impaired": float(imp.mean()),
                "mean_unimpaired": float(unimp.mean()),
                "mann_whitney_u": float(u),
                "p_value": float(p),
                "rank_biserial": float(rank_biserial),
                "median_diff": float(imp.median() - unimp.median()),
            }
            row.update(bootstrap_effect(imp_arr, unimp_arr, rng))
            rows.append(row)

    results = pd.DataFrame(rows)

    # FDR-BH correction within each cutoff family (6 metrics per cutoff)
    results["p_value_fdr_bh"] = np.nan
    results["significant_fdr_0_05"] = False
    for cutoff in CUTOFFS:
        m = (results["cutoff"] == cutoff) & results["p_value"].notna()
        if m.any():
            reject, p_fdr, _, _ = multipletests(
                results.loc[m, "p_value"].to_numpy(), alpha=0.05, method="fdr_bh"
            )
            results.loc[m, "p_value_fdr_bh"] = p_fdr
            results.loc[m, "significant_fdr_0_05"] = reject

    return results


def main() -> None:
    results = run()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(OUT_PATH, index=False)

    pd.set_option("display.width", 240)
    pd.set_option("display.max_columns", 40)
    print(f"Wrote: {OUT_PATH}")
    print(f"Bootstrap: {N_BOOT} stratified resamples, seed={SEED}, 95% percentile CIs\n")

    def _fmt_ci(r):
        return f"[{r['rank_biserial_ci_low']:+.4f}, {r['rank_biserial_ci_high']:+.4f}]"

    for cutoff in CUTOFFS:
        block = results[results["cutoff"] == cutoff].copy()
        n_imp = int(block["n_impaired"].iloc[0])
        n_unimp = int(block["n_unimpaired"].iloc[0])
        print(f"===== Cutoff {cutoff}: impaired (MoCA < {cutoff}) n={n_imp} vs unimpaired n={n_unimp} =====")
        block["rank_biserial_95ci"] = block.apply(_fmt_ci, axis=1)
        show = block[
            [
                "metric",
                "median_impaired",
                "median_unimpaired",
                "p_value",
                "p_value_fdr_bh",
                "significant_fdr_0_05",
                "rank_biserial",
                "rank_biserial_95ci",
                "rank_biserial_ci_excludes_0",
            ]
        ]
        print(show.to_string(index=False))
        print()


if __name__ == "__main__":
    main()
