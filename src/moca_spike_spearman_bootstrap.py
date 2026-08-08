"""
Spearman correlation between MoCA score and the three person-level glucose-spike
aggregate features, each with a bootstrap 95% confidence interval.

Features tested (from full_person_aggregates.csv):
  - avg_spike_duration_minutes
  - avg_spikes_per_day
  - avg_cgm_per_spike_mg_dl

For each feature we report:
  - n (complete pairs after dropping NaNs)
  - Spearman rho and its (asymptotic) p-value
  - bootstrap 95% CI for rho (BCa, resampling person rows with replacement)
  - FDR (Benjamini-Hochberg) adjusted p-value across the three features

Run from repo root:
    .venv/bin/python src/moca_spike_spearman_bootstrap.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path.home() / "aireadi"
DATA_PATH = ROOT / "outputs" / "full_person_aggregates.csv"
OUT_PATH = ROOT / "outputs" / "moca_spike_spearman_bootstrap.csv"

TARGET = "moca_score"
FEATURES = [
    "avg_spike_duration_minutes",
    "avg_spikes_per_day",
    "avg_cgm_per_spike_mg_dl",
]

N_BOOT = 10_000
SEED = 42
CI_LEVEL = 0.95


def _spearman_stat(x: np.ndarray, y: np.ndarray) -> float:
    return spearmanr(x, y).correlation


def main() -> None:
    df = pd.read_csv(DATA_PATH)

    rows = []
    raw_pvals = []

    for col in FEATURES:
        tmp = df[[TARGET, col]].copy()
        tmp[TARGET] = pd.to_numeric(tmp[TARGET], errors="coerce")
        tmp[col] = pd.to_numeric(tmp[col], errors="coerce")
        tmp = tmp.dropna()

        n = len(tmp)
        if n < 3:
            rows.append(
                {
                    "feature": col,
                    "n": n,
                    "spearman_rho": float("nan"),
                    "p_value": float("nan"),
                    "ci_low": float("nan"),
                    "ci_high": float("nan"),
                }
            )
            raw_pvals.append(float("nan"))
            continue

        target = tmp[TARGET].to_numpy()
        feat = tmp[col].to_numpy()

        res = spearmanr(target, feat)
        rho = float(res.correlation)
        pval = float(res.pvalue)

        boot = bootstrap(
            (target, feat),
            _spearman_stat,
            n_resamples=N_BOOT,
            confidence_level=CI_LEVEL,
            method="BCa",
            paired=True,
            random_state=SEED,
        )

        rows.append(
            {
                "feature": col,
                "n": n,
                "spearman_rho": rho,
                "p_value": pval,
                "ci_low": float(boot.confidence_interval.low),
                "ci_high": float(boot.confidence_interval.high),
            }
        )
        raw_pvals.append(pval)

    results = pd.DataFrame(rows)

    valid_mask = results["p_value"].notna()
    pvals_valid = results.loc[valid_mask, "p_value"].to_numpy()
    if len(pvals_valid) > 0:
        reject, pvals_fdr, _, _ = multipletests(pvals_valid, alpha=0.05, method="fdr_bh")
        results.loc[valid_mask, "p_value_fdr_bh"] = pvals_fdr
        results.loc[valid_mask, "significant_fdr_0_05"] = reject
    else:
        results["p_value_fdr_bh"] = float("nan")
        results["significant_fdr_0_05"] = False

    results = results.sort_values("p_value", na_position="last").reset_index(drop=True)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(OUT_PATH, index=False)

    print(f"Wrote: {OUT_PATH}")
    print(f"(n_boot={N_BOOT}, method=BCa, ci_level={CI_LEVEL}, seed={SEED})")
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
