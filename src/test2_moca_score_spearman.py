from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path.home() / "aireadi"
FEATURES_PATH = ROOT / "outputs" / "eda_advanced" / "tables" / "master_patient_features.csv"
OUT_PATH = ROOT / "outputs" / "moca_score_spearman_results.csv"

TARGET = "moca_score"
PREDICTORS = [
    "mean_glucose",
    "night_mean_hr",
    "night_movement",
    "min_spo2",
    "day_hr",
    "total_movement",
    "meta_age",
]


def main() -> None:
    df = pd.read_csv(FEATURES_PATH)
    df[TARGET] = pd.to_numeric(df[TARGET], errors="coerce")

    rows = []
    raw_pvals = []

    for col in PREDICTORS:
        tmp = df[[TARGET, col]].copy()
        tmp[TARGET] = pd.to_numeric(tmp[TARGET], errors="coerce")
        tmp[col] = pd.to_numeric(tmp[col], errors="coerce")
        tmp = tmp.dropna()

        if len(tmp) < 3:
            rows.append(
                {
                    "predictor": col,
                    "n": len(tmp),
                    "spearman_rho": float("nan"),
                    "p_value": float("nan"),
                }
            )
            raw_pvals.append(float("nan"))
            continue

        rho, pval = spearmanr(tmp[TARGET], tmp[col], nan_policy="omit")
        rows.append(
            {
                "predictor": col,
                "n": len(tmp),
                "spearman_rho": float(rho),
                "p_value": float(pval),
            }
        )
        raw_pvals.append(float(pval))

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
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
