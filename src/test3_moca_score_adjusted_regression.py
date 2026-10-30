from pathlib import Path

import pandas as pd
import statsmodels.formula.api as smf


ROOT = Path(__file__).resolve().parents[1]
FEATURES_PATH = ROOT / "outputs" / "eda_advanced" / "tables" / "master_patient_features.csv"
SPIKE_PATH = ROOT / "outputs" / "glycemic_phenotype" / "spike_aggregates_by_person.csv"

COEF_OUT_PATH = ROOT / "outputs" / "cognition" / "moca_score_adjusted_regression_coefficients.csv"
MODEL_OUT_PATH = ROOT / "outputs" / "cognition" / "moca_score_adjusted_regression_model_summary.csv"


def main() -> None:
    features = pd.read_csv(FEATURES_PATH)
    spikes = pd.read_csv(SPIKE_PATH)

    features["participant_id"] = features["participant_id"].astype(str).str.strip()
    spikes["person_id"] = spikes["person_id"].astype(str).str.strip()

    merged = features.merge(
        spikes,
        left_on="participant_id",
        right_on="person_id",
        how="inner",
        validate="one_to_one",
    )

    numeric_cols = [
        "moca_score",
        "avg_spikes_per_day",
        "avg_spike_duration_minutes",
        "avg_cgm_per_spike_mg_dl",
        "meta_age",
    ]
    for col in numeric_cols:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")

    model_df = merged[
        [
            "moca_score",
            "avg_spikes_per_day",
            "avg_spike_duration_minutes",
            "avg_cgm_per_spike_mg_dl",
            "meta_age",
            "meta_study_group",
        ]
    ].dropna()

    formula = (
        "moca_score ~ avg_spikes_per_day + avg_spike_duration_minutes + "
        "avg_cgm_per_spike_mg_dl + meta_age + C(meta_study_group)"
    )
    model = smf.ols(formula=formula, data=model_df).fit(cov_type="HC3")

    coef = model.summary2().tables[1].reset_index().rename(columns={"index": "term"})
    coef.to_csv(COEF_OUT_PATH, index=False)

    model_summary = pd.DataFrame(
        [
            {
                "n_rows_after_dropna": int(len(model_df)),
                "r_squared": float(model.rsquared),
                "adj_r_squared": float(model.rsquared_adj),
                "f_statistic": float(model.fvalue) if model.fvalue is not None else float("nan"),
                "f_pvalue": float(model.f_pvalue) if model.f_pvalue is not None else float("nan"),
                "aic": float(model.aic),
                "bic": float(model.bic),
            }
        ]
    )
    model_summary.to_csv(MODEL_OUT_PATH, index=False)

    print(f"Wrote: {COEF_OUT_PATH}")
    print(f"Wrote: {MODEL_OUT_PATH}")
    print(f"N used: {len(model_df)}")
    print(f"R^2: {model.rsquared:.4f} | Adj R^2: {model.rsquared_adj:.4f}")


if __name__ == "__main__":
    main()
