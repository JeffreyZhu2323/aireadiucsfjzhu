"""
Severity-adjusted re-estimation of every significant survey-glycemia correlation.

src/survey_glycemic_spearman_bootstrap.py reports crude bivariate Spearman
correlations. Diabetes severity (meta_study_group) is a common cause of both
sides of many of those pairs -- it drives glycemia directly, and it drives
survey responses like diabetes distress and clinician-prescribed self-care --
so a crude correlation can be substantially severity mixing rather than a
within-person association. This script re-estimates each one holding severity
constant.

Workflow, which is the standard stratified-confounding analysis:

  1. Re-run the Spearman test INSIDE each severity stratum, so the four
     estimates are each computed on participants of comparable disease status.

  2. Test homogeneity across the four strata with Cochran's Q on the Fisher-z
     scale. This asks whether the stratum estimates are mutually consistent.
     If they are not, severity is an effect modifier and no single pooled
     number describes the cohort -- the strata must be reported separately.

  3. Pool the stratum estimates into one severity-adjusted rho: Fisher
     z-transform each, average with inverse-variance weights (n_k - 3), then
     back-transform. Pooling is what recovers the statistical power lost by
     splitting the sample, and it is the number that should replace the crude
     estimate.

  4. Judge confounding by how far the estimate MOVED (crude vs pooled), not by
     how many strata stayed individually significant. Splitting ~2100 people
     into four groups costs enough power that stratum p-values drop out of
     significance even when the underlying association is identical in every
     group, so counting significant strata systematically mistakes power loss
     for confounding. The conventional threshold is a >=10% change in the
     estimate.

Because the pooled estimates are a fresh family of tests, Benjamini-Hochberg is
applied to the pooled p-values within each target, matching how the original
script corrected its crude p-values.

Run from repo root:
    .venv/bin/python src/survey_glycemic_severity_stratified_pooling.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"
OUT_PATH = ROOT / "outputs" / "survey_glycemia" / "survey_glycemic_severity_stratified_pooling.csv"

TARGETS = {
    "hba1c_percent": ROOT / "outputs" / "survey_glycemia" / "survey_hba1c_spearman_bootstrap.csv",
    "mean_cgm_overall_mg_dl": ROOT / "outputs" / "survey_glycemia" / "survey_mean_cgm_spearman_bootstrap.csv",
}

STRATUM_COL = "meta_study_group"
# Ordered least to most severe; also fixes column order in the output.
STRATA = [
    "healthy",
    "pre_diabetes_lifestyle_controlled",
    "oral_medication_and_or_non_insulin_injectable_medication_controlled",
    "insulin_dependent",
]
SHORT = {
    "healthy": "healthy",
    "pre_diabetes_lifestyle_controlled": "pre_dm",
    "oral_medication_and_or_non_insulin_injectable_medication_controlled": "oral_med",
    "insulin_dependent": "insulin",
}

ALPHA = 0.05
# Conventional epidemiological threshold for declaring meaningful confounding.
CONFOUND_CHANGE = 0.10
# Below this a stratum cannot support a stable rank correlation.
MIN_STRATUM_N = 20


def stratified_pool(df: pd.DataFrame, feature: str, target: str) -> dict:
    """Re-run Spearman within each severity stratum, then pool on the Fisher-z scale."""
    pair = df[[feature, target, STRATUM_COL]].apply(
        lambda s: pd.to_numeric(s, errors="coerce") if s.name != STRATUM_COL else s
    ).dropna()

    crude = spearmanr(pair[feature], pair[target])
    row = {
        "feature": feature,
        "target": target,
        "n_total": len(pair),
        "rho_crude": float(crude.correlation),
        "p_crude": float(crude.pvalue),
    }

    zs, ws = [], []
    for stratum in STRATA:
        sub = pair[pair[STRATUM_COL] == stratum]
        key = SHORT[stratum]
        # A stratum with no variation in the feature (e.g. a binary item nobody
        # in this group endorsed) yields nan rho; it contributes nothing.
        if len(sub) < MIN_STRATUM_N or sub[feature].nunique() < 2:
            row[f"n_{key}"] = len(sub)
            row[f"rho_{key}"] = np.nan
            continue
        r = spearmanr(sub[feature], sub[target]).correlation
        row[f"n_{key}"] = len(sub)
        row[f"rho_{key}"] = float(r)
        if np.isfinite(r) and abs(r) < 1.0:
            zs.append(np.arctanh(r))
            ws.append(len(sub) - 3)

    row["k_strata"] = len(zs)
    if len(zs) < 2:
        for c in ("rho_pooled", "ci_low", "ci_high", "p_pooled", "pct_change",
                  "Q", "df_Q", "p_homogeneity", "I2"):
            row[c] = np.nan
        row["homogeneous"] = "n/a"
        row["confounded"] = "n/a"
        return row

    zs, ws = np.array(zs), np.array(ws, dtype=float)
    z_bar = float((ws * zs).sum() / ws.sum())
    se = float(1.0 / np.sqrt(ws.sum()))

    row["rho_pooled"] = float(np.tanh(z_bar))
    row["ci_low"] = float(np.tanh(z_bar - 1.96 * se))
    row["ci_high"] = float(np.tanh(z_bar + 1.96 * se))
    row["p_pooled"] = float(2 * norm.sf(abs(z_bar / se)))

    # Signed magnitude change; the sign of rho is not meaningful here, only size.
    row["pct_change"] = float(
        (abs(row["rho_pooled"]) - abs(row["rho_crude"])) / abs(row["rho_crude"]) * 100
    )

    Q = float((ws * (zs - z_bar) ** 2).sum())
    df_Q = len(zs) - 1
    row["Q"] = Q
    row["df_Q"] = df_Q
    row["p_homogeneity"] = float(chi2.sf(Q, df_Q))
    row["I2"] = float(max(0.0, (Q - df_Q) / Q) * 100) if Q > 0 else 0.0
    row["homogeneous"] = "\u2713" if row["p_homogeneity"] >= ALPHA else "\u2717"
    row["confounded"] = "\u2713" if abs(row["pct_change"]) >= CONFOUND_CHANGE * 100 else "\u2717"
    return row


def main() -> None:
    df = pd.read_csv(DATA_PATH)
    assert STRATUM_COL in df.columns
    print("Severity strata (all participants):")
    print(df[STRATUM_COL].value_counts().reindex(STRATA).to_string())

    all_rows = []
    for target, results_path in TARGETS.items():
        prior = pd.read_csv(results_path)
        sig = prior[prior["significant_fdr_0_05"]]
        print(f"\n{'=' * 110}")
        print(f"TARGET: {target}  -- re-estimating {len(sig)} crude-significant features "
              f"within {len(STRATA)} severity strata")
        print("=" * 110)

        rows = [stratified_pool(df, f, target) for f in sig["feature"]]
        res = pd.DataFrame(rows)

        valid = res["p_pooled"].notna()
        res["p_pooled_fdr_bh"] = np.nan
        res["pooled_significant_fdr_0_05"] = False
        if valid.any():
            reject, p_fdr, _, _ = multipletests(
                res.loc[valid, "p_pooled"].to_numpy(), alpha=ALPHA, method="fdr_bh"
            )
            res.loc[valid, "p_pooled_fdr_bh"] = p_fdr
            res.loc[valid, "pooled_significant_fdr_0_05"] = reject
        res["survives"] = np.where(res["pooled_significant_fdr_0_05"], "\u2713", "\u2717")

        res = res.reindex(res["rho_crude"].abs().sort_values(ascending=False).index)
        show = res.copy()
        for c, w in [("rho_crude", "crude"), ("rho_pooled", "pooled")]:
            show[w] = show[c].map(lambda v: f"{v:+.3f}" if np.isfinite(v) else "n/a")
        show["chg%"] = show["pct_change"].map(lambda v: f"{v:+.0f}" if np.isfinite(v) else "n/a")
        show["I2%"] = show["I2"].map(lambda v: f"{v:.0f}" if np.isfinite(v) else "n/a")
        show["pQ"] = show["p_homogeneity"].map(lambda v: f"{v:.3f}" if np.isfinite(v) else "n/a")
        show["pFDR"] = show["p_pooled_fdr_bh"].map(lambda v: f"{v:.3f}" if np.isfinite(v) else "n/a")
        print(
            show[["feature", "n_total", "crude", "pooled", "chg%", "pQ", "I2%",
                  "homogeneous", "confounded", "pFDR", "survives"]].to_string(index=False)
        )

        n = len(res)
        hom = int((res["homogeneous"] == "\u2713").sum())
        conf = int((res["confounded"] == "\u2713").sum())
        surv = int(res["pooled_significant_fdr_0_05"].sum())
        print(f"\n  homogeneous across strata : {hom}/{n}")
        print(f"  effect-modified (Q p<0.05): {n - hom}/{n}  <- do NOT pool these")
        print(f"  confounded (>={int(CONFOUND_CHANGE * 100)}% change) : {conf}/{n}")
        print(f"  pooled estimate survives FDR: {surv}/{n}")
        print(f"  median magnitude change   : {res['pct_change'].median():+.1f}%")
        all_rows.append(res)

    out = pd.concat(all_rows, ignore_index=True)
    out.to_csv(OUT_PATH, index=False)
    print(f"\nWrote: {OUT_PATH}  ({len(out)} rows)")


if __name__ == "__main__":
    main()
