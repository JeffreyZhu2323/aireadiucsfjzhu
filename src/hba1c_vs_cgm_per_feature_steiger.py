"""
Per-feature test of whether a survey feature correlates more strongly with lab
HbA1c than with mean CGM.

src/hba1c_vs_cgm_magnitude_comparison.py answers this at the aggregate level:
across the six survey instruments, correlations with mean CGM are weaker by
about 0.014 in |rho|. That says nothing about any individual feature, so claims
like "education attenuates by half under CGM" are unsupported by it. This
script supplies the per-feature evidence.

The two correlations cannot be compared with independent tests. They share the
survey feature, they are measured on the same participants, and HbA1c and mean
CGM correlate with each other at rho ~ 0.68. That is the overlapping dependent
correlation case, and the standard test is Williams' t (Steiger 1980, eq. 14).

Two things are reported for every feature:

  1. Williams' t on the signed difference rho(feat, HbA1c) - rho(feat, CGM),
     with df = n - 3 and Benjamini-Hochberg across the 44 features. This is the
     named significance test. Caveat: Williams' t is derived for Pearson
     correlations. Applying it to Spearman is common practice but remains an
     approximation, which is why (2) exists.

  2. A paired bootstrap CI for the difference in |rho|, resampling participants
     with replacement. This is distribution-free, inherits the dependency
     between the two correlations automatically, and is the number to trust if
     it disagrees with Williams.

Both correlations are recomputed on the COMMON sample -- participants with a
non-null feature, HbA1c and mean CGM. The published tables use pairwise
deletion (HbA1c n=2182, CGM n=2245), and a paired comparison across different
participant sets is not a paired comparison. The original pairwise estimates
are carried through as rho_*_pairwise so the shift is visible.

Sign handling: the null is rho_hba1c == rho_cgm, so the test runs on signed
correlations. For the handful of features where the two correlations have
opposite signs, the signed difference is not a magnitude difference; those are
flagged signs_agree = False and the |rho| columns should not be read as
attenuation.

Caveat carried from src/survey_glycemic_severity_stratified_pooling.py: these
are CRUDE correlations, and roughly 76% of the median crude rank covariance is
between-severity-stratum. This asks which target shows the larger apparent
association, not whether either is causal.

Run from repo root:
    .venv/bin/python src/hba1c_vs_cgm_per_feature_steiger.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, t as tdist
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"
HBA1C_TABLE = ROOT / "outputs" / "survey_glycemia" / "survey_hba1c_spearman_bootstrap.csv"
CGM_TABLE = ROOT / "outputs" / "survey_glycemia" / "survey_mean_cgm_spearman_bootstrap.csv"
OUT_PATH = ROOT / "outputs" / "survey_glycemia" / "hba1c_vs_cgm_per_feature_steiger.csv"

HBA1C = "hba1c_percent"
CGM = "mean_cgm_overall_mg_dl"

N_BOOT = 10_000
BOOT_CHUNK = 500
SEED = 42
ALPHA = 0.05
MIN_N = 30


def spearman_axis(x: np.ndarray, y: np.ndarray, axis: int = -1) -> np.ndarray:
    """Spearman rho along `axis`, vectorised over bootstrap resamples.

    Ranks are recomputed within each resample because resampling with
    replacement creates ties that change the ranks.
    """
    rx = rankdata(x, axis=axis).astype(float)
    ry = rankdata(y, axis=axis).astype(float)
    rx -= rx.mean(axis=axis, keepdims=True)
    ry -= ry.mean(axis=axis, keepdims=True)
    num = (rx * ry).sum(axis=axis)
    den = np.sqrt((rx**2).sum(axis=axis) * (ry**2).sum(axis=axis))
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def williams_t(r12: float, r13: float, r23: float, n: int) -> tuple[float, float]:
    """Williams' t for two dependent correlations sharing variable 1.

    r12 = rho(feature, HbA1c), r13 = rho(feature, CGM), r23 = rho(HbA1c, CGM).
    """
    det = 1 - r12**2 - r13**2 - r23**2 + 2 * r12 * r13 * r23
    rbar = (r12 + r13) / 2
    denom = 2 * ((n - 1) / (n - 3)) * det + rbar**2 * (1 - r23) ** 3
    if denom <= 0:
        return float("nan"), float("nan")
    tval = (r12 - r13) * np.sqrt(((n - 1) * (1 + r23)) / denom)
    return float(tval), float(2 * tdist.sf(abs(tval), n - 3))


def bootstrap_diff(
    feat: np.ndarray, h: np.ndarray, c: np.ndarray, rng: np.random.Generator
) -> tuple[float, float]:
    """Percentile CI for |rho(feat, HbA1c)| - |rho(feat, CGM)|.

    Participants are resampled as whole rows, so the two correlations stay
    coupled exactly as they are in the data. Chunked to keep the resample
    matrices small enough to rank in memory.
    """
    n = len(feat)
    diffs = np.empty(N_BOOT)
    done = 0
    while done < N_BOOT:
        size = min(BOOT_CHUNK, N_BOOT - done)
        idx = rng.integers(0, n, size=(size, n))
        rh = spearman_axis(feat[idx], h[idx])
        rc = spearman_axis(feat[idx], c[idx])
        diffs[done : done + size] = np.abs(rh) - np.abs(rc)
        done += size
    finite = diffs[np.isfinite(diffs)]
    if finite.size < N_BOOT // 2:
        return float("nan"), float("nan")
    return float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))


def main() -> None:
    df = pd.read_csv(DATA_PATH, low_memory=False)
    h_tab = pd.read_csv(HBA1C_TABLE)
    c_tab = pd.read_csv(CGM_TABLE)
    meta = h_tab[["feature", "survey", "direction"]].merge(
        c_tab[["feature", "spearman_rho", "significant_fdr_0_05"]].rename(
            columns={
                "spearman_rho": "rho_cgm_pairwise",
                "significant_fdr_0_05": "sig_cgm_crude",
            }
        ),
        on="feature",
    )
    meta = meta.merge(
        h_tab[["feature", "spearman_rho", "significant_fdr_0_05"]].rename(
            columns={
                "spearman_rho": "rho_hba1c_pairwise",
                "significant_fdr_0_05": "sig_hba1c_crude",
            }
        ),
        on="feature",
    )

    rng = np.random.default_rng(SEED)
    rows = []
    for i, m in enumerate(meta.itertuples(index=False), start=1):
        col = m.feature
        sub = df[[col, HBA1C, CGM]].apply(pd.to_numeric, errors="coerce").dropna()
        n = len(sub)
        print(f"  [{i:2d}/{len(meta)}] {col:<34} n={n:<5}", end="", flush=True)

        if n < MIN_N or sub[col].nunique() < 2:
            print(" skipped (degenerate)")
            continue

        feat = sub[col].to_numpy(float)
        hv = sub[HBA1C].to_numpy(float)
        cv = sub[CGM].to_numpy(float)

        r12 = float(spearmanr(feat, hv).correlation)
        r13 = float(spearmanr(feat, cv).correlation)
        r23 = float(spearmanr(hv, cv).correlation)
        tval, pval = williams_t(r12, r13, r23, n)
        lo, hi = bootstrap_diff(feat, hv, cv, rng)

        rows.append(
            {
                "feature": col,
                "survey": m.survey,
                "direction": m.direction,
                "n_common": n,
                "rho_hba1c": r12,
                "rho_cgm": r13,
                "rho_hba1c_cgm": r23,
                "rho_hba1c_pairwise": m.rho_hba1c_pairwise,
                "rho_cgm_pairwise": m.rho_cgm_pairwise,
                "sig_hba1c_crude": m.sig_hba1c_crude,
                "sig_cgm_crude": m.sig_cgm_crude,
                "signs_agree": bool(np.sign(r12) == np.sign(r13)),
                "diff_signed": r12 - r13,
                "diff_abs": abs(r12) - abs(r13),
                "williams_t": tval,
                "p_value": pval,
                "boot_ci_low": lo,
                "boot_ci_high": hi,
                "boot_excludes_zero": bool(np.isfinite(lo) and lo * hi > 0),
            }
        )
        print(f" rho_h={r12:+.3f} rho_c={r13:+.3f} t={tval:+.2f} p={pval:.3g}")

    res = pd.DataFrame(rows)
    valid = res["p_value"].notna()
    res["p_value_fdr_bh"] = np.nan
    res["significant_fdr_0_05"] = False
    if valid.any():
        reject, p_fdr, _, _ = multipletests(
            res.loc[valid, "p_value"].to_numpy(), alpha=ALPHA, method="fdr_bh"
        )
        res.loc[valid, "p_value_fdr_bh"] = p_fdr
        res.loc[valid, "significant_fdr_0_05"] = reject
    res = res.sort_values("p_value").reset_index(drop=True)
    res.to_csv(OUT_PATH, index=False)

    n_sig = int(res["significant_fdr_0_05"].sum())
    print("\n" + "=" * 100)
    print(f"WILLIAMS' t: {n_sig}/{len(res)} features differ between HbA1c and mean CGM (FDR < {ALPHA})")
    print("=" * 100)
    show = res[res["significant_fdr_0_05"]] if n_sig else res.head(10)
    cols = [
        "feature", "survey", "n_common", "rho_hba1c", "rho_cgm",
        "diff_abs", "williams_t", "p_value_fdr_bh", "boot_ci_low",
        "boot_ci_high", "boot_excludes_zero",
    ]
    print(show[cols].to_string(index=False, float_format=lambda v: f"{v:+.4f}"))

    disagree = res[res["significant_fdr_0_05"] != res["boot_excludes_zero"]]
    print(f"\nWilliams and bootstrap disagree on {len(disagree)} features"
          f"{': ' + ', '.join(disagree['feature']) if len(disagree) else ''}")
    flipped = res[~res["signs_agree"]]
    print(f"Opposite-sign features (diff_abs not an attenuation): {len(flipped)}"
          f"{': ' + ', '.join(flipped['feature']) if len(flipped) else ''}")
    print(f"\nWrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
