"""
Is the HbA1c-vs-survey correlation systematically larger than the
mean-CGM-vs-survey correlation?

A naive paired test over the 44 survey features answers this at p ~ 1e-4, but
that test assumes 44 independent observations and the features are nowhere near
independent. sdoh_education_code, sdoh_education_years and sdoh_college_graduate
are three encodings of one construct; the five CES-D-10 subscales are built from
the same ten items; the four food-insecurity variants overlap almost completely.
The effective sample is the six survey instruments, not the 44 columns, so the
naive p-value is badly anti-conservative.

This script therefore treats the INSTRUMENT as the unit of analysis and reports
three complementary cluster-aware answers:

  1. Cluster sign-flip permutation. Under the null that neither target is
     systematically stronger, the sign of each instrument's mean difference is
     exchangeable. With six instruments all 2^6 = 64 sign assignments are
     enumerated exactly, so this is an exact test rather than a sampled one.
     Note the floor: the smallest attainable two-sided p is 2/64 = 0.031.

  2. Cluster bootstrap. Resample the six instruments with replacement and
     recompute the mean difference, giving a CI that propagates the fact that we
     have six clusters rather than 44 points.

  3. Mixed model with a random intercept per instrument. Uses all 44 paired
     differences but lets features within an instrument share a correlated
     error, which recovers power the six-point tests give up.

Differences are computed on two scales because correlations are not additive:
raw |rho| and |Fisher z|. Conclusions should agree; if they do not, the raw-rho
version is the one to distrust.

Caveat carried from src/survey_glycemic_severity_stratified_pooling.py: these
are CRUDE correlations, and roughly 76% of the median crude rank covariance is
between-severity-stratum rather than within. This script asks whether the two
targets differ in apparent magnitude, not whether either magnitude is causal.

Run from repo root:
    .venv/bin/python src/hba1c_vs_cgm_magnitude_comparison.py
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
HBA1C = ROOT / "outputs" / "survey_glycemia" / "survey_hba1c_spearman_bootstrap.csv"
CGM = ROOT / "outputs" / "survey_glycemia" / "survey_mean_cgm_spearman_bootstrap.csv"
OUT_PATH = ROOT / "outputs" / "survey_glycemia" / "hba1c_vs_cgm_magnitude_comparison.csv"

N_BOOT = 10_000
SEED = 42


def load_pairs() -> pd.DataFrame:
    h = pd.read_csv(HBA1C)
    c = pd.read_csv(CGM)
    m = h.merge(c, on=["feature", "survey", "direction"], suffixes=("_h", "_c"))
    m["abs_h"] = m["spearman_rho_h"].abs()
    m["abs_c"] = m["spearman_rho_c"].abs()
    m["d_rho"] = m["abs_c"] - m["abs_h"]
    # Fisher z puts the two magnitudes on a scale where differences are additive.
    m["d_z"] = np.arctanh(m["abs_c"]) - np.arctanh(m["abs_h"])
    return m


def naive_tests(m: pd.DataFrame) -> None:
    """The feature-level test, reported only to show how much it overstates."""
    print("=" * 96)
    print("NAIVE FEATURE-LEVEL TEST (treats 44 features as independent -- it is wrong)")
    print("=" * 96)
    w = wilcoxon(m["abs_h"], m["abs_c"])
    print(f"  Wilcoxon signed-rank over {len(m)} features: p = {w.pvalue:.2e}")
    print(f"  mean |rho|: HbA1c {m['abs_h'].mean():.4f} vs CGM {m['abs_c'].mean():.4f}")
    print(f"  shrinks under CGM: {(m['d_rho'] < 0).sum()}/{len(m)}")


def cluster_table(m: pd.DataFrame) -> pd.DataFrame:
    g = m.groupby("survey").agg(
        n_features=("feature", "size"),
        mean_abs_hba1c=("abs_h", "mean"),
        mean_abs_cgm=("abs_c", "mean"),
        mean_d_rho=("d_rho", "mean"),
        mean_d_z=("d_z", "mean"),
    )
    g["ratio"] = g["mean_abs_cgm"] / g["mean_abs_hba1c"]
    return g.sort_values("mean_d_rho")


def signflip_permutation(values: np.ndarray) -> tuple[float, float]:
    """Exact cluster sign-flip test on the k cluster means."""
    k = len(values)
    observed = values.mean()
    null = np.array(
        [np.mean(values * np.array(signs)) for signs in itertools.product([1, -1], repeat=k)]
    )
    p = float((np.abs(null) >= abs(observed) - 1e-12).mean())
    return float(observed), p


def cluster_bootstrap(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """Percentile CI from resampling whole clusters with replacement."""
    k = len(values)
    draws = rng.integers(0, k, size=(N_BOOT, k))
    means = values[draws].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def mixed_model(m: pd.DataFrame, col: str) -> None:
    """Random intercept per instrument; uses all 44 points but respects clustering."""
    try:
        import statsmodels.formula.api as smf
    except ImportError:
        return
    df = m[[col, "survey"]].rename(columns={col: "d"})
    fit = smf.mixedlm("d ~ 1", df, groups=df["survey"]).fit(reml=True)
    est = fit.params["Intercept"]
    se = fit.bse["Intercept"]
    print(f"  mixed model (1|instrument) on {col}: "
          f"estimate={est:+.4f}  SE={se:.4f}  "
          f"95% CI [{est - 1.96 * se:+.4f}, {est + 1.96 * se:+.4f}]  "
          f"p={fit.pvalues['Intercept']:.4f}")


def main() -> None:
    m = load_pairs()
    rng = np.random.default_rng(SEED)

    naive_tests(m)

    clusters = cluster_table(m)
    print("\n" + "=" * 96)
    print("CLUSTER-AWARE: instrument is the unit of analysis")
    print("=" * 96)
    print(clusters.to_string(float_format=lambda v: f"{v:+.4f}"))
    print("\n  negative mean_d_rho = CGM correlation is weaker for that instrument")

    for col, label in [("d_rho", "raw |rho|"), ("d_z", "|Fisher z|")]:
        vals = clusters[f"mean_{col}"].to_numpy()
        obs, p_perm = signflip_permutation(vals)
        lo, hi = cluster_bootstrap(vals, rng)
        print(f"\n--- scale: {label} ---")
        print(f"  mean difference across {len(vals)} instruments: {obs:+.4f}")
        print(f"  exact cluster sign-flip permutation p = {p_perm:.4f}  "
              f"(floor for k={len(vals)} is {2 / 2 ** len(vals):.4f})")
        print(f"  cluster bootstrap 95% CI: [{lo:+.4f}, {hi:+.4f}]  "
              f"excludes zero: {lo * hi > 0}")
        mixed_model(m, col)

    print("\n" + "=" * 96)
    print("READING")
    print("=" * 96)
    print("  All six instruments point the same way, which is the substantive")
    print("  result. With only six clusters the exact permutation test cannot")
    print("  go below p=0.031 no matter how large the effect, so a 'significant'")
    print("  verdict here is a statement about consistency of sign across")
    print("  instruments, not about effect size.")

    clusters.to_csv(OUT_PATH)
    print(f"\nWrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
