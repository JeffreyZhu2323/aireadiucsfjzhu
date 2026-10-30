"""
Visual summary of sleep interval overlaps (the ~2495 rows from sleep_overlap_deep_dive).

Reads sleep_overlap_pair_summary.csv (run sleep_overlap_deep_dive.py first).
Same chart style idea as plot_sleep_interval_summary.py: bar charts + KPI text.

Output: sleep_overlap_pairs_chart.png in outputs folder.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from config import OUTDIR


def main() -> None:
    pair_path = OUTDIR / "archive" / "exploratory_2026-05-01" / "legacy_outputs" / "sleep_overlap_pair_summary.csv"
    rows_path = OUTDIR / "archive" / "exploratory_2026-05-01" / "legacy_outputs" / "sleep_overlap_interval_rows.csv"
    if not pair_path.exists() or not rows_path.exists():
        raise FileNotFoundError(
            f"Missing {pair_path} or {rows_path}. Run: python3 src/sleep_overlap_deep_dive.py"
        )

    sig = pd.read_csv(pair_path)
    detail = pd.read_csv(rows_path)
    n_ov = len(detail)

    bucket_counts = detail["pair_bucket"].value_counts().reset_index()
    bucket_counts.columns = ["pair_bucket", "n_rows"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), layout="constrained")
    fig.suptitle(
        f"Sleep interval overlaps (n = {n_ov:,}): prior vs current stage (by start order)",
        fontsize=12,
    )

    # Left: pair_signature horizontal bars (sorted by count)
    sig_plot = sig.sort_values("n_rows", ascending=True)
    y = range(len(sig_plot))
    axes[0].barh(y, sig_plot["n_rows"], color="#4472c4", edgecolor="white", linewidth=0.5)
    axes[0].set_yticks(list(y))
    axes[0].set_yticklabels(sig_plot["pair_signature"].astype(str), fontsize=8)
    axes[0].set_xlabel("Count")
    axes[0].set_title("Overlap pair (prev|curr)")
    xmax = float(sig_plot["n_rows"].max()) if len(sig_plot) else 1.0
    for yi, (_, row) in enumerate(sig_plot.iterrows()):
        axes[0].text(
            row["n_rows"] + xmax * 0.01,
            yi,
            f"{row['pct_of_overlaps']:.1f}%",
            va="center",
            fontsize=7,
        )

    # Right: pair_bucket
    b = bucket_counts.sort_values("n_rows", ascending=True)
    yb = range(len(b))
    axes[1].barh(yb, b["n_rows"], color="#70ad47", edgecolor="white", linewidth=0.5)
    axes[1].set_yticks(list(yb))
    axes[1].set_yticklabels(b["pair_bucket"].astype(str), fontsize=8)
    axes[1].set_xlabel("Count")
    axes[1].set_title("Pair bucket (heuristic stage order: awake→light→deep→rem)")
    xmax_b = float(b["n_rows"].max()) if len(b) else 1.0
    for yi, (_, row) in enumerate(b.iterrows()):
        pct = 100.0 * row["n_rows"] / n_ov if n_ov else 0.0
        axes[1].text(row["n_rows"] + xmax_b * 0.01, yi, f"{pct:.1f}%", va="center", fontsize=7)

    note = (
        "Overlap rule matches cadence_analysis.summarize_interval_cadence: "
        "within person, sort by start; row i overlaps if start_i < end_{i-1}."
    )
    fig.text(0.5, 0.01, note, ha="center", va="bottom", fontsize=8, style="italic")

    out = OUTDIR / "figures" / "sleep_overlap_pairs_chart.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
