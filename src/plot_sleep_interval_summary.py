"""
Build the same style of summary figure as for physical activity intervals:
  - duration distribution (binned minutes)
  - gap distribution
  - overlap KPIs (n intervals, any overlap %, big overlap %)

Reads sleep_interval_summary.csv (from cadence_analysis.main()).
Run from repo root: python3 src/plot_sleep_interval_summary.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from config import OUTDIR


def main() -> None:
    path = OUTDIR / "sleep_interval_summary.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run: python3 src/cadence_analysis.py (or main) to generate it."
        )

    df = pd.read_csv(path)
    dur = df[df["kind"] == "duration"].reset_index(drop=True)
    gap = df[df["kind"] == "gap"].reset_index(drop=True)
    ov = df[df["kind"] == "overlap"].copy()

    fig, axes = plt.subplots(2, 1, figsize=(9, 5.5), layout="constrained")
    fig.suptitle("Sleep stage intervals — duration & gaps (same bins as cadence_analysis)", fontsize=12)

    x = range(len(dur))
    c0 = dur["count"].fillna(0)
    ymax0 = float(c0.max()) if len(c0) else 0.0
    axes[0].bar(x, c0, color="#4472c4", edgecolor="white", linewidth=0.5)
    axes[0].set_xticks(list(x))
    axes[0].set_xticklabels(dur["bin"].astype(str), rotation=0, fontsize=9)
    axes[0].set_ylabel("Count")
    axes[0].set_title("Interval duration")
    for xi in range(len(dur)):
        row = dur.iloc[xi]
        h = float(row["count"]) if pd.notna(row["count"]) else 0.0
        pct = float(row["pct"]) if pd.notna(row["pct"]) else 0.0
        axes[0].text(xi, h + ymax0 * 0.02, f"{pct:.1f}%", ha="center", va="bottom", fontsize=7)

    xg = range(len(gap))
    cg = gap["count"].fillna(0)
    ymax1 = float(cg.max()) if len(cg) else 0.0
    axes[1].bar(xg, cg, color="#70ad47", edgecolor="white", linewidth=0.5)
    axes[1].set_xticks(list(xg))
    axes[1].set_xticklabels(gap["bin"].astype(str), rotation=0, fontsize=9)
    axes[1].set_ylabel("Count")
    axes[1].set_title("Gap to previous interval (same person, by start order)")
    for xi in range(len(gap)):
        row = gap.iloc[xi]
        h = float(row["count"]) if pd.notna(row["count"]) else 0.0
        pct = float(row["pct"]) if pd.notna(row["pct"]) else 0.0
        axes[1].text(xi, h + ymax1 * 0.02, f"{pct:.1f}%", ha="center", va="bottom", fontsize=7)

    # Overlap text box
    ov_map = dict(zip(ov["metric"], ov["value"]))
    lines = [
        "Overlap (vs immediately prior interval by start time)",
        f"  n_intervals:              {ov_map.get('n_intervals', float('nan')):,.0f}",
        f"  n_with_any_overlap:      {ov_map.get('n_with_any_overlap', float('nan')):,.0f}",
        f"  pct_with_any_overlap:    {ov_map.get('pct_with_any_overlap', float('nan')):.4f}%",
        f"  n_with_big_overlap:      {ov_map.get('n_with_big_overlap', float('nan')):,.0f}",
        f"  pct_with_big_overlap:    {ov_map.get('pct_with_big_overlap', float('nan')):.4f}%",
    ]
    fig.text(
        0.5,
        0.02,
        "\n".join(lines),
        ha="center",
        va="bottom",
        fontsize=9,
        family="monospace",
        transform=fig.transFigure,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.35),
    )

    out = OUTDIR / "sleep_interval_summary_chart.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
