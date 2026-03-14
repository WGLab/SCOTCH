#!/usr/bin/env python3
"""Generate benchmark comparison plots from collected metrics.

Usage:
    python3 05_plot_results.py benchmark_metrics.tsv
"""

import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

READ_ORDER = {"1M": 1, "5M": 5, "10M": 10, "50M": 50}


def load_data(path):
    df = pd.read_csv(path, sep="\t")
    df["read_millions"] = df["read_label"].map(READ_ORDER)
    return df


def highest_single_label(df):
    single = df[df["bam_mode"] == "single"]
    if single.empty:
        return None
    labels = single["read_label"].dropna().unique().tolist()
    labels.sort(key=lambda label: READ_ORDER.get(label, -1))
    return labels[-1] if labels else None


def single_xticks(df):
    single = df[df["bam_mode"] == "single"]["read_millions"].dropna()
    if single.empty:
        return sorted(READ_ORDER.values())
    return sorted(single.unique())


def plot_scalability_wall(df, ax):
    """Wall time vs read count (single BAM, totals)."""
    single = df[(df["bam_mode"] == "single") & (df["step"] == "total")]
    colors = {"scotch": "#1f77b4", "isoquant": "#ff7f0e"}
    markers = {"scotch": "o", "isoquant": "s"}
    labels = {"scotch": "SCOTCH", "isoquant": "IsoQuant"}

    for tool in ["scotch", "isoquant"]:
        sub = single[single["tool"] == tool].sort_values("read_millions")
        if sub.empty:
            continue
        ax.plot(sub["read_millions"], sub["wall_seconds"] / 60,
                marker=markers[tool], color=colors[tool], label=labels[tool])

    ax.set_xlabel("Reads (millions)")
    ax.set_ylabel("Wall time (minutes)")
    ax.set_title("Runtime Scalability")
    ax.legend()
    ax.set_xticks(single_xticks(df))


def plot_scalability_mem(df, ax):
    """Peak memory vs read count (single BAM, totals)."""
    single = df[(df["bam_mode"] == "single") & (df["step"] == "total")]
    colors = {"scotch": "#1f77b4", "isoquant": "#ff7f0e"}

    for tool in ["scotch", "isoquant"]:
        sub = single[single["tool"] == tool].sort_values("read_millions")
        if sub.empty:
            continue
        label = "SCOTCH" if tool == "scotch" else "IsoQuant"
        ax.plot(sub["read_millions"], sub["peak_rss_gb"],
                marker="s", color=colors[tool], label=label)

    ax.set_xlabel("Reads (millions)")
    ax.set_ylabel("Peak memory (GB)")
    ax.set_title("Memory Scalability")
    ax.legend()
    ax.set_xticks(single_xticks(df))


def plot_cpu_hours(df, ax):
    """CPU-hours vs read count (single BAM, totals)."""
    single = df[(df["bam_mode"] == "single") & (df["step"] == "total")]
    colors = {"scotch": "#1f77b4", "isoquant": "#ff7f0e"}

    for tool in ["scotch", "isoquant"]:
        sub = single[single["tool"] == tool].sort_values("read_millions")
        if sub.empty:
            continue
        ax.plot(sub["read_millions"], sub["cpu_hours"],
                marker="^", color=colors[tool], label="SCOTCH" if tool == "scotch" else "IsoQuant")

    ax.set_xlabel("Reads (millions)")
    ax.set_ylabel("CPU-hours")
    ax.set_title("CPU-Hours Scalability")
    ax.legend()
    ax.set_xticks(single_xticks(df))


def plot_scotch_breakdown(df, ax):
    """SCOTCH per-step breakdown for the largest single-sample input."""
    focus_label = highest_single_label(df)
    sub = df[
        (df["tool"] == "scotch") &
        (df["bam_mode"] == "single") &
        (df["read_label"] == focus_label) &
        (~df["step"].isin(["total"]))
    ]

    if sub.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    step_order = ["annotation", "compatible", "summary", "count"]
    sub = sub.set_index("step").reindex(step_order).dropna(subset=["wall_seconds"])
    all_colors = {"annotation": "#1f77b4", "compatible": "#ff7f0e", "summary": "#2ca02c", "count": "#d62728"}
    colors = [all_colors[s] for s in sub.index]

    ax.barh(range(len(sub)), sub["wall_seconds"] / 60, color=colors, alpha=0.8)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([s.capitalize() for s in sub.index])
    ax.set_xlabel("Wall time (minutes)")
    ax.set_title(f"SCOTCH Step Breakdown ({focus_label})" if focus_label else "SCOTCH Step Breakdown")
    ax.invert_yaxis()


def plot_multisample(df, ax):
    """Compare the largest single-sample benchmark with the multi-sample condition."""
    totals = df[df["step"] == "total"]
    focus_single = highest_single_label(df)

    data = []
    for tool in ["scotch", "isoquant"]:
        single = totals[
            (totals["tool"] == tool) &
            (totals["bam_mode"] == "single") &
            (totals["read_label"] == focus_single)
        ]
        multi = totals[
            (totals["tool"] == tool) &
            (totals["bam_mode"] == "multi")
        ]
        if not single.empty:
            data.append({
                "label": f"{'SCOTCH' if tool == 'scotch' else 'IsoQuant'}\n(single {focus_single})",
                "wall_min": single.iloc[0]["wall_seconds"] / 60,
            })
        if not multi.empty:
            data.append({
                "label": f"{'SCOTCH' if tool == 'scotch' else 'IsoQuant'}\n(multi {multi.iloc[0]['read_label']})",
                "wall_min": multi.iloc[0]["wall_seconds"] / 60,
            })

    if not data:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    labels = [d["label"] for d in data]
    vals = [d["wall_min"] for d in data]
    colors = ["#1f77b4", "#1f77b4", "#ff7f0e", "#ff7f0e"][:len(data)]
    alphas = [0.8, 0.5, 0.8, 0.5][:len(data)]

    bars = ax.bar(range(len(data)), vals, color=colors, alpha=0.8)
    for i, b in enumerate(bars):
        b.set_alpha(alphas[i])
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Wall time (minutes)")
    ax.set_title("Single vs Multi-Sample")


def main():
    parser = argparse.ArgumentParser(description="Generate benchmark plots")
    parser.add_argument("metrics_tsv", help="Path to benchmark_metrics.tsv")
    parser.add_argument(
        "-o",
        "--output-prefix",
        default=None,
        help="Output path prefix (default: <metrics dir>/benchmark_results)",
    )
    args = parser.parse_args()
    df = load_data(args.metrics_tsv)
    output_prefix = args.output_prefix or str(Path(args.metrics_tsv).resolve().parent / "benchmark_results")

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    plot_scalability_wall(df, axes[0, 0])
    plot_scalability_mem(df, axes[0, 1])
    plot_cpu_hours(df, axes[0, 2])
    plot_scotch_breakdown(df, axes[1, 0])
    plot_multisample(df, axes[1, 1])

    # Panel F: Summary table as text
    ax = axes[1, 2]
    ax.axis("off")
    totals = df[(df["step"] == "total") & (df["bam_mode"] == "single")]
    if not totals.empty:
        table_data = []
        for _, row in totals.iterrows():
            tool = "SCOTCH" if row["tool"] == "scotch" else "IsoQuant"
            table_data.append([
                tool, row["read_label"],
                f'{row["wall_seconds"] / 60:.1f}',
                f'{row.get("cpu_hours", 0):.2f}',
                f'{row.get("peak_rss_gb", 0):.1f}'
            ])
        table = ax.table(
            cellText=table_data,
            colLabels=["Tool", "Reads", "Wall (min)", "CPU-hrs", "Mem (GB)"],
            loc="center", cellLoc="center"
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.5)
    ax.set_title("Summary", fontsize=11)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.pdf", dpi=150, bbox_inches="tight")
    plt.savefig(f"{output_prefix}.png", dpi=150, bbox_inches="tight")
    print(f"Saved {output_prefix}.pdf / .png")


if __name__ == "__main__":
    main()
