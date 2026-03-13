#!/usr/bin/env python3
"""Collect benchmark metrics from /usr/bin/time -v output files into a single TSV.

Handles SCOTCH (per-step + array tasks) and IsoQuant (single job).
For SCOTCH compatible matrix array: wall time = max across tasks, CPU-hours = sum.

Usage:
    python3 04_collect_metrics.py <results_dir> [-o benchmark_metrics.tsv]
"""

import os
import glob
import csv
import argparse


def parse_time_v(filepath):
    """Parse GNU time -v output file."""
    metrics = {}
    if not os.path.exists(filepath):
        return metrics
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if "wall clock" in line:
                val = line.split(": ")[-1]
                metrics["wall_seconds"] = hms_to_seconds(val)
            elif "Maximum resident set size" in line:
                kb = int(line.split(": ")[-1])
                metrics["peak_rss_gb"] = round(kb / 1024 / 1024, 2)
            elif "User time" in line:
                metrics["user_time_s"] = float(line.split(": ")[-1])
            elif "System time" in line:
                metrics["sys_time_s"] = float(line.split(": ")[-1])
            elif "Exit status" in line:
                metrics["exit_status"] = int(line.split(": ")[-1])
    return metrics


def hms_to_seconds(s):
    """Convert h:mm:ss or m:ss.ss to seconds."""
    parts = s.split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    return float(parts[0])


def cpu_hours(metrics):
    """Compute CPU-hours from user + system time."""
    user = metrics.get("user_time_s", 0)
    sys_t = metrics.get("sys_time_s", 0)
    return round((user + sys_t) / 3600, 4)


def collect_scotch(rep_path, bam_mode, read_label):
    """Collect metrics for SCOTCH (array mode with dependent jobs)."""
    rows = []
    total_wall = 0
    total_cpu_hours = 0
    peak_mem = 0

    # Non-array steps: annotation, summary, count
    for step in ["annotation", "summary", "count"]:
        tf = os.path.join(rep_path, f"time_{step}.txt")
        m = parse_time_v(tf)
        if not m:
            continue
        ch = cpu_hours(m)
        total_wall += m.get("wall_seconds", 0)
        total_cpu_hours += ch
        peak_mem = max(peak_mem, m.get("peak_rss_gb", 0))
        rows.append({
            "tool": "scotch", "step": step,
            "bam_mode": bam_mode, "read_label": read_label,
            "cpu_hours": ch, **m
        })

    # Compatible matrix array tasks
    compat_files = sorted(glob.glob(os.path.join(rep_path, "time_compatible_*.txt")))
    if compat_files:
        max_wall = 0
        sum_cpu_hours = 0
        max_mem = 0
        for cf in compat_files:
            m = parse_time_v(cf)
            max_wall = max(max_wall, m.get("wall_seconds", 0))
            sum_cpu_hours += cpu_hours(m)
            max_mem = max(max_mem, m.get("peak_rss_gb", 0))

        # Wall time = slowest task (parallel execution)
        total_wall += max_wall
        total_cpu_hours += sum_cpu_hours
        peak_mem = max(peak_mem, max_mem)
        rows.append({
            "tool": "scotch", "step": "compatible",
            "bam_mode": bam_mode, "read_label": read_label,
            "wall_seconds": max_wall, "peak_rss_gb": max_mem,
            "cpu_hours": sum_cpu_hours,
            "n_array_tasks": len(compat_files)
        })

    # Total row
    rows.append({
        "tool": "scotch", "step": "total",
        "bam_mode": bam_mode, "read_label": read_label,
        "wall_seconds": total_wall, "peak_rss_gb": peak_mem,
        "cpu_hours": total_cpu_hours
    })
    return rows


def collect_isoquant(rep_path, bam_mode, read_label):
    """Collect metrics for IsoQuant (dedup + isoquant steps).

    IsoQuant total = dedup time + isoquant run time.
    This is fair because SCOTCH's annotation step includes the same dedup logic internally.
    """
    rows = []
    total_wall = 0
    total_cpu_hours = 0
    peak_mem = 0

    # Per-step metrics: dedup and isoquant
    for step in ["dedup", "isoquant"]:
        tf = os.path.join(rep_path, f"time_{step}.txt")
        m = parse_time_v(tf)
        if not m:
            continue
        ch = cpu_hours(m)
        total_wall += m.get("wall_seconds", 0)
        total_cpu_hours += ch
        peak_mem = max(peak_mem, m.get("peak_rss_gb", 0))
        rows.append({
            "tool": "isoquant", "step": step,
            "bam_mode": bam_mode, "read_label": read_label,
            "cpu_hours": ch, **m
        })

    if rows:
        rows.append({
            "tool": "isoquant", "step": "total",
            "bam_mode": bam_mode, "read_label": read_label,
            "wall_seconds": total_wall, "peak_rss_gb": peak_mem,
            "cpu_hours": total_cpu_hours
        })

    return rows


def collect_all(results_dir):
    """Walk results directory and collect all metrics."""
    rows = []

    for tool_name, collector in [("scotch", collect_scotch), ("isoquant", collect_isoquant)]:
        for bam_mode in ["single", "multi"]:
            mode_dir = os.path.join(results_dir, tool_name, bam_mode)
            if not os.path.isdir(mode_dir):
                continue
            for reads_dir in sorted(os.listdir(mode_dir)):
                if not reads_dir.startswith("reads_"):
                    continue
                read_label = reads_dir.replace("reads_", "")
                reads_path = os.path.join(mode_dir, reads_dir)
                if not os.path.isdir(reads_path):
                    continue
                rows.extend(collector(reads_path, bam_mode, read_label))

    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect benchmark metrics into TSV")
    parser.add_argument("results_dir", help="Path to benchmark/results/")
    parser.add_argument("-o", "--output", default="benchmark_metrics.tsv")
    args = parser.parse_args()

    rows = collect_all(args.results_dir)
    if rows:
        fieldnames = [
            "tool", "step", "bam_mode", "read_label",
            "wall_seconds", "peak_rss_gb", "user_time_s", "sys_time_s",
            "cpu_hours", "exit_status", "n_array_tasks"
        ]
        with open(args.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"Wrote {len(rows)} rows to {args.output}")
    else:
        print("No metrics found.")
