# SCOTCH vs IsoQuant Computational Benchmark Report

Platform: 10x-ONT, single-cell long-read RNA-seq.
SCOTCH compatible matrix: 10 parallel array jobs; wall time = max (slowest task).
IsoQuant: 10 threads.

---

## 1. Summary Comparison

### 1a. Wall-Clock Time (hours)

| Tool | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| SCOTCH | 1.14 | 3.94 | — | 26.52 | 14.01 |
| IsoQuant | 0.80 | 1.16 | — | 11.52 | 2.48 |

### 1b. CPU-Hours

| Tool | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| SCOTCH | 5.95 | 17.15 | — | 92.26 | 52.77 |
| IsoQuant | 4.07 | 6.06 | — | 63.41 | 12.09 |

### 1c. Peak Memory (GB)

| Tool | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| SCOTCH | 2.15 | 6.50 | — | 38.65 | 6.85 |
| IsoQuant | 1.43 | 1.60 | — | 10.49 | 4.73 |

---

## 2. SCOTCH Step Breakdown

### 2a. Wall-Clock Time per Step (hours)

| Step | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| Annotation | 0.02 | 0.24 | — | 3.22 | 0.42 |
| Compatible | 0.93 | 3.39 | — | 22.26 | 12.43 |
| Summary | 0.03 | 0.06 | — | 0.27 | 0.29 |
| Count | 0.15 | 0.25 | — | 0.77 | 0.87 |
| **Total** | **1.14** | **3.94** | **—** | **26.52** | **14.01** |

### 2b. CPU-Hours per Step

| Step | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| Annotation | 0.10 | 0.56 | — | 1.21 | 0.43 |
| Compatible | 5.71 | 16.33 | — | 89.99 | 51.22 |
| Summary | 0.02 | 0.05 | — | 0.35 | 0.38 |
| Count | 0.12 | 0.21 | — | 0.71 | 0.74 |
| **Total** | **5.95** | **17.15** | **—** | **92.26** | **52.77** |

### 2c. Peak Memory per Step (GB)

| Step | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| Annotation | 0.85 | 3.34 | — | 25.55 | 3.25 |
| Compatible | 1.99 | 2.00 | — | 6.00 | 2.07 |
| Summary | 2.15 | 6.50 | — | 38.65 | 6.85 |
| Count | 0.86 | 2.39 | — | 10.74 | 4.14 |

---

## 3. SCOTCH Practical Resource Guide

Based on benchmark results (10 array jobs, 10x-ONT, hg38).

### 3a. Memory Recommendations per Step

Memory scales primarily with read count. Annotation and summary load full gene structures into memory; compatible matrix per-task memory stays low regardless of input size.

| Step | 1M | 5M | 50M | Rule of Thumb |
|---|---|---|---|---|
| Annotation | 1 GB | 4 GB | 26 GB | ~0.5 GB per 1M reads |
| Compatible (per task) | 2 GB | 2 GB | 6 GB | 2-8 GB flat (depends on gene complexity, not read count) |
| Summary | 3 GB | 7 GB | 39 GB | ~0.8 GB per 1M reads |
| Count | 1 GB | 3 GB | 11 GB | ~0.2 GB per 1M reads |

**Recommended SLURM `--mem` (with ~2x safety margin):**

| Reads | Annotation | Compatible (per task) | Summary | Count |
|---|---|---|---|---|
| 1-5M | 8G | 4G | 16G | 8G |
| 5-15M | 20G | 8G | 36G | 12G |
| 15-50M | 50G | 12G | 80G | 24G |
| >50M | 80G | 16G | 120G | 40G |

### 3b. CPU Recommendations

| Step | Parallelism Model | Recommended CPUs |
|---|---|---|
| Annotation | Internal multiprocessing (`--workers`) | 8-16 |
| Compatible | SLURM job array (1 CPU per task) | 1 per task |
| Summary | Single-threaded merge | 1-2 |
| Count | Internal joblib (`--workers`) | 4-8 |

### 3c. Job Array Sizing for Compatible Matrix

The compatible matrix step is the bottleneck (82-90% of wall time). All benchmarks above used **10 array jobs**. The wall time equals the slowest task, so more array jobs = smaller chunks = lower max wall time.

**Current results with 10 array jobs:**

| Reads | Compatible Wall (h) | Compatible CPU-hours | Avg per task (h) | Speedup potential |
|---|---|---|---|---|
| 1M | 0.93 | 5.71 | 0.57 | Moderate — tasks already short |
| 5M | 3.39 | 16.33 | 1.63 | Good |
| 50M | 22.26 | 89.99 | 9.00 | High — large per-task variance expected |

**Theoretical scaling:** With N array jobs, wall time ~ CPU-hours / N (ideal) but is bounded by the slowest task. Diminishing returns occur when:
- Per-task overhead (loading BAM index, annotation pickle) dominates compute time
- Gene count is unevenly distributed across metagenes (some tasks get heavy genes)
- N approaches the number of metagenes (~15k for hg38)

**Rough guidance:**

| Reads | Array Jobs | Expected Wall Reduction vs 10 jobs |
|---|---|---|
| 1-5M | 10-20 | Marginal (tasks already <1h) |
| 5-15M | 20-30 | ~2x faster compatible wall |
| 15-50M | 30-50 | ~3-5x faster compatible wall |
| >50M | 50-100 | Significant; per-task mem stays low |

**Note:** More array jobs consume the same total CPU-hours — the savings are in wall time only. Cluster scheduling overhead may increase with very large arrays (>100 tasks).

### 3d. Proposed Experiment: Array Job Scaling

To quantify the actual speedup curve, we propose running the 15M setting with varying array sizes:

| Experiment | `--total_jobs` | Expected compatible wall (h) | Notes |
|---|---|---|---|
| baseline | 10 | ~12 | Current default |
| 2x | 20 | ~6 | Doubling should nearly halve |
| 3x | 30 | ~4 | Standard HPC sweet spot |
| 5x | 50 | ~2.5 | Diminishing returns expected |

This will produce a speedup curve showing where the practical limit is, which informs the recommended default for users.

---

## 4. Key Observations

- **Compatible matrix dominates SCOTCH runtime**: 82-90% of wall time and 95-98% of CPU-hours across all settings.
- **SCOTCH wall time scales roughly linearly** with read count (1M->50M ~23x for ~50x reads).
- **IsoQuant is faster in wall time** due to multi-threaded execution within a single process, but **CPU-hours gap narrows at scale** (50M: 92 vs 63 CPU-hours).
- **SCOTCH memory is higher** primarily in annotation and summary steps, which load full gene structures. Compatible matrix per-task memory stays low (~2 GB) regardless of input size.
- **Multi-sample (3x5M)**: SCOTCH wall time is ~3.5x its single 5M run (processes 3 BAMs with shared annotation); IsoQuant is ~2.1x.
- **Increasing array jobs is the most effective lever** for reducing SCOTCH wall time, with no additional memory cost per task.
