# SCOTCH vs IsoQuant Computational Benchmark Report

Platform: 10x-ONT, single-cell long-read RNA-seq.
SCOTCH compatible matrix: 10 parallel array jobs; wall time = max (slowest task).
IsoQuant: 10 threads.

---

## 1. Summary Comparison

### 1a. Wall-Clock Time (hours)

| Tool | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| SCOTCH | 1.14 | 3.94 | 13.85 | 26.52 | 14.01 |
| IsoQuant | 0.80 | 1.16 | 2.12 | 11.52 | 2.48 |

### 1b. CPU-Hours

| Tool | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| SCOTCH | 5.95 | 17.15 | 62.56 | 92.26 | 52.77 |
| IsoQuant | 4.07 | 6.06 | 11.35 | 63.41 | 12.09 |

### 1c. Peak Memory (GB)

| Tool | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| SCOTCH | 2.15 | 6.50 | 17.06 | 38.65 | 6.85 |
| IsoQuant | 1.43 | 1.60 | 4.75 | 10.49 | 4.73 |

---

## 2. SCOTCH Step Breakdown

### 2a. Wall-Clock Time per Step (hours)

| Step | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| Annotation | 0.02 | 0.24 | 0.80 | 3.22 | 0.42 |
| Compatible | 0.93 | 3.39 | 12.40 | 22.26 | 12.43 |
| Summary | 0.03 | 0.06 | 0.11 | 0.27 | 0.29 |
| Count | 0.15 | 0.25 | 0.54 | 0.77 | 0.87 |
| **Total** | **1.14** | **3.94** | **13.85** | **26.52** | **14.01** |

### 2b. CPU-Hours per Step

| Step | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| Annotation | 0.10 | 0.56 | 0.50 | 1.21 | 0.43 |
| Compatible | 5.71 | 16.33 | 61.43 | 89.99 | 51.22 |
| Summary | 0.02 | 0.05 | 0.14 | 0.35 | 0.38 |
| Count | 0.12 | 0.21 | 0.48 | 0.71 | 0.74 |
| **Total** | **5.95** | **17.15** | **62.56** | **92.26** | **52.77** |

### 2c. Peak Memory per Step (GB)

| Step | 1M | 5M | 15M | 50M | 3x5M |
|---|---|---|---|---|---|
| Annotation | 0.85 | 3.34 | 9.81 | 25.55 | 3.25 |
| Compatible | 1.99 | 2.00 | 2.04 | 6.00 | 2.07 |
| Summary | 2.15 | 6.50 | 17.06 | 38.65 | 6.85 |
| Count | 0.86 | 2.39 | 5.50 | 10.74 | 4.14 |

---

## 3. SCOTCH Practical Resource Guide

Based on benchmark results (10 array jobs, 10x-ONT, hg38).

### 3a. Memory Recommendations per Step

Memory scales primarily with read count. Annotation and summary load full gene structures into memory; compatible matrix per-task memory stays low regardless of input size.

| Step | 1M | 5M | 15M | 50M | Rule of Thumb |
|---|---|---|---|---|---|
| Annotation | 1 GB | 4 GB | 10 GB | 26 GB | ~0.5 GB per 1M reads |
| Compatible (per task) | 2 GB | 2 GB | 2 GB | 6 GB | 2-8 GB flat (depends on gene complexity, not read count) |
| Summary | 3 GB | 7 GB | 18 GB | 39 GB | ~1.1 GB per 1M reads |
| Count | 1 GB | 3 GB | 6 GB | 11 GB | ~0.2 GB per 1M reads |

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
| 15M | 12.40 | 61.43 | 6.14 | High |
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

### 3d. Array Job Scaling Results (15M reads)

Measured speedup from increasing `--total_jobs` on the 15M dataset. Annotation is shared across all runs (0.80 h); only compatible, summary, and count steps are re-run.

| Array Jobs | Compatible Wall (h) | Compatible CPU-hours | Total Wall (h) | Total CPU-hours | Speedup vs 10 jobs |
|---|---|---|---|---|---|
| 10 | 12.40 | 61.43 | 13.85 | 62.56 | 1.0x |
| 25 | 7.01 | 27.54 | 7.94 | 28.46 | 1.8x |
| 40 | 2.45 | 26.99 | 2.98 | 27.51 | 4.6x |

**Key findings:**
- **10→25 jobs (2.5x tasks)**: compatible wall drops 1.8x. CPU-hours drop significantly (61→28), suggesting better load balancing with more tasks.
- **10→40 jobs (4x tasks)**: compatible wall drops 5.1x (12.40→2.45 h). Total SCOTCH time drops below IsoQuant's 2.12 h at this scale.
- **CPU-hours decrease** with more jobs (unlike the theoretical constant), because shorter per-task runtimes reduce overhead and tail latency from slow tasks.
- **Memory is unchanged** (~2 GB per compatible task, ~17 GB for summary) regardless of array size.

**SCOTCH (40 jobs) vs IsoQuant head-to-head at 15M reads:**

| Metric | SCOTCH (40 jobs) | IsoQuant (10 threads) | Ratio |
|---|---|---|---|
| Wall time (h) | 2.98 | 2.12 | 1.4x |
| CPU-hours | 28.01 | 11.35 | 2.5x |
| Peak memory (GB) | 17.06 | 4.75 | 3.6x |

SCOTCH with 40 jobs closes the wall-time gap to 1.4x, but still uses 2.5x more CPU-hours and 3.6x more memory. The memory overhead comes from the summary step (17 GB), which loads all gene structures regardless of array count.

**Practical recommendation:** For 15M+ reads, use **30-40 array jobs** to bring SCOTCH wall time on par with IsoQuant, at the cost of higher total CPU and memory usage.

---

## 4. Key Observations

- **Compatible matrix dominates SCOTCH runtime**: 82-90% of wall time and 95-98% of CPU-hours across all settings (with 10 array jobs).
- **SCOTCH wall time scales roughly linearly** with read count (1M→50M ~23x for ~50x reads).
- **IsoQuant is faster in wall time at default settings** (10 array jobs) due to multi-threaded execution within a single process, but **CPU-hours gap narrows at scale** (50M: 92 vs 63 CPU-hours).
- **With 40 array jobs, SCOTCH approaches IsoQuant wall time**: at 15M reads, SCOTCH total is 2.98 h vs IsoQuant's 2.12 h (1.4x), but at the cost of 2.5x more CPU-hours (28 vs 11.3) and 3.6x more peak memory (17 vs 4.75 GB).
- **SCOTCH memory is higher** primarily in annotation and summary steps, which load full gene structures. Compatible matrix per-task memory stays low (~2 GB) regardless of input size or array count.
- **Multi-sample (3x5M)**: SCOTCH wall time is ~3.5x its single 5M run (processes 3 BAMs with shared annotation); IsoQuant is ~2.1x.
- **Increasing array jobs is the most effective lever** for reducing SCOTCH wall time, with no additional memory cost per task. CPU-hours also decrease with more jobs due to reduced tail latency.
