# SCOTCH vs IsoQuant Benchmarking Plan

## 1. Experiment Design

### 1.1 Overview

| Axis | Levels | Purpose |
|------|--------|---------|
| **Tool** | SCOTCH, IsoQuant | Head-to-head |
| **Read count** | 5M, 10M, 15M | Scalability curve (cell count reported post hoc) |
| **Sample mode** | Single (1 sample), Multi (3 samples x 5M = 15M total) | Multi-sample overhead |
| **Platform** | 10x-ONT | Single platform |

### 1.2 Experiment Matrix

| # | Tool | Read level | Mode | Description |
|---|------|-----------|------|-------------|
| 1 | SCOTCH | 5M | single | 1 BAM, 5M reads |
| 2 | SCOTCH | 10M | single | 1 BAM, 10M reads |
| 3 | SCOTCH | 15M | single | 1 BAM, 15M reads |
| 4 | SCOTCH | 15M | multi | 3 BAMs x 5M reads |
| 5 | IsoQuant | 5M | single | 1 BAM, 5M reads |
| 6 | IsoQuant | 10M | single | 1 BAM, 10M reads |
| 7 | IsoQuant | 15M | single | 1 BAM, 15M reads |
| 8 | IsoQuant | 15M | multi | 3 BAMs x 5M reads |

**Total: 8 experiments.** No replicates.

### 1.3 Parallelism Settings

| | SCOTCH | IsoQuant |
|---|---|---|
| **Annotation** | 1 job, 10 CPUs (`--workers 10`) | N/A (built-in) |
| **Compatible matrix** | 10 SLURM array jobs, 1 CPU each | N/A (built-in) |
| **Summary** | 1 job, 1 CPU | N/A (built-in) |
| **Count matrix** | 1 job, 10 CPUs (`--workers 10`) | N/A (built-in) |
| **IsoQuant** | N/A | 1 job, 10 threads (`--threads 10`) |

### 1.4 Dataset Preparation

**Random read subsampling** from real 10x-ONT BAMs using `samtools view -s SEED.FRACTION`.

| Level | Target reads | Used in |
|-------|-------------|---------|
| 5M | 5,000,000 | Single-sample scalability; multi-sample (per sample) |
| 10M | 10,000,000 | Single-sample scalability |
| 15M | 15,000,000 | Single-sample scalability |

For multi-sample: 3 independent BAMs subsampled to 5M reads each (different seeds), totaling 15M reads -- directly comparable to the 15M single-sample run.

### 1.5 Shared Resources

Both tools use the same:
- Reference genome FASTA (GRCh38)
- Gene annotation GTF (GENCODE or 10x genes.gtf)
- Pre-built annotation pickle for SCOTCH (`geneStructureInformation.pkl`)
- Identical subsampled BAM files

### 1.6 IsoQuant Multi-Sample: Barcode Collision Handling

**Problem:** IsoQuant's `--read_group tag:CB` merges identical CB tags from different BAMs into one cell. There is no built-in way to create composite `{sample}_{barcode}` identifiers.

**Solution:** Before running IsoQuant multi-sample, preprocess each BAM to prefix the CB tag with a sample identifier:
- `ACGTACGT-1` in sample 1 becomes `S1_ACGTACGT-1`
- `ACGTACGT-1` in sample 2 becomes `S2_ACGTACGT-1`

This is done by `00b_prefix_barcodes.py`, which rewrites the CB tag in each BAM. The same prefixed BAMs are also used for SCOTCH multi-sample to keep inputs identical.

After prefixing, IsoQuant can safely use `--read_group tag:CB` and each sample's cells are distinct.

### 1.7 UMI Deduplication: Fair Preprocessing

SCOTCH internally deduplicates reads during its annotation step: for each CB+UMI (`CB` + `UB` tags for 10x-ONT) combination within a sample, only the **longest read** is kept for downstream analysis. Shorter reads sharing the same CB+UMI are discarded.

IsoQuant does not perform this deduplication by default. To ensure a fair comparison, we preprocess IsoQuant's input BAMs with `00c_dedup_longest_read.py`, which applies the same logic:
1. Group all reads by `CB_UB` (cell barcode + UMI)
2. Keep only the longest read per group
3. Write a filtered BAM

**The dedup time is measured with `/usr/bin/time -v` and counted toward IsoQuant's total runtime**, since SCOTCH's annotation step includes this cost internally.

### 1.8 Timing: Real Run Time Only

All timing uses `/usr/bin/time -v` **inside** the SLURM job script, which measures actual process execution time (not SLURM queue wait time). SLURM's `sacct Elapsed` includes queue wait and is NOT used for benchmarking.

For SCOTCH array mode, the compatible matrix step has 10 parallel tasks. We record:
- **Wall time** = annotation + max(array task wall times) + summary + count
- **CPU-hours** = annotation CPU + sum(all array task CPU) + summary CPU + count CPU

---

## 2. Metrics to Record

| Metric | How to measure | Unit |
|--------|---------------|------|
| **Wall-clock time** | `/usr/bin/time -v` -> "Elapsed (wall clock) time" | seconds |
| **CPU time (user+sys)** | `/usr/bin/time -v` -> "User time" + "System time" | seconds |
| **Peak RSS memory** | `/usr/bin/time -v` -> "Maximum resident set size" (kB) | GB |
| **CPU-hours** | (user + sys time) / 3600 | hours |
| **Disk I/O written** | `du -sh <output_dir>` | GB |
| **Exit status** | `/usr/bin/time -v` -> "Exit status" | 0/non-0 |
| **Number of cells** | Unique CB tags in subsampled BAM | count |
| **Number of genes/isoforms** | From output count matrix | count |

### Per-step timing

**SCOTCH** (4 dependent SLURM jobs):

| Step | SLURM jobs | CPUs |
|------|-----------|------|
| Annotation (includes UMI dedup internally) | 1 job | 10 |
| Compatible matrix | 10 array jobs | 1 each |
| Summary | 1 job | 1 |
| Count matrix | 1 job | 10 |

**IsoQuant** (1 SLURM job, 2 timed steps):

| Step | Purpose | CPUs |
|------|---------|------|
| Dedup (`00c_dedup_longest_read.py`) | Keep longest read per CB+UMI (matches SCOTCH) | 1 |
| IsoQuant run | Transcript discovery + quantification | 10 |

IsoQuant total = dedup time + IsoQuant run time.

---

## 3. Scripts

| File | Purpose |
|------|---------|
| **`run_experiment.sh`** | **Main entry point. Edit config at top, then `bash run_experiment.sh`** |
| `00_prepare_data.sh` | Random read subsampling via `samtools view -s` |
| `00b_prefix_barcodes.py` | Prefix CB tags with sample ID for multi-sample |
| `00c_dedup_longest_read.py` | Keep longest read per CB+UMI (IsoQuant preprocessing) |
| `01_run_scotch.sh` | SCOTCH standalone runner (used by `run_experiment.sh` logic) |
| `02_run_isoquant.sh` | IsoQuant standalone runner (used by `run_experiment.sh` logic) |
| `04_collect_metrics.py` | Parse `/usr/bin/time -v` logs into TSV |
| `05_plot_results.py` | Generate figures |

---

## 4. Expected Output Tables

### Table 1: Runtime Scalability (single BAM)

| Reads | Cells | SCOTCH wall (min) | SCOTCH CPU-hrs | SCOTCH peak mem (GB) | IsoQuant wall (min) | IsoQuant CPU-hrs | IsoQuant peak mem (GB) |
|------:|------:|-------------------:|---------------:|---------------------:|--------------------:|-----------------:|----------------------:|
| 5M | | | | | | | |
| 10M | | | | | | | |
| 15M | | | | | | | |

### Table 2: SCOTCH Per-Step Breakdown (single BAM, 15M reads)

| Step | Wall time (min) | Peak RSS (GB) | CPU time (min) | CPUs |
|------|----------------:|--------------:|---------------:|-----:|
| Annotation | | | | 10 |
| Compatible matrix (max of 10 tasks) | | | | 1 x 10 |
| Summary | | | | 1 |
| Count matrix | | | | 10 |
| **Total** | | | | |

### Table 3: Multi-Sample Performance (3 samples x 5M reads = 15M total)

| Tool | Wall time (min) | CPU-hours | Peak mem (GB) |
|------|----------------:|----------:|--------------:|
| SCOTCH | | | |
| IsoQuant | | | |

### Table 4: Output Comparison (single BAM, 15M reads)

| Metric | SCOTCH | IsoQuant |
|--------|-------:|---------:|
| Cells detected | | |
| Genes detected | | |
| Known isoforms detected | | |
| Novel isoforms detected | | |
| Output disk size (GB) | | |

---

## 5. Directory Structure

```
benchmark/
├── BENCHMARK_PLAN.md                # This file
├── 00_prepare_data.sh               # Random read subsampling
├── 00b_prefix_barcodes.py           # Prefix CB tags for multi-sample
├── 01_run_scotch.sh                 # SCOTCH runner (array mode)
├── 02_run_isoquant.sh               # IsoQuant runner
├── 03_submit_all.sh                 # Master SLURM launcher
├── 04_collect_metrics.py            # Parse time/memory logs -> TSV
├── 05_plot_results.py               # Generate figures
├── data/                            # Subsampled BAMs
│   ├── reads_5M.bam
│   ├── reads_10M.bam
│   ├── reads_15M.bam
│   ├── multi_s1_5M.bam              # Multi-sample (original)
│   ├── multi_s2_5M.bam
│   ├── multi_s3_5M.bam
│   ├── multi_s1_5M.prefixed.bam     # Multi-sample (CB prefixed)
│   ├── multi_s2_5M.prefixed.bam
│   └── multi_s3_5M.prefixed.bam
├── results/
│   ├── scotch/
│   │   ├── single/reads_5M/
│   │   ├── single/reads_10M/
│   │   ├── single/reads_15M/
│   │   └── multi/reads_15M/
│   └── isoquant/
│       ├── single/reads_5M/
│       ├── single/reads_10M/
│       ├── single/reads_15M/
│       └── multi/reads_15M/
├── logs/                            # SLURM stdout/stderr
├── benchmark_metrics.tsv            # Collected metrics
└── benchmark_results.pdf            # Final figure
```

## 6. How to Run

```bash
cd benchmark/

# 1. Edit the USER CONFIG section at the top of run_experiment.sh:
#    - Set all paths (SCOTCH_DIR, REF_FASTA, REF_GTF, etc.)
#    - Set INPUT_BAM_DIR to the folder containing your 10x-ONT BAM(s)
#    - Adjust SLURM resource settings (memory, CPUs, time limits) if needed
#    - Email is pre-set to cyranvvv@hotmail.com

# 2. Run:
bash run_experiment.sh

# This automatically:
#   Phase 0: Subsamples BAMs + prefixes barcodes for multi-sample
#   Phase 1-4: Submits all 8 experiments with afterok dependencies
#   Phase 5: Collects metrics + generates plots after all experiments finish
#
# Monitor: squeue -u $(whoami)
# You'll get email on job FAIL or END.
```

## 7. Practical Notes

- **Timing accuracy**: `/usr/bin/time -v` measures real process time, not SLURM queue wait. This is critical for accurate benchmarking.
- **Fair CPU allocation**: SCOTCH gets 10 CPUs for annotation/count, 10 x 1 CPU for compatible matrix. IsoQuant gets 10 threads. Both use comparable resources.
- **Multi-sample barcode prefixing**: Required for IsoQuant to distinguish cells across samples. Applied to all tools for consistency.
- **Summary step**: SCOTCH requires a summary step after compatible matrix to merge per-job outputs before generating count matrices.
- **Cell count reporting**: After subsampling, count unique CB tags: `samtools view subset.bam | awk '{for(i=12;i<=NF;i++) if($i~/^CB:Z:/) print substr($i,6)}' | sort -u | wc -l`.
