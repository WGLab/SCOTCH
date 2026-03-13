#!/bin/bash
# =============================================================================
# SCOTCH vs IsoQuant Benchmark — Master Experiment Runner
# =============================================================================
#
# This script runs all 8 benchmark experiments with proper SLURM dependency
# chains. Each experiment only starts after its prerequisites complete
# successfully (afterok).
#
# Usage:
#   bash run_experiment.sh
#
# Workflow:
#   Phase 0: Data preparation (subsampling + barcode prefixing + dedup)
#   Phase 1: SCOTCH single-sample (5M, 10M, 15M) — run in parallel
#   Phase 2: SCOTCH multi-sample (3 x 5M)
#   Phase 3: IsoQuant single-sample (5M, 10M, 15M) — run in parallel
#   Phase 4: IsoQuant multi-sample (3 x 5M)
#   Phase 5: Collect metrics (after all experiments finish)
#
# All SLURM jobs send email on FAIL and END to the configured address.
# =============================================================================

set -euo pipefail

# =============================================================================
# USER CONFIG — Edit these before running
# =============================================================================

# --- Email for SLURM notifications (FAIL + END) ---
EMAIL="cyranvvv@hotmail.com"

# --- Paths ---
SCOTCH_DIR=/path/to/SCOTCH                          # Root of SCOTCH repo
BENCHMARK_DIR=/path/to/benchmark                    # This benchmark/ directory
REF_FASTA=/path/to/genome.fa                        # Reference genome FASTA
REF_GTF=/path/to/genes.gtf                          # Gene annotation GTF
REF_PKL=/path/to/geneStructureInformation.pkl       # SCOTCH annotation pickle
INPUT_BAM_DIR=/path/to/bam_folder                   # Directory with input BAM(s)
DATA_DIR=/path/to/benchmark/data                    # Where subsampled BAMs go
RESULTS_BASE=/path/to/benchmark/results             # Where results go

# --- SLURM settings ---
PARTITION=default                                   # SLURM partition name
ACCOUNT=""                                          # SLURM account (leave empty if N/A)

# --- Resource settings: SCOTCH ---
SCOTCH_ANNOT_CPUS=10                                # CPUs for annotation step
SCOTCH_ANNOT_MEM=300G                               # Memory for annotation step
SCOTCH_ANNOT_TIME=12:00:00                          # Time limit for annotation

SCOTCH_COMPAT_NJOBS=10                              # Number of array jobs for compatible matrix
SCOTCH_COMPAT_MEM=100G                              # Memory per array task
SCOTCH_COMPAT_TIME=12:00:00                         # Time limit per array task

SCOTCH_SUMMARY_MEM=32G                              # Memory for summary step
SCOTCH_SUMMARY_TIME=02:00:00                        # Time limit for summary

SCOTCH_COUNT_CPUS=10                                # CPUs for count matrix step
SCOTCH_COUNT_MEM=200G                               # Memory for count matrix
SCOTCH_COUNT_TIME=06:00:00                          # Time limit for count matrix

# --- Resource settings: IsoQuant ---
ISOQUANT_CPUS=10                                    # Threads for IsoQuant
ISOQUANT_MEM=200G                                   # Memory for IsoQuant
ISOQUANT_TIME=24:00:00                              # Time limit for IsoQuant

# --- Resource settings: Data preparation ---
PREP_CPUS=4                                         # CPUs for subsampling
PREP_MEM=32G                                        # Memory for subsampling
PREP_TIME=06:00:00                                  # Time limit for subsampling

# =============================================================================
# END USER CONFIG — Do not edit below unless you know what you're doing
# =============================================================================

# Build common sbatch flags
SBATCH_COMMON="--partition=${PARTITION}"
if [ -n "${ACCOUNT}" ]; then
    SBATCH_COMMON="${SBATCH_COMMON} --account=${ACCOUNT}"
fi
SBATCH_COMMON="${SBATCH_COMMON} --mail-user=${EMAIL} --mail-type=FAIL,END"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p logs "${DATA_DIR}" "${RESULTS_BASE}"

echo "============================================================"
echo "SCOTCH vs IsoQuant Benchmark — Submitting all experiments"
echo "============================================================"
echo "Email:   ${EMAIL}"
echo "Data:    ${DATA_DIR}"
echo "Results: ${RESULTS_BASE}"
echo ""

# -------------------------------------------------------
# Helper: submit SCOTCH experiment (4 dependent jobs)
# Usage: submit_scotch <read_label> <mode> [dep_jobid]
# Returns: final job ID (count matrix) in SCOTCH_LAST_JOB
# -------------------------------------------------------
submit_scotch() {
    local READ_LABEL=$1
    local MODE=$2
    local DEP_FLAG=""
    if [ "${3:-}" != "" ]; then
        DEP_FLAG="--dependency=afterok:$3"
    fi

    local OUT_BASE="${RESULTS_BASE}/scotch/${MODE}/reads_${READ_LABEL}"
    mkdir -p "${OUT_BASE}" logs

    # Build BAM and target arguments
    local BAM_ARGS TARGET_ARGS
    if [ "${MODE}" == "single" ]; then
        BAM_ARGS="--bam ${DATA_DIR}/reads_${READ_LABEL}.bam"
        TARGET_ARGS="--target ${OUT_BASE}/sample1"
        mkdir -p "${OUT_BASE}/sample1"
    else
        BAM_ARGS="--bam ${DATA_DIR}/multi_s1_5M.prefixed.bam ${DATA_DIR}/multi_s2_5M.prefixed.bam ${DATA_DIR}/multi_s3_5M.prefixed.bam"
        TARGET_ARGS="--target ${OUT_BASE}/sample1 ${OUT_BASE}/sample2 ${OUT_BASE}/sample3"
        mkdir -p "${OUT_BASE}/sample1" "${OUT_BASE}/sample2" "${OUT_BASE}/sample3"
    fi

    # Job 1: Annotation
    local JOB1
    JOB1=$(sbatch --parsable ${SBATCH_COMMON} ${DEP_FLAG} \
      --job-name=scotch_annot_${READ_LABEL}_${MODE} \
      --cpus-per-task=${SCOTCH_ANNOT_CPUS} --mem=${SCOTCH_ANNOT_MEM} --time=${SCOTCH_ANNOT_TIME} \
      --output=logs/scotch_annot_${READ_LABEL}_${MODE}_%j.out \
      --error=logs/scotch_annot_${READ_LABEL}_${MODE}_%j.err \
      --wrap="set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_annotation.txt \
        python3 ${SCOTCH_DIR}/src/main_preprocessing.py \
          --task annotation --platform 10x-ont \
          ${TARGET_ARGS} ${BAM_ARGS} \
          --reference ${REF_GTF} --reference_pkl ${REF_PKL} \
          --update_gtf --workers ${SCOTCH_ANNOT_CPUS} \
        2>&1 | tee ${OUT_BASE}/annotation.log")
    echo "    Annotation: ${JOB1}"

    # Job 2: Compatible matrix (array)
    local JOB2
    JOB2=$(sbatch --parsable ${SBATCH_COMMON} \
      --dependency=afterok:${JOB1} \
      --job-name=scotch_compat_${READ_LABEL}_${MODE} \
      --array=0-$((SCOTCH_COMPAT_NJOBS - 1)) \
      --cpus-per-task=1 --mem=${SCOTCH_COMPAT_MEM} --time=${SCOTCH_COMPAT_TIME} \
      --output=logs/scotch_compat_${READ_LABEL}_${MODE}_%A_%a.out \
      --error=logs/scotch_compat_${READ_LABEL}_${MODE}_%A_%a.err \
      --wrap="set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_compatible_\${SLURM_ARRAY_TASK_ID}.txt \
        python3 ${SCOTCH_DIR}/src/main_preprocessing.py \
          --task 'compatible matrix' --platform 10x-ont \
          ${TARGET_ARGS} ${BAM_ARGS} \
          --reference ${REF_GTF} \
          --job_index \${SLURM_ARRAY_TASK_ID} --total_jobs ${SCOTCH_COMPAT_NJOBS} \
        2>&1 | tee ${OUT_BASE}/compatible_\${SLURM_ARRAY_TASK_ID}.log")
    echo "    Compatible matrix array: ${JOB2}"

    # Job 3: Summary
    local JOB3
    JOB3=$(sbatch --parsable ${SBATCH_COMMON} \
      --dependency=afterok:${JOB2} \
      --job-name=scotch_summary_${READ_LABEL}_${MODE} \
      --cpus-per-task=1 --mem=${SCOTCH_SUMMARY_MEM} --time=${SCOTCH_SUMMARY_TIME} \
      --output=logs/scotch_summary_${READ_LABEL}_${MODE}_%j.out \
      --error=logs/scotch_summary_${READ_LABEL}_${MODE}_%j.err \
      --wrap="set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_summary.txt \
        python3 ${SCOTCH_DIR}/src/main_preprocessing.py \
          --task summary ${TARGET_ARGS} \
        2>&1 | tee ${OUT_BASE}/summary.log")
    echo "    Summary: ${JOB3}"

    # Job 4: Count matrix
    local JOB4
    JOB4=$(sbatch --parsable ${SBATCH_COMMON} \
      --dependency=afterok:${JOB3} \
      --job-name=scotch_count_${READ_LABEL}_${MODE} \
      --cpus-per-task=${SCOTCH_COUNT_CPUS} --mem=${SCOTCH_COUNT_MEM} --time=${SCOTCH_COUNT_TIME} \
      --output=logs/scotch_count_${READ_LABEL}_${MODE}_%j.out \
      --error=logs/scotch_count_${READ_LABEL}_${MODE}_%j.err \
      --wrap="set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_count.txt \
        python3 ${SCOTCH_DIR}/src/main_preprocessing.py \
          --task 'count matrix' --platform 10x-ont \
          ${TARGET_ARGS} --workers ${SCOTCH_COUNT_CPUS} --group_novel \
        2>&1 | tee ${OUT_BASE}/count.log && du -sh ${OUT_BASE}/ > ${OUT_BASE}/disk_usage.txt")
    echo "    Count matrix: ${JOB4}"

    echo -e "annotation\t${JOB1}\ncompatible\t${JOB2}\nsummary\t${JOB3}\ncount\t${JOB4}" \
      > "${OUT_BASE}/slurm_job_ids.txt"

    SCOTCH_LAST_JOB=${JOB4}
}

# -------------------------------------------------------
# Helper: submit IsoQuant experiment (1 job with dedup)
# Usage: submit_isoquant <read_label> <mode> [dep_jobid]
# Returns: job ID in ISOQUANT_LAST_JOB
# -------------------------------------------------------
submit_isoquant() {
    local READ_LABEL=$1
    local MODE=$2
    local DEP_FLAG=""
    if [ "${3:-}" != "" ]; then
        DEP_FLAG="--dependency=afterok:$3"
    fi

    local OUT_BASE="${RESULTS_BASE}/isoquant/${MODE}/reads_${READ_LABEL}"
    local DEDUP_DIR="${OUT_BASE}/dedup_bams"
    mkdir -p "${OUT_BASE}" "${DEDUP_DIR}" logs

    # Build the full command inline (dedup + isoquant)
    local DEDUP_CMD ISOQUANT_CMD INPUT_ARG LABELS_ARG
    if [ "${MODE}" == "single" ]; then
        DEDUP_CMD="python3 ${BENCHMARK_DIR}/00c_dedup_longest_read.py \
          ${DATA_DIR}/reads_${READ_LABEL}.bam ${DEDUP_DIR}/reads_${READ_LABEL}.dedup.bam --cb CB --umi UB"
        INPUT_ARG="--bam ${DEDUP_DIR}/reads_${READ_LABEL}.dedup.bam"
        LABELS_ARG="--labels sample1"
    else
        DEDUP_CMD="python3 ${BENCHMARK_DIR}/00c_dedup_longest_read.py \
          ${DATA_DIR}/multi_s1_5M.prefixed.bam ${DEDUP_DIR}/multi_s1_5M.dedup.bam --cb CB --umi UB && \
        python3 ${BENCHMARK_DIR}/00c_dedup_longest_read.py \
          ${DATA_DIR}/multi_s2_5M.prefixed.bam ${DEDUP_DIR}/multi_s2_5M.dedup.bam --cb CB --umi UB && \
        python3 ${BENCHMARK_DIR}/00c_dedup_longest_read.py \
          ${DATA_DIR}/multi_s3_5M.prefixed.bam ${DEDUP_DIR}/multi_s3_5M.dedup.bam --cb CB --umi UB"
        INPUT_ARG="--bam ${DEDUP_DIR}/multi_s1_5M.dedup.bam ${DEDUP_DIR}/multi_s2_5M.dedup.bam ${DEDUP_DIR}/multi_s3_5M.dedup.bam"
        LABELS_ARG="--labels S1 S2 S3"
    fi

    local JOBID
    JOBID=$(sbatch --parsable ${SBATCH_COMMON} ${DEP_FLAG} \
      --job-name=isoquant_${READ_LABEL}_${MODE} \
      --cpus-per-task=${ISOQUANT_CPUS} --mem=${ISOQUANT_MEM} --time=${ISOQUANT_TIME} \
      --output=logs/isoquant_${READ_LABEL}_${MODE}_%j.out \
      --error=logs/isoquant_${READ_LABEL}_${MODE}_%j.err \
      --wrap="set -euo pipefail; \
        /usr/bin/time -v -o ${OUT_BASE}/time_dedup.txt bash -c '${DEDUP_CMD}' 2>&1 | tee ${OUT_BASE}/dedup.log && \
        /usr/bin/time -v -o ${OUT_BASE}/time_isoquant.txt \
          isoquant.py \
            --reference ${REF_FASTA} --genedb ${REF_GTF} --complete_genedb \
            ${INPUT_ARG} --data_type nanopore --barcoded_bam \
            --barcode_tag CB --umi_tag UB ${LABELS_ARG} \
            --read_group tag:CB --threads ${ISOQUANT_CPUS} \
            -o ${OUT_BASE}/output \
          2>&1 | tee ${OUT_BASE}/isoquant.log && \
        du -sh ${OUT_BASE}/ > ${OUT_BASE}/disk_usage.txt")
    echo "    IsoQuant: ${JOBID}"

    ISOQUANT_LAST_JOB=${JOBID}
}

# =============================================================================
# PHASE 0: Data Preparation
# =============================================================================
echo "--- Phase 0: Data Preparation ---"

# Step 0a: Subsample BAMs
echo "  Submitting: subsample BAMs..."
PREP_JOB=$(sbatch --parsable ${SBATCH_COMMON} \
  --job-name=bench_prep \
  --cpus-per-task=${PREP_CPUS} --mem=${PREP_MEM} --time=${PREP_TIME} \
  --output=logs/prep_%j.out --error=logs/prep_%j.err \
  --wrap="bash ${SCRIPT_DIR}/00_prepare_data.sh ${INPUT_BAM_DIR} ${DATA_DIR}")
echo "    Subsample: ${PREP_JOB}"

# Step 0b: Prefix barcodes for multi-sample
echo "  Submitting: prefix barcodes..."
PREFIX_JOB=$(sbatch --parsable ${SBATCH_COMMON} \
  --dependency=afterok:${PREP_JOB} \
  --job-name=bench_prefix \
  --cpus-per-task=1 --mem=16G --time=04:00:00 \
  --output=logs/prefix_%j.out --error=logs/prefix_%j.err \
  --wrap="python3 ${SCRIPT_DIR}/00b_prefix_barcodes.py ${DATA_DIR}")
echo "    Prefix barcodes: ${PREFIX_JOB}"

# =============================================================================
# PHASE 1: SCOTCH Single-Sample (5M, 10M, 15M) — parallel
# =============================================================================
echo ""
echo "--- Phase 1: SCOTCH Single-Sample ---"

ALL_SCOTCH_JOBS=()
for LABEL in 5M 10M 15M; do
    echo "  SCOTCH single ${LABEL}:"
    submit_scotch "${LABEL}" single "${PREFIX_JOB}"
    ALL_SCOTCH_JOBS+=("${SCOTCH_LAST_JOB}")
done

# =============================================================================
# PHASE 2: SCOTCH Multi-Sample (3 x 5M)
# =============================================================================
echo ""
echo "--- Phase 2: SCOTCH Multi-Sample ---"
echo "  SCOTCH multi 15M:"
submit_scotch 15M multi "${PREFIX_JOB}"
ALL_SCOTCH_JOBS+=("${SCOTCH_LAST_JOB}")

# =============================================================================
# PHASE 3: IsoQuant Single-Sample (5M, 10M, 15M) — parallel
# =============================================================================
echo ""
echo "--- Phase 3: IsoQuant Single-Sample ---"

ALL_ISOQUANT_JOBS=()
for LABEL in 5M 10M 15M; do
    echo "  IsoQuant single ${LABEL}:"
    submit_isoquant "${LABEL}" single "${PREFIX_JOB}"
    ALL_ISOQUANT_JOBS+=("${ISOQUANT_LAST_JOB}")
done

# =============================================================================
# PHASE 4: IsoQuant Multi-Sample (3 x 5M)
# =============================================================================
echo ""
echo "--- Phase 4: IsoQuant Multi-Sample ---"
echo "  IsoQuant multi 15M:"
submit_isoquant 15M multi "${PREFIX_JOB}"
ALL_ISOQUANT_JOBS+=("${ISOQUANT_LAST_JOB}")

# =============================================================================
# PHASE 5: Collect Metrics (after ALL experiments finish)
# =============================================================================
echo ""
echo "--- Phase 5: Metrics Collection ---"

# Build dependency string: afterok:job1:job2:job3:...
ALL_JOBS=("${ALL_SCOTCH_JOBS[@]}" "${ALL_ISOQUANT_JOBS[@]}")
DEP_STR=$(IFS=:; echo "${ALL_JOBS[*]}")

COLLECT_JOB=$(sbatch --parsable ${SBATCH_COMMON} \
  --dependency=afterok:${DEP_STR} \
  --job-name=bench_collect \
  --cpus-per-task=1 --mem=8G --time=00:30:00 \
  --output=logs/collect_%j.out --error=logs/collect_%j.err \
  --wrap="python3 ${SCRIPT_DIR}/04_collect_metrics.py ${RESULTS_BASE} -o ${RESULTS_BASE}/benchmark_metrics.tsv && \
    python3 ${SCRIPT_DIR}/05_plot_results.py ${RESULTS_BASE}/benchmark_metrics.tsv && \
    echo 'Benchmark complete. Results in ${RESULTS_BASE}/benchmark_metrics.tsv'")
echo "  Collect metrics: ${COLLECT_JOB}"

echo ""
echo "============================================================"
echo "All jobs submitted! Dependency chain:"
echo "  Prep (${PREP_JOB}) -> Prefix (${PREFIX_JOB})"
echo "  -> SCOTCH 5M/10M/15M/multi + IsoQuant 5M/10M/15M/multi (parallel)"
echo "  -> Collect metrics (${COLLECT_JOB})"
echo ""
echo "Monitor: squeue -u \$(whoami)"
echo "Email notifications on FAIL/END -> ${EMAIL}"
echo "============================================================"
