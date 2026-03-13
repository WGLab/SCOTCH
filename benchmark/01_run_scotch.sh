#!/bin/bash
# SCOTCH Benchmark: Submits 4 dependent SLURM jobs (annotation -> compatible array -> summary -> count).
#
# Usage: bash 01_run_scotch.sh <read_label> <mode>
#   read_label: "5M", "10M", "15M"
#   mode:       "single" or "multi"
#
# Parallelism:
#   Annotation:        10 CPUs, 1 job
#   Compatible matrix: 1 CPU x 10 array jobs
#   Summary:           1 CPU, 1 job
#   Count matrix:      10 CPUs, 1 job
#
# Timing: /usr/bin/time -v is used inside each job to measure real execution time
# (not SLURM queue wait time).

set -euo pipefail

READ_LABEL=$1
MODE=${2:-single}
TOTAL_JOBS=10
WORKERS=10

# ===== USER CONFIG (edit these paths) =====
SCOTCH_DIR=/path/to/SCOTCH
REF_GTF=/path/to/genes.gtf
REF_PKL=/path/to/geneStructureInformation.pkl
DATA_DIR=/path/to/benchmark/data
RESULTS_BASE=/path/to/benchmark/results
PARTITION=default          # SLURM partition name
# ===========================================

OUT_BASE=${RESULTS_BASE}/scotch/${MODE}/reads_${READ_LABEL}

mkdir -p ${OUT_BASE} logs

# Build BAM and target arguments
if [ "${MODE}" == "single" ]; then
    BAM_ARGS="--bam ${DATA_DIR}/reads_${READ_LABEL}.bam"
    TARGET_ARGS="--target ${OUT_BASE}/sample1"
    mkdir -p ${OUT_BASE}/sample1
else
    # Multi-sample: 3 prefixed BAMs x 5M reads
    BAM_ARGS="--bam ${DATA_DIR}/multi_s1_5M.prefixed.bam ${DATA_DIR}/multi_s2_5M.prefixed.bam ${DATA_DIR}/multi_s3_5M.prefixed.bam"
    TARGET_ARGS="--target ${OUT_BASE}/sample1 ${OUT_BASE}/sample2 ${OUT_BASE}/sample3"
    mkdir -p ${OUT_BASE}/sample1 ${OUT_BASE}/sample2 ${OUT_BASE}/sample3
fi

echo "========================================"
echo "SCOTCH Benchmark (array mode)"
echo "Read level: ${READ_LABEL}, Mode: ${MODE}"
echo "Output: ${OUT_BASE}"
echo "========================================"

# --- Job 1: Annotation (10 CPUs) ---
JOB1=$(sbatch --parsable \
  --partition=${PARTITION} \
  --job-name=scotch_annot_${READ_LABEL}_${MODE} \
  --cpus-per-task=${WORKERS} --mem=300G --time=12:00:00 \
  --output=logs/scotch_annot_${READ_LABEL}_${MODE}_%j.out \
  --error=logs/scotch_annot_${READ_LABEL}_${MODE}_%j.err \
  --wrap="set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_annotation.txt \
    python3 ${SCOTCH_DIR}/src/main_preprocessing.py \
      --task annotation --platform 10x-ont \
      ${TARGET_ARGS} ${BAM_ARGS} \
      --reference ${REF_GTF} --reference_pkl ${REF_PKL} \
      --update_gtf --workers ${WORKERS} \
    2>&1 | tee ${OUT_BASE}/annotation.log")
echo "Submitted annotation: ${JOB1}"

# --- Job 2: Compatible matrix (10 array tasks x 1 CPU) ---
JOB2=$(sbatch --parsable \
  --dependency=afterok:${JOB1} \
  --partition=${PARTITION} \
  --job-name=scotch_compat_${READ_LABEL}_${MODE} \
  --array=0-$((TOTAL_JOBS - 1)) \
  --cpus-per-task=1 --mem=100G --time=12:00:00 \
  --output=logs/scotch_compat_${READ_LABEL}_${MODE}_%A_%a.out \
  --error=logs/scotch_compat_${READ_LABEL}_${MODE}_%A_%a.err \
  --wrap="set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_compatible_\${SLURM_ARRAY_TASK_ID}.txt \
    python3 ${SCOTCH_DIR}/src/main_preprocessing.py \
      --task 'compatible matrix' --platform 10x-ont \
      ${TARGET_ARGS} ${BAM_ARGS} \
      --reference ${REF_GTF} \
      --job_index \${SLURM_ARRAY_TASK_ID} --total_jobs ${TOTAL_JOBS} \
    2>&1 | tee ${OUT_BASE}/compatible_\${SLURM_ARRAY_TASK_ID}.log")
echo "Submitted compatible matrix array: ${JOB2}"

# --- Job 3: Summary (1 CPU) ---
JOB3=$(sbatch --parsable \
  --dependency=afterok:${JOB2} \
  --partition=${PARTITION} \
  --job-name=scotch_summary_${READ_LABEL}_${MODE} \
  --cpus-per-task=1 --mem=32G --time=02:00:00 \
  --output=logs/scotch_summary_${READ_LABEL}_${MODE}_%j.out \
  --error=logs/scotch_summary_${READ_LABEL}_${MODE}_%j.err \
  --wrap="set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_summary.txt \
    python3 ${SCOTCH_DIR}/src/main_preprocessing.py \
      --task summary ${TARGET_ARGS} \
    2>&1 | tee ${OUT_BASE}/summary.log")
echo "Submitted summary: ${JOB3}"

# --- Job 4: Count matrix (10 CPUs) ---
JOB4=$(sbatch --parsable \
  --dependency=afterok:${JOB3} \
  --partition=${PARTITION} \
  --job-name=scotch_count_${READ_LABEL}_${MODE} \
  --cpus-per-task=${WORKERS} --mem=200G --time=06:00:00 \
  --output=logs/scotch_count_${READ_LABEL}_${MODE}_%j.out \
  --error=logs/scotch_count_${READ_LABEL}_${MODE}_%j.err \
  --wrap="set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_count.txt \
    python3 ${SCOTCH_DIR}/src/main_preprocessing.py \
      --task 'count matrix' --platform 10x-ont \
      ${TARGET_ARGS} --workers ${WORKERS} --group_novel \
    2>&1 | tee ${OUT_BASE}/count.log && du -sh ${OUT_BASE}/ > ${OUT_BASE}/disk_usage.txt")
echo "Submitted count matrix: ${JOB4}"

# Save job IDs for tracking
echo -e "annotation\t${JOB1}\ncompatible\t${JOB2}\nsummary\t${JOB3}\ncount\t${JOB4}" \
  > ${OUT_BASE}/slurm_job_ids.txt
echo "Job chain submitted. IDs saved to ${OUT_BASE}/slurm_job_ids.txt"
