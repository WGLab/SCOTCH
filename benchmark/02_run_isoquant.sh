#!/bin/bash
#SBATCH --job-name=bench_isoquant
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/isoquant_%j.out
#SBATCH --error=logs/isoquant_%j.err

# IsoQuant Benchmark: Single SLURM job with 10 threads.
#
# Includes a preprocessing step to deduplicate reads by keeping only the longest
# read per CB+UB (cell barcode + UMI) combination. This matches SCOTCH's internal
# deduplication during its annotation step, ensuring a fair comparison.
# The dedup time is counted toward IsoQuant's total runtime.
#
# Usage: sbatch 02_run_isoquant.sh <read_label> <mode>
#   read_label: "5M", "10M", "15M"
#   mode:       "single" or "multi"
#
# Timing: /usr/bin/time -v measures real execution time (not queue wait).

set -euo pipefail

READ_LABEL=$1
MODE=${2:-single}
THREADS=10

# ===== USER CONFIG (edit these paths) =====
BENCHMARK_DIR=/path/to/benchmark
REF_FASTA=/path/to/genome.fa
REF_GTF=/path/to/genes.gtf
DATA_DIR=/path/to/benchmark/data
RESULTS_BASE=/path/to/benchmark/results
# ===========================================

OUT_BASE=${RESULTS_BASE}/isoquant/${MODE}/reads_${READ_LABEL}
DEDUP_DIR=${OUT_BASE}/dedup_bams

mkdir -p ${OUT_BASE} ${DEDUP_DIR} logs

echo "========================================"
echo "IsoQuant Benchmark"
echo "Read level: ${READ_LABEL}, Mode: ${MODE}"
echo "Output: ${OUT_BASE}"
echo "========================================"

# -------------------------------------------------------
# Step 1: Deduplicate BAMs (keep longest read per CB+UB)
# This mirrors SCOTCH's annotation-step deduplication.
# Time is counted toward IsoQuant total.
# -------------------------------------------------------
echo "[$(date)] Step 1: Deduplicating BAMs (longest read per CB+UB)..."

if [ "${MODE}" == "single" ]; then
    ORIG_BAM=${DATA_DIR}/reads_${READ_LABEL}.bam
    DEDUP_BAM=${DEDUP_DIR}/reads_${READ_LABEL}.dedup.bam

    /usr/bin/time -v -o ${OUT_BASE}/time_dedup.txt \
      python3 ${BENCHMARK_DIR}/00c_dedup_longest_read.py \
        ${ORIG_BAM} ${DEDUP_BAM} --cb CB --umi UB \
      2>&1 | tee ${OUT_BASE}/dedup.log

    INPUT="--bam ${DEDUP_BAM}"
    LABELS="--labels sample1"
else
    # Multi-sample: dedup each prefixed BAM separately
    /usr/bin/time -v -o ${OUT_BASE}/time_dedup.txt bash -c "
      python3 ${BENCHMARK_DIR}/00c_dedup_longest_read.py \
        ${DATA_DIR}/multi_s1_5M.prefixed.bam ${DEDUP_DIR}/multi_s1_5M.dedup.bam --cb CB --umi UB && \
      python3 ${BENCHMARK_DIR}/00c_dedup_longest_read.py \
        ${DATA_DIR}/multi_s2_5M.prefixed.bam ${DEDUP_DIR}/multi_s2_5M.dedup.bam --cb CB --umi UB && \
      python3 ${BENCHMARK_DIR}/00c_dedup_longest_read.py \
        ${DATA_DIR}/multi_s3_5M.prefixed.bam ${DEDUP_DIR}/multi_s3_5M.dedup.bam --cb CB --umi UB
    " 2>&1 | tee ${OUT_BASE}/dedup.log

    INPUT="--bam ${DEDUP_DIR}/multi_s1_5M.dedup.bam ${DEDUP_DIR}/multi_s2_5M.dedup.bam ${DEDUP_DIR}/multi_s3_5M.dedup.bam"
    LABELS="--labels S1 S2 S3"
fi

echo "[$(date)] Dedup complete."

# -------------------------------------------------------
# Step 2: Run IsoQuant on deduplicated BAMs
# -------------------------------------------------------
echo "[$(date)] Step 2: Running IsoQuant..."

/usr/bin/time -v -o ${OUT_BASE}/time_isoquant.txt \
  isoquant.py \
    --reference ${REF_FASTA} \
    --genedb ${REF_GTF} \
    --complete_genedb \
    ${INPUT} \
    --data_type nanopore \
    --barcoded_bam \
    --barcode_tag CB --umi_tag UB \
    ${LABELS} \
    --read_group tag:CB \
    --threads ${THREADS} \
    -o ${OUT_BASE}/output \
  2>&1 | tee ${OUT_BASE}/isoquant.log

# Record output size
du -sh ${OUT_BASE}/ > ${OUT_BASE}/disk_usage.txt

echo "[$(date)] IsoQuant benchmark complete."
