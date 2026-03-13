#!/bin/bash
# Subsample BAM files for benchmarking via random read sampling.
# Resource settings (CPU, mem, time) are controlled by run_experiment.sh.
#
# Usage:
#   bash 00_prepare_data.sh <input_bam_dir> <output_dir>
#
# Input: directory containing one or more BAM files.
#   If multiple BAMs exist (e.g. per-chromosome), they are merged first.
#
# Creates:
#   Single-sample: reads_5M.bam, reads_10M.bam, reads_15M.bam
#   Multi-sample:  multi_s1_5M.bam, multi_s2_5M.bam, multi_s3_5M.bam

set -euo pipefail

INPUT_DIR=$1
OUTDIR=$2

mkdir -p ${OUTDIR} logs

# Find BAM files in input directory
BAMS=($(ls ${INPUT_DIR}/*.bam 2>/dev/null))
if [ ${#BAMS[@]} -eq 0 ]; then
    echo "ERROR: No BAM files found in ${INPUT_DIR}"
    exit 1
fi
echo "Found ${#BAMS[@]} BAM file(s) in ${INPUT_DIR}"

# --- Merge all BAMs if multiple exist (e.g. per-chromosome split) ---
MERGED_BAM="${OUTDIR}/merged_input.bam"
if [ ${#BAMS[@]} -eq 1 ]; then
    echo "Single BAM found, using directly."
    MERGED_BAM="${BAMS[0]}"
else
    if [ -f "${MERGED_BAM}" ]; then
        echo "Merged BAM already exists: ${MERGED_BAM}"
    else
        echo "Merging ${#BAMS[@]} BAM files into ${MERGED_BAM}..."
        samtools merge --threads 4 "${MERGED_BAM}" "${BAMS[@]}"
        samtools index "${MERGED_BAM}"
        echo "Merge complete."
    fi
fi

TOTAL_READS=$(samtools view -c -F 0x904 "${MERGED_BAM}")
echo "Total primary reads in merged BAM: ${TOTAL_READS}"

# Function: subsample a BAM to a target read count
subsample_bam() {
    local bam=$1
    local target_reads=$2
    local seed=$3
    local outbam=$4

    local total_reads
    total_reads=$(samtools view -c -F 0x904 ${bam})

    if [ ${total_reads} -le ${target_reads} ]; then
        echo "WARNING: BAM has fewer reads (${total_reads}) than target (${target_reads}). Copying as-is."
        cp ${bam} ${outbam}
    else
        local fraction
        fraction=$(python3 -c "print(f'{min(${target_reads}/${total_reads}, 0.999999):.6f}')")
        echo "Subsampling: target=${target_reads}, fraction=${fraction}, seed=${seed}"
        samtools view -b -s ${seed}.${fraction#0.} --threads 4 ${bam} > ${outbam}
    fi
    samtools index ${outbam}

    # Report stats
    local actual_reads ncells
    actual_reads=$(samtools view -c -F 0x904 ${outbam})
    ncells=$(samtools view ${outbam} | awk '{for(i=12;i<=NF;i++) if($i~/^CB:Z:/) print substr($i,6)}' | sort -u | wc -l)
    echo "Result: ${actual_reads} reads, ${ncells} cells -> ${outbam}"
    echo -e "${outbam}\t${actual_reads}\t${ncells}" >> ${OUTDIR}/subsample_stats.tsv
}

echo -e "bam_file\tread_count\tcell_count" > ${OUTDIR}/subsample_stats.tsv

# --- Single-sample subsampling (from merged BAM) ---
echo ""
echo "=== Single-sample subsampling from merged BAM (${TOTAL_READS} reads) ==="

for TARGET in 5000000 10000000 15000000; do
    LABEL=$((TARGET / 1000000))M
    OUTBAM=${OUTDIR}/reads_${LABEL}.bam
    if [ -f "${OUTBAM}" ]; then
        echo "Skipping ${OUTBAM} (already exists)"
        continue
    fi
    echo "--- Subsampling ${LABEL} ---"
    subsample_bam "${MERGED_BAM}" ${TARGET} 42 ${OUTBAM}
done

# --- Multi-sample subsampling (3 samples x 5M reads, different seeds) ---
echo ""
echo "=== Multi-sample subsampling (3 x 5M reads from merged BAM, different seeds) ==="

for i in 1 2 3; do
    OUTBAM=${OUTDIR}/multi_s${i}_5M.bam
    if [ -f "${OUTBAM}" ]; then
        echo "Skipping ${OUTBAM} (already exists)"
        continue
    fi
    echo "--- Multi-sample S${i} ---"
    subsample_bam "${MERGED_BAM}" 5000000 $((42 + i)) ${OUTBAM}
done

echo ""
echo "Done. Stats: ${OUTDIR}/subsample_stats.tsv"
