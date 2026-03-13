#!/bin/bash
#SBATCH --job-name=bench_prep
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/prep_%j.out
#SBATCH --error=logs/prep_%j.err

# Subsample BAM files for benchmarking via random read sampling.
#
# Usage:
#   sbatch 00_prepare_data.sh <input_bam_dir> <output_dir>
#
# Input: directory containing one or more BAM files (uses the first found).
# For multi-sample: uses the first 3 BAMs found in the directory.
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

# Function: subsample a BAM to a target read count
subsample_bam() {
    local bam=$1
    local target_reads=$2
    local seed=$3
    local outbam=$4

    echo "Counting reads in ${bam}..."
    local total_reads
    total_reads=$(samtools view -c -F 0x904 ${bam})
    echo "Total primary reads: ${total_reads}"

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

# --- Single-sample subsampling (from first BAM) ---
PRIMARY_BAM=${BAMS[0]}
echo ""
echo "=== Single-sample subsampling from ${PRIMARY_BAM} ==="

for TARGET in 5000000 10000000 15000000; do
    LABEL=$((TARGET / 1000000))M
    OUTBAM=${OUTDIR}/reads_${LABEL}.bam
    if [ -f "${OUTBAM}" ]; then
        echo "Skipping ${OUTBAM} (already exists)"
        continue
    fi
    subsample_bam ${PRIMARY_BAM} ${TARGET} 42 ${OUTBAM}
done

# --- Multi-sample subsampling (3 samples x 5M reads) ---
echo ""
echo "=== Multi-sample subsampling (3 x 5M reads) ==="

if [ ${#BAMS[@]} -ge 3 ]; then
    # Use first 3 BAMs as separate samples
    SAMPLE_BAMS=("${BAMS[0]}" "${BAMS[1]}" "${BAMS[2]}")
else
    # Use the same BAM with different seeds
    echo "NOTE: Only ${#BAMS[@]} BAM(s) found. Using same BAM with different seeds for multi-sample."
    SAMPLE_BAMS=("${BAMS[0]}" "${BAMS[0]}" "${BAMS[0]}")
fi

for i in 1 2 3; do
    idx=$((i - 1))
    OUTBAM=${OUTDIR}/multi_s${i}_5M.bam
    if [ -f "${OUTBAM}" ]; then
        echo "Skipping ${OUTBAM} (already exists)"
        continue
    fi
    subsample_bam ${SAMPLE_BAMS[$idx]} 5000000 $((42 + i)) ${OUTBAM}
done

echo ""
echo "Done. Stats: ${OUTDIR}/subsample_stats.tsv"
