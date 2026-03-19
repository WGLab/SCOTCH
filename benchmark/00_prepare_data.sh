#!/bin/bash
# Subsample BAM files for benchmarking via random read sampling.
# Resource settings (CPU, mem, time) are controlled by run_experiment.sh.
#
# Usage:
#   bash 00_prepare_data.sh <input_bam_dir> <output_dir>
#
# Input: directory containing one or more BAM files.
#   Multiple input BAMs are treated as shards of the same sample
#   (for example per-chromosome BAMs). Each shard is counted separately,
#   subsampled with the same fraction, then concatenated with samtools cat.
#
# Creates:
#   Single-sample: reads_1M.bam, reads_5M.bam, reads_15M.bam, reads_50M.bam
#   Multi-sample:  multi_s1_5M.bam, multi_s2_5M.bam, multi_s3_5M.bam

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: bash 00_prepare_data.sh <input_bam_dir> <output_dir>"
    exit 1
fi

INPUT_DIR=$1
OUTDIR=$2
THREADS=${THREADS:-4}

mkdir -p "${OUTDIR}" logs

mapfile -t BAMS < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.bam" | sort -V)
if [ ${#BAMS[@]} -eq 0 ]; then
    echo "ERROR: No BAM files found in ${INPUT_DIR}"
    exit 1
fi
echo "Found ${#BAMS[@]} BAM file(s) in ${INPUT_DIR}"

declare -a BAM_PRIMARY_READS=()
TOTAL_PRIMARY_READS=0

echo ""
echo "=== Counting primary reads per BAM shard ==="
for bam in "${BAMS[@]}"; do
    count=$(samtools view -c -F 0x904 "${bam}")
    BAM_PRIMARY_READS+=("${count}")
    TOTAL_PRIMARY_READS=$((TOTAL_PRIMARY_READS + count))
    echo "  ${bam}: ${count}"
done
echo "Total primary reads across all BAMs: ${TOTAL_PRIMARY_READS}"

count_cells() {
    local bam=$1
    samtools view -F 0x904 "${bam}" | \
        awk '{for(i=12;i<=NF;i++) if($i~/^CB:Z:/) print substr($i,6)}' | \
        sort -u | wc -l
}

subsample_concat_bams() {
    local target_reads=$1
    local seed=$2
    local outbam=$3
    local fraction
    local tmpdir=""
    local -a shard_outputs=()

    cleanup_tmpdir() {
        if [ -n "${tmpdir}" ] && [ -d "${tmpdir}" ]; then
            rm -rf "${tmpdir}"
        fi
    }

    if [ "${TOTAL_PRIMARY_READS}" -le "${target_reads}" ]; then
        echo "WARNING: Total reads (${TOTAL_PRIMARY_READS}) <= target (${target_reads}). Concatenating full BAM set."
        if [ ${#BAMS[@]} -eq 1 ]; then
            cp "${BAMS[0]}" "${outbam}"
        else
            samtools cat --threads "${THREADS}" -o "${outbam}" "${BAMS[@]}"
        fi
    else
        fraction=$(python3 - "${target_reads}" "${TOTAL_PRIMARY_READS}" <<'PY'
import sys
target = int(sys.argv[1])
total = int(sys.argv[2])
fraction = min(target / total, 0.999999)
print(f"{fraction:.6f}")
PY
)
        echo "Subsampling: target=${target_reads}, fraction=${fraction}, seed=${seed}"

        tmpdir=$(mktemp -d "${OUTDIR}/tmp.subsample.XXXXXX")
        trap cleanup_tmpdir RETURN

        for i in "${!BAMS[@]}"; do
            bam=${BAMS[$i]}
            bam_reads=${BAM_PRIMARY_READS[$i]}
            shard_out="${tmpdir}/$(basename "${bam%.bam}").subsampled.bam"

            if [ "${bam_reads}" -eq 0 ]; then
                echo "  Skipping empty shard: ${bam}"
                continue
            fi

            samtools view -b --threads "${THREADS}" -s "${seed}.${fraction#0.}" \
                -o "${shard_out}" "${bam}"
            shard_outputs+=("${shard_out}")
        done

        if [ ${#shard_outputs[@]} -eq 0 ]; then
            echo "ERROR: No shard outputs were produced for ${outbam}"
            exit 1
        elif [ ${#shard_outputs[@]} -eq 1 ]; then
            mv "${shard_outputs[0]}" "${outbam}"
        else
            samtools cat --threads "${THREADS}" -o "${outbam}" "${shard_outputs[@]}"
        fi

        cleanup_tmpdir
        trap - RETURN
    fi

    samtools index "${outbam}"

    actual_reads=$(samtools view -c -F 0x904 "${outbam}")
    ncells=$(count_cells "${outbam}")
    echo "Result: ${actual_reads} reads, ${ncells} cells -> ${outbam}"
    echo -e "${outbam}\t${actual_reads}\t${ncells}" >> "${OUTDIR}/subsample_stats.tsv"
}

echo -e "bam_file\tread_count\tcell_count" > "${OUTDIR}/subsample_stats.tsv"

echo ""
echo "=== Single-sample subsampling from ${TOTAL_PRIMARY_READS} total reads ==="

for TARGET in 1000000 5000000 15000000 50000000; do
    LABEL=$((TARGET / 1000000))M
    OUTBAM="${OUTDIR}/reads_${LABEL}.bam"
    if [ -f "${OUTBAM}" ]; then
        echo "Skipping ${OUTBAM} (already exists)"
        continue
    fi
    echo "--- Subsampling ${LABEL} ---"
    subsample_concat_bams "${TARGET}" 42 "${OUTBAM}"
done

echo ""
echo "=== Multi-sample subsampling (3 x 5M reads, different seeds) ==="

for i in 1 2 3; do
    OUTBAM="${OUTDIR}/multi_s${i}_5M.bam"
    if [ -f "${OUTBAM}" ]; then
        echo "Skipping ${OUTBAM} (already exists)"
        continue
    fi
    echo "--- Multi-sample S${i} ---"
    subsample_concat_bams 5000000 $((42 + i)) "${OUTBAM}"
done

echo ""
echo "Done. Stats: ${OUTDIR}/subsample_stats.tsv"
