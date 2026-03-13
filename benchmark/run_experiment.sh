#!/bin/bash
# =============================================================================
# SCOTCH vs IsoQuant Benchmark — Master Experiment Runner
# =============================================================================
#
# Run each phase independently. Phases 1-4 can run in any order, whenever
# you have resources available. Timing is measured per-job with /usr/bin/time
# so running phases separately does NOT affect the final metrics.
#
# Usage:
#   bash run_experiment.sh <phase>
#
# Phases:
#   0  — Data preparation (subsampling + barcode prefixing)
#   1  — SCOTCH single-sample (5M, 10M, 15M)
#   2  — SCOTCH multi-sample (3 x 5M)
#   3  — IsoQuant single-sample (5M, 10M, 15M)
#   4  — IsoQuant multi-sample (3 x 5M)
#   5  — Collect metrics + plot (checks all Phase 1-4 outputs exist first)
#   all — Run phases 0-4 with SLURM dependency chains (old behavior)
#
# Examples:
#   bash run_experiment.sh 0        # Prep data first
#   bash run_experiment.sh 1        # Then run SCOTCH single whenever ready
#   bash run_experiment.sh 3        # Run IsoQuant single in parallel with SCOTCH
#   bash run_experiment.sh 5        # Collect results after everything finishes
#   bash run_experiment.sh all      # Submit everything with dependency chains
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
SCOTCH_DIR="/home/xu3/SCOTCH"    # Root of SCOTCH codes
BENCHMARK_DIR="/mnt/isilon/wang_lab/karen/scotch/benchmark/computation"    # This benchmark/ directory
REF_FASTA="/scr1/users/xu3/singlecell/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa"    # Reference genome FASTA
REF_GTF="/scr1/users/xu3/singlecell/ref/refdata-gex-GRCh38-2020-A/genes.gtf"           # Gene annotation GTF
REF_PKL="/home/xu3/SCOTCH/data/geneStructureInformation.pkl"    # SCOTCH annotation pickle
INPUT_BAM_DIR="/mnt/isilon/wang_lab/xinya/projects/single_cell_pipeline/CAG_SingleCell/sample7-R10-allpass-v4/wf-single-cell-v1-output-sample7R10-allpass-ed1/reseq/bams"                # Directory with input BAM(s)
DATA_DIR="/mnt/isilon/wang_lab/karen/scotch/benchmark/computation/data"                    # Where subsampled BAMs go
RESULTS_BASE="/mnt/isilon/wang_lab/karen/scotch/benchmark/computation/results"             # Where results go

# --- Resource settings: SCOTCH ---
SCOTCH_ANNOT_CPUS=10                                # CPUs for annotation step
SCOTCH_ANNOT_MEM=300G                               # Memory for annotation step
SCOTCH_ANNOT_TIME=1-00:00:00                        # Time limit for annotation

SCOTCH_COMPAT_NJOBS=10                              # Number of array jobs for compatible matrix
SCOTCH_COMPAT_MEM=100G                              # Memory per array task
SCOTCH_COMPAT_TIME=1-00:00:00                       # Time limit per array task

SCOTCH_SUMMARY_MEM=80G                              # Memory for summary step
SCOTCH_SUMMARY_TIME=04:00:00                        # Time limit for summary

SCOTCH_COUNT_CPUS=1                                 # CPUs for count matrix step
SCOTCH_COUNT_MEM=200G                               # Memory for count matrix
SCOTCH_COUNT_TIME=1-00:00:00                        # Time limit for count matrix

# --- Resource settings: IsoQuant ---
ISOQUANT_CPUS=10                                    # Threads for IsoQuant
ISOQUANT_MEM=300G                                   # Memory for IsoQuant
ISOQUANT_TIME=5-00:00:00                            # Time limit for IsoQuant

# --- Resource settings: Data preparation ---
PREP_CPUS=4                                         # CPUs for subsampling
PREP_MEM=32G                                        # Memory for subsampling
PREP_TIME=06:00:00                                  # Time limit for subsampling

# =============================================================================
# END USER CONFIG — Do not edit below unless you know what you're doing
# =============================================================================

# Parse phase argument
PHASE="${1:-}"
if [ -z "${PHASE}" ]; then
    echo "Usage: bash run_experiment.sh <phase>"
    echo ""
    echo "Phases:"
    echo "  0    Data preparation (subsampling + barcode prefixing)"
    echo "  1    SCOTCH single-sample (5M, 10M, 15M)"
    echo "  2    SCOTCH multi-sample (3 x 5M)"
    echo "  3    IsoQuant single-sample (5M, 10M, 15M)"
    echo "  4    IsoQuant multi-sample (3 x 5M)"
    echo "  5    Collect metrics + plot (checks outputs exist first)"
    echo "  all  Run phases 0-4 with SLURM dependency chains"
    exit 1
fi

# Build common sbatch flags
SBATCH_COMMON="--mail-user=${EMAIL} --mail-type=FAIL,END"

# Conda activation prefix — prepended to every --wrap command
# Adjust the conda.sh path if your conda is installed elsewhere
CONDA_INIT="source \$(conda info --base)/etc/profile.d/conda.sh && conda activate singlecell &&"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p logs "${DATA_DIR}" "${RESULTS_BASE}"

echo "============================================================"
echo "SCOTCH vs IsoQuant Benchmark — Phase: ${PHASE}"
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
      --wrap="${CONDA_INIT} set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_annotation.txt \
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
      --wrap="${CONDA_INIT} set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_compatible_\${SLURM_ARRAY_TASK_ID}.txt \
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
      --wrap="${CONDA_INIT} set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_summary.txt \
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
      --wrap="${CONDA_INIT} set -o pipefail; /usr/bin/time -v -o ${OUT_BASE}/time_count.txt \
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
      --wrap="${CONDA_INIT} set -euo pipefail; \
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

# -------------------------------------------------------
# Helper: check if all required time files exist for Phase 5
# Returns 0 if all exist, 1 if any missing
# -------------------------------------------------------
check_phase5_ready() {
    local MISSING=()

    # SCOTCH single-sample: 5M, 10M, 15M
    for LABEL in 5M 10M 15M; do
        local BASE="${RESULTS_BASE}/scotch/single/reads_${LABEL}"
        for STEP in annotation summary count; do
            [ -f "${BASE}/time_${STEP}.txt" ] || MISSING+=("scotch/single/${LABEL}/time_${STEP}.txt")
        done
        # Check at least one compatible array task
        ls "${BASE}"/time_compatible_*.txt &>/dev/null || MISSING+=("scotch/single/${LABEL}/time_compatible_*.txt")
    done

    # SCOTCH multi-sample: 15M
    local BASE="${RESULTS_BASE}/scotch/multi/reads_15M"
    for STEP in annotation summary count; do
        [ -f "${BASE}/time_${STEP}.txt" ] || MISSING+=("scotch/multi/15M/time_${STEP}.txt")
    done
    ls "${BASE}"/time_compatible_*.txt &>/dev/null || MISSING+=("scotch/multi/15M/time_compatible_*.txt")

    # IsoQuant single-sample: 5M, 10M, 15M
    for LABEL in 5M 10M 15M; do
        local BASE="${RESULTS_BASE}/isoquant/single/reads_${LABEL}"
        for STEP in dedup isoquant; do
            [ -f "${BASE}/time_${STEP}.txt" ] || MISSING+=("isoquant/single/${LABEL}/time_${STEP}.txt")
        done
    done

    # IsoQuant multi-sample: 15M
    local BASE="${RESULTS_BASE}/isoquant/multi/reads_15M"
    for STEP in dedup isoquant; do
        [ -f "${BASE}/time_${STEP}.txt" ] || MISSING+=("isoquant/multi/15M/time_${STEP}.txt")
    done

    if [ ${#MISSING[@]} -gt 0 ]; then
        echo ""
        echo "WARNING: The following time files are missing (${#MISSING[@]} files):"
        for m in "${MISSING[@]}"; do
            echo "  - ${m}"
        done
        echo ""
        echo "Run the corresponding phases first before collecting metrics."
        return 1
    fi
    return 0
}

# =============================================================================
# PHASE DISPATCH
# =============================================================================

run_phase_0() {
    echo "--- Phase 0: Data Preparation ---"

    # Step 0a: Subsample BAMs
    echo "  Submitting: subsample BAMs..."
    PREP_JOB=$(sbatch --parsable ${SBATCH_COMMON} \
      --job-name=bench_prep \
      --cpus-per-task=${PREP_CPUS} --mem=${PREP_MEM} --time=${PREP_TIME} \
      --output=logs/prep_%j.out --error=logs/prep_%j.err \
      --wrap="${CONDA_INIT} bash ${SCRIPT_DIR}/00_prepare_data.sh ${INPUT_BAM_DIR} ${DATA_DIR}")
    echo "    Subsample: ${PREP_JOB}"

    # Step 0b: Prefix barcodes for multi-sample
    echo "  Submitting: prefix barcodes..."
    PREFIX_JOB=$(sbatch --parsable ${SBATCH_COMMON} \
      --dependency=afterok:${PREP_JOB} \
      --job-name=bench_prefix \
      --cpus-per-task=1 --mem=16G --time=04:00:00 \
      --output=logs/prefix_%j.out --error=logs/prefix_%j.err \
      --wrap="${CONDA_INIT} python3 ${SCRIPT_DIR}/00b_prefix_barcodes.py ${DATA_DIR}")
    echo "    Prefix barcodes: ${PREFIX_JOB}"
}

run_phase_1() {
    echo "--- Phase 1: SCOTCH Single-Sample ---"
    for LABEL in 5M 10M 15M; do
        echo "  SCOTCH single ${LABEL}:"
        submit_scotch "${LABEL}" single
    done
}

run_phase_2() {
    echo "--- Phase 2: SCOTCH Multi-Sample ---"
    echo "  SCOTCH multi 15M:"
    submit_scotch 15M multi
}

run_phase_3() {
    echo "--- Phase 3: IsoQuant Single-Sample ---"
    for LABEL in 5M 10M 15M; do
        echo "  IsoQuant single ${LABEL}:"
        submit_isoquant "${LABEL}" single
    done
}

run_phase_4() {
    echo "--- Phase 4: IsoQuant Multi-Sample ---"
    echo "  IsoQuant multi 15M:"
    submit_isoquant 15M multi
}

run_phase_5() {
    echo "--- Phase 5: Collect Metrics ---"
    echo "  Checking if all experiment outputs exist..."

    if ! check_phase5_ready; then
        echo ""
        read -p "Continue anyway with partial results? [y/N] " -n 1 -r
        echo ""
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Aborted."
            exit 1
        fi
    else
        echo "  All time files found."
    fi

    echo "  Collecting metrics..."
    python3 "${SCRIPT_DIR}/04_collect_metrics.py" "${RESULTS_BASE}" -o "${RESULTS_BASE}/benchmark_metrics.tsv"
    python3 "${SCRIPT_DIR}/05_plot_results.py" "${RESULTS_BASE}/benchmark_metrics.tsv"
    echo "  Done. Results in ${RESULTS_BASE}/benchmark_metrics.tsv"
}

run_all() {
    echo "--- Running all phases with SLURM dependency chains ---"
    echo ""

    # Phase 0
    run_phase_0
    local LAST_PREP_JOB="${PREFIX_JOB}"

    # Phases 1-4: all depend on Phase 0, run in parallel
    ALL_FINAL_JOBS=()

    echo ""
    echo "--- Phase 1: SCOTCH Single-Sample ---"
    for LABEL in 5M 10M 15M; do
        echo "  SCOTCH single ${LABEL}:"
        submit_scotch "${LABEL}" single "${LAST_PREP_JOB}"
        ALL_FINAL_JOBS+=("${SCOTCH_LAST_JOB}")
    done

    echo ""
    echo "--- Phase 2: SCOTCH Multi-Sample ---"
    echo "  SCOTCH multi 15M:"
    submit_scotch 15M multi "${LAST_PREP_JOB}"
    ALL_FINAL_JOBS+=("${SCOTCH_LAST_JOB}")

    echo ""
    echo "--- Phase 3: IsoQuant Single-Sample ---"
    for LABEL in 5M 10M 15M; do
        echo "  IsoQuant single ${LABEL}:"
        submit_isoquant "${LABEL}" single "${LAST_PREP_JOB}"
        ALL_FINAL_JOBS+=("${ISOQUANT_LAST_JOB}")
    done

    echo ""
    echo "--- Phase 4: IsoQuant Multi-Sample ---"
    echo "  IsoQuant multi 15M:"
    submit_isoquant 15M multi "${LAST_PREP_JOB}"
    ALL_FINAL_JOBS+=("${ISOQUANT_LAST_JOB}")

    # Phase 5: after all experiments
    echo ""
    echo "--- Phase 5: Metrics Collection ---"
    DEP_STR=$(IFS=:; echo "${ALL_FINAL_JOBS[*]}")
    COLLECT_JOB=$(sbatch --parsable ${SBATCH_COMMON} \
      --dependency=afterok:${DEP_STR} \
      --job-name=bench_collect \
      --cpus-per-task=1 --mem=8G --time=00:30:00 \
      --output=logs/collect_%j.out --error=logs/collect_%j.err \
      --wrap="${CONDA_INIT} python3 ${SCRIPT_DIR}/04_collect_metrics.py ${RESULTS_BASE} -o ${RESULTS_BASE}/benchmark_metrics.tsv && \
        python3 ${SCRIPT_DIR}/05_plot_results.py ${RESULTS_BASE}/benchmark_metrics.tsv && \
        echo 'Benchmark complete. Results in ${RESULTS_BASE}/benchmark_metrics.tsv'")
    echo "  Collect metrics: ${COLLECT_JOB}"
}

# Dispatch
case "${PHASE}" in
    0) run_phase_0 ;;
    1) run_phase_1 ;;
    2) run_phase_2 ;;
    3) run_phase_3 ;;
    4) run_phase_4 ;;
    5) run_phase_5 ;;
    all) run_all ;;
    *)
        echo "ERROR: Unknown phase '${PHASE}'"
        echo "Valid phases: 0, 1, 2, 3, 4, 5, all"
        exit 1
        ;;
esac

echo ""
echo "============================================================"
echo "Phase ${PHASE} submitted."
echo "Monitor: squeue -u \$(whoami)"
echo "Email notifications on FAIL/END -> ${EMAIL}"
echo "============================================================"
