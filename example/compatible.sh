#! /bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --array=0-29

path="path/to/SCOTCH"

python3 ${path}/src/main_preprocessing.py \
--task 'compatible matrix' --platform '10x-ont' \
--target path/to/sample1 path/to/sample2 \
--bam path/to/bam/file1 path/to/bam/file2 \
--job_index ${SLURM_ARRAY_TASK_ID} \
--reference path/to/reference/genes.gtf \
--total_jobs 30
