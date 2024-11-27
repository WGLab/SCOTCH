#! /bin/bash -l
#SBATCH --cpus-per-task=30
#SBATCH --mem=300G

path="path/to/SCOTCH"

python3 ${path}/src/main_preprocessing.py \
--task annotation --platform '10x-ont' \
--target path/to/sample1 path/to/sample2 \
--bam path/to/bam/file1 path/to/bam/file2 \
--reference_pkl ${path}/data/geneStructureInformation.pkl \
--reference path/to/reference/genes.gtf \
--update_gtf --workers 30








