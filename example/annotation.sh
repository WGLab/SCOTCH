#! /bin/bash -l
#SBATCH --cpus-per-task=30
#SBATCH --mem=300G

path="path/to/SCOTCH"

python3 ${path}/src/main_preprocessing.py \
--task annotation --platform '10x' \
--target path/to/sample1 path/to/sample2 \
--bam path/to/bam/file1 path/to/bam/file2 \
--reference ${path}/data/geneStructureInformation.pkl \
--update_gtf --workers 30








