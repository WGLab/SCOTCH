#! /bin/bash -l
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G

path="path/to/SCOTCH"

python3 ${path}/src/main_preprocessing.py \
--task 'count matrix' --platform '10x' \
--target path/to/sample1 path/to/sample2 \
--novel_read_n 20 --workers 100