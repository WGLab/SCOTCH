#! /bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G

path="path/to/SCOTCH"

python3 ${path}/src/main_preprocessing.py \
--task 'count matrix' --platform '10x-ont' \
--target path/to/sample1 path/to/sample2 \
--workers 8 --save_mtx --save_csv --group_novel
