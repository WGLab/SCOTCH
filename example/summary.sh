#! /bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

path="path/to/SCOTCH"

python3 ${path}/src/main_preprocessing.py \
--task 'summary' \
--target path/to/sample1 path/to/sample2