#! /bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=parse
#SBATCH --time=7-00:00:00   # HH/MM/SS
#SBATCH --mem=320G   # memory requested, units available: K,M,G,T
#SBATCH --output='/home/xu3/slurm_sc/out/%j_parse.out'
#SBATCH --error='/home/xu3/slurm_sc/err/%j_parse.err'
#SBATCH --array=0
source ~/.bashrc
conda activate scats


code_path='/path/to/onp_create_reads.py' #onp_create_reads.py is provided by Parse
ref='/path/to/reference/genomes'
output_folder='/path/to/output/folder'
input_folder='/path/to/input/fastqgz_file'
build="GRCh38"

cd $code_path
bash parse.sh -o $output_folder -f $input_folder --workers 32 -r $ref --build "GRCh38"


cd $ref
IDX=GRCh38.pre.mmi
BED=genome_reference.pre.bed
FQ=${output_folder}/400bps-out/process/barcode_head.fastq.gz
minimap2 --MD -a -u f -x splice -t 16 --junc-bed "$BED" "$IDX" "$FQ" > ${output_folder}/400bps-out/process/aligned.sam

fq1=$(find "$output_folder" -name "*R1_cat.fastq.gz" -print -quit)
split-pipe --fq1 $fq1 --mode all --chemistry v1 \
    --nthreads 32 --genome_dir $ref/$build \
    --parfile ${code_path}/parfile.txt --output_dir ${output_folder}/400bps-out/




