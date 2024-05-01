#!/bin/bash

source ~/.bashrc
conda activate singlecell

showHelp() {
    echo "Usage: ./parse.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help           Show this help message"
    echo "Mandatory arguments:"
    echo "  -o, --output_dir     Output directory"
    echo "  -f, --fastq      path to fastq folder"
    echo "  -r, --ref      path to reference folder file"
    echo "Optional arguments:"
    echo "  --workers      number of workers"
    echo "  --build      genome build of reference"
}
options=$(getopt -l "help,output_dir:,fastq:,workers:,ref:,build:," \
          -o "ho:f:r:" -a -- "$@")
[ $? -ne 0 ] && exit 1
eval set -- "$options"

while true; do
    case "$1" in
        -h | --help)
            showHelp
            exit 0
            ;;
        -o | --output_dir)
            output_dir="$2"
            shift 2
            ;;
        -f | --fastq)
            fastq="$2"
            shift 2
            ;;
        -r | --ref)
            ref="$2"
            shift 2
            ;;
        --workers)
            workers="$2"
            shift 2
            ;;
        --build)
            build="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Invalid option: $1"
            showHelp
            exit 1
            ;;
    esac
done

SCRIPT="onp_create_reads.py"
SCRIPT_PATH=$(pwd)
MAIN=$output_dir
FQ=$fastq
if [ -z "$workers" ]; then
    WORKER=1
else
    WORKER=$workers
fi

if [ -z "$build" ]; then
    BUILD='GRCh38'
else
    BUILD=$build
fi

echo "processing fastq"
# Iterate over each fastq file in the FQ directory
for file in "$FQ"/*.fastq.gz; do
    # Run the Python script in parallel for each file
    parallel -j "$WORKER" python3 "$SCRIPT" \
    --out_dir "$MAIN" \
    --fastq "{}" \
    --l1dist 6 \
    --l2dist 4 ::: "$file"
done 2>&1 #| tee "$(date "+%Y%m%d_%H%M")_log.txt"

# Concatenate results
cd "$MAIN"
echo "concatenating resutls"
cat 260bps_parsed_R1_*.fastq > 260bps_parsed_R1_cat.fastq
cat 260bps_parsed_R2_*.fastq > 260bps_parsed_R2_cat.fastq

#compress files
echo "compressing results"
pigz -p 8 260bps_parsed_R1_cat.fastq
pigz -p 8 260bps_parsed_R2_cat.fastq

#this step need to install parse pipeline software from https://support.parsebiosciences.com/
echo "make reference"
cd "$ref"
fasta=$(find "$ref" -name "*.fa.gz" -print -quit)
genes=$(find "$ref" -name "*.gtf.gz" -print -quit)
split-pipe --mode mkref --genome_name $build --fasta $fasta --genes $genes --output_dir $ref/$build

echo "create reference files for minimap2"
# Create bed file from GTF
gtf_ref=$(find "$ref" -name "*.gtf" -print -quit)
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id\"\";"; }' $gtf_ref | convert2bed --input=gtf --output=bed > genome_reference.bed


# Add prefix to fasta
fa_ref=$(find "$ref" -name "*.fa" -print -quit)
cat "$fa_ref" | sed "s/^>/>${BUILD}_/g" > genome_reference.pre.dna.fa
pigz -k genome_reference.pre.dna.fa
# Add prefix to bed reference
sed "s/^/${BUILD}_/g" genome_reference.bed > genome_reference.pre.bed
# Create index
FA="genome_reference.pre.dna.fa.gz"
minimap2 -x map-ont -d ${BUILD}.pre.mmi "$FA"


echo "Run split-pipe"
fq1=$(find "$MAIN" -name "*R1_cat.fastq.gz" -print -quit)
split-pipe --fq1 $fq1 --mode pre --one_step --chemistry v1 --nthreads $workers --genome_dir $ref/$build \
--parfile ${SCRIPT_PATH}/parfile.txt --output_dir ${MAIN}/400bps-out


IDX=${BUILD}.pre.mmi
BED=genome_reference.pre.bed
FQ=${MAIN}/400bps-out/process/barcode_head.fastq.gz
minimap2 --MD -a -u f -x splice -t 16 --junc-bed "$BED" "$IDX" "$FQ" > ${MAIN}/400bps-out/process/aligned.sam


# Get reads input
grep "reads_valid_barcode" pre_align_stats.csv \
| cut -d"," -f2 \
| sed 's/^/reads_align_input,/' >> pre_align_stats.csv

# Unique reads
grep -v '^@' aligned.sam \
| cut -f5 \
| grep "60" \
| wc -l \
| sed 's/^/reads_align_unique,/' >> pre_align_stats.csv

# Multimappers
grep -v '^@' aligned.sam \
| cut -f5 \
| grep -v "60" \
| wc -l \
| sed 's/^/reads_align_multimap,/' >> pre_align_stats.csv

# Check to see if the number of SAM entries matches the sum of multimappers and unique reads.
grep -v '^@' aligned.sam | wc -l
cut -d"," -f2 pre_align_stats.csv | tail -2 | paste -sd+ - | bc


samtools view -Sb ${MAIN}/400bps-out/process/aligned.sam | \
samtools sort -o ${MAIN}/400bps-out/process/barcode_headAligned_sorted.bam
samtools index ${MAIN}/400bps-out/process/barcode_headAligned_sorted.bam


split-pipe --mode post --chemistry v1 \
    --nthreads $workers --genome_dir $ref/$build \
    --parfile ${SCRIPT_PATH}/parfile.txt --output_dir ${MAIN}/400bps-out/