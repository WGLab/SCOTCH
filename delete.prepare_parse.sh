#!/bin/bash

source ~/.bashrc
conda activate singlecell

##TODO: add parse options
showHelp() {
    echo "Usage: ./prepare.sh [OPTIONS]"
    echo "Options:"
    echo "  -h, --help           Show this help message"
    echo "Mandatory arguments:"
    echo "  -t, --task TASK      Task to run (nanopore, annotation, matrix, or count)"
    echo "  -o, --output_dir     Output directory"
    echo "annotation task-specific options:"
    echo "  --reference      Gene annotation file in gtf format"
    echo "  --bam          bam file path, need to be sorted, filtered, and indexed before running"
    echo "Parse task-specific options:"
    echo "  --fastq          fastq file path"
    echo "  --kit_name           multiome, 3prime, 5prime"
    echo "  --kit_version        kit version"
    echo "  --reference_nanopore path to reference directory"
    echo "matrix task-specific options:"
    echo "  --bam          bam file path, need to be sorted, filtered, and indexed before running"
    echo "  --match           threshold to decide reads map to exons"
    echo "  --workers       number of processes per job"
    echo "  --jobs    number of jobs"
    echo "  --job_index   index of jobs"
    echo "  --cover_existing   whether to overwirte existing compatible matrix"
    echo "count task-specific options:"
    echo "  --novel_read_n          read number to support a novel isoform"
    echo "  --sample_name          sample name"
    echo "  --key_words             key words for parse sublibraries"
}


options=$(getopt -l "help,task:,output_dir:,reference:,fastq:,kit_name:,kit_version:,reference_nanopore:,bam:,match:,workers:,job_index:,jobs:,cover_existing,novel_read_n:,sample_name:,key_words:," \
          -o "ht:o:" -a -- "$@")
[ $? -ne 0 ] && exit 1
eval set -- "$options"

while true
do
    case $1 in
        -h|--help)
            showHelp
            exit 0
            ;;
        -o|--output_dir)
            output_dir=$2
            echo "set output directory as: "$output_dir
            shift
            ;;
        -t|--task)
            task=$2
            echo "running task of: "$task
            shift
            ;;
        --reference)
            reference=$2
            echo "use reference annotation file at: "$reference
            shift
            ;;
        --reference_nanopore)
            reference_nanopore=$2
            echo "use reference annotation file at: "$reference_nanopore
            shift
            ;;
        --fastq)
            fastq=$2
            shift
            ;;
        --kit_name)
            kit_name=$2
            shift
            ;;
        --kit_version)
            kit_version=$2
            shift
            ;;
        --bam)
            bam=$2
            shift
            ;;
        --match)
            match=$2
            shift
            ;;
        --workers)
            workers=$2
            shift
            ;;
        --job_index)
            job_index=$2
            shift
            ;;
        --jobs)
            jobs=$2
            shift
            ;;
        --cover_existing)
            cover_existing=$2
            shift
            ;;
        --novel_read_n)
            novel_read_n=$2
            shift
            ;;
        --key_words)
            key_words=$2
            shift
            ;;
        --sample_name)
            sample_name=$2
            shift
            ;;
        --)
            shift
            break;;
    esac
    shift
done

if [ -z "$task" ]; then
    echo "Task is missing." >&2
    exit 1
fi

if [ -z "$output_dir" ]; then
    echo "output_dir is missing." >&2
    exit 1
fi

if [ "$task" = "annotation" ]; then
    # Check if task specific options are provided
    if [ -z "$reference" ]||[ -z "$bam" ]; then
        echo "annotation-specific option is missing" >&2
        exit 1
    fi
fi

if [ "$task" = "nanopore" ]; then
    if [ -z "$reference_nanopore" ]||[ -z "$fastq" ]||[ -z "$kit_name" ]||[ -z "$kit_version" ]; then
        echo "nanopore-specific options are missing" >&2
        exit 1
    fi
fi

if [ "$task" = "matrix" ]; then
    if [ -z "$bam" ]||[ -z "$match" ]||[ -z "$job_index" ]||[ -z "$jobs" ]; then
        echo "matrix-specific options are missing" >&2
        exit 1
    fi
fi

if [ "$task" = "count" ]; then
    if [ -z "$novel_read_n" ]; then
        echo "count-specific options are missing" >&2
        exit 1
    fi
fi



#----------------task annotation: prepare annotation file-----------------#
if [ "$task" = "annotation" ]; then
    #Get the directory path of the shell script
    current_path="$(dirname "$(readlink -f "$0")")"
    echo "preparing annotation file"
    mkdir -p $output_dir/reference
    #extract annotation information
    python3 $current_path/src/main_preprocessing_parse.py --ref $reference --bam $bam --workers $workers --target $output_dir --task annotation
    files=$(ls $output_dir/reference)
    echo "task done, $files are in $output_dir/reference"
fi


#----------------task matrix: prepare compatible matrix -----------------#
if [ "$task" = "matrix" ]; then
    current_path="$(dirname "$(readlink -f "$0")")"
    mkdir -p $output_dir/compatible_matrix
    if [ -z "$cover_existing" ] || [ "$cover_existing" -eq 1 ]; then
  # If cover_existing is null or equals 1, run with --cover_existing
    python3 "$current_path/src/main_preprocessing_parse.py" --target "$output_dir" --task matrix \
    --bam "$bam" --match "$match" --geneinfo "$output_dir/reference/geneStructureInformation.pkl" \
    --job_index "$job_index" --total_jobs "$jobs" --cover_existing
    else
      # Otherwise, run with --cover_existing_false
      python3 "$current_path/src/main_preprocessing_parse.py" --target "$output_dir" --task matrix \
      --bam "$bam" --match "$match" --geneinfo "$output_dir/reference/geneStructureInformation.pkl" \
      --job_index "$job_index" --total_jobs "$jobs" --cover_existing_false
    fi
fi

#----------------task count: prepare count matrix -----------------#
if [ "$task" = "count" ]; then
    current_path="$(dirname "$(readlink -f "$0")")"
    mkdir -p $output_dir/count_matrix

    if [ -z "$workers" ]; then
      workers=8
    fi
    python3 $current_path/src/main_preprocessing_parse.py --target $output_dir --task count --novel_read_n $novel_read_n --workers $workers --sample_name $sample_name --key_words $key_words
fi


#----------------nanopore workflow for alignment and tagging bam files --------------#
if [ "$task" = "nanopore" ]; then
  echo "running nanopore pipeline for tagged bam file"
  ~/nextflow run epi2me-labs/wf-singlecell \
  -w $output_dir -profile singularity \
  --matrix_min_genes 1 --fastq $fastq \
  --ref_genome_dir $reference_nanopore --out_dir $output_dir \
  --kit_name $kit_name --kit_version $kit_version
fi
#multiome, v1; 3prime, v3











