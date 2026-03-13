# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment Setup

```bash
conda env create --name SCOTCH --file scotch.yml
conda activate SCOTCH
```

No `setup.py` or `pyproject.toml` — conda is the sole package manager (Python 3.10.12).

## Running the Pipeline

The single entry point is `src/main_preprocessing.py`. Tasks must run in order:

```bash
# Step 1: Annotation — parse GTF + BAM, build gene structure pickle
python3 src/main_preprocessing.py \
  --task annotation --platform 10x-ont \
  --target <output_dir> --bam <file.bam> \
  --reference <genes.gtf> --reference_pkl data/geneStructureInformation.pkl \
  --update_gtf --workers 30

# Step 2: Compatible matrix — map reads to isoforms (parallelizable via job arrays)
python3 src/main_preprocessing.py \
  --task 'compatible matrix' \
  --target <output_dir> --bam <file.bam> \
  --reference <genes.gtf> \
  --job_index 0 --total_jobs 30    # omit for single-job run

# Step 3: Summary — merge per-job outputs
python3 src/main_preprocessing.py --task summary --target <output_dir>

# Step 4: Count matrix — generate AnnData/MTX/CSV outputs
python3 src/main_preprocessing.py \
  --task 'count matrix' --platform 10x-ont \
  --target <output_dir> --workers 8 --group_novel
```

**Gene subset mode** (for testing or targeted analysis):
```bash
--gene_subset BRCA1 TP53       # pass gene names directly
--gene_subset /path/genes.txt  # or a file with one gene per line
```

There are no automated tests. Manual validation uses the SLURM example scripts in `example/`.

## Architecture

### Python Preprocessing Pipeline (`src/`)

Four main classes, each mapping to a pipeline step:

| Class | File | Step |
|---|---|---|
| `Annotator` | `annotation.py` | Step 1 — build gene structure from GTF + BAM |
| `ReadMapper` | `compatible.py` | Step 2 — map reads to isoforms, detect novel isoforms |
| `summarise_annotation` / `summarise_auxillary` | `compatible.py` | Step 3 — merge outputs |
| `CountMatrix` | `count_matrix.py` | Step 4 — generate count matrices |

`main_preprocessing.py` is a thin CLI orchestrator (argparse) that instantiates and calls these classes. `preprocessing.py` contains shared utilities (`load_pickle`, `merge_exons`, etc.) used across all steps. `reference.py` handles GTF parsing using bedtools.

### Core Data Structures

**`metageneStructureInformation`** — central annotation dict:
```python
{
  metagene_name: [
    [geneInfo, exonInfo, isoformInfo],  # one entry per gene in the metagene
    ...
  ]
}
# geneInfo keys: geneName, geneID, geneChr, geneStrand, geneStart, geneEnd
# isoformInfo: { isoform_name: [exon_indices] }
```

Compatible matrix CSVs are named `{geneName}_{geneID}.csv`.

### R Statistical Pipeline (`R/`)

`DTU.R` implements Differential Transcript Usage via a Dirichlet-multinomial model:
- `scotch_gene()` — gene-level differential expression
- `scotch_transcript()` — transcript-level DTU with LRT (`LRT_test()`)

Example count matrices for the R pipeline are in `data/`.

### Parallelization Strategy

- **Compatible matrix step**: designed for SLURM job arrays — `--job_index` / `--total_jobs` split genes across array tasks (see `example/compatible.sh`)
- **Count matrix step**: uses joblib internally via `--workers`

### Platform Differences

Three platforms are supported (`--platform`): `10x-ont`, `10x-pacbio`, `parse-ont`. Parse Biosciences requires upstream FASTQ→BAM conversion; see `example/parse_preprocessing/`.
