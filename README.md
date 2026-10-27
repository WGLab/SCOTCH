# SCOTCH

**Isoform-level analysis of long-read RNA sequencing.** SCOTCH maps genome-aligned reads to isoforms, produces expression count matrices, and provides an R package for comparing gene expression and transcript usage between cell populations.

## What can SCOTCH do?

- Map reads to isoforms, including novel isoforms.
- Generate gene- and isoform-level count matrices.
- Test differential gene expression, differential transcript usage (DTU), and isoform switching.
- Process both bulk and single-cell long-read RNA-seq data from Oxford Nanopore and PacBio.

## Quick start

**Prefer to run each step yourself? Follow the [step-by-step tutorial](docs/preprocessing.md), with commands and parameter explanations for all four steps.**

The wrapper below automatically runs all four preprocessing steps.

### Install

With Conda installed:

```bash
git clone https://github.com/WGLab/SCOTCH.git
cd SCOTCH
conda env create --name SCOTCH --file scotch.yml
conda activate SCOTCH
```

### Run

Start with a **genome-aligned, sorted and indexed BAM**. Single-cell BAMs must contain cell barcode and UMI tags from upstream preprocessing. Use a GTF and genome FASTA that match the BAM's genome assembly and chromosome names. `--target` is the output directory.

Copy the configuration template. Set `PLATFORM`, `BULK`, and `RUNNER`, fill in the input/output paths, and check the barcode and computing-resource settings for your run. Optional analysis defaults are at the bottom; every setting has an explanatory comment:

```bash
cp config_template.sh my_run.sh
# Edit my_run.sh: data type, inputs, barcodes, and local/SLURM resources
bash scotch.sh my_run.sh
```

The wrapper runs **annotation → compatible matrix → summary → count matrix**. Set `RUNNER="local"` for local CPU execution or `RUNNER="slurm"` to submit jobs with dependencies and a compatible-matrix job array. Use `bash scotch.sh my_run.sh --dry-run` to check your config and preview commands.

Intermediate files are retained by default. Set `KEEP_INTERMEDIATE=false` to clean reusable intermediates after successful completion; see the [cleanup details](docs/configuration.md#keeping-or-removing-intermediate-files) before disabling retention.

Set `PLATFORM="10x-pacbio"` for PacBio or `PLATFORM="parse-ont"` for Parse + ONT. For bulk data, set `BULK=true`. See the [configuration guide](docs/configuration.md) for CPU and SLURM settings, or the [preprocessing guide](docs/preprocessing.md) for individual steps.

## Outputs

| Output | Location | Use |
| --- | --- | --- |
| Gene and isoform count matrices | `count_matrix/` | Expression analysis and import into Seurat or other tools |
| Updated annotation GTF | `reference/SCOTCH_updated_annotation.gtf` | Known and newly identified isoform structures |
| Read–isoform mapping TSV | `auxiliary/all_read_isoform_exon_mapping.tsv` | Inspect read assignments and exon mappings |

Matrices are saved as `.mtx` files with cell and feature names in accompanying `.pickle` files. Add `--save_csv` for CSV output. For single-cell data, rows are cells and columns are genes or isoforms; transpose the matrix before importing into Seurat.

For Parse runs, matrices and read mappings are under `samples/<sample>/` within the output directory; the shared GTF remains in `reference/`.

## DTU analysis

After preprocessing, use the **SCOTCH R package** to compare gene expression and relative transcript usage between cell populations. See the [R installation and analysis example](docs/statistical_analysis.md); example count matrices are in [`data/`](data).

The Python preprocessing command generates the counts and annotations. Run the R analysis separately with your cell-population assignments.

## Documentation

- [Step-by-step tutorial: inputs, annotation modes, all four steps, and gene subsets](docs/preprocessing.md)
- [Configuration: local CPU execution and SLURM job arrays](docs/configuration.md)
- [Statistical analysis: R installation, example code, and result columns](docs/statistical_analysis.md)
- [SLURM example scripts](example)
- [Report a problem](https://github.com/WGLab/SCOTCH/issues) with your command and the log from the failing step.

## Citation

Xu Z, Qu HQ, Chan J, Mu S, Kao C, Hakonarson H, Wang K. [SCOTCH: isoform-level characterization of gene expression through long-read single-cell RNA sequencing](https://www.nature.com/articles/s41467-026-72665-5). Nat Commun. 2026 May 7. doi: 10.1038/s41467-026-72665-5. PMID: 42098110.
