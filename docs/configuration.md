# Running from a configuration file

[Back to README](../README.md)

Activate the SCOTCH Conda environment, copy [`config_template.sh`](../config_template.sh), and edit your copy:

```bash
cp config_template.sh my_run.sh
bash scotch.sh my_run.sh --dry-run
bash scotch.sh my_run.sh
```

The template separates required run choices and paths, conditional barcode/SLURM settings, computing resources to review, and optional analysis defaults at the bottom. Having a preset value does not mean a setting is appropriate for your data. Every setting has an inline comment. Configs are sourced as Bash code, so use your own trusted config. Relative paths resolve against the config file's directory. Commands can be launched from any working directory using the path to `scotch.sh`.

## Required choices

Set `PLATFORM` (`10x-ont`, `10x-pacbio`, or `parse-ont`), `BULK` (`true` or `false`), and `RUNNER` (`local` or `slurm`) explicitly. These fields are blank in the template, and the wrapper rejects a config that leaves them blank.

For single-cell 10x-ont or 10x-pacbio data, explicitly set `BARCODE_CELL` and `BARCODE_UMI` to the tags used by your upstream workflow. Common examples are `CB`/`UB` for ONT and `CB`/`XM` or `XC`/`XM` for PacBio. The wrapper does not guess or inspect these tags for you. All samples in one run use the same tag names.

For bulk data, barcode fields can be empty and are ignored. Parse mode also ignores these fields: the current Parse parser takes the cell barcode and UMI from read names and the sample identifier from the `pS` tag. See [Parse preprocessing](../example/parse_preprocessing/README.md) for the expected input preparation.

## Required paths

| Setting | Meaning |
| --- | --- |
| `BAM_PATHS` | Bash array of genome-aligned BAM files, one entry per sample. BAMs must be sorted and indexed. |
| `TARGET_DIRS` | Bash array of output directories, one per BAM input in matching order. |
| `REFERENCE_GTF` | GTF matching the genome assembly and chromosome names in the BAM. |
| `REFERENCE_FASTA` | Matching genome FASTA, with its `.fai` index. |

For two samples, for example:

```bash
BAM_PATHS=("/data/sample1.bam" "/data/sample2.bam")
TARGET_DIRS=("/results/sample1" "/results/sample2")
```

Multiple entries represent different samples and must share a platform, reference, and barcode tag convention. For multiple sequencing BAMs from the same sample, merge them into one sorted, indexed BAM before running. The underlying CLI also accepts chromosome-split BAM directories with specific filenames; these are not a generic way to pool arbitrary BAMs (see the individual-step implementation).

Parse uses one BAM input/target; its sample outputs are created under `samples/<sample>/`. This wrapper runs the reference-based workflow. For annotation-free workflows or additional advanced CLI options, use the [individual steps](preprocessing.md).

## Local CPU execution

Set `RUNNER="local"` and `CPUS` to the total CPU budget. Annotation, summary, and count each default to that many workers. Compatible mapping divides the metagenes into `COMPATIBLE_JOBS` partitions and executes multiple processes, with concurrency limited by both `MAX_CONCURRENT_JOBS` and `CPUS / COMPATIBLE_CPUS`.

For example, `CPUS=8`, `COMPATIBLE_JOBS=32`, and `COMPATIBLE_CPUS=1` runs up to eight partitions at once until all 32 finish. All partitions must succeed before summary runs. Each compatible process loads its own annotation data, so lower `MAX_CONCURRENT_JOBS` if memory is limited. SLURM memory/time settings do not impose limits in local mode.

`ANNOTATION_CPUS`, `SUMMARY_CPUS`, and `COUNT_CPUS` can override the per-stage workers; their default value of `0` means use `CPUS`. Compatible mapping processes each partition's genes sequentially, so increasing `COMPATIBLE_JOBS` and concurrency is how to parallelize that stage. Numerical-library threads are limited to one per worker to avoid nested thread oversubscription.

## SLURM execution

When you select `RUNNER="slurm"`, review the entire SLURM section, including partition, account, CPU, memory, and time requests. Set `PARTITION`, `ACCOUNT`, and `QOS` according to your cluster's requirements. A blank partition is usable only if the cluster has a default partition; account and QOS can be blank when your cluster does not require them explicitly. The wrapper submits:

```text
annotation job
    -> compatible array (one task per partition)
        -> summary job
            -> count job
```

Each transition uses an `afterok` dependency. Summary waits for the entire array to succeed. Set CPU, memory, and time limits separately for each stage in the template. `COMPATIBLE_MEM` is the memory per array task, not the entire array, and `MAX_CONCURRENT_JOBS` limits active array tasks. The template's resource requests are starting values to adjust for your dataset and cluster.

Activate the SCOTCH environment before submitting. The resolved Python executable is saved with the run; the repository, inputs, environment, and output directories must be accessible from compute nodes. `PYTHON` can explicitly select another SCOTCH Python executable. The wrapper returns after submission; check job completion with SLURM.

## Logs and failures

Each invocation creates `TARGET_DIRS[0]/.scotch_runs/<timestamp>/` containing resolved settings (`config.json`) and logs. SLURM runs also store the submitted scripts and job IDs (`jobs.json`). Editing your config after submission does not change the saved settings for queued jobs.

Local execution stops on failure and terminates its remaining active partitions. SLURM dependencies prevent later stages from running after a failed predecessor. Inspect the failed stage's log before retrying. If submission itself fails partway through, the wrapper reports the jobs already submitted so you can inspect or cancel them with SLURM.

Existing annotation/checkpoint behavior comes from the underlying SCOTCH steps. Use fresh target directories when changing the input data or annotation settings; do not run two wrappers against the same targets simultaneously. If resuming compatible jobs with existing checkpoints, keep the partition count unchanged. The wrapper does not automatically decide which stages to skip; to rerun just one stage, use the [preprocessing commands](preprocessing.md).
