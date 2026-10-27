# Copy this file to my_run.sh, edit it, then run: bash scotch.sh my_run.sh
# This is a Bash configuration. Paths may be absolute or relative to this file.

# 1. REQUIRED CHOICES: set these for your data and compute environment.
PLATFORM=""                # REQUIRED: 10x-ont, 10x-pacbio, or parse-ont.
BULK=""                    # REQUIRED: true for bulk RNA-seq; false for single-cell RNA-seq.
RUNNER=""                  # REQUIRED: local for this machine, or slurm for cluster submission.

# 2. REQUIRED PATHS: one BAM and output directory per sample, in matching order.
BAM_PATHS=("/path/to/sample.sorted.bam") # Genome-aligned, sorted, indexed BAM file(s).
TARGET_DIRS=("/path/to/output/sample1") # One output directory per entry in BAM_PATHS.
REFERENCE_GTF="/path/to/genes.gtf"       # Gene annotation matching the BAM genome assembly and contig names.
REFERENCE_FASTA="/path/to/genome.fa"     # Matching genome FASTA; prepare its .fai index before running.
# Multiple samples (all must use the same platform and barcode tag names):
# BAM_PATHS=("/data/sample1.bam" "/data/sample2.bam")
# TARGET_DIRS=("/results/sample1" "/results/sample2")
# For multiple BAMs belonging to ONE sample, merge, sort, and index first.
# Parse requires ONE BAM input/target; samples are identified by the pS tag.

# 3. SINGLE-CELL BARCODE SETTINGS: required for 10x-ont / 10x-pacbio when BULK=false.
BARCODE_CELL=""            # Set the actual cell-barcode tag in your BAM (e.g. CB or XC); no automatic guess.
BARCODE_UMI=""             # Set the actual UMI tag in your BAM (e.g. UB or XM); no automatic guess.
# Common combinations: 10x-ont CB/UB; PacBio CB/XM or XC/XM, depending on preprocessing.
# BULK=true: leave both empty; barcode tags are not used.
# parse-ont: leave both empty; SCOTCH reads cell/UMI from read names and sample from pS.
BUILD=""                   # Set if Parse contigs have a genome prefix, e.g. GRCh38 for GRCh38_chr1.

# 4. EXECUTION RESOURCES: review these for the machine or cluster you will use.
CPUS=8                     # Total local CPU budget; default CPUs for non-array SLURM jobs.
COMPATIBLE_JOBS=32         # Number of gene partitions; one process locally or array task per partition.
MAX_CONCURRENT_JOBS=8      # Maximum simultaneous compatible jobs; local mode also respects CPUS.
ANNOTATION_CPUS=0          # CPUs for annotation; 0 uses CPUS.
COMPATIBLE_CPUS=1          # CPUs per compatible job; genes within a partition run sequentially.
SUMMARY_CPUS=0             # CPUs for read-mapping summary; 0 uses CPUS.
COUNT_CPUS=0               # CPUs for count-matrix generation; 0 uses CPUS.
PYTHON=""                  # Optional SCOTCH Python executable; empty uses the wrapper's Python.

# 5. SLURM SETTINGS: review this entire section when RUNNER=slurm; ignored for local.
PARTITION=""               # Set your cluster partition; empty works only if a default partition is configured.
ACCOUNT=""                 # Set your allocation/account if required by your cluster.
QOS=""                     # Set your quality of service if required by your cluster.
ANNOTATION_MEM="64G"        # Total memory requested for the annotation job.
COMPATIBLE_MEM="16G"        # Total memory requested for EACH compatible array task.
SUMMARY_MEM="16G"           # Total memory requested for the summary job.
COUNT_MEM="32G"             # Total memory requested for the count-matrix job.
ANNOTATION_TIME="24:00:00"  # Annotation wall-time limit (hours:minutes:seconds).
COMPATIBLE_TIME="24:00:00"  # Wall-time limit for EACH compatible array task.
SUMMARY_TIME="04:00:00"     # Summary wall-time limit.
COUNT_TIME="12:00:00"       # Count-matrix wall-time limit.

# 6. OPTIONAL ANALYSIS DEFAULTS: change these only when you want different behavior.
UPDATE_GTF=true            # Refine exon annotations using reads; false uses the supplied exon annotation.
REFERENCE_PKL=""           # Optional matching SCOTCH annotation pickle; empty builds it from REFERENCE_GTF.
SAVE_MEM_ANN=false         # Process BAM annotation in chunks to reduce memory use.
SAVE_MEM=true              # Use on-disk read lookups during isoform mapping to reduce memory use.

# Isoform discovery and output
NOVEL_DISCOVERY_MAX_READS=10000 # Maximum reads used for novel discovery per gene; 0 removes the cap.
NOVEL_DISCOVERY_MAX_ROUNDS=50   # Maximum novel-discovery re-clustering rounds; 0 removes the cap.
NOVEL_READ_N=0                 # Minimum supporting reads for retaining a novel isoform.
NOVEL_READ_PCT=0               # Minimum fraction of gene coverage for retaining a novel isoform (0-1).
GROUP_NOVEL=true               # Group novel isoforms potentially caused by read truncation.
SAVE_CSV=false                 # Also write CSV count matrices; MTX output is always enabled here.
