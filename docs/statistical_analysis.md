# Statistical analysis

[Back to README](../README.md)

The R package compares two cell populations using gene and transcript count matrices. `scotch_gene()` tests differential gene expression; `scotch_transcript()` tests differential transcript usage and reports isoform-switching statistics. Matrices have cells as rows and genes or transcripts as columns.

### Installation

In R, run below codes to install SCOTCH statistical pipeline.

```
if (!requireNamespace("devtools", quietly = TRUE))
install.packages("devtools")
library("devtools")
install_github("WGLab/SCOTCH")
```

### DTU analysis
Below is sample codes for DTU analysis. Example count matrices are in [`data/`](../data). Run the example from that directory after installing the R package.

```
library(SCOTCH)
library(dplyr)
library(stringr)

#----read gene-level count matrix-----#
sample8_CD4_gene=as.matrix(read.csv('gene_count_1.csv',row.names = 'X'))
sample8_CD8_gene=as.matrix(read.csv('gene_count_2.csv',row.names = 'X'))

#----read transcript-level count matrix-----#
sample8_CD4_transcript=as.matrix(read.csv("transcript_count_1.csv",row.names = 'X'))
gene_transcript_CD4_df = data.frame(genes=str_remove(colnames(sample8_CD4_transcript),"\\.(ENST|novel|uncategorized).+"),
                                    transcripts=colnames(sample8_CD4_transcript))

sample8_CD8_transcript=as.matrix(read.csv("transcript_count_2.csv",row.names = 'X'))
gene_transcript_CD8_df = data.frame(genes=str_remove(colnames(sample8_CD8_transcript),"\\.(ENST|novel|uncategorized).+"),
                                    transcripts=colnames(sample8_CD8_transcript))

#----gene-level analysis-----#
df_gene = scotch_gene(sample8_CD4_gene, sample8_CD8_gene, epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#----transcript-level analysis-----#
df_transcript = scotch_transcript(gene_transcript_CD4_df,gene_transcript_CD8_df,
                                  sample8_CD4_transcript, sample8_CD8_transcript, ncores=10)
df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))
```

SCOTCH will output the following statistics in results:
```
genes: gene name
pct1: percent of cells expression the gene in cell population 1
pct2: percent of cells expression the gene in cell population 2
logFC: fold change in log scale
p_gene: p value on the gene level
pct_diff: difference between pct1 and pct2
p_gene_adj: adjusted p value of differential gene expression
isoforms: isoform name
alpha1: estimated parameter of the dirichlet distribution for cell population 1
alpha2: estimated parameter of the dirichlet distribution for cell population 2
TU1: relative transcript usage of cell population 1
TU2: relative transcript usage of cell population 2
TU_diff: difference of TU1 and TU2
TU_var1: variance of TU1
TU_var2: variance of TU2
dispersion1: dispersion parameter for cell population 1
dispersion2: dispersion parameter for cell population 2
isoform_switch: whether there is isoform switching event
isoform_switch_ES: effect size of isoform switching event
p_DTU_gene: p value for differential transcript usage analysis for the whole gene
p_transcript: p value for differential transcript usage analysis for the specific transcript
p_transcript_adj: adjusted p value for differential transcript usage analysis for the specific transcript
p_DTU_gene_adj: adjusted p value for differential transcript usage analysis for the whole gene
```
