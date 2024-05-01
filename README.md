# SCOTCH
Single-Cell Omics for Transcriptome CHaracterization (SCOTCH): isoform-level characterization of gene expression through long-read single-cell RNA sequencing

## Background
Recent development involving long-read single-cell transcriptome sequencing (lr-scRNA-Seq) represents a significant leap forward in single-cell genomics. With the recent introduction of R10 flowcells by Oxford Nanopore, computational methods should now shift focus on harnessing the unique benefits of long reads to analyze transcriptome complexity. In this context, we introduce a comprehensive suite of computational methods named Single-Cell Omics for Transcriptome CHaracterization (SCOTCH). Our method is compatible with the single-cell library preparation platform from both 10X Genomics and Parse Biosciences, facilitating the analysis of special cell populations, such as neurons, hepatocytes and developing cardiomyocytes. 

## Installation
### Step1: install sources
To avoid package version conflicts, we strongly recommand to use conda to set up the environment. If you don't have conda installed, you can run codes below in linux to install.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
bash Miniconda3-latest-Linux-x86_64.sh
```
After installing conda, SCOTCH sources can be downloaded:
git clone https://github.com/WGLab/SCOTCH.git

```
cd SCOTCH
conda env create --name SCOTCH --file SCOTCH.yml
conda activate SCOTCH
```
