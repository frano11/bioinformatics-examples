# Repository Overview

This repository contains curated examples of my bioinformatics work, spanning bulk RNA-seq and spatial transcriptomics analyses. Below is a brief description of each file. More files and annotations will be added shortly.

# Mouse Bulk RNA-seq Pipeline

## 1. Quality Control and Trimming of `fastq.gz` Files.

### Processes
- MD5 checksum validation
- FastQC on raw and trimmed reads
- Trim Galore (adapter and quality trimming)
- MultiQC report aggregation

**Nextflow Pipeline:** [`main_fastqc_multiqc_trimming.nf`](main_fastqc_multiqc_trimming.nf)

---

## 2. Alignment and Gene Count Quantification.

### Processes
- Alignment using HISAT2
- Identification of duplicates  
- Feature counts
- QC metrics (rRNA counts, strandness, raw duplicates)
- MultiQC

Go to: `main_analysis_prebuilt_HISAT2_only_FeatCounts_multi_vM25.nf`


## R pipeline for count matrix generation and data processing.

### Processes
- Generation of count matrix .csv from raw read counts files
- Data modification
- PCA
- DESeq2:
	- MA and Volcano plots
	- DEG analysis
- GSEA

Go to: `bulk_RNAseq_L_M_count_matrix_DESeq2.R`


# Mouse Spatial transcriptomics Pipeline

## R pipeline of 10X Genomics Visium dataset from mouse brain slice

### Processes
- Load spatial data and QC analysis
- Data pre-processing
- Dimensional reduction, clustering and visualizations
- Summary table of spots identifiers and spot amount
- DEG analysis
- Brain region annotation
- Flextable of highly expressed genes per brain region

Go to: `260225_Spatial_Mouse_Brain_Coronal_10x.qmd`


