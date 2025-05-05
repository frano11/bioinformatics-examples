# Mouse Bulk RNAseq Pipeline

## NextFlow pipeline for quality control and trimming of fastq.gz files
from mouse samples.

### Processes
- MD5 checks
- FastQC (raw/trimmed)
- Trim Galore
- MultiQC reports

Go to: `main_fastqc_multiqc_trimming.nf`


## NextFlow pipeline for trimmed fq.gz files to raw read counts from mouse samples.

### Processes
- HISAT2 alignment
- Identification of duplicates  
- Feature counts
- QC metrics (rRNA counts, strandness, raw duplicates)
- MultiQC

Go to: `main_analysis_prebuilt_HISAT2_only_FeatCounts_multi_vM25.nf`


## R pipeline from count matrix generation and data processing using DESeq2 from mouse samples.

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

R pipeline from 10X Genomics Visium dataset from mouse brain slice

## Processes
- Load spatial data and QC analysis
- Data pre-processing
- Dimensional reduction, clustering and visualizations
- Summary table of spots identifiers and spot amount
- DEG analysis
- Brain region annotation
- Flextable of highly expressed genes per brain region

Go to: `260225_Spatial_Mouse_Brain_Coronal_10x.qmd`


