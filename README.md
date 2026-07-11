# Complement-Evasion HNSCC

Single-cell and bulk RNA-seq analysis of tumor-intrinsic complement regulatory gene 
expression in head and neck squamous cell carcinoma (HNSCC).

## Overview

This project characterizes a five-gene complement-evasion signature (CD55, CD46, CD59, 
CFH, SERPING1) across malignant cell states (single-cell RNA-seq) and bulk tumors (TCGA), 
testing its relationship to immune infiltration, survival, and tumor cell differentiation state.

**Preprint:** [link to be added once live on bioRxiv]

## Data Sources

- **Single-cell RNA-seq:** GEO accession [GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322) 
  (Puram et al., 2017), accessed via Bioconductor ExperimentHub (record EH5419)
- **Bulk RNA-seq and clinical data:** TCGA-HNSC, accessed via [TCGAbiolinks](https://bioconductor.org/packages/TCGAbiolinks/) 
  and the NCI Genomic Data Commons

Raw data is not included in this repository; scripts re-download it directly from the 
sources above.

## Pipeline / Scripts

| Script | Description |
|---|---|
| `01_scRNA_QC_clustering.R` | Load GSE103322, QC, batch-effect correction, malignant cell clustering |
| `02_complement_scoring.R` | Complement-evasion gene panel scoring |
| `03_marker_genes_enrichment.R` | Cluster marker genes and GO enrichment |
| `04_TCGA_bulk_analysis.R` | TCGA-HNSC download, complement/immune scoring |
| `05_survival_analysis.R` | Kaplan-Meier and Cox regression |

## Requirements

R (≥4.5), with Bioconductor packages: `ExperimentHub`, `SingleCellExperiment`, `scran`, 
`scater`, `bluster`, `org.Hs.eg.db`, `TCGAbiolinks`, `clusterProfiler`, `survival`, `survminer`.

## Author

Alfina Mariam Sanoj 
