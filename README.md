# Complement-Evasion HNSCC

Single-cell and bulk RNA-seq analysis of tumor-intrinsic complement regulatory gene 
expression in head and neck squamous cell carcinoma (HNSCC).

## Overview

This project characterizes a five-gene complement-evasion signature (CD55, CD46, CD59, 
CFH, SERPING1) across malignant cell states (single-cell RNA-seq) and bulk tumors (TCGA), 
testing its relationship to immune infiltration, survival, and tumor cell differentiation state.

**Preprint:** [link to be added once live on bioRxiv]

## Pipeline Overview

```mermaid
flowchart TD
    A[GSE103322 scRNA-seq<br/>5,902 cells, 18 patients] --> B[QC + malignant cell ID<br/>2,215 cells]
    B --> C[Batch-effect correction<br/>per-patient centering]
    C --> D[Clustering<br/>7 malignant cell states]
    D --> E[Complement-evasion scoring<br/>CD55/CD46/CD59/CFH/SERPING1]
    D --> F[Marker genes + GO enrichment<br/>high vs. low clusters]
    E --> G[TCGA-HNSC bulk data<br/>520 tumors]
    G --> H[Survival analysis<br/>Kaplan-Meier + Cox]
    G --> I[Immune infiltration correlation]
```

## Key Results

**Figure 1 — Complement-evasion score differs significantly across malignant cell states**
![Figure 1](figures/Figure1_complement_evasion_by_cluster.png)

**Figure 2 — No independent survival association**
![Figure 2](figures/Figure2_survival_curves.png)

**Figure 3 — High complement-evasion cluster shows squamous differentiation + antimicrobial defense enrichment**
![Figure 3](figures/Figure3_GO_enrichment.png)

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
| `Scripts/01_scRNA_QC_clustering.R` | Load GSE103322, QC, batch-effect correction, malignant cell clustering |
| `Scripts/02_complement_scoring.R` | Complement-evasion gene panel scoring |
| `Scripts/03_marker_genes_enrichment.R` | Cluster marker genes and GO enrichment |
| `Scripts/04_TCGA_bulk_analysis.R` | TCGA-HNSC download, complement/immune scoring |
| `Scripts/05_survival_analysis.R` | Kaplan-Meier and Cox regression |

## Requirements

R (≥4.5), with Bioconductor packages: `ExperimentHub`, `SingleCellExperiment`, `scran`, 
`scater`, `bluster`, `org.Hs.eg.db`, `TCGAbiolinks`, `clusterProfiler`, `survival`, `survminer`.

## Author

Alfina Mariam Sanoj 
