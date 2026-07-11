## 04_TCGA_bulk_analysis.R
## Download TCGA-HNSC bulk RNA-seq and clinical data, and score tumors for the
## complement-evasion signature and an independent cytotoxic/T-cell infiltration signature.

library(TCGAbiolinks)
library(org.Hs.eg.db)

## --- Query and download TCGA-HNSC RNA-seq data ---
query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
tcga_data <- GDCprepare(query)

## --- Match complement-evasion genes to Ensembl IDs (TCGA uses Ensembl, with version
## suffixes stripped for matching) ---
tcga_gene_ids <- sub("\\..*", "", rownames(tcga_data))

gene_map2 <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = c("CD55", "CD46", "CD59", "CFH", "SERPING1"),
                                    keytype = "SYMBOL",
                                    columns = "ENSEMBL")
complement_ensembl <- gene_map2$ENSEMBL

## --- Restrict to primary tumor samples ---
table(colData(tcga_data)$shortLetterCode)
tumor_only <- tcga_data[, colData(tcga_data)$shortLetterCode == "TP"]  # 520 tumor samples

## --- Complement-evasion score per tumor (log2(TPM+1), same scoring method as scRNA data) ---
gene_rows <- match(complement_ensembl, tcga_gene_ids)
tcga_complement_expr <- assay(tumor_only, "tpm_unstrand")[gene_rows, ]
rownames(tcga_complement_expr) <- gene_map2$SYMBOL

tcga_complement_log <- log2(tcga_complement_expr + 1)
tcga_complement_score <- colMeans(tcga_complement_log)

summary(tcga_complement_score)
hist(tcga_complement_score, breaks = 30,
     main = "Complement-Evasion Score Across TCGA-HNSC Patients",
     xlab = "Mean expression (log2 TPM+1)")
# Supplementary Figure S1 in manuscript

## --- Independent cytotoxic/T-cell infiltration score (validation + Section 3.4) ---
immune_genes <- c("CD3D", "CD3E", "CD8A", "GZMB", "PRF1", "NKG7")
gene_map3 <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = immune_genes,
                                    keytype = "SYMBOL",
                                    columns = "ENSEMBL")
immune_ensembl <- gene_map3$ENSEMBL

immune_rows <- match(immune_ensembl, tcga_gene_ids)
tcga_immune_expr <- assay(tumor_only, "tpm_unstrand")[immune_rows, ]
rownames(tcga_immune_expr) <- gene_map3$SYMBOL

tcga_immune_log <- log2(tcga_immune_expr + 1)
tcga_immune_score <- colMeans(tcga_immune_log)

## Objects used downstream: tumor_only, tcga_complement_score, tcga_immune_score
