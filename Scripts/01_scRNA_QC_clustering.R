## 01_scRNA_QC_clustering.R
## Load GSE103322, identify malignant cells, correct patient-of-origin batch effect,
## and cluster malignant cells into transcriptional states.

library(ExperimentHub)
library(scran)
library(scater)
library(bluster)
library(irlba)

## --- Load data ---
eh <- ExperimentHub()
query(eh, "GSE103322")
sce <- eh[["EH5419"]]
sce

## --- Identify malignant cells ---
table(colData(sce)$classified..as.cancer.cell)
malignant <- sce[, colData(sce)$classified..as.cancer.cell == "1"]
dim(malignant)

## --- Fix cell-barcode naming inconsistency (HNSCC_17 vs HNSCC17; exclude "combo" cells
## whose true patient identity cannot be resolved) ---
malignant_names <- colnames(malignant)
malignant_is_combo <- grepl("combo", malignant_names)
malignant_normalized <- sub("^HNSCC_([0-9]+)_", "HNSCC\\1_", malignant_names)
patient_id_fixed <- sub("_.*", "", malignant_normalized)
patient_id_fixed[malignant_is_combo] <- NA

keep <- !is.na(patient_id_fixed)
length(unique(patient_id_fixed[keep]))  # 18 uniquely identifiable patients

## --- Highly variable genes ---
logcounts(malignant) <- assay(malignant, "TPM")
gene_var <- modelGeneVar(malignant)
top_genes <- getTopHVGs(gene_var, n = 2000)

## --- Batch correction: center each gene within each patient before PCA/clustering.
## (Raw clustering without this step produced clusters that matched patient identity,
## not biology - a known confound in single-patient tumor scRNA-seq.) ---
expr_mat2 <- assay(malignant, "TPM")[top_genes, keep]
patient_id2 <- patient_id_fixed[keep]

relative_expr2 <- expr_mat2
for (p in unique(patient_id2)) {
  cells_in_patient <- which(patient_id2 == p)
  if (length(cells_in_patient) > 1) {
    relative_expr2[, cells_in_patient] <- expr_mat2[, cells_in_patient, drop = FALSE] -
      rowMeans(expr_mat2[, cells_in_patient, drop = FALSE])
  } else {
    relative_expr2[, cells_in_patient] <- 0
  }
}

## --- PCA and graph-based clustering on the batch-corrected matrix ---
pca_relative2 <- prcomp_irlba(t(relative_expr2), n = 20)
pca_scores2 <- pca_relative2$x
rownames(pca_scores2) <- colnames(expr_mat2)

clusters_v3 <- clusterRows(pca_scores2, NNGraphParam())
table(clusters_v3)                     # cluster sizes
table(clusters_v3, patient_id2)        # confirm clusters mix across patients (not patient-driven)

## Objects used downstream: malignant, keep, patient_id2, clusters_v3, top_genes
