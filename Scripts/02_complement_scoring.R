## 02_complement_scoring.R
## Score malignant cells for the five-gene complement-evasion signature, test whether
## the score differs across malignant cell states, and correlate with T-cell infiltration
## in the single-cell cohort.

library(org.Hs.eg.db)

## --- Complement-evasion gene panel (Entrez IDs matching GSE103322 annotation) ---
complement_genes <- c("1604", "4179", "966", "3075", "710")  # CD55, CD46, CD59, CFH, SERPING1
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = complement_genes,
                                   keytype = "ENTREZID",
                                   columns = "SYMBOL")

## --- Score each malignant cell (mean log-normalized expression across the panel) ---
complement_expr2 <- assay(malignant, "TPM")[complement_genes, keep]
rownames(complement_expr2) <- gene_map$SYMBOL[match(rownames(complement_expr2), gene_map$ENTREZID)]
complement_score2 <- colMeans(complement_expr2)

cluster_score_df2 <- data.frame(
  cluster = clusters_v3,
  complement_score = complement_score2,
  patient = patient_id2
)

## --- Does complement-evasion score differ across malignant cell states? ---
aggregate(complement_score ~ cluster, data = cluster_score_df2, FUN = mean)
kruskal.test(complement_score ~ cluster, data = cluster_score_df2)

boxplot(complement_score ~ cluster, data = cluster_score_df2,
        main = "Complement-Evasion Score by Tumor Cell State",
        xlab = "Cluster", ylab = "Complement-Evasion Score",
        col = "lightblue")
# Figure 1 in manuscript

## --- Per-patient T-cell infiltration fraction (single-cell cohort) ---
all_names <- colnames(sce)
all_is_combo <- grepl("combo", all_names)
all_normalized <- sub("^HNSCC_([0-9]+)_", "HNSCC\\1_", all_names)
all_patient_fixed <- sub("_.*", "", all_normalized)
all_patient_fixed[all_is_combo] <- NA

immune_df2 <- data.frame(
  patient = all_patient_fixed,
  cell_type = colData(sce)$non.cancer.cell.type
)
immune_df2 <- immune_df2[!is.na(immune_df2$patient), ]

tcell_counts2 <- table(immune_df2$patient[immune_df2$cell_type == "T cell"])
total_counts2 <- table(immune_df2$patient)

tcell_fraction2 <- data.frame(
  patient = names(total_counts2),
  tcell_frac = as.numeric(tcell_counts2[names(total_counts2)]) / as.numeric(total_counts2)
)
tcell_fraction2[is.na(tcell_fraction2)] <- 0

## --- Correlate patient-level complement-evasion score with T-cell fraction ---
patient_complement2 <- aggregate(complement_score ~ patient, data = cluster_score_df2, FUN = mean)
merged_df2 <- merge(patient_complement2, tcell_fraction2, by = "patient")
cor.test(merged_df2$complement_score, merged_df2$tcell_frac, method = "spearman")

## Objects used downstream: cluster_score_df2, complement_score2
