## 03_marker_genes_enrichment.R
## Identify marker genes distinguishing the highest and lowest complement-evasion
## malignant cell clusters, and run GO Biological Process enrichment on the
## high-complement-evasion cluster's markers.

library(scran)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

malignant_clean <- malignant[, keep]

## --- Strict marker gene test: a gene must be significantly up-regulated in the target
## cluster relative to ALL other clusters (pval.type = "all"), not just one (the
## default "any"), to avoid genes being called markers of multiple clusters at once. ---
markers_strict <- findMarkers(malignant_clean, groups = clusters_v3,
                               test.type = "wilcox",
                               pval.type = "all",
                               direction = "up")

## Cluster 6 = highest mean complement-evasion score; cluster 1 = lowest
## (see aggregate() output in 02_complement_scoring.R)
top_cluster6_strict <- markers_strict[["6"]]
top_cluster1_strict <- markers_strict[["1"]]

ids_cluster6 <- rownames(top_cluster6_strict)[1:15]
ids_cluster1 <- rownames(top_cluster1_strict)[1:15]

genes_cluster6 <- AnnotationDbi::select(org.Hs.eg.db, keys = ids_cluster6,
                                         keytype = "ENTREZID", columns = "SYMBOL")
genes_cluster1 <- AnnotationDbi::select(org.Hs.eg.db, keys = ids_cluster1,
                                         keytype = "ENTREZID", columns = "SYMBOL")
genes_cluster6
genes_cluster1

## --- GO (Biological Process) enrichment on top 100 markers of cluster 6 ---
ids_cluster6_enrich <- rownames(top_cluster6_strict)[1:100]

ego <- enrichGO(gene = ids_cluster6_enrich,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)

ego_df <- as.data.frame(ego)
ego_df[1:9, c("Description", "pvalue", "p.adjust", "Count")]

barplot(ego, showCategory = 9) +
  ggtitle("GO Biological Process Enrichment: High Complement-Evasion Cluster (6)")
ggsave("GO_enrichment_barplot.pdf", width = 8, height = 5)
# Figure 3 in manuscript
