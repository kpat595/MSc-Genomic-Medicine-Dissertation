#!/usr/bin/env Rscript

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
})

obj <- readRDS("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/rds/03_merged_sct_pca.rds")

pca_dims <- dim(Embeddings(obj, "pca"))
pc_sd <- Stdev(obj[["pca"]])
pc_var <- pc_sd^2

pca_summary <- data.frame(
  PC = paste0("PC_", seq_along(pc_sd)),
  stdev = pc_sd,
  percent_variance = 100 * pc_var / sum(pc_var),
  cumulative_percent_variance = 100 * cumsum(pc_var) / sum(pc_var)
)

write.table(pca_summary, "/well/bsg/projects/CGG0M1_Patel/seurat_analysis/tables/03_pca_variance_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

loadings <- Loadings(obj[["pca"]])
top_pc_genes <- do.call(rbind, lapply(seq_len(ncol(loadings)), function(pc) {
  x <- sort(abs(loadings[, pc]), decreasing = TRUE)
  data.frame(PC = paste0("PC_", pc), gene = names(head(x, 20)), abs_loading = as.numeric(head(x, 20)))
}))

write.table(top_pc_genes, "/well/bsg/projects/CGG0M1_Patel/seurat_analysis/tables/03_pca_top_loading_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("PCA embedding dimensions cells x PCs:\n")
print(pca_dims)
cat("Saved PCA summaries.\n")
