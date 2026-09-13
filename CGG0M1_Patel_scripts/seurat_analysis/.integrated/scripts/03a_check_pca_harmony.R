#!/usr/bin/env Rscript

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
})

integrated_dir <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis/.integrated"
tables_dir <- file.path(integrated_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(integrated_dir, "rds", "03_merged_sct_pca.rds"))

harmony_dims <- dim(Embeddings(obj, "harmony"))
harmony_sd <- Stdev(obj[["harmony"]])
harmony_var <- harmony_sd^2

harmony_summary <- data.frame(
  PC = paste0("harmony_", seq_along(harmony_sd)),
  stdev = harmony_sd,
  percent_variance = 100 * harmony_var / sum(harmony_var),
  cumulative_percent_variance = 100 * cumsum(harmony_var) / sum(harmony_var)
)

write.table(harmony_summary, file.path(tables_dir, "03_pca_variance_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

loadings <- Loadings(obj[["harmony"]])
top_pc_genes <- do.call(rbind, lapply(seq_len(ncol(loadings)), function(pc) {
  x <- sort(abs(loadings[, pc]), decreasing = TRUE)
  data.frame(PC = paste0("harmony_", pc), gene = names(head(x, 20)), abs_loading = as.numeric(head(x, 20)))
}))

write.table(top_pc_genes, file.path(tables_dir, "03_pca_top_loading_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Harmony embedding dimensions cells x PCs:\n")
print(harmony_dims)
cat("Saved harmony PCA summaries.\n")
