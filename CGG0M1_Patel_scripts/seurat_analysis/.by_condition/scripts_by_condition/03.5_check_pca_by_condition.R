#!/usr/bin/env Rscript

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
})

conditions <- c("NR1", "CASPR2", "Control", "NoTreated")
args <- commandArgs(trailingOnly = TRUE)
target_condition <- conditions[as.integer(args[1])]
if (is.na(target_condition)) stop("Usage: 03.5_check_pca_by_condition.R <1-4>, got array index ", args[1])

by_condition_dir <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition"
rds_dir <- file.path(by_condition_dir, "rds_by_condition", target_condition)
tables_dir <- file.path(by_condition_dir, "tables_by_condition", target_condition)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(rds_dir, "03_merged_sct_pca.rds"))

pca_dims <- dim(Embeddings(obj, "pca"))
pc_sd <- Stdev(obj[["pca"]])
pc_var <- pc_sd^2

pca_summary <- data.frame(
  PC = paste0("PC_", seq_along(pc_sd)),
  stdev = pc_sd,
  percent_variance = 100 * pc_var / sum(pc_var),
  cumulative_percent_variance = 100 * cumsum(pc_var) / sum(pc_var)
)

write.table(pca_summary, file.path(tables_dir, "03_pca_variance_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

loadings <- Loadings(obj[["pca"]])
top_pc_genes <- do.call(rbind, lapply(seq_len(ncol(loadings)), function(pc) {
  x <- sort(abs(loadings[, pc]), decreasing = TRUE)
  data.frame(PC = paste0("PC_", pc), gene = names(head(x, 20)), abs_loading = as.numeric(head(x, 20)))
}))

write.table(top_pc_genes, file.path(tables_dir, "03_pca_top_loading_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Condition:", target_condition, "\n")
cat("PCA embedding dimensions cells x PCs:\n")
print(pca_dims)
cat("Saved PCA summaries.\n")
