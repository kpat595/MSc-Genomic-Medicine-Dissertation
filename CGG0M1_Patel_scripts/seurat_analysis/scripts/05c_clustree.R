#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(clustree)
  library(ggplot2)
})

## ---- CONFIG ------------------------------------------------------------
project_dir <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
input_rds   <- file.path(project_dir, "rds", "04_clustered_sct.rds")
plot_dir    <- file.path(project_dir, "plots", "05_stability")

resolutions_to_test <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0)
res_label <- function(res) gsub("\\.", "_", formatC(res, format = "f", digits = 1))
cluster_col_for <- function(res) paste0("cluster_sct_res_", res_label(res))

## ---- Build the resolution-flow tree straight from existing Louvain cluster
## assignments (Seurat::FindClusters default algorithm = 1) -- no density
## estimation or merge-order chaining involved, so this sidesteps the RSL
## chaining problem in the old clustertree-package approach entirely.
obj <- readRDS(input_rds)

res_cols <- vapply(resolutions_to_test, cluster_col_for, character(1))
stopifnot(all(res_cols %in% colnames(obj@meta.data)))

md <- obj@meta.data[, res_cols, drop = FALSE]
colnames(md) <- paste0("cluster_sct_res_", formatC(resolutions_to_test, format = "f", digits = 1))

p <- clustree(md, prefix = "cluster_sct_res_")
ggsave(file.path(plot_dir, "05c_clustree.png"), p, width = 12, height = 14, units = "in", dpi = 300)

## SC3 stability score: how consistently each cluster's members stay together
## across neighbouring resolutions (1 = perfectly stable, 0 = highly unstable).
p_stab <- clustree(md, prefix = "cluster_sct_res_", node_colour = "sc3_stability")
ggsave(file.path(plot_dir, "05c_clustree_stability.png"), p_stab, width = 12, height = 14, units = "in", dpi = 300)

message("Done. Output in ", plot_dir)
