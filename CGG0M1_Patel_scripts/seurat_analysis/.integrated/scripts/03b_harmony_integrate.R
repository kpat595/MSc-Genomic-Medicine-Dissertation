#!/usr/bin/env Rscript

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

DEBUG <- FALSE
out_suffix <- if (DEBUG) "_debug" else ""

project_root <- "/well/bsg/projects/CGG0M1_Patel"
analysis_dir <- file.path(project_root, "seurat_analysis")
integrated_dir <- file.path(analysis_dir, ".integrated")

input_rds_dir <- file.path(analysis_dir, "rds")
rds_dir <- file.path(integrated_dir, "rds")
plots_dir <- file.path(integrated_dir, "plots", "03_merge_sct_pca")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

merged <- readRDS(file.path(input_rds_dir, paste0("03_merged_sct_pca", out_suffix, ".rds")))

npcs_use <- ncol(Embeddings(merged, "pca"))
merged <- RunHarmony(merged, group.by.vars = "pool", reduction.use = "pca", reduction.save = "harmony", dims.use = 1:npcs_use, verbose = FALSE)

saveRDS(merged, file.path(rds_dir, paste0("03_merged_sct_pca", out_suffix, ".rds")))

ggsave(file.path(plots_dir, "merged_harmony_elbowplot.png"), ElbowPlot(merged, ndims = npcs_use, reduction = "harmony"), width = 8, height = 5)
ggsave(file.path(plots_dir, "merged_harmony_by_condition.png"), DimPlot(merged, reduction = "harmony", group.by = "condition", raster = TRUE), width = 7, height = 5)
ggsave(file.path(plots_dir, "merged_harmony_by_pool.png"), DimPlot(merged, reduction = "harmony", group.by = "pool", raster = TRUE), width = 7, height = 5)

message("Saved: ", file.path(rds_dir, paste0("03_merged_sct_pca", out_suffix, ".rds")))
