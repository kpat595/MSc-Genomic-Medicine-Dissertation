#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(sceasy)
  library(reticulate)
})

root      <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
ref_dir   <- file.path(root, "reference")
rds_dir   <- file.path(root, "rds")
max_cells_per_label <- 500
set.seed(1234)

use_condaenv("/well/bsg/projects/resources/shared_conda/conda_environs/cellhint", required = TRUE)
anndata <- import("anndata", convert = FALSE)
np      <- import("numpy", convert = FALSE)

subsample_h5ad <- function(in_path, out_path, label_col, symbol_col, use_raw, backed) {
  ad <- anndata$read_h5ad(in_path, backed = if (backed) "r" else NULL)
  labels <- as.character(py_to_r(ad$obs[[label_col]]))
  keep0 <- unlist(lapply(split(seq_along(labels) - 1L, labels), function(idx) {
    sample(idx, min(length(idx), max_cells_per_label))
  }))

  sub <- ad[np$array(as.integer(keep0)), ]
  sub <- if (backed) sub$to_memory() else sub$copy()
  if (use_raw) sub$X <- sub$raw$X
  sub$var_names <- make.unique(as.character(py_to_r(sub$var[[symbol_col]])))
  sub$write(out_path)

  rm(ad, sub); gc()
  invisible(out_path)
}

## HNOCA
hnoca_small <- file.path(ref_dir, "ref_hnoca_subsampled.h5ad")
subsample_h5ad(file.path(ref_dir, "ref_hnoca.h5ad"), hnoca_small,
                label_col = "annot_level_2", symbol_col = "feature_name",
                use_raw = TRUE, backed = FALSE)
hnoca_seurat <- sceasy::convertFormat(hnoca_small, from = "anndata", to = "seurat",
                                       main_layer = "counts",
                                       outFile = file.path(rds_dir, "ref_hnoca.rds"))
message("ref_hnoca.rds: ", ncol(hnoca_seurat), " cells, ",
        length(unique(hnoca_seurat$annot_level_2)), " annot_level_2 labels")
rm(hnoca_seurat); gc()

## Siletti (neurons + non-neuronal)
sil_neu_small <- file.path(ref_dir, "ref_siletti_neurons_subsampled.h5ad")
sil_non_small <- file.path(ref_dir, "ref_siletti_nonneuronal_subsampled.h5ad")
subsample_h5ad(file.path(ref_dir, "ref_siletti_neurons.h5ad"), sil_neu_small,
                label_col = "supercluster_term", symbol_col = "feature_name",
                use_raw = FALSE, backed = TRUE)
subsample_h5ad(file.path(ref_dir, "ref_siletti_nonneuronal.h5ad"), sil_non_small,
                label_col = "supercluster_term", symbol_col = "feature_name",
                use_raw = FALSE, backed = TRUE)

sil_neu <- anndata$read_h5ad(sil_neu_small)
sil_non <- anndata$read_h5ad(sil_non_small)
sil_combined <- anndata$concat(list(sil_neu, sil_non), join = "outer")
siletti_small <- file.path(ref_dir, "ref_siletti_subsampled.h5ad")
sil_combined$write(siletti_small)
rm(sil_neu, sil_non, sil_combined); gc()

siletti_seurat <- sceasy::convertFormat(siletti_small, from = "anndata", to = "seurat", main_layer = "counts")
siletti_sce <- as.SingleCellExperiment(siletti_seurat)
saveRDS(siletti_sce, file.path(rds_dir, "ref_siletti.rds"))
message("ref_siletti.rds: ", ncol(siletti_sce), " cells, ",
        length(unique(siletti_sce$supercluster_term)), " supercluster_term labels")

## Braun 2023 (developing human brain)
braun_small <- file.path(ref_dir, "ref_braun_subsampled.h5ad")
subsample_h5ad(file.path(ref_dir, "ref_braun.h5ad"), braun_small,
                label_col = "CellClass", symbol_col = "feature_name",
                use_raw = FALSE, backed = TRUE)
braun_seurat <- sceasy::convertFormat(braun_small, from = "anndata", to = "seurat",
                                      main_layer = "counts",
                                      outFile = file.path(rds_dir, "ref_braun.rds"))
message("ref_braun.rds: ", ncol(braun_seurat), " cells, ",
        length(unique(braun_seurat$CellClass)), " CellClass labels")
rm(braun_seurat); gc()

message("\nDone.")
