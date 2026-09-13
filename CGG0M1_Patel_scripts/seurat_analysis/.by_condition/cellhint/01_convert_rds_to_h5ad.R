#!/usr/bin/env Rscript
# Convert per-condition clustered Seurat objects to h5ad for CellHint.
# Requires: sceasy, reticulate (with anndata python package importable)

suppressPackageStartupMessages({library(Seurat); library(sceasy)})

# sceasy 0.0.7's internal seurat2anndata() calls Seurat::GetAssayData(..., slot=)
# which is defunct as of SeuratObject 5.0.0 (use `layer=` instead). Patch the
# function in-place for this session rather than touching the installed package.
seurat2anndata_fixed <- function(obj, outFile = NULL, assay = "RNA", main_layer = "data",
                                  transfer_layers = NULL, drop_single_values = TRUE) {
  main_layer <- match.arg(main_layer, c("data", "counts", "scale.data"))
  transfer_layers <- transfer_layers[transfer_layers %in% c("data", "counts", "scale.data")]
  transfer_layers <- transfer_layers[transfer_layers != main_layer]

  if (compareVersion(as.character(obj@version), "3.0.0") < 0) {
    obj <- Seurat::UpdateSeuratObject(object = obj)
  }

  X <- Seurat::GetAssayData(object = obj, assay = assay, layer = main_layer)
  obs <- sceasy:::.regularise_df(obj@meta.data, drop_single_values = drop_single_values)
  var <- sceasy:::.regularise_df(Seurat::GetAssay(obj, assay = assay)[[]],
                                  drop_single_values = drop_single_values)

  obsm <- NULL
  reductions <- names(obj@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(reductions, function(name) as.matrix(Seurat::Embeddings(obj, reduction = name)),
                    simplify = FALSE)
    names(obsm) <- paste0("X_", tolower(names(obj@reductions)))
  }

  layers <- list()
  for (layer in transfer_layers) {
    mat <- Seurat::GetAssayData(object = obj, assay = assay, layer = layer)
    if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
  }

  anndata <- reticulate::import("anndata", convert = FALSE)
  adata <- anndata$AnnData(X = Matrix::t(X), obs = obs, var = var, obsm = obsm, layers = layers)
  if (!is.null(outFile)) adata$write(outFile, compression = "gzip")
  adata
}
assignInNamespace("seurat2anndata", seurat2anndata_fixed, ns = "sceasy")

conditions <- c("NR1", "CASPR2", "Control", "NoTreated")
args <- commandArgs(trailingOnly = TRUE)
target_condition <- conditions[as.integer(args[1])]
if (is.na(target_condition)) stop("Usage: 01_convert_rds_to_h5ad.R <1-4>, got array index ", args[1])

resolution <- "0_4"
cluster_col <- paste0("cluster_sct_res_", resolution)

root <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
input_rds <- file.path(root, ".by_condition", "rds_by_condition", target_condition, "04_clustered_sct.rds")
out_dir <- file.path(root, ".by_condition", "cellhint", "h5ads")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output_h5ad <- file.path(out_dir, paste0("01_", target_condition, ".h5ad"))

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")

if (!cluster_col %in% colnames(obj@meta.data)) {
  stop("Cluster column ", cluster_col, " not found in ", input_rds)
}

obj$condition <- target_condition
obj$cellhint_cluster <- as.character(obj@meta.data[[cluster_col]])

obj <- DietSeurat(obj, assays = "RNA")

sceasy::convertFormat(
  obj,
  from = "seurat",
  to = "anndata",
  assay = "RNA",
  main_layer = "counts",
  outFile = output_h5ad
)

cat("Wrote", output_h5ad, "\n")
