#!/usr/bin/env Rscript

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))

suppressPackageStartupMessages({
  library(future)
  library(Seurat)
  library(ggplot2)
  library(glmGamPoi)
})

stopifnot(requireNamespace("glmGamPoi", quietly = TRUE))
options(future.globals.maxSize = Inf)
future::plan("sequential")

DEBUG <- FALSE
debug_cells_per_sample <- 500
out_suffix <- if (DEBUG) "_debug" else ""

conditions <- c("NR1", "CASPR2", "Control", "NoTreated")

project_root <- "/well/bsg/projects/CGG0M1_Patel"
analysis_dir <- file.path(project_root, "seurat_analysis")
by_condition_dir <- file.path(analysis_dir, ".by_condition")
manifest_path <- file.path(analysis_dir, "tables", "sample_manifest.tsv")
filtered_rds_dir <- file.path(analysis_dir, "rds", "02_filtered")

manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
sample_ids <- manifest$sample_id
if (DEBUG) sample_ids <- head(sample_ids, 2)

pools <- unique(manifest$pool[manifest$sample_id %in% sample_ids])

sample_objs <- lapply(sample_ids, function(x) readRDS(file.path(filtered_rds_dir, paste0("02_filtered_", x, ".rds"))))
names(sample_objs) <- sample_ids
if (DEBUG) {
  set.seed(1)
  sample_objs <- lapply(sample_objs, function(x) subset(x, cells = sample(colnames(x), min(debug_cells_per_sample, ncol(x)))))
}

# merge within pool (all conditions in that pool together) so SCTransform is fit per pool
pool_objs <- lapply(pools, function(p) {
  ids <- manifest$sample_id[manifest$pool == p & manifest$sample_id %in% sample_ids]
  objs_p <- sample_objs[ids]
  if (length(objs_p) == 1) return(RenameCells(objs_p[[1]], add.cell.id = ids))
  merge(objs_p[[1]], y = objs_p[-1], add.cell.ids = ids, project = p)
})
names(pool_objs) <- pools
pool_objs <- lapply(pool_objs, JoinLayers, assay = "RNA")
rm(sample_objs)
gc()

pool_objs <- lapply(pool_objs, function(x) {
  DefaultAssay(x) <- "RNA"
  SCTransform(x, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", conserve.memory = TRUE, verbose = FALSE)
})

filter_summary_all_path <- file.path(analysis_dir, "tables", "01_filter_summary_all.tsv")
if (!file.exists(filter_summary_all_path)) stop("Missing combined filter summary: ", filter_summary_all_path)

# split the pool-normalized data out by condition and finish the rest of the pipeline per condition
for (target_condition in conditions) {

  rds_dir <- file.path(by_condition_dir, "rds_by_condition", target_condition)
  tables_dir <- file.path(by_condition_dir, "tables_by_condition", target_condition)
  plots_dir <- file.path(by_condition_dir, "plots_by_condition", target_condition, if (DEBUG) "03_merge_sct_pca_debug" else "03_merge_sct_pca")

  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  cond_objs <- lapply(pool_objs, function(x) subset(x, subset = condition == target_condition))
  cond_objs <- Filter(function(x) ncol(x) > 0, cond_objs)
  if (length(cond_objs) < 1) stop("No pools contain cells for condition: ", target_condition)

  merged <- if (length(cond_objs) == 1) cond_objs[[1]] else merge(cond_objs[[1]], y = cond_objs[-1], project = target_condition, merge.data = TRUE)
  rm(cond_objs)

  DefaultAssay(merged) <- "SCT"
  merged <- FindVariableFeatures(merged, assay = "SCT", selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  pca_features <- VariableFeatures(merged)
  if (length(pca_features) < 2) stop("Too few SCT variable features before PCA: ", length(pca_features))

  # VariableFeatures() already returns genes ranked most- to least-variable
  feature_table <- data.frame(gene = pca_features, rank = seq_along(pca_features))
  write.table(feature_table, file.path(tables_dir, paste0("02_sct_variable_features_", target_condition, out_suffix, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  missing_features <- setdiff(pca_features, rownames(merged[["SCT"]]@scale.data))
  if (length(missing_features) > 0) merged <- GetResidual(object = merged, features = missing_features, assay = "SCT", umi.assay = "RNA", verbose = FALSE)

  pca_features <- intersect(pca_features, rownames(merged[["SCT"]]@scale.data))
  if (length(pca_features) < 2) stop("Too few merged SCT residual features for PCA: ", length(pca_features))

  npcs_use <- min(50, length(pca_features) - 1, ncol(merged) - 1)
  merged <- RunPCA(merged, features = pca_features, npcs = npcs_use, verbose = FALSE)

  merged <- NormalizeData(merged, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = FALSE)

  saveRDS(merged, file.path(rds_dir, paste0("03_merged_sct_pca", out_suffix, ".rds")))

  rna_vln <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_id", pt.size = 0, ncol = 3)
  ggsave(file.path(plots_dir, "merged_postfilter_vln_rna_qc_by_sample.png"), rna_vln, width = 14, height = 5)

  top_features <- head(feature_table, 30)
  top_features$plot_score <- rev(seq_len(nrow(top_features)))
  vf_plot <- ggplot(top_features, aes(x = reorder(gene, plot_score), y = plot_score)) + geom_col() + coord_flip() + labs(x = "Gene", y = "Variable feature rank (higher = more variable)", title = paste("Condition:", target_condition))
  ggsave(file.path(plots_dir, "merged_sct_variable_features.png"), vf_plot, width = 8, height = 8)

  ggsave(file.path(plots_dir, "merged_pca_elbowplot.png"), ElbowPlot(merged, ndims = npcs_use), width = 8, height = 5)
  ggsave(file.path(plots_dir, "merged_pca_by_pool.png"), DimPlot(merged, reduction = "pca", group.by = "pool", raster = TRUE), width = 7, height = 5)

  adt_clr <- VlnPlot(merged, assay = "ADT", features = rownames(merged[["ADT"]]), group.by = "pool", slot = "data", pt.size = 0, ncol = 3)
  ggsave(file.path(plots_dir, "merged_vln_adt_clr_by_pool.png"), adt_clr, width = 12, height = 5)

  message("Condition: ", target_condition)
  message("Saved: ", file.path(rds_dir, paste0("03_merged_sct_pca", out_suffix, ".rds")))
  message("Saved: ", file.path(tables_dir, paste0("02_sct_variable_features_", target_condition, out_suffix, ".tsv")))

  rm(merged)
  gc()
}

message("Checked: ", filter_summary_all_path)
message("Done: all conditions processed from a single pool-level SCTransform pass.")