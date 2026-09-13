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

project_root <- "/well/bsg/projects/CGG0M1_Patel"
analysis_dir <- file.path(project_root, "seurat_analysis")
manifest_path <- file.path(analysis_dir, "tables", "sample_manifest.tsv")
rds_dir <- file.path(analysis_dir, "rds")
tables_dir <- file.path(analysis_dir, "tables")
plots_dir <- file.path(analysis_dir, "plots", if (DEBUG) "03_merge_sct_pca_debug" else "03_merge_sct_pca")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
sample_ids <- manifest$sample_id
if (DEBUG) sample_ids <- head(sample_ids, 2)

pools <- unique(manifest$pool[manifest$sample_id %in% sample_ids])

sample_objs <- lapply(sample_ids, function(x) readRDS(file.path(rds_dir, "02_filtered", paste0("02_filtered_", x, ".rds"))))
names(sample_objs) <- sample_ids
if (DEBUG) {
  set.seed(1)
  sample_objs <- lapply(sample_objs, function(x) subset(x, cells = sample(colnames(x), min(debug_cells_per_sample, ncol(x)))))
}

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
  SCTransform(x, vst.flavor = "v2", method = "glmGamPoi", variable.features.n = 3000, vars.to.regress = "percent.mt", conserve.memory = TRUE, verbose = FALSE)
})

vf_list <- lapply(pool_objs, VariableFeatures)
feature_scores <- sort(table(unlist(vf_list)), decreasing = TRUE)

pca_features <- names(feature_scores)[seq_len(min(3000, length(feature_scores)))]
if (length(pca_features) < 2) stop("Too few consensus SCT variable features before PCA: ", length(pca_features))

feature_table <- data.frame(gene = names(feature_scores), score = as.numeric(feature_scores), selected_for_pca = names(feature_scores) %in% pca_features)
write.table(feature_table, file.path(tables_dir, paste0("02_sct_consensus_variable_features", out_suffix, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

pool_objs <- lapply(pool_objs, function(x) {
  DefaultAssay(x) <- "SCT"
  missing_features <- setdiff(pca_features, rownames(x[["SCT"]]@scale.data))
  if (length(missing_features) > 0) x <- GetResidual(object = x, features = missing_features, assay = "SCT", umi.assay = "RNA", verbose = FALSE)
  x
})

merged <- merge(pool_objs[[1]], y = pool_objs[-1], project = "CGG0M1_Patel", merge.data = TRUE)
rm(pool_objs)
gc()

DefaultAssay(merged) <- "SCT"
pca_features <- intersect(pca_features, rownames(merged[["SCT"]]@scale.data))
if (length(pca_features) < 2) stop("Too few merged SCT residual features for PCA: ", length(pca_features))

npcs_use <- min(50, length(pca_features) - 1, ncol(merged) - 1)
merged <- RunPCA(merged, features = pca_features, npcs = npcs_use, verbose = FALSE)

merged <- NormalizeData(merged, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = FALSE)

saveRDS(merged, file.path(rds_dir, paste0("03_merged_sct_pca", out_suffix, ".rds")))

filter_summary_all_path <- file.path(tables_dir, "01_filter_summary_all.tsv")
if (!file.exists(filter_summary_all_path)) stop("Missing combined filter summary: ", filter_summary_all_path)

rna_vln <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_id", pt.size = 0, ncol = 3)
ggsave(file.path(plots_dir, "merged_postfilter_vln_rna_qc_by_sample.png"), rna_vln, width = 14, height = 5)

top_features <- head(feature_table, 30)
vf_plot <- ggplot(top_features, aes(x = reorder(gene, score), y = score)) + geom_col() + coord_flip() + labs(x = "Gene", y = "Consensus SCT variable-feature score")
ggsave(file.path(plots_dir, "merged_sct_consensus_variable_features.png"), vf_plot, width = 8, height = 8)

ggsave(file.path(plots_dir, "merged_pca_elbowplot.png"), ElbowPlot(merged, ndims = npcs_use), width = 8, height = 5)
ggsave(file.path(plots_dir, "merged_pca_by_condition.png"), DimPlot(merged, reduction = "pca", group.by = "condition", raster = TRUE), width = 7, height = 5)
ggsave(file.path(plots_dir, "merged_pca_by_pool.png"), DimPlot(merged, reduction = "pca", group.by = "pool", raster = TRUE), width = 7, height = 5)

adt_clr <- VlnPlot(merged, assay = "ADT", features = rownames(merged[["ADT"]]), group.by = "condition", slot = "data", pt.size = 0, ncol = 3)
ggsave(file.path(plots_dir, "merged_vln_adt_clr_by_condition.png"), adt_clr, width = 12, height = 5)

message("Saved: ", file.path(rds_dir, paste0("03_merged_sct_pca", out_suffix, ".rds")))
message("Checked: ", filter_summary_all_path)
message("Saved: ", file.path(tables_dir, paste0("02_sct_consensus_variable_features", out_suffix, ".tsv")))
