#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(Seurat); library(ggplot2)})

set.seed(1234)
options(future.globals.maxSize = 300 * 1024^3)

conditions <- c("NR1", "CASPR2", "Control", "NoTreated")
args <- commandArgs(trailingOnly = TRUE)
target_condition <- conditions[as.integer(args[1])]
if (is.na(target_condition)) stop("Usage: 04_cluster_umap_by_condition.R <1-4>, got array index ", args[1])

root <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition"
rds_dir <- file.path(root, "rds_by_condition", target_condition)
input_rds <- file.path(rds_dir, "03_merged_sct_pca.rds")
output_rds <- file.path(rds_dir, "04_clustered_sct.rds")
plot_dir <- file.path(root, "plots_by_condition", target_condition, "04_clustering")
table_dir <- file.path(root, "tables_by_condition", target_condition)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

dims_use <- 1:10
resolutions <- c(0.3, 0.4, 0.5, 0.6, 1.0)
res_labels <- formatC(resolutions, format = "f", digits = 1)
pipeline <- "SCT"

obj <- readRDS(input_rds)
DefaultAssay(obj) <- pipeline

obj <- FindNeighbors(obj, reduction = "pca", dims = dims_use, verbose = FALSE)

cluster_cols <- setNames(
  paste0("cluster_", tolower(pipeline), "_res_", gsub("\\.", "_", res_labels)),
  res_labels
)

for (i in seq_along(resolutions)) {
  obj <- FindClusters(obj, resolution = resolutions[i], verbose = FALSE)
  obj[[cluster_cols[[i]]]] <- as.character(Idents(obj))
}

obj <- RunUMAP(obj, reduction = "pca", dims = dims_use, verbose = FALSE)
Idents(obj) <- cluster_cols[["0.3"]]
saveRDS(obj, output_rds)

save_plot <- function(p, name, width = 8, height = 6) ggsave(filename = file.path(plot_dir, paste0(name, ".png")), plot = p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)

for (r in names(cluster_cols)) {
  p <- DimPlot(obj, reduction = "umap", group.by = cluster_cols[[r]], label = TRUE, raster = TRUE) +
  NoLegend() +
  ggtitle(paste0("UMAP clusters, ", target_condition, ", resolution ", r))
  save_plot(p, paste0("04_umap_clusters_res_", gsub("\\.", "_", r)))
}

for (g in c("pool", "sample_id")) {
  p <- DimPlot(obj, reduction = "umap", group.by = g, raster = TRUE) +
    ggtitle(paste0("UMAP by ", g, " (", target_condition, ")"))
  save_plot(p, paste0("04_umap_by_", g), width = if (g == "sample_id") 10 else 8, height = if (g == "sample_id") 7 else 6)
}

count_clusters <- function(col, res) {
  x <- as.data.frame(table(cluster = obj[[col]][, 1]))
  names(x) <- c("cluster", "n_cells")
  cbind(resolution = res, x)
}

count_by_group <- function(col, res, group) {
  x <- as.data.frame(table(cluster = obj[[col]][, 1], group = obj[[group]][, 1]))
  names(x) <- c("cluster", group, "n_cells")
  x <- x[x$n_cells > 0, ]
  cbind(resolution = res, x)
}

write.table(
  do.call(rbind, Map(count_clusters, cluster_cols, names(cluster_cols))),
  file.path(table_dir, "04_cluster_cell_counts_by_resolution.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

for (g in c("pool", "sample_id")) {
  out_name <- paste0("04_cluster_", if (g == "sample_id") "sample" else g, "_counts.tsv")
  write.table(
    do.call(rbind, Map(function(col, res) count_by_group(col, res, g), cluster_cols, names(cluster_cols))),
    file.path(table_dir, out_name),
    sep = "\t", quote = FALSE,
    row.names = FALSE
  )
}

capture.output(sessionInfo(), file = file.path(table_dir, "04_cluster_umap_session_info.txt"))
