#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

project_root <- "/well/bsg/projects/CGG0M1_Patel"
analysis_dir <- file.path(project_root, "seurat_analysis")
manifest_path <- file.path(analysis_dir, "tables", "sample_manifest.tsv")
rds_dir <- file.path(analysis_dir, "rds")
tables_dir <- file.path(analysis_dir, "tables")
plots_dir <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis/plots/02_subset"

manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
sample_info <- manifest[as.integer(commandArgs(trailingOnly = TRUE)[1]), ]
sample_id <- sample_info$sample_id

rds_in <- file.path(rds_dir, "00_import", paste0("00_individual_qc_", sample_id, ".rds"))
rds_out_dir <- file.path(rds_dir, "02_filtered")
dir.create(rds_out_dir, recursive = TRUE, showWarnings = FALSE)
rds_out <- file.path(rds_out_dir, paste0("02_filtered_", sample_id, ".rds"))

obj <- readRDS(rds_in)
obj$sample_id <- sample_id; obj$pool <- sample_info$pool; obj$condition <- sample_info$condition

n_before <- ncol(obj)

obj_filt <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < quantile(nFeature_RNA, 0.90, na.rm = TRUE) & percent.mt < 10)

n_after <- ncol(obj_filt)

saveRDS(obj_filt, rds_out)

summary_df <- data.frame(
  sample_id = sample_id,
  pool = sample_info$pool,
  condition = sample_info$condition,
  n_cells_before = n_before,
  n_cells_after = n_after,
  percent_retained = round(100 * n_after / n_before, 2),
  median_nFeature_RNA_before = median(obj$nFeature_RNA, na.rm = TRUE),
  median_nFeature_RNA_after = median(obj_filt$nFeature_RNA, na.rm = TRUE),
  median_nCount_RNA_before = median(obj$nCount_RNA, na.rm = TRUE),
  median_nCount_RNA_after = median(obj_filt$nCount_RNA, na.rm = TRUE),
  median_percent_mt_before = median(obj$percent.mt, na.rm = TRUE),
  median_percent_mt_after = median(obj_filt$percent.mt, na.rm = TRUE),
  rds_path = rds_out
)

summary_path <- file.path(tables_dir, paste0("01_filter_summary", sample_id, ".tsv"))
write.table(summary_df, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

rna_vln <- VlnPlot(obj_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_id", pt.size = 0.03, ncol = 3)
ggsave(file.path(plots_dir, paste0(sample_id, "_postfilter_vln_rna_qc.png")), rna_vln, width = 12, height = 4)

message("Finished: ", sample_id)
message("Cells retained: ", n_after, " / ", n_before)
message("Saved: ", rds_out)
message("Summary: ", summary_path)