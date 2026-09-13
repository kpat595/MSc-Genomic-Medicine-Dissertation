.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(hdf5r)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
sample_index <- as.integer(args[1])

project_dir <- "/well/bsg/projects/CGG0M1_Patel"
manifest_path <- file.path(project_dir, "seurat_analysis/tables/sample_manifest.tsv")

rds_dir <- file.path(project_dir, "seurat_analysis/rds/00_import")
plot_dir <- file.path(project_dir, "seurat_analysis/plots")
table_dir <- file.path(project_dir, "seurat_analysis/tables")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
sample_row <- manifest[sample_index, , drop = FALSE]
sample_id <- sample_row$sample_id

cat("Processing ", sample_id, "\n", sep = "")

read_cellbender_h5 <- function(h5_path) {
  h5 <- H5File$new(h5_path, mode = "r")
  on.exit(h5$close_all())

  mat <- h5[["matrix"]]
  features <- mat[["features"]]

  counts <- new(
    "dgCMatrix",
    x = as.numeric(mat[["data"]][]),
    i = as.integer(mat[["indices"]][]),
    p = as.integer(mat[["indptr"]][]),
    Dim = as.integer(mat[["shape"]][]),
    Dimnames = list(
      make.unique(as.character(features[["name"]][])),
      as.character(mat[["barcodes"]][])
    )
  )

  feature_type <- as.character(features[["feature_type"]][])

  list(
    RNA = counts[feature_type == "Gene Expression", , drop = FALSE],
    ADT = counts[feature_type == "Antibody Capture", , drop = FALSE]
  )
}

counts <- read_cellbender_h5(sample_row$h5_path)

adt_rename <- c(
  "ADT_Abcam_10X5002" = "NR1-mAb",
  "ADT_Abcam_10X5005" = "CASPR2-mAb",
  "ADT_Abcam_10X5006" = "A33-control-mAb",
  "Abcam_10X5002" = "NR1-mAb",
  "Abcam_10X5005" = "CASPR2-mAb",
  "Abcam_10X5006" = "A33-control-mAb"
)

adt_counts <- counts$ADT
new_adt_names <- rownames(adt_counts)
new_adt_names[new_adt_names %in% names(adt_rename)] <- adt_rename[new_adt_names[new_adt_names %in% names(adt_rename)]]
rownames(adt_counts) <- new_adt_names

obj <- CreateSeuratObject(counts = counts$RNA, assay = "RNA", project = sample_id, min.cells = 3, min.features = 200)

adt_counts <- adt_counts[, colnames(obj), drop = FALSE]
obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)

obj$sample_id <- sample_id
obj$pool <- sample_row$pool
obj$condition <- sample_row$condition
obj$condition_class <- sample_row$condition_class
obj$antibody_target <- sample_row$antibody_target
obj$expected_adt_feature_original <- sample_row$expected_adt_feature
obj$replicate_type <- sample_row$replicate_type
obj$h5_path <- sample_row$h5_path
obj$orig.ident <- sample_id

expected_dash <- c(
  "NR1_mAb" = "NR1-mAb",
  "CASPR2_mAb" = "CASPR2-mAb",
  "A33_control_mAb" = "A33-control-mAb",
  "None" = "None"
)

obj$expected_adt_feature <- ifelse(
  sample_row$expected_adt_feature %in% names(expected_dash),
  expected_dash[sample_row$expected_adt_feature],
  sample_row$expected_adt_feature
)

DefaultAssay(obj) <- "RNA"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")

adt_count_matrix <- GetAssayData(obj, assay = "ADT", layer = "counts")
obj$nCount_ADT <- Matrix::colSums(adt_count_matrix)
obj$nFeature_ADT <- Matrix::colSums(adt_count_matrix > 0)

obj <- NormalizeData(
  obj,
  assay = "ADT",
  normalization.method = "CLR",
  margin = 2,
  verbose = FALSE
)

DefaultAssay(obj) <- "RNA"
Idents(obj) <- obj$sample_id

rds_path <- file.path(rds_dir, paste0("00_individual_qc_", sample_id, ".rds"))
saveRDS(obj, rds_path)

qc_summary <- data.frame(
  sample_id = sample_id,
  pool = sample_row$pool,
  condition = sample_row$condition,
  n_cells = ncol(obj),
  n_rna_features = nrow(obj[["RNA"]]),
  n_adt_features = nrow(obj[["ADT"]]),
  median_nFeature_RNA = median(obj$nFeature_RNA),
  median_nCount_RNA = median(obj$nCount_RNA),
  median_percent_mt = median(obj$percent.mt, na.rm = TRUE),
  median_nFeature_ADT = median(obj$nFeature_ADT),
  median_nCount_ADT = median(obj$nCount_ADT),
  q01_nFeature_RNA = as.numeric(quantile(obj$nFeature_RNA, 0.01)),
  q99_nFeature_RNA = as.numeric(quantile(obj$nFeature_RNA, 0.99)),
  q99_percent_mt = as.numeric(quantile(obj$percent.mt, 0.99, na.rm = TRUE)),
  adt_features = paste(rownames(obj[["ADT"]]), collapse = ","),
  rds_path = rds_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.table(
  qc_summary,
  file = file.path(table_dir, paste0("qc_summary_", sample_id, ".tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (file.exists(sample_row$cellbender_metrics_path)) {
  cb <- read.csv(sample_row$cellbender_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
  cb <- cbind(data.frame(sample_id = sample_id, pool = sample_row$pool, condition = sample_row$condition), cb)

  write.table(
    cb,
    file = file.path(table_dir, paste0("cellbender_metrics_", sample_id, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.005)
ggsave(file.path(plot_dir, paste0(sample_id, "_seurat_vln_rna_qc.png")), p1, width = 12, height = 5)

p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
p_scatter <- p2 + p3 + plot_annotation(caption = sample_id) & theme(plot.caption = element_text(hjust = 0.5, size = 10))
ggsave(file.path(plot_dir, paste0(sample_id, "_seurat_scatter_rna_qc_combined.png")), p_scatter, width = 12, height = 5)

adt_features <- c("NR1-mAb", "CASPR2-mAb", "A33-control-mAb")
p4 <- VlnPlot(obj, features = adt_features, assay = "ADT", ncol = 3, pt.size = 0.05)
ggsave(file.path(plot_dir, paste0(sample_id, "_seurat_vln_adt_raw.png")), p4, width = 12, height = 5)

cat("Saved ", rds_path, "\n", sep = "")