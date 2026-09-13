#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(scclusteval)
  library(tidyverse)
  library(parallel)
})

options(future.globals.maxSize = 8 * 1024^3)

## ---- CONFIG ------------------------------------------------------------
conditions <- c("NR1", "CASPR2", "Control", "NoTreated")
args <- commandArgs(trailingOnly = TRUE)
target_condition <- conditions[as.integer(args[1])]
if (is.na(target_condition)) stop("Usage: 05_bootstrap_cluster_stability_by_condition.R <1-4>, got array index ", args[1])

by_condition_dir <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition"
rds_dir     <- file.path(by_condition_dir, "rds_by_condition", target_condition)
input_rds   <- file.path(rds_dir, "04_clustered_sct.rds")

out_dir   <- file.path(rds_dir, "05_stability")
plot_dir  <- file.path(by_condition_dir, "plots_by_condition", target_condition, "05_stability")
table_dir <- file.path(by_condition_dir, "tables_by_condition", target_condition)

resolutions_to_test <- c(0.3, 0.4, 0.5, 0.6)  # resolution not yet chosen; comparing these

dims_use <- 1:10
k_param  <- 20
n_start  <- 10

n_boot         <- 30
subsample_rate <- 0.8

## Zumel & Mount (2014) rule of thumb cited in Tang et al. 2021:
## median/mean Jaccard < 0.6 unstable, 0.6-0.75 weak pattern, > 0.85 highly stable
jaccard_cutoff_median  <- 0.6
jaccard_cutoff_percent <- 0.8
percent_cutoff         <- 0.8

## Each worker reprocesses SCTransform+PCA+FindNeighbors on an 80% subsample
n_cores <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")), 4L)
seed    <- 20260707

dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Load + resolve cluster columns ------------------------------------
message("Condition: ", target_condition)
message("Loading: ", input_rds)
obj <- readRDS(input_rds)

stopifnot("RNA" %in% names(obj@assays))   # SCTransform is re-run from raw RNA counts each round
stopifnot("pool" %in% colnames(obj@meta.data))

vars_to_regress <- "percent.mt"
stopifnot(vars_to_regress %in% colnames(obj@meta.data))

DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = "RNA", layers = "counts")

## matches this project's column naming (see 06_cluster_annotation_by_condition.R::find_cluster_col())
res_label <- function(res) gsub("\\.", "_", formatC(res, format = "f", digits = 1))
cluster_col_for <- function(res) paste0("cluster_sct_res_", res_label(res))

ident_cols <- setNames(sapply(resolutions_to_test, cluster_col_for), resolutions_to_test)
missing <- ident_cols[!ident_cols %in% colnames(obj@meta.data)]
if (length(missing) > 0) stop("Missing cluster columns in object: ", paste(missing, collapse = ", "))

## ---- Bootstrap one round ------------------------------------------------
reprocess_subsample <- function(object, resolution) {
  pool_ids  <- unique(object$pool)
  pool_list <- lapply(pool_ids, function(p) object[, object$pool == p])
  names(pool_list) <- pool_ids

  pool_list <- lapply(pool_list, function(x) {
    DefaultAssay(x) <- "RNA"
    SCTransform(x, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = vars_to_regress,
                new.assay.name = "SCT", conserve.memory = TRUE, verbose = FALSE)
  })

  merged <- if (length(pool_list) == 1) pool_list[[1]] else merge(pool_list[[1]], y = pool_list[-1], merge.data = TRUE)
  DefaultAssay(merged) <- "SCT"
  merged <- FindVariableFeatures(merged, assay = "SCT", selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  pca_features <- VariableFeatures(merged)

  missing_features <- setdiff(pca_features, rownames(merged[["SCT"]]@scale.data))
  if (length(missing_features) > 0) merged <- GetResidual(merged, features = missing_features, assay = "SCT", umi.assay = "RNA", verbose = FALSE)
  pca_features <- intersect(pca_features, rownames(merged[["SCT"]]@scale.data))

  merged <- RunPCA(merged, features = pca_features, npcs = max(dims_use), verbose = FALSE)
  merged <- FindNeighbors(merged, dims = dims_use, k.param = k_param, reduction = "pca", verbose = FALSE)
  FindClusters(merged, resolution = resolution, n.start = n_start, verbose = FALSE)
}

run_one_bootstrap <- function(i, resolution) {
  ## same seed per round across resolutions -> identical subsampled cells,
  ## so resolutions are compared on the same subsamples
  tryCatch({
    sub_obj <- RandomSubsetData(obj, rate = subsample_rate, random.subset.seed = seed + i)
    original_ident <- Idents(sub_obj)
    sub_obj <- reprocess_subsample(sub_obj, resolution)
    list(round = i, original_ident = original_ident, recluster_ident = Idents(sub_obj))
  }, error = function(e) {
    message("round ", i, " (res ", resolution, ") failed: ", conditionMessage(e))
    NULL
  })
}

## ---- Loop over resolutions ----------------------------------------------
all_summaries <- list()

for (resolution in resolutions_to_test) {
  res_out   <- res_label(resolution)
  ident_col <- ident_cols[[as.character(resolution)]]
  message("\n== condition ", target_condition, ", resolution ", resolution, " (", ident_col, ") ==")

  Idents(obj) <- obj@meta.data[[ident_col]]
  original_ident_full <- Idents(obj)

  boot_results <- parallel::mclapply(seq_len(n_boot), function(i) run_one_bootstrap(i, resolution),
                                      mc.cores = n_cores)
  boot_results <- Filter(function(x) !is.null(x) && !inherits(x, "try-error"), boot_results)
  if (length(boot_results) < n_boot) {
    warning(res_out, ": only ", length(boot_results), "/", n_boot, " rounds succeeded")
  }
  saveRDS(boot_results, file.path(out_dir, sprintf("bootstrap_results_res_%s.rds", res_out)))

  if (length(boot_results) == 0) {
    warning(res_out, ": all bootstrap rounds failed (likely OOM -- check n_cores/--mem); skipping resolution")
    next
  }

  idents1 <- lapply(boot_results, `[[`, "original_ident")
  idents2 <- lapply(boot_results, `[[`, "recluster_ident")

  stability_median <- AssignStableCluster(idents1, idents2, method = "jaccard_median",
                                           jaccard_cutoff = jaccard_cutoff_median)
  ## paper notes the Jaccard distribution can be bimodal, where % > cutoff summarizes better than the median
  stability_percent <- AssignStableCluster(idents1, idents2, method = "jaccard_percent",
                                            jaccard_cutoff = jaccard_cutoff_percent, percent_cutoff = percent_cutoff)

  summary_tbl <- tibble(
    resolution              = resolution,
    cluster                 = names(stability_median$stable_cluster),
    n_cells                 = as.integer(table(original_ident_full)[names(stability_median$stable_cluster)]),
    jaccard_median          = apply(stability_median$jaccardIndex, 2, median),
    jaccard_mean            = apply(stability_median$jaccardIndex, 2, mean),
    stable_by_median_gt_0.6 = stability_median$stable_cluster,
    pct_rounds_gt_0.8       = stability_percent$stable_index,
    stable_by_percent_rule  = stability_percent$stable_cluster
  )
  write_tsv(summary_tbl, file.path(table_dir, sprintf("05_bootstrap_stability_summary_res_%s.tsv", res_out)))
  all_summaries[[res_out]] <- summary_tbl

  message(sprintf("res %s: %d/%d stable clusters (median rule); %.1f%% of cells in stable clusters", res_out,
                   stability_median$number_of_stable_cluster, length(stability_median$stable_cluster),
                   CalculatePercentCellInStable(original_ident_full, stability_median$stable_cluster) * 100))

  p_rain <- JaccardRainCloudPlot(idents1, idents2) +
    geom_hline(yintercept = c(0.6, 0.75, 0.85), linetype = 2) +
    ggtitle(sprintf("Cluster stability (Jaccard), %s, resolution %s, n_boot=%d, rate=%.2f",
                    target_condition, res_out, length(boot_results), subsample_rate)) +
    xlab(paste0("cluster id (", ident_col, ")"))
  ggsave(file.path(plot_dir, sprintf("05_jaccard_raincloud_res_%s.png", res_out)), p_rain, width = 10, height = 6, dpi = 300)
}

## ---- Cross-resolution comparison ----------------------
combined <- bind_rows(all_summaries)
if (nrow(combined) == 0) {
  message("\nNo resolution produced any successful bootstrap rounds -- nothing to compare. See per-round errors above.")
} else {
  write_tsv(combined, file.path(table_dir, "05_bootstrap_stability_summary_all_resolutions.tsv"))

  comparison_tbl <- combined %>%
    group_by(resolution) %>%
    summarise(n_clusters = n(),
              n_stable_median = sum(stable_by_median_gt_0.6),
              pct_cells_stable_median = sum(n_cells[stable_by_median_gt_0.6]) / sum(n_cells) * 100,
              .groups = "drop")
  write_tsv(comparison_tbl, file.path(table_dir, "05_bootstrap_stability_resolution_comparison.tsv"))

  p_compare <- comparison_tbl %>%
    pivot_longer(c(n_clusters, n_stable_median), names_to = "category", values_to = "n") %>%
    ggplot(aes(factor(resolution), n, color = category, group = category)) +
    geom_point(size = 2) + geom_line() +
    labs(x = "resolution", y = "number of clusters", color = NULL,
         title = paste("Total vs. stable clusters across resolutions,", target_condition)) +
    theme_classic()
  ggsave(file.path(plot_dir, "05_resolution_comparison.png"), p_compare, width = 7, height = 5, dpi = 300)
}

message("\nDone. Condition: ", target_condition, ". Outputs in ", out_dir, ", ", plot_dir, ", ", table_dir)
