#!/usr/bin/env Rscript

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

project_dir <- "/well/bsg/projects/CGG0M1_Patel"
manifest_path <- file.path(project_dir, "seurat_analysis/tables/sample_manifest.tsv")

rds_dir <- file.path(project_dir, "seurat_analysis/rds/00_import")
plot_dir <- file.path(project_dir, "seurat_analysis/plots/01.5_qc_conditions")
table_dir <- file.path(project_dir, "seurat_analysis/tables")

manifest <- read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

obj_list <- setNames(vector("list", nrow(manifest)), manifest$sample_id)
qc_list <- lapply(seq_len(nrow(manifest)), function(i) {
  sample_id <- manifest$sample_id[i]
  rds_path <- file.path(rds_dir, paste0("00_individual_qc_", sample_id, ".rds"))
  obj <- readRDS(rds_path)
  obj_list[[sample_id]] <<- obj

  data.frame(
    sample_id = sample_id,
    pool = manifest$pool[i],
    condition = manifest$condition[i],
    nFeature_RNA = obj$nFeature_RNA,
    nCount_RNA = obj$nCount_RNA,
    percent.mt = obj$percent.mt,
    stringsAsFactors = FALSE
  )
})

qc_df <- do.call(rbind, qc_list)
qc_df$condition <- factor(qc_df$condition, levels = c("NoTreated", "Control", "NR1", "CASPR2"))

condition_levels <- c("NoTreated", "Control", "NR1", "CASPR2")
pool_levels <- sort(unique(manifest$pool))
manifest$condition <- factor(manifest$condition, levels = condition_levels)
manifest$pool <- factor(manifest$pool, levels = pool_levels)

grid_condition_levels <- c("NR1", "CASPR2", "Control", "NoTreated")
grid_order <- order(factor(manifest$condition, levels = grid_condition_levels), manifest$pool)

cat("Loaded QC metrics for ", length(unique(qc_df$sample_id)), " samples across ", nlevels(qc_df$condition), " conditions\n", sep = "")

condition_summary <- do.call(rbind, lapply(split(qc_df, qc_df$condition), function(d) {
  data.frame(
    condition = d$condition[1],
    n_cells = nrow(d),
    median_nFeature_RNA = median(d$nFeature_RNA),
    median_nCount_RNA = median(d$nCount_RNA),
    median_percent_mt = median(d$percent.mt, na.rm = TRUE),
    q01_nFeature_RNA = as.numeric(quantile(d$nFeature_RNA, 0.01)),
    q99_nFeature_RNA = as.numeric(quantile(d$nFeature_RNA, 0.99)),
    q99_percent_mt = as.numeric(quantile(d$percent.mt, 0.99, na.rm = TRUE))
  )
}))

write.table(
  condition_summary,
  file = file.path(table_dir, "01.5_qc_summary_by_condition.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

VlnPlot_metric <- function(metric) {
  ggplot(qc_df, aes(x = condition, y = .data[[metric]], fill = condition)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_jitter(width = 0.15, size = 0.05, alpha = 0.1) +
    theme_bw() +
    NoLegend() +
    labs(x = NULL, y = metric)
}

p_vln <- VlnPlot_metric("nFeature_RNA") + VlnPlot_metric("nCount_RNA") + VlnPlot_metric("percent.mt")
ggsave(file.path(plot_dir, "01.5_vln_rna_qc_by_condition.png"), p_vln, width = 12, height = 5)

p_scatter1 <- ggplot(qc_df, aes(x = nCount_RNA, y = percent.mt, color = condition)) +
  geom_point(size = 0.05, alpha = 0.1) +
  theme_bw() +
  facet_wrap(~condition)

p_scatter2 <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA, color = condition)) +
  geom_point(size = 0.05, alpha = 0.1) +
  theme_bw() +
  facet_wrap(~condition)

ggsave(file.path(plot_dir, "01.5_scatter_count_vs_mt_by_condition.png"), p_scatter1, width = 10, height = 8)
ggsave(file.path(plot_dir, "01.5_scatter_count_vs_feature_by_condition.png"), p_scatter2, width = 10, height = 8)

p_pool <- ggplot(qc_df, aes(x = condition, y = nCount_RNA, fill = pool)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  labs(x = NULL, y = "nCount_RNA")
ggsave(file.path(plot_dir, "01.5_boxplot_ncount_by_condition_and_pool.png"), p_pool, width = 8, height = 5)

qc_plot_dir <- file.path(project_dir, "seurat_analysis/plots/01_seurat_qc")

n_pools <- length(pool_levels)
n_conditions <- length(condition_levels)

get_row_samples <- function(cond) manifest$sample_id[as.character(manifest$condition) == cond]

build_grid <- function(plot_fn, ncol, title, xlim_by_row = NULL, ylim_by_row = NULL, strip_x_identity = FALSE, strip_y_full = FALSE, strip_x_nonlast_row = FALSE) {
  n_rows <- ceiling(length(grid_order) / ncol)
  panels <- lapply(seq_along(grid_order), function(pos) {
    i <- grid_order[pos]
    sample_id <- manifest$sample_id[i]
    condition <- as.character(manifest$condition[i])
    is_first_col <- ((pos - 1) %% ncol) == 0
    is_last_row <- ((pos - 1) %/% ncol) == (n_rows - 1)

    p <- plot_fn(obj_list[[sample_id]], condition) +
      ggtitle(gsub("_", " ", sample_id)) +
      theme(
        plot.title = element_text(size = 9),
        plot.margin = margin(t = 12, b = 12, l = 5, r = 5)
      )

    if (strip_x_identity) {
      p <- p + xlab(NULL) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }

    if (strip_x_nonlast_row && !is_last_row) {
      p <- p + xlab(NULL) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }

    row_xlim <- if (!is.null(xlim_by_row)) xlim_by_row[[condition]] else NULL
    row_ylim <- if (!is.null(ylim_by_row)) ylim_by_row[[condition]] else NULL
    if (!is.null(row_xlim) || !is.null(row_ylim)) {
      p <- p + coord_cartesian(xlim = row_xlim, ylim = row_ylim)
    }
    if (!is_first_col) {
      p <- p + theme(axis.title.y = element_blank())
      if (!strip_x_identity) p <- p + theme(axis.title.x = element_blank())
      if (strip_y_full) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    p
  })
  wrap_plots(panels, ncol = ncol) + plot_annotation(title = title)
}

adt_vln_plot <- function(obj, feature) {
  df <- FetchData(obj, vars = feature, layer = "data")
  colnames(df) <- "value"
  ggplot(df, aes(x = "", y = value)) +
    geom_violin(fill = "#F8766D", trim = TRUE, scale = "width") +
    geom_jitter(width = 0.3, size = 0.02, alpha = 0.3) +
    theme_bw() +
    labs(x = NULL, y = "Expression Level")
}

row_range <- function(values) {
  setNames(
    lapply(grid_condition_levels, function(cond) range(values[qc_df$sample_id %in% get_row_samples(cond)])),
    grid_condition_levels
  )
}

row_rna_range <- row_range(qc_df$nFeature_RNA)
row_scatter_x_range <- row_range(qc_df$nCount_RNA)
row_scatter_y_range <- row_range(qc_df$percent.mt)

row_adt_range <- setNames(lapply(grid_condition_levels, function(cond) {
  ids <- get_row_samples(cond)
  vals <- unlist(lapply(ids, function(id) GetAssayData(obj_list[[id]], assay = "ADT", layer = "data")["CASPR2-mAb", ]))
  range(vals)
}), grid_condition_levels)

metric_labels <- c(
  nFeature_RNA = "Number of detected genes",
  nCount_RNA   = "RNA count",
  percent.mt   = "Percentage of mitochondrial reads"
)
vln_condition_order <- c("NR1", "CASPR2", "Control", "NoTreated")
metric_range <- setNames(lapply(names(metric_labels), function(m) range(qc_df[[m]], na.rm = TRUE)), names(metric_labels))

vln_panels <- list()
for (metric in names(metric_labels)) {
  for (pool in pool_levels) {
    df_panel <- qc_df[as.character(qc_df$pool) == as.character(pool), ]
    df_panel$condition <- factor(as.character(df_panel$condition), levels = vln_condition_order)

    p <- ggplot(df_panel, aes(x = condition, y = .data[[metric]], fill = condition)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_jitter(width = 0.15, size = 0.05, alpha = 0.15) +
      coord_cartesian(ylim = metric_range[[metric]]) +
      theme_bw() +
      NoLegend() +
      labs(x = NULL, y = if (pool == pool_levels[1]) metric_labels[[metric]] else NULL)

    if (metric != tail(names(metric_labels), 1)) {
      p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    if (metric == names(metric_labels)[1]) {
      p <- p + ggtitle(as.character(pool)) + theme(plot.title = element_text(size = 10, hjust = 0.5))
    }
    vln_panels[[length(vln_panels) + 1]] <- p
  }
}

p_vln_grid <- wrap_plots(vln_panels, ncol = n_pools)
ggsave(file.path(qc_plot_dir, "01_vln_rna_qc_grid.png"), p_vln_grid, width = 4 * n_pools, height = 4 * length(metric_labels))

p_vln_adt_grid <- build_grid(
  function(obj, condition) adt_vln_plot(obj, "CASPR2-mAb"),
  ncol = n_pools,
  title = "ADT CLR (CASPR2-mAb) by sample",
  ylim_by_row = row_adt_range,
  strip_x_identity = TRUE
)
ggsave(file.path(qc_plot_dir, "01_vln_adt_clr_grid.png"), p_vln_adt_grid, width = 4 * n_pools, height = 4 * n_conditions)

condition_colors <- setNames(scales::hue_pal()(length(vln_condition_order)), vln_condition_order)

p_scatter_grid <- build_grid(
  function(obj, condition) {
    p <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1) + NoLegend()
    point_idx <- which(sapply(p$layers, function(l) inherits(l$geom, "GeomPoint")))
    for (idx in point_idx) {
      p$layers[[idx]]$mapping <- NULL
      p$layers[[idx]]$aes_params$colour <- condition_colors[[condition]]
    }
    p
  },
  ncol = n_pools,
  title = NULL,
  xlim_by_row = row_scatter_x_range,
  ylim_by_row = row_scatter_y_range,
  strip_x_nonlast_row = TRUE
)
p_scatter_grid <- p_scatter_grid & theme(axis.title = element_blank())

y_axis_label <- wrap_elements(grid::textGrob("Percentage of mitochondrial reads per sample", rot = 90, gp = grid::gpar(fontsize = 12)))
x_axis_label <- wrap_elements(grid::textGrob("RNA count per sample", gp = grid::gpar(fontsize = 12)))

p_scatter_final <- y_axis_label + p_scatter_grid + x_axis_label +
  plot_layout(design = "AB\n#C", widths = c(0.03, 1), heights = c(1, 0.03))
ggsave(file.path(qc_plot_dir, "01_scatter_rna_qc_grid.png"), p_scatter_final, width = 4 * n_pools, height = 4 * n_conditions)

cat("Saved condition-level QC summary to ", file.path(table_dir, "01.5_qc_summary_by_condition.tsv"), "\n", sep = "")
cat("Saved plots to ", plot_dir, "\n", sep = "")
cat("Saved multi-panel QC grids to ", qc_plot_dir, "\n", sep = "")
