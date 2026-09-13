suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(tidyverse)
  library(UCell); library(clusterProfiler); library(pheatmap)
})

options(future.globals.maxSize = 8 * 1024^3) 
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

root <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
input_rds <- file.path(root, "rds", "04_clustered_sct.rds")
output_rds <- file.path(root, "rds", "06_annotated_sct.rds")
plot_dir <- file.path(root, "plots", "06_marker_evidence")
table_dir <- file.path(root, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

hnoca_ref_rds <- file.path(root, "rds", "ref_hnoca.rds")
siletti_ref_rds <- file.path(root, "rds", "ref_siletti.rds")
braun_ref_rds <- file.path(root, "rds", "ref_braun.rds")
atlas_max_ref_cells_per_label <- 500
## annot_level_2 = HNOCA's harmonized label; `cell_type` is per-dataset, unrecoded.
hnoca_label_col <- "annot_level_2"
## CellClass = Braun's own native classification (12 categories), cleaner than the
## CxG-standardized cell_type field.
braun_label_col <- "CellClass"

if (!file.exists(hnoca_ref_rds))   stop("Missing ", hnoca_ref_rds)
if (!file.exists(siletti_ref_rds)) stop("Missing ", siletti_ref_rds)
if (!file.exists(braun_ref_rds))   stop("Missing ", braun_ref_rds)

## Global (all-condition) clustering annotated at this resolution only; by-condition
## resolution is chosen separately (see by_condition/ scripts).
final_resolution <- 0.4

## 06a signatures: hand-curated marker gene sets, sourced from HNOCA
## (Fig1_HNOCA_establishment/manual_annotation.ipynb marker_genes_dict), Siletti
## (linnarsson-lab/auto-annotation-ah, Human_adult/Class/*.md), and Braun (Table S2
## per-cluster PoolEnriched genes, linnarsson-lab/developing-human-brain) -- the last
## added off_target_immune/off_target_erythrocyte and reinforced several existing sets.
## Source of truth is tables/06b_marker_gene_signatures.tsv -- edit genes there, not here.
signature_table_file <- file.path(table_dir, "06b_marker_gene_signatures.tsv")
if (!file.exists(signature_table_file)) stop("Missing ", signature_table_file)
signature_table <- read_tsv(signature_table_file, show_col_types = FALSE)
signatures <- setNames(
  lapply(strsplit(signature_table$genes, ",\\s*"), trimws),
  signature_table$id
)
signature_labels <- setNames(str_wrap(signature_table$cell_type, width = 14), signature_table$id
)

## CellHint ran per-condition (by_condition/cellhint/02_cellhint_harmonize.py); join by
## barcode, stripping the "-<condition>" suffix ad.concat(index_unique="-") appended.
cellhint_file <- file.path(root, "by_condition", "cellhint", "harmonisation_output", "02_harmonisation_reannotation.csv")
if (!file.exists(cellhint_file)) stop("Missing ", cellhint_file)
cellhint_reannotation <- read.csv(cellhint_file, stringsAsFactors = FALSE) %>%
  mutate(barcode = sub("-[^-]+$", "", ID))

## CellTypist ran per-condition (by_condition/cellhint/06_celltypist_predict.py); unlike CellHint's
## output, barcodes were never index_unique-suffixed, so no stripping needed before joining.
celltypist_file <- file.path(root, "by_condition", "cellhint", "celltypist_output", "06_celltypist_predictions.csv")
if (!file.exists(celltypist_file)) stop("Missing ", celltypist_file)
celltypist_predictions <- read.csv(celltypist_file, stringsAsFactors = FALSE)

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "SCT"
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)
obj$cellhint_group <- setNames(cellhint_reannotation$group, cellhint_reannotation$barcode)[colnames(obj)]
obj$celltypist_hnoca_label <- setNames(celltypist_predictions$celltypist_hnoca_label, celltypist_predictions$barcode)[colnames(obj)]
obj$celltypist_braun_label <- setNames(celltypist_predictions$celltypist_braun_label, celltypist_predictions$barcode)[colnames(obj)]

## expressed_universe below is the enricher() background (06c); min_signature_size drops
## sets under 4 genes after the detected-gene filter -- below that, UCell/enrichment scores
## get noisy (single-cell dropout on 1-3 genes isn't averaged out). signatures_detected keeps
## everything with >=1 detected gene, for the DotPlot only: automated composite scoring
## shouldn't trust a 1-gene signature, but a human eyeballing raw expression still can.
signatures_sourced <- signatures  # pre-filter, for signature_gene_table below
gene_detected_cells <- Matrix::rowSums(GetAssayData(obj, assay = "SCT", layer = "counts") > 0)
signature_sizes_before <- lengths(signatures)
signatures_detected <- lapply(signatures, function(genes) genes[genes %in% names(gene_detected_cells)[gene_detected_cells >= 10]])
signatures_detected <- signatures_detected[lengths(signatures_detected) > 0]
min_signature_size <- 4
signatures <- signatures_detected[lengths(signatures_detected) >= min_signature_size]
message("Signature sizes before/after filtering (genes detected in >=10 cells):")
print(data.frame(signature = names(signature_sizes_before), n_before = signature_sizes_before,
                  n_after = lengths(signatures)[names(signature_sizes_before)]))
expressed_universe <- names(gene_detected_cells)[gene_detected_cells >= 10]
write.table(data.frame(gene = expressed_universe), file.path(plot_dir, "06_expressed_gene_universe.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
term2gene <- stack(signatures)[, c(2, 1)]
colnames(term2gene) <- c("term", "gene")

## Every sourced marker gene per cell type (06a), and whether it made it into the
## UCell/enrichment signatures actually used (detected in data + set above size threshold).
signature_gene_table <- bind_rows(lapply(names(signatures_sourced), function(ct) {
  genes <- signatures_sourced[[ct]]
  tibble(cell_type = ct, gene = genes,
         detected_in_data = genes %in% names(gene_detected_cells)[gene_detected_cells >= 10],
         signature_used = ct %in% names(signatures))
}))
write_tsv(signature_gene_table, file.path(plot_dir, "06a_signature_genes.tsv"))

res_label <- function(res) formatC(res, format = "f", digits = 1)
res_file  <- function(res) gsub("\\.", "_", res_label(res))

sort_clusters <- function(x) {
  x <- unique(as.character(x))
  if (all(grepl("^[0-9]+$", x))) x[order(as.integer(x))] else sort(x)
}

find_cluster_col <- function(res) {
  candidates <- c(
    paste0("cluster_sct_res_", res_file(res)),
    paste0("SCT_snn_res.", res_label(res)),
    paste0("SCT_snn_res.", sub("\\.0$", "", res_label(res)))
  )
  hit <- candidates[candidates %in% colnames(obj@meta.data)]
  if (length(hit) == 0) stop("No cluster column found for resolution ", res_label(res))
  hit[1]
}

save_plot <- function(p, out_dir, name, width = 14, height = 7) {
  ggsave(filename = file.path(out_dir, paste0("06_", name, ".png")), plot = p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
}

heatmap_png <- function(mat, title, file, scale = "none") {
  pheatmap(mat, cluster_rows = nrow(mat) >= 2, cluster_cols = ncol(mat) >= 2, scale = scale, main = title, filename = file)
}

## DotPlot(features = <named list>) errors on a gene appearing under two group labels
## ("factor level ... is duplicated"); drop it from later signatures, keeping it under
## the first one listed. UCell/enrichment don't need this -- they handle shared genes fine.
dedupe_across_signatures <- function(signature_list) {
  seen <- character(0)
  lapply(signature_list, function(genes) {
    keep <- genes[!genes %in% seen]
    seen <<- c(seen, keep)
    keep
  })
}

top_markers <- function(markers, n = 10) {
  empty <- data.frame(cluster = character(0), gene = character(0))
  if (nrow(markers) == 0) return(empty)
  fc_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
  markers <- markers[markers$p_val_adj < 0.05 & markers[[fc_col]] > 0, ]
  if (nrow(markers) == 0) return(empty)
  do.call(rbind, lapply(split(markers, markers$cluster), function(x) head(x[order(x$p_val_adj, -x[[fc_col]]), ], n)))
}

## Dominant atlas label within a cluster (top label + its fraction of cells)
dominant_label <- function(idx, labels, scores = NULL) {
  tab <- sort(table(labels[idx]), decreasing = TRUE)
  top <- names(tab)[1]
  fraction <- unname(tab[1]) / sum(idx)
  median_score <- if (is.null(scores)) NA_real_ else median(scores[idx & labels == top], na.rm = TRUE)
  list(top = top, fraction = fraction, median_score = median_score)
}

## CellHint/CellTypist group x cluster cross-tab
crosstab_top <- function(clusters, labels, prefix, filename) {
  tab <- table(cluster = clusters, group = labels, useNA = "ifany")
  write.table(as.data.frame.matrix(tab), filename, sep = "\t", quote = FALSE)
  as.data.frame(tab) %>%
    group_by(cluster) %>%
    mutate(fraction = Freq / sum(Freq)) %>%
    slice_max(Freq, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(cluster = as.character(cluster),
              !!paste0(prefix, "_top_group") := as.character(group),
              !!paste0(prefix, "_fraction") := fraction)
}

## Atlas label transfer (Seurat anchors)
transfer_atlas_labels <- function(ref_rds, label_col, prefix, query) {
  ref <- readRDS(ref_rds)
  if (ncol(ref) > atlas_max_ref_cells_per_label * length(unique(ref[[label_col]]))) {
    keep <- unlist(lapply(split(colnames(ref), ref[[label_col]]), function(cells) {
      sample(cells, min(length(cells), atlas_max_ref_cells_per_label))
    }))
    ref <- subset(ref, cells = keep)
  }
  ref <- SCTransform(ref, verbose = FALSE)
  ref <- RunPCA(ref, npcs = 30, verbose = FALSE)
  ref <- RunUMAP(ref, dims = 1:30, return.model = TRUE, verbose = FALSE)

  anchors <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", dims = 1:30, reference.reduction = "pca")
  mapped <- MapQuery(anchorset = anchors, query = query, reference = ref,
                      refdata = setNames(list(label_col), prefix), reference.reduction = "pca")
  mapped[[paste0(prefix, "_label")]] <- mapped[[paste0("predicted.", prefix)]]
  mapped[[paste0(prefix, "_score")]] <- mapped[[paste0("predicted.", prefix, ".score")]]
  mapped
}

obj <- transfer_atlas_labels(hnoca_ref_rds, hnoca_label_col, "hnoca", obj)
obj <- transfer_atlas_labels(braun_ref_rds, braun_label_col, "braun", obj)
gc()

## Siletti 2023 broad class (SingleR)
siletti_ref <- readRDS(siletti_ref_rds)
label_col <- intersect(c("supercluster_term", "cell_type", "class"), colnames(SummarizedExperiment::colData(siletti_ref)))[1]
if (is.na(label_col)) stop("Siletti reference missing cell_type/class colData")
pred <- SingleR::SingleR(test = GetAssayData(obj, assay = "SCT", layer = "data"), ref = siletti_ref, labels = siletti_ref[[label_col]],
                          BPPARAM = BiocParallel::MulticoreParam(n_cores))
obj$siletti_class <- pred$labels
rm(siletti_ref, pred); gc()

## 06b -- UCell signature scoring: per-cell, ~0-centered, resolution-independent
obj <- AddModuleScore_UCell(obj, features = signatures, assay = "SCT", name = "_UCell")
ucell_cols <- paste0(names(signatures), "_UCell")
ucell_plots <- lapply(ucell_cols, function(col) {
  FeaturePlot(obj, features = col, reduction = "umap", raster = TRUE, order = TRUE) +
    scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0) +
    ggtitle(sub("_UCell$", "", col))
})
n_ucell <- length(ucell_plots)
ncol_ucell <- ceiling(sqrt(n_ucell))
nrow_ucell <- ceiling(n_ucell / ncol_ucell)
p <- wrap_plots(ucell_plots, ncol = ncol_ucell, nrow = nrow_ucell)
save_plot(p, plot_dir, "ucell_signatures_grid", width = 7 * ncol_ucell, height = 6 * nrow_ucell)

for (res in final_resolution) {
  res_out <- res_file(res)
  res_dir <- file.path(plot_dir, paste0("res_", res_out))
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

  cluster_col <- find_cluster_col(res)
  clusters <- as.character(obj@meta.data[[cluster_col]])
  Idents(obj) <- factor(clusters, levels = sort_clusters(clusters))

  p <- DimPlot(obj, reduction = "umap", group.by = cluster_col, label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() + ggtitle(paste("UMAP clusters, resolution", res_label(res)))
  save_plot(p, res_dir, paste0("umap_clusters_res_", res_out), width = 8, height = 6)

  ## Marker dotplot, genes grouped by signature.
  dotplot_features <- dedupe_across_signatures(signatures_detected)
  names(dotplot_features) <- signature_labels[names(dotplot_features)]
  p <- DotPlot(obj, features = dotplot_features, group.by = cluster_col, assay = "SCT") + RotatedAxis() +
    ggtitle(paste("Cell-type marker DotPlot, resolution", res_label(res))) +
    theme(axis.text.x = element_text(size = 7), strip.text = element_text(size = 8))
  save_plot(p, res_dir, paste0("celltype_marker_dotplot_res_", res_out), width = 20, height = 8)

  markers <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE, test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, max.cells.per.ident = 2000, random.seed = 1234, verbose = TRUE)
  if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)

  write.table(markers, file.path(res_dir, paste0("06_findallmarkers_res_", res_out, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  ## Significant, positive markers (p_val_adj < 0.05, avg_log2FC > 0) drive top10 and top50 labels
  top50 <- top_markers(markers, n = 50)
  sig_pos <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0, ]

  cluster_levels <- levels(Idents(obj))
  top10 <- vapply(cluster_levels, function(cl) {
    genes <- top50$gene[top50$cluster == cl][1:10]
    paste(genes[!is.na(genes)], collapse = ", ")
  }, character(1))
  marker_evidence_top50 <- vapply(cluster_levels, function(cl) {
    paste(top50$gene[top50$cluster == cl], collapse = ", ")
  }, character(1))
  marker_evidence_genes_available <- vapply(cluster_levels, function(cl) {
    sum(sig_pos$cluster == cl)
  }, integer(1))

  agreement <- do.call(rbind, lapply(cluster_levels, function(cl) {
    idx <- clusters == cl
    h <- dominant_label(idx, obj$hnoca_label, obj$hnoca_score)
    b <- dominant_label(idx, obj$braun_label, obj$braun_score)
    s <- dominant_label(idx, obj$siletti_class)
    data.frame(
      cluster = cl,
      hnoca_top_label = h$top, hnoca_fraction = h$fraction, hnoca_median_score = h$median_score,
      braun_top_label = b$top, braun_fraction = b$fraction, braun_median_score = b$median_score,
      siletti_top_label = s$top, siletti_fraction = s$fraction
    )
  }))
  reconciliation <- agreement %>%
    mutate(
      top10 = unname(top10[as.character(cluster)]),
      marker_evidence_genes_available = unname(marker_evidence_genes_available[as.character(cluster)]),
      marker_evidence_top50 = unname(marker_evidence_top50[as.character(cluster)]),
      final_cell_type = "",
      notes = ""
    )

  ## 06b (cont.) -- cluster x signature mean UCell score; top signature + margin to 2nd
  ucell_mat <- sapply(ucell_cols, function(col) tapply(obj@meta.data[[col]], clusters, mean))
  colnames(ucell_mat) <- sub("_UCell$", "", ucell_cols)
  ucell_mat <- ucell_mat[sort_clusters(rownames(ucell_mat)), , drop = FALSE]
  write.table(as.data.frame(ucell_mat) %>% rownames_to_column("cluster"),
              file.path(res_dir, paste0("06b_ucell_cluster_means_res_", res_out, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  ucell_top <- data.frame(
    cluster = rownames(ucell_mat),
    ucell_top = colnames(ucell_mat)[apply(ucell_mat, 1, which.max)],
    ucell_runnerup = colnames(ucell_mat)[apply(ucell_mat, 1, function(x) order(x, decreasing = TRUE)[2])],
    ucell_margin = apply(ucell_mat, 1, function(x) sort(x, decreasing = TRUE)[1] - sort(x, decreasing = TRUE)[2])
  )

  ## 06c -- over-representation of each cluster's top-200 markers against `signatures`
  top200 <- top_markers(markers, n = 200)
  overrep <- compareCluster(gene ~ cluster, data = top200[, c("cluster", "gene")], fun = "enricher",
                             TERM2GENE = term2gene, universe = expressed_universe,
                             pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                             minGSSize = 1, maxGSSize = 500)
  overrep_df <- as.data.frame(overrep)
  write.table(overrep_df, file.path(res_dir, paste0("06c_overrep_res_", res_out, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  if (nrow(overrep_df) > 0) {
    overrep_top <- overrep_df %>% group_by(Cluster) %>% slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
      ungroup() %>% transmute(cluster = as.character(Cluster), overrep_top = ID, overrep_padj = p.adjust)
    save_plot(dotplot(overrep, x = "Cluster") + ggtitle(paste("Over-representation, resolution", res_label(res))),
              res_dir, paste0("overrep_dotplot_res_", res_out), width = 10, height = 8)
  } else {
    message("No significant over-representation for any cluster at resolution ", res_label(res), " -- skipping overrep dotplot")
    overrep_top <- tibble(cluster = character(0), overrep_top = character(0), overrep_padj = numeric(0))
  }

  ## 06d -- CellHint group x cluster cross-tab
  cellhint_top <- crosstab_top(clusters, obj$cellhint_group, "cellhint",
                                file.path(res_dir, paste0("06d_cellhint_crosstab_res_", res_out, ".tsv")))

  ## 06f -- CellTypist (HNOCA- and Braun-trained) group x cluster cross-tabs
  celltypist_hnoca_top <- crosstab_top(clusters, obj$celltypist_hnoca_label, "celltypist_hnoca",
                                        file.path(res_dir, paste0("06f_celltypist_hnoca_crosstab_res_", res_out, ".tsv")))
  celltypist_braun_top <- crosstab_top(clusters, obj$celltypist_braun_label, "celltypist_braun",
                                        file.path(res_dir, paste0("06f_celltypist_braun_crosstab_res_", res_out, ".tsv")))

  ## 06e -- hierarchical cluster x signature/marker heatmaps (visualization aid, not a new result)
  heatmap_png(ucell_mat, paste("Cluster x signature UCell score, resolution", res_label(res)),
              file.path(res_dir, paste0("06e_heatmap_signatures_res_", res_out, ".png")))

  top_marker_genes <- unique(na.omit(unlist(lapply(cluster_levels, function(cl) head(top50$gene[top50$cluster == cl], 5)))))
  expr <- GetAssayData(obj, assay = "SCT", layer = "data")[top_marker_genes, , drop = FALSE]
  marker_mat <- t(sapply(cluster_levels, function(cl) Matrix::rowMeans(expr[, clusters == cl, drop = FALSE])))
  heatmap_png(marker_mat, paste("Cluster x top-marker expression, resolution", res_label(res)),
              file.path(res_dir, paste0("06e_heatmap_markers_res_", res_out, ".png")), scale = "column")

  ## Jaccard cluster stability: 100-round bootstrap
  boot_100_file <- file.path(table_dir, paste0("05_bootstrap_stability_summary_res_", res_out, "_100.tsv"))
  boot_all_file <- file.path(table_dir, "05_bootstrap_stability_summary_all_resolutions.tsv")
  jaccard_stability <- if (file.exists(boot_100_file)) {
    read_tsv(boot_100_file, show_col_types = FALSE) %>% mutate(cluster = as.character(cluster)) %>%
      select(cluster, jaccard_median, stable_by_median_gt_0.6)
  } else if (file.exists(boot_all_file)) {
    read_tsv(boot_all_file, show_col_types = FALSE) %>% filter(resolution == res) %>%
      mutate(cluster = as.character(cluster)) %>% select(cluster, jaccard_median, stable_by_median_gt_0.6)
  } else {
    tibble(cluster = character(0), jaccard_median = numeric(0), stable_by_median_gt_0.6 = logical(0))
  }
  ## Zumel & Mount (2014) / Tang et al. 2021 tiers, matching the 05 raincloud plot hlines.
  jaccard_stability <- jaccard_stability %>%
    mutate(jaccard_stability_tier = cut(jaccard_median, breaks = c(-Inf, 0.6, 0.75, 0.85, Inf),
                                         labels = c("unstable", "weak", "stable", "highly stable")))

 
  reconciliation <- reconciliation %>%
    left_join(ucell_top, by = "cluster") %>%
    left_join(overrep_top, by = "cluster") %>%
    left_join(cellhint_top, by = "cluster") %>%
    left_join(celltypist_hnoca_top, by = "cluster") %>%
    left_join(celltypist_braun_top, by = "cluster") %>%
    left_join(jaccard_stability, by = "cluster") %>%
    mutate(methods_agree = ucell_top == overrep_top) %>%
    select(
      cluster,
      marker_evidence_genes_available,
      jaccard_median, stable_by_median_gt_0.6, jaccard_stability_tier,
      marker_evidence_top50,
      overrep_top, overrep_padj,
      hnoca_top_label, hnoca_fraction, hnoca_median_score,
      braun_top_label, braun_fraction, braun_median_score,
      siletti_top_label, siletti_fraction,
      ucell_top, ucell_runnerup, ucell_margin,
      cellhint_top_group, cellhint_fraction,
      celltypist_hnoca_top_group, celltypist_hnoca_fraction,
      celltypist_braun_top_group, celltypist_braun_fraction,
      top10, methods_agree,
      final_cell_type,
      notes
    )

  ## Refreshes evidence columns every run; preserves hand-filled final_cell_type/notes.
  annotation_file <- file.path(res_dir, paste0("06_marker_evidence_template_res_", res_out, ".tsv"))
  if (file.exists(annotation_file)) {
    prior <- read.delim(annotation_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
      transmute(cluster = as.character(cluster), final_cell_type, notes)
    reconciliation <- reconciliation %>% select(-final_cell_type, -notes) %>% left_join(prior, by = "cluster")
    reconciliation$final_cell_type[is.na(reconciliation$final_cell_type)] <- ""
    reconciliation$notes[is.na(reconciliation$notes)] <- ""
  }
  write.table(reconciliation, annotation_file, sep = "\t", quote = FALSE, row.names = FALSE)

  if (any(nzchar(trimws(reconciliation$final_cell_type)))) {
    map <- setNames(reconciliation$final_cell_type, reconciliation$cluster)
    cell_type <- unname(map[clusters])
    cell_type[is.na(cell_type) | cell_type == ""] <- paste0("Cluster ", clusters[is.na(cell_type) | cell_type == ""])
    celltype_col <- paste0("celltype_res_", res_out)
    obj[[celltype_col]] <- cell_type

    p <- DimPlot(obj, reduction = "umap", group.by = celltype_col, label = TRUE, repel = TRUE, raster = TRUE) + NoLegend() + ggtitle(paste("Annotated UMAP, resolution", res_label(res)))
    save_plot(p, res_dir, paste0("annotated_umap_res_", res_out), width = 9, height = 7)
  }
}

## ---- Atlas purity check at the final resolution ----

## Per-cluster dominant atlas label + purity. hnoca/siletti_label_shared_by_n_clusters > 1
## flags a merge candidate: two clusters both dominated by the same atlas label.
final_clusters <- as.character(obj@meta.data[[find_cluster_col(final_resolution)]])
purity_all <- bind_rows(lapply(sort_clusters(final_clusters), function(cl) {
  idx <- final_clusters == cl
  h <- dominant_label(idx, obj$hnoca_label)
  b <- dominant_label(idx, obj$braun_label)
  s <- dominant_label(idx, obj$siletti_class)
  tibble(cluster = cl, n_cells = sum(idx),
         hnoca_top_label = h$top, hnoca_purity = h$fraction,
         braun_top_label = b$top, braun_purity = b$fraction,
         siletti_top_label = s$top, siletti_purity = s$fraction)
})) %>%
  group_by(hnoca_top_label) %>% mutate(hnoca_label_shared_by_n_clusters = n()) %>% ungroup() %>%
  group_by(braun_top_label) %>% mutate(braun_label_shared_by_n_clusters = n()) %>% ungroup() %>%
  group_by(siletti_top_label) %>% mutate(siletti_label_shared_by_n_clusters = n()) %>% ungroup()
write_tsv(purity_all, file.path(table_dir, "06_atlas_purity.tsv"))
message("\n== Atlas purity, resolution ", final_resolution, " ==")
print(purity_all)

boot_file <- file.path(table_dir, "05_bootstrap_stability_summary_all_resolutions.tsv")
if (file.exists(boot_file)) {
  boot <- read_tsv(boot_file, show_col_types = FALSE) %>%
    filter(resolution == final_resolution) %>% mutate(cluster = as.character(cluster)) %>% select(-resolution)
  write_tsv(left_join(purity_all, boot, by = "cluster"), file.path(table_dir, "06_atlas_and_bootstrap_combined.tsv"))
}

label_umap_specs <- list(
  hnoca_label = "HNOCA transferred labels",
  braun_label = "Braun transferred labels",
  siletti_class = "Siletti transferred class",
  cellhint_group = "CellHint group",
  celltypist_hnoca_label = "CellTypist (HNOCA-trained) label",
  celltypist_braun_label = "CellTypist (Braun-trained) label"
)
label_umap_plots <- lapply(names(label_umap_specs), function(col) {
  DimPlot(obj, reduction = "umap", group.by = col, label = TRUE, repel = TRUE, raster = TRUE) +
    ggtitle(label_umap_specs[[col]])
})
save_plot(wrap_plots(label_umap_plots, ncol = 3, nrow = 2), plot_dir, "label_umaps_panel", width = 33, height = 14)

saveRDS(obj, output_rds)

capture.output(sessionInfo(), file = file.path(plot_dir, "06_marker_evidence_session_info.txt"))