## Cluster-vs-annotation companion figure (review: shared colour rule between
## the res-0.4 cluster UMAP and the final annotation UMAP, so a cluster's fate
## after merging/renaming can be traced by colour). Regenerates both individual
## panels with matching formatting and combines them as a tagged two-panel figure.

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(patchwork)
  library(dplyr); library(colorspace)
})

sa <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
cluster_plot_dir <- file.path(sa, "plots", "04_clustering")
annot_plot_dir <- file.path(sa, "plots", "07_annotation")
combined_plot_dir <- cluster_plot_dir

cluster_col <- "cluster_sct_res_0_4"

obj <- readRDS(file.path(sa, "rds", "07_annotated_sct.rds"))

umap_df <- FetchData(obj, vars = c("umap_1", "umap_2", cluster_col, "final_cell_type"))
names(umap_df)[3] <- "cluster"
umap_df$cluster <- factor(as.character(umap_df$cluster), levels = as.character(sort(as.numeric(unique(umap_df$cluster)))))
umap_df$final_cell_type <- factor(umap_df$final_cell_type)

## --- shared colour rule: one hue per final annotation, clusters that merge
## into it get evenly spaced lightness steps within that hue -----------------
clusters_by_label <- umap_df %>%
  distinct(cluster, final_cell_type) %>%
  arrange(final_cell_type, as.numeric(as.character(cluster)))

label_levels <- levels(umap_df$final_cell_type)

## default-ggplot-style palette: the two large, related radial-glia populations
## get salmon/rose shades (cluster 0 -- the single biggest cluster -- is
## already this salmon in the unannotated res-0.4 plot); the rest keep the
## same clean, moderately saturated hues as that original palette.
label_colors <- c(
  "Possible radial glia"        = "#F8766D",  # salmon/coral (matches cluster 0)
  "Radial glia"                 = "#B23A48",  # deeper rose -- same family, related population
  "Unresolved ribosomal C3"     = "#00BF7D",  # spring green
  "Maturing inhibitory"         = "#E64B35",  # vivid red
  "Non-telencephalic"           = "#7C6BC4",  # purple
  "Unresolved neuron C8"        = "#D9A404",  # gold/amber
  "Probable maturing inhibitory" = "#00B0F6", # sky blue
  "Maturing excitatory"         = "#FF9F1C",  # warm orange
  "Possible immature astrocyte" = "#5E7C4A",  # olive green
  "Cycling progenitor"          = "#00BFC4"   # teal cyan
)
stopifnot(all(label_levels %in% names(label_colors)))
base_hues <- label_colors[label_levels]

cluster_cols <- unlist(lapply(label_levels, function(lab) {
  ids <- as.character(clusters_by_label$cluster[clusters_by_label$final_cell_type == lab])
  n <- length(ids)
  shade <- if (n == 1) 0 else seq(-0.25, 0.3, length.out = n)
  setNames(lighten(base_hues[[lab]], shade, method = "relative"), ids)
}))
cluster_cols <- cluster_cols[levels(umap_df$cluster)]

## --- panel A: clusters, coloured by shared rule, boxed number labels -------
n_per_cluster <- umap_df %>% count(cluster, name = "n")
umap_df <- left_join(umap_df, n_per_cluster, by = "cluster")
cluster_label_levels <- n_per_cluster %>%
  arrange(as.numeric(as.character(cluster))) %>%
  mutate(lab = paste0(cluster, "  (n=", format(n, big.mark = ","), ")")) %>%
  pull(lab)
umap_df$cluster_label <- factor(paste0(umap_df$cluster, "  (n=", format(umap_df$n, big.mark = ","), ")"), levels = cluster_label_levels)
names(cluster_cols) <- cluster_label_levels[match(names(cluster_cols), sub("  \\(n=.*\\)", "", cluster_label_levels))]

centroids_cluster <- umap_df %>% group_by(cluster) %>% summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

p_clusters <- ggplot(umap_df, aes(umap_1, umap_2, color = cluster_label)) +
  geom_point(size = 0.3, alpha = 0.6) +
  geom_label_repel(
    data = centroids_cluster, mapping = aes(x = umap_1, y = umap_2, label = cluster),
    inherit.aes = FALSE, color = "black", fill = "white", alpha = 0.85,
    label.size = 0.3, fontface = "bold", size = 4,
    segment.color = "grey30", max.overlaps = Inf, seed = 1234
  ) +
  scale_color_manual(values = cluster_cols) +
  guides(color = guide_legend(title = "Cluster (n cells)", override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 8))

ggsave(file.path(cluster_plot_dir, "04_umap_clusters_res_0_4_annotated.png"),
       p_clusters, width = 10, height = 6, dpi = 300, bg = "white")

## --- panel B: final annotation, same boxed-label style, no title -----------
n_per_label <- umap_df %>% count(final_cell_type, name = "n")
label_label_levels <- n_per_label %>%
  arrange(desc(n)) %>%
  mutate(lab = paste0(final_cell_type, "  (n=", format(n, big.mark = ","), ")")) %>%
  pull(lab)
umap_df <- left_join(umap_df, n_per_label, by = "final_cell_type")
umap_df$label_label <- factor(paste0(umap_df$final_cell_type, "  (n=", format(umap_df$n.y, big.mark = ","), ")"), levels = label_label_levels)
label_cols <- base_hues
names(label_cols) <- label_label_levels[match(names(label_cols), sub("  \\(n=.*\\)", "", label_label_levels))]

centroids_label <- umap_df %>% group_by(final_cell_type) %>% summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

p_annotation <- ggplot(umap_df, aes(umap_1, umap_2, color = label_label)) +
  geom_point(size = 0.3, alpha = 0.6) +
  geom_label_repel(
    data = centroids_label, mapping = aes(x = umap_1, y = umap_2, label = final_cell_type),
    inherit.aes = FALSE, color = "black", fill = "white", alpha = 0.85,
    label.size = 0.3, fontface = "bold", size = 3.6,
    segment.color = "grey30", max.overlaps = Inf, seed = 1234
  ) +
  scale_color_manual(values = label_cols) +
  guides(color = guide_legend(title = "Annotation (n cells)", override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 8))

ggsave(file.path(annot_plot_dir, "07_annotated_umap.png"),
       p_annotation, width = 11, height = 6, dpi = 300, bg = "white")

## --- combined tagged panel figure ------------------------------------------
combined <- p_clusters + p_annotation + plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
ggsave(file.path(combined_plot_dir, "08_umap_clusters_vs_annotation.png"),
       combined, width = 20, height = 6.5, dpi = 300, bg = "white")

message("Written: 04_umap_clusters_res_0_4_annotated.png, 07_annotated_umap.png, 08_umap_clusters_vs_annotation.png")
