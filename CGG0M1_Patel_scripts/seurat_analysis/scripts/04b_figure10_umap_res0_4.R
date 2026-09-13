## Figure 10 fix (review G6): legible cluster labels, legend with per-cluster n,
## no in-plot title (caption carries it), and a condition-coloured companion panel
## showing condition overlap in embedding space (supports R11).

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(patchwork); library(dplyr)
})

sa <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
plot_dir <- file.path(sa, "plots", "04_clustering")
cluster_col <- "cluster_sct_res_0_4"

obj <- readRDS(file.path(sa, "rds", "04_clustered_sct.rds"))

umap_df <- FetchData(obj, vars = c("umap_1", "umap_2", cluster_col, "condition"))
names(umap_df)[3] <- "cluster"
umap_df$cluster <- factor(as.character(umap_df$cluster), levels = as.character(sort(as.numeric(unique(umap_df$cluster)))))

n_per_cluster <- umap_df %>% count(cluster, name = "n")
label_levels <- n_per_cluster %>%
  arrange(as.numeric(as.character(cluster))) %>%
  mutate(lab = paste0(cluster, "  (n=", format(n, big.mark = ","), ")")) %>%
  pull(lab)
umap_df <- left_join(umap_df, n_per_cluster, by = "cluster")
umap_df$cluster_label <- factor(paste0(umap_df$cluster, "  (n=", format(umap_df$n, big.mark = ","), ")"), levels = label_levels)

centroids <- umap_df %>% group_by(cluster) %>% summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

p_clusters <- ggplot(umap_df, aes(umap_1, umap_2, color = cluster_label)) +
  geom_point(size = 0.3, alpha = 0.6) +
  geom_label_repel(
    data = centroids, mapping = aes(x = umap_1, y = umap_2, label = cluster),
    inherit.aes = FALSE, color = "black", fill = "white", alpha = 0.85,
    label.size = 0.3, fontface = "bold", size = 4,
    segment.color = "grey30", max.overlaps = Inf, seed = 1234
  ) +
  guides(color = guide_legend(title = "Cluster (n cells)", override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 8))

ggsave(file.path(plot_dir, "04_umap_clusters_res_0_4_annotated.png"),
       p_clusters, width = 10, height = 6, dpi = 300, bg = "white")

p_condition <- ggplot(umap_df, aes(umap_1, umap_2, color = condition)) +
  geom_point(size = 0.3, alpha = 0.6) +
  guides(color = guide_legend(title = "Condition", override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_classic() +
  theme(legend.position = "right")

combined <- p_clusters + p_condition + plot_layout(widths = c(1.15, 1))
ggsave(file.path(plot_dir, "04_umap_clusters_res_0_4_with_condition.png"),
       combined, width = 16, height = 6, dpi = 300, bg = "white")

message("Figure 10 fix written: 04_umap_clusters_res_0_4_annotated.png, 04_umap_clusters_res_0_4_with_condition.png")
