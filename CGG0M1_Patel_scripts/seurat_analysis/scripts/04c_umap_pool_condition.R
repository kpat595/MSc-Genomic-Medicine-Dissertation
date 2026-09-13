## Two-panel UMAP: (a) pool, (b) condition, side by side.
## No in-plot titles (panel tags + legend titles + figure legend carry the info,
## per common journal practice - avoids redundant/duplicated text in the figure).

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})

sa <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
plot_dir <- file.path(sa, "plots", "04_clustering")

obj <- readRDS(file.path(sa, "rds", "04_clustered_sct.rds"))

umap_df <- FetchData(obj, vars = c("umap_1", "umap_2", "pool", "condition"))

message("Pool counts:"); print(table(umap_df$pool))
message("Condition counts:"); print(table(umap_df$condition))

## Shuffle plotting order so no single category systematically overplots
## the others in dense regions (rows were originally grouped by sample).
set.seed(1234)
umap_df <- umap_df[sample(nrow(umap_df)), ]

base_theme <- theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.spacing.y = unit(1, "pt"),
    legend.box.spacing = unit(4, "pt"),
    legend.margin = margin(l = 0)
  )

p_pool <- ggplot(umap_df, aes(umap_1, umap_2, color = pool)) +
  geom_point(size = 0.15, alpha = 0.7) +
  guides(color = guide_legend(title = "Pool", override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  base_theme +
  theme(plot.margin = margin(t = 5, r = 14, b = 5, l = 5))

p_condition <- ggplot(umap_df, aes(umap_1, umap_2, color = condition)) +
  geom_point(size = 0.15, alpha = 0.7) +
  guides(color = guide_legend(title = "Condition", override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  base_theme +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 14))

combined <- (p_pool | p_condition) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

ggsave(file.path(plot_dir, "04_umap_pool_condition_2panel.png"),
       combined, width = 13, height = 5.5, dpi = 300, bg = "white")

message("Written 04_umap_pool_condition_2panel.png")
