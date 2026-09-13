## Three-panel batch/condition diagnostic figure.
## (A) UMAP coloured by pool, (B) UMAP coloured by condition -- both on the
## non-integrated embedding used throughout the main pipeline, so a reader can
## see directly whether pools or conditions form separate islands (and, for
## (B), that the effect of interest isn't just a pool/batch split).
## (C) per-cell kNN pool-mixing score (inverse Simpson) for the non-integrated
## vs. Harmony embeddings, against the permuted-label null and the theoretical
## maximum for 3 pools -- quantifies what (A) shows by eye and documents why
## Harmony integration was not required for the main pipeline.

.libPaths(c("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs", .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})

sa <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
plot_dir <- file.path(sa, "plots", "09_pool_condition_mixing")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

## ---- panels A/B: non-integrated UMAP by pool and condition ----------------
umap_cache <- file.path(sa, "tables", "09_umap_pool_condition.tsv")
if (file.exists(umap_cache)) {
  umap_df <- read.delim(umap_cache, stringsAsFactors = FALSE)
} else {
  obj <- readRDS(file.path(sa, "rds", "04_clustered_sct.rds"))
  umap_df <- FetchData(obj, vars = c("umap_1", "umap_2", "pool", "condition"))
  rm(obj); gc()

  ## random plotting order so no single category systematically overplots
  ## the others in dense regions (406,380 cells; rows are grouped by sample)
  set.seed(1234)
  umap_df <- umap_df[sample(nrow(umap_df)), ]
  write.table(umap_df, umap_cache, sep = "\t", quote = FALSE, row.names = FALSE)
}

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
  base_theme

p_condition <- ggplot(umap_df, aes(umap_1, umap_2, color = condition)) +
  geom_point(size = 0.15, alpha = 0.7) +
  guides(color = guide_legend(title = "Condition", override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  base_theme

## ---- panel C: kNN pool-mixing score, non-integrated vs. Harmony -----------
percell <- read.delim(file.path(sa, "tables", "05c_pool_mixing_pca_vs_harmony_percell.tsv"), stringsAsFactors = FALSE)
percell$embedding <- factor(percell$embedding, levels = c("Non-integrated", "Harmony"))

mixing_summary <- read.delim(file.path(sa, "tables", "05c_pool_mixing_pca_vs_harmony.tsv"), stringsAsFactors = FALSE)
null_value <- unique(round(mixing_summary$null_median_mean, 2))[1]
n_pools <- length(unique(umap_df$pool))

p_mix <- ggplot(percell, aes(embedding, inverse_simpson, fill = embedding)) +
  geom_violin(trim = TRUE, colour = NA, alpha = 0.8, linewidth = 0) +
  stat_summary(fun = median, geom = "point", size = 2, shape = 21, fill = "white", colour = "black") +
  geom_hline(aes(yintercept = null_value, linetype = paste0("Permuted-label null (", null_value, ")")),
             colour = "red", linewidth = 0.6) +
  geom_hline(aes(yintercept = n_pools, linetype = paste0("Theoretical maximum (", n_pools, " pools)")),
             colour = "black", linewidth = 0.6) +
  scale_linetype_manual(name = NULL, values = setNames(c("dashed", "dotted"),
                         c(paste0("Permuted-label null (", null_value, ")"), paste0("Theoretical maximum (", n_pools, " pools)")))) +
  scale_fill_manual(values = c("Non-integrated" = "#118AB2", "Harmony" = "#EF476F"), guide = "none") +
  guides(linetype = guide_legend(nrow = 1, override.aes = list(colour = c("red", "black")))) +
  labs(x = NULL, y = "kNN pool-mixing score\n(inverse Simpson)") +
  theme_classic() +
  theme(legend.position = "bottom", legend.text = element_text(size = 9))

## ---- combine ----------------------------------------------------------
combined <- (p_pool | p_condition | p_mix) +
  plot_layout(widths = c(1, 1, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(file.path(plot_dir, "09_umap_pool_condition_mixing_3panel.png"),
       combined, width = 19, height = 5.8, dpi = 300, bg = "white")

message("Written 09_umap_pool_condition_mixing_3panel.png")
message("Null=", null_value, ", n_pools=", n_pools,
        ", medians: ", paste(mixing_summary$embedding, round(mixing_summary$observed_median_score, 2), collapse = "; "))
