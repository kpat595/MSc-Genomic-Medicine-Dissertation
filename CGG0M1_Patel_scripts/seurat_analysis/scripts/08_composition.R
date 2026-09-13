suppressPackageStartupMessages({library(Seurat); library(speckle); library(tidyverse)})

root <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
input_rds <- file.path(root, "rds", "07_annotated_sct.rds")
table_dir <- file.path(root, "tables")
plot_dir <- file.path(root, "plots", "08_composition")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)

prop_test <- propeller(clusters = obj$final_cell_type, sample = obj$sample_id, group = obj$condition)
write_tsv(rownames_to_column(prop_test, "cell_type"), file.path(table_dir, "08_propeller_results.tsv"))

counts <- obj@meta.data %>% count(condition, pool, final_cell_type)

props_by_sample <- counts %>% group_by(condition, pool) %>% mutate(proportion = n / sum(n))
write_tsv(props_by_sample, file.path(table_dir, "08_proportions_by_sample.tsv"))
p <- ggplot(props_by_sample, aes(x = interaction(pool, condition), y = proportion, fill = final_cell_type)) +
  geom_col() + RotatedAxis() + ggtitle("Cell-type proportions per sample (pool x condition)")
ggsave(file.path(plot_dir, "08_proportions_stacked_bar_by_sample.png"), p, width = 10, height = 6, dpi = 300, bg = "white")

props_by_condition <- counts %>% group_by(condition, final_cell_type) %>% summarise(n = sum(n), .groups = "drop_last") %>%
  mutate(proportion = n / sum(n))
write_tsv(props_by_condition, file.path(table_dir, "08_proportions_by_condition.tsv"))
p <- ggplot(props_by_condition, aes(x = condition, y = proportion, fill = final_cell_type)) +
  geom_col() + ggtitle("Cell type proportions by condition") + labs(fill = "Annotated cell type")
ggsave(file.path(plot_dir, "08_proportions_stacked_bar.png"), p, width = 8, height = 6, dpi = 300, bg = "white")

cell_type_order <- props_by_sample %>%
  group_by(final_cell_type) %>% summarise(m = mean(proportion)) %>%
  arrange(desc(m)) %>% pull(final_cell_type)

props_by_sample_plot <- props_by_sample %>%
  mutate(
    final_cell_type = factor(final_cell_type, levels = cell_type_order),
    condition = factor(condition, levels = c("Control", "NoTreated", "CASPR2", "NR1"))
  )

ink <- "#0b0b0b"
muted <- "#898781"
grid_col <- "#e1e0d9"
baseline <- "#c3c2b7"

# match the stacked-bar's default ggplot hue palette, keyed by the same
# alphabetical factor order ggplot used there
cell_type_colors <- setNames(scales::hue_pal()(length(levels(factor(counts$final_cell_type)))),
                              levels(factor(counts$final_cell_type)))

p <- ggplot(props_by_sample_plot, aes(x = condition, y = proportion, color = final_cell_type)) +
  geom_point(size = 2.6,
             position = position_jitter(width = 0.08, height = 0)) +
  scale_color_manual(values = cell_type_colors, guide = "none") +
  facet_wrap(~final_cell_type, scales = "free_y", ncol = 4) +
  labs(x = NULL, y = "Proportion of cells per sample") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = grid_col, linewidth = 0.3),
    axis.line = element_line(color = baseline, linewidth = 0.3),
    axis.ticks = element_line(color = baseline, linewidth = 0.3),
    axis.text = element_text(color = muted),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(color = ink),
    strip.text = element_text(color = ink, face = "bold", size = 9),
    strip.background = element_blank(),
    panel.spacing = unit(1, "lines")
  )
ggsave(file.path(plot_dir, "08_proportions_by_sample_faceted.png"), p, width = 11, height = 8, dpi = 300, bg = "white")
