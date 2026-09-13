#!/usr/bin/env Rscript

## Justifies not using Harmony in the main pipeline: cluster concordance
## between branches (ARI, majority vote), and kNN pool-mixing score for both
## embeddings vs. a shared permuted-label null.

suppressPackageStartupMessages({
  library(Seurat)
  library(RANN)
  library(ggplot2)
})

set.seed(1234)

root <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
plot_dir  <- file.path(root, ".integrated", "plots", "05c_harmony_vs_nonintegrated")
table_dir <- file.path(root, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
wt <- function(x, name) write.table(x, file.path(table_dir, name), sep = "\t", quote = FALSE, row.names = FALSE)

res <- 0.4 
cluster_col <- paste0("cluster_sct_res_", gsub("\\.", "_", formatC(res, format = "f", digits = 1)))
dims_use <- 1:15
k_neighbors <- 30
n_perm <- 100

load_branch <- function(rds, reduction) {
  obj <- readRDS(rds)
  out <- list(cluster = setNames(as.character(obj@meta.data[[cluster_col]]), colnames(obj)),
              emb = Embeddings(obj, reduction)[, dims_use])
  if (reduction == "pca") out$pool <- setNames(as.character(obj@meta.data[["pool"]]), colnames(obj))
  rm(obj); gc()
  out
}

pca <- load_branch(file.path(root, "rds", "04_clustered_sct.rds"), "pca")
harmony <- load_branch(file.path(root, ".integrated", "rds", "04_clustered_sct.rds"), "harmony")

cells <- intersect(names(pca$cluster), names(harmony$cluster))
cluster_pca <- pca$cluster[cells]; cluster_harmony <- harmony$cluster[cells]; pool <- pca$pool[cells]
emb_pca <- pca$emb[cells, ]; emb_harmony <- harmony$emb[cells, ]

## ---- cluster concordance: non-integrated vs. harmony ----------------------

choose2 <- function(x) x * (x - 1) / 2
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sa <- sum(choose2(rowSums(tab))); sb <- sum(choose2(colSums(tab)))
  expected <- sa * sb / choose2(n)
  (sum(choose2(tab)) - expected) / ((sa + sb) / 2 - expected)
}

tab_concord <- table(pca = cluster_pca, harmony = cluster_harmony)
ari_val <- ari(cluster_pca, cluster_harmony)
majority_agreement <- sum(apply(tab_concord, 1, max)) / sum(tab_concord)
wt(data.frame(resolution = res, n_clusters_pca = length(unique(cluster_pca)),
              n_clusters_harmony = length(unique(cluster_harmony)), adjusted_rand_index = ari_val,
              majority_vote_agreement = majority_agreement), "05c_cluster_concordance_summary.tsv")

p_conf <- ggplot(as.data.frame(tab_concord), aes(pca, harmony, fill = Freq)) +
  geom_tile() + scale_fill_viridis_c() +
  labs(title = paste0("Cluster concordance, res ", res),
       subtitle = paste0("ARI = ", round(ari_val, 3), ", majority vote = ", round(majority_agreement * 100, 1), "%")) +
  theme_minimal()
ggsave(file.path(plot_dir, "05c_cluster_concordance_heatmap.png"), p_conf, width = 7, height = 6)

## ---- kNN pool-mixing score: pca vs. harmony, shared permutation null ------

inv_simpson <- function(idx, labels) {
  vapply(seq_len(nrow(idx)), function(i) 1 / sum((table(labels[idx[i, ]]) / ncol(idx))^2), numeric(1))
}
perm_labels <- replicate(n_perm, sample(pool), simplify = FALSE)
n_pools <- length(unique(pool))

mixing_for <- function(emb, name) {
  nn <- RANN::nn2(emb, k = k_neighbors + 1)$nn.idx
  percell <- inv_simpson(nn, pool)
  observed <- median(percell)
  null <- vapply(perm_labels, function(p) median(inv_simpson(nn, p)), numeric(1))
  list(summary = data.frame(embedding = name, observed_median_score = observed, null_median_mean = mean(null),
                            null_median_sd = sd(null), fraction_null_le_observed = mean(null <= observed)),
       percell = data.frame(embedding = name, inverse_simpson = percell))
}

res_pca <- mixing_for(emb_pca, "pca")
res_harmony <- mixing_for(emb_harmony, "harmony")

mixing <- rbind(res_pca$summary, res_harmony$summary)
mixing$delta_vs_pca <- mixing$observed_median_score - mixing$observed_median_score[1]
wt(mixing, "05c_pool_mixing_pca_vs_harmony.tsv")

percell <- rbind(res_pca$percell, res_harmony$percell)
percell$embedding <- factor(percell$embedding, levels = c("pca", "harmony"),
                             labels = c("Non-integrated", "Harmony"))
wt(percell, "05c_pool_mixing_pca_vs_harmony_percell.tsv")

null_line <- mean(mixing$null_median_mean)

p_mix <- ggplot(percell, aes(embedding, inverse_simpson, fill = embedding)) +
  geom_violin(trim = TRUE, colour = NA, alpha = 0.8) +
  geom_hline(aes(yintercept = null_line, linetype = "Permuted-label null"), colour = "red", linewidth = 0.6) +
  geom_hline(aes(yintercept = n_pools, linetype = "Theoretical maximum (n pools)"), colour = "black", linewidth = 0.6) +
  scale_linetype_manual(name = NULL, values = c("Permuted-label null" = "dashed", "Theoretical maximum (n pools)" = "dotted")) +
  scale_fill_manual(values = c("Non-integrated" = "grey70", "Harmony" = "grey40"), guide = "none") +
  labs(x = "Embedding", y = "kNN pool-mixing score (inverse Simpson)") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file.path(plot_dir, "05c_pool_mixing_pca_vs_harmony.png"), p_mix, width = 6.5, height = 5.5)

message("ARI=", round(ari_val, 3), ", majority vote=", round(majority_agreement * 100, 1), "%")
message("Mixing score: pca=", round(mixing$observed_median_score[1], 3), ", harmony=", round(mixing$observed_median_score[2], 3),
        ", delta=", round(mixing$delta_vs_pca[2], 3))
