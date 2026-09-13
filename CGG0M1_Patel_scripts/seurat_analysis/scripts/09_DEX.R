suppressPackageStartupMessages({
  library(Seurat); library(DESeq2); library(tidyverse); library(patchwork)
  library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr)
})

root <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
input_rds <- file.path(root, "rds", "07_annotated_sct.rds")
table_dir <- file.path(root, "tables")
plot_dir <- file.path(root, "plots", "09_DEX")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

control_level <- "Control"
contrasts <- list(
  NR1_vs_control         = c("NR1", control_level),
  CASPR2_vs_control      = c("CASPR2", control_level),
  NR1_vs_CASPR2          = c("NR1", "CASPR2"),
  control_vs_NoTreated   = c(control_level, "NoTreated")
)
nmdar_genes <- c("GRIN1", "GRIN2A", "GRIN2B", "GRIN3A", "CNTNAP2", "NPAS4")  # CNTNAP2 = CASPR2; NPAS4 = a priori activity-dependent IEG (4h exposure window)

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

## Pseudobulk one cell type: counts summed per (condition, pool) pseudo-sample, then DESeq2.
run_de <- function(cell_type) {
  sub <- subset(obj, final_cell_type == cell_type)
  pb <- AggregateExpression(sub, assays = "RNA", slot = "counts", group.by = c("condition", "pool"))$RNA
  meta <- tibble(sample = colnames(pb)) %>% separate(sample, into = c("condition", "pool"), sep = "_")

  dds <- DESeqDataSetFromMatrix(countData = pb, colData = meta, design = ~ pool + condition)
  dds <- DESeq(dds, quiet = TRUE)

  map_dfr(names(contrasts), function(name) {
    lvl <- contrasts[[name]]
    if (!all(lvl %in% meta$condition)) return(NULL)
    as.data.frame(results(dds, contrast = c("condition", lvl[1], lvl[2]))) %>%
      rownames_to_column("gene") %>%
      mutate(cell_type = cell_type, contrast = name, .before = 1)
  })
}

cell_types <- unique(obj$final_cell_type)
de_results <- map_dfr(cell_types, function(ct) {
  tryCatch(run_de(ct), error = function(e) {
    message(ct, ": DESeq2 failed (", conditionMessage(e), ") -- skipping")
    NULL
  })
})
write_tsv(de_results, file.path(table_dir, "09_de_results.tsv"))

## Top hits: significant (padj < 0.05) and large-effect (|log2FC| > 1)
top_hits <- de_results %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1,
         !grepl("^(Possible|Probable|Unresolved)", cell_type, ignore.case = TRUE)) %>%
  arrange(padj, desc(abs(log2FoldChange)))
write_tsv(top_hits, file.path(table_dir, "09_de_top_hits_confirmed_celltypes.tsv"))

## Targeted readout: NMDAR subunits + CNTNAP2, mean expression / pct expressed by
## condition x cell type (descriptive only -- cell-level, not the pseudobulk DE above).
genes_present <- intersect(nmdar_genes, rownames(obj))
expr <- FetchData(obj, vars = c(genes_present, "condition", "final_cell_type"))
targeted <- expr %>%
  pivot_longer(all_of(genes_present), names_to = "gene", values_to = "expr") %>%
  group_by(gene, final_cell_type, condition) %>%
  summarise(mean_expr = mean(expr), pct_expressed = mean(expr > 0), .groups = "drop")
write_tsv(targeted, file.path(table_dir, "09_nmdar_cntnap2_targeted.tsv"))

p <- ggplot(targeted, aes(x = condition, y = gene, size = pct_expressed, colour = mean_expr)) +
  geom_point() + facet_wrap(~final_cell_type) +
  scale_colour_gradient(low = "grey85", high = "#D04A32") +
  ggtitle("NMDAR subunits + CNTNAP2 (CASPR2) + NPAS4 by condition and cell type")
ggsave(file.path(plot_dir, "09_nmdar_cntnap2_dotplot.png"), p, width = 12, height = 8, dpi = 300, bg = "white")

## Pathway enrichment: GSEA on the full ranked gene list, per cell type x contrast.
## Two independent gene set collections: GO:BP (gseGO) and MSigDB Hallmark (GSEA), as a
## cross-check -- Hallmark's 50 curated, non-redundant sets test whether GO:BP's fine-grained,
## overlapping terms reflect the same underlying biology.
symbol_to_entrez <- function(symbols) bitr(unique(symbols), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

de_entrez <- de_results %>% inner_join(symbol_to_entrez(.$gene), by = c("gene" = "SYMBOL"))

rank_genes <- function(de) {
  de %>% filter(!is.na(log2FoldChange)) %>% distinct(ENTREZID, .keep_all = TRUE) %>%
    arrange(desc(log2FoldChange)) %>% {setNames(.$log2FoldChange, .$ENTREZID)}
}

run_gsea_go <- function(de) {
  as.data.frame(gseGO(geneList = rank_genes(de), OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05))
}

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
run_gsea_hallmark <- function(de) {
  as.data.frame(GSEA(geneList = rank_genes(de), TERM2GENE = hallmark_sets, pvalueCutoff = 0.05)) %>%
    mutate(Description = str_to_sentence(str_replace_all(str_remove(ID, "^HALLMARK_"), "_", " ")))
}

by_group <- split(de_entrez, paste(de_entrez$cell_type, de_entrez$contrast, sep = "::"))
gsea_go_all <- bind_rows(map(by_group, run_gsea_go), .id = "cell_type_contrast")
gsea_hallmark_all <- bind_rows(map(by_group, run_gsea_hallmark), .id = "cell_type_contrast")
write_tsv(gsea_go_all, file.path(table_dir, "09_condition_gsea_go.tsv"))
write_tsv(gsea_hallmark_all, file.path(table_dir, "09_condition_gsea_hallmark.tsv"))

## GSEA dotplot -- x-axis grouped by contrast (fixed order), cell types in a consistent order
## within each contrast; y-axis = top 5 enriched terms per cell_type::contrast group by
## p.adjust, colored by NES sign/direction.
contrast_labels <- c(
  NR1_vs_control        = "NR1 vs. Control",
  CASPR2_vs_control     = "CASPR2 vs. Control",
  NR1_vs_CASPR2         = "NR1 vs. CASPR2",
  control_vs_NoTreated  = "Control vs. NoTreated"
)

make_gsea_dotplot <- function(gsea_df, out_file) {
  plot_df <- gsea_df %>%
    separate(cell_type_contrast, into = c("cell_type", "contrast"), sep = "::") %>%
    mutate(contrast = factor(contrast_labels[contrast], levels = contrast_labels),
           cell_type = factor(cell_type, levels = sort(unique(cell_type)))) %>%
    group_by(cell_type, contrast) %>% slice_min(p.adjust, n = 5, with_ties = FALSE) %>% ungroup() %>%
    mutate(Description = fct_reorder(Description, p.adjust, .desc = TRUE))

  dotplot <- ggplot(plot_df, aes(x = cell_type, y = Description, size = abs(NES), color = NES)) +
    geom_point() +
    facet_wrap(~contrast, nrow = 1) +
    scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0) +
    xlab("Cell type contrast") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 13))
  ggsave(out_file, dotplot, width = 30, height = 14, dpi = 300, bg = "white")
}

if (nrow(gsea_go_all) > 0) make_gsea_dotplot(gsea_go_all, file.path(plot_dir, "09_condition_gsea_go_dotplot.png"))
if (nrow(gsea_hallmark_all) > 0) make_gsea_dotplot(gsea_hallmark_all, file.path(plot_dir, "09_condition_gsea_hallmark_dotplot.png"))
