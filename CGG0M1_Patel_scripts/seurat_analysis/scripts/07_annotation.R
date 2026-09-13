suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(tidyverse)})

root <- "/well/bsg/projects/CGG0M1_Patel/seurat_analysis"
input_rds <- file.path(root, "rds", "06_annotated_sct.rds")
output_rds <- file.path(root, "rds", "07_annotated_sct.rds")
plot_dir <- file.path(root, "plots", "07_annotation")
table_dir <- file.path(root, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

final_resolution <- 0.4
res_out <- gsub("\\.", "_", formatC(final_resolution, format = "f", digits = 1))
res_dir <- file.path(root, "plots", "06_marker_evidence", paste0("res_", res_out))
cluster_col <- paste0("cluster_sct_res_", res_out)

annotation_file <- file.path(res_dir, paste0("06_marker_evidence_template_res_", res_out, ".tsv"))
if (!file.exists(annotation_file)) stop("Missing ", annotation_file)
annotation <- read.delim(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
if (any(!nzchar(trimws(annotation$final_cell_type)))) {
  stop("final_cell_type is blank for cluster(s) ", paste(annotation$cluster[!nzchar(trimws(annotation$final_cell_type))], collapse = ", "),
       " in ", annotation_file, " -- fill these in by hand before committing labels.")
}

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "SCT"
clusters <- as.character(obj@meta.data[[cluster_col]])
obj$final_cell_type <- unname(setNames(annotation$final_cell_type, annotation$cluster)[clusters])

p <- DimPlot(obj, reduction = "umap", group.by = "final_cell_type", label = TRUE, repel = TRUE, raster = TRUE) +
  NoLegend() + ggtitle(paste("Final annotation, resolution", formatC(final_resolution, format = "f", digits = 1)))
ggsave(file.path(plot_dir, "07_annotated_umap.png"), p, width = 9, height = 7, dpi = 300, bg = "white")

canonical_markers <- c(
  "SOX2", "HES1",       # radial glia
  "MKI67",              # cycling
  "EOMES",               # IPC
  "SLC17A7", "NEUROD6", # excitatory
  "GAD1", "GAD2",        # inhibitory
  "AQP4", "GFAP",        # astrocyte
  "OLIG2", "PDGFRA"      # OPC
)
canonical_markers <- canonical_markers[canonical_markers %in% rownames(obj)]
p <- DotPlot(obj, features = canonical_markers, group.by = "final_cell_type", assay = "SCT") + RotatedAxis() +
  ggtitle("Canonical marker sanity check by final cell type")
ggsave(file.path(plot_dir, "07_canonical_marker_dotplot.png"), p, width = 10, height = 7, dpi = 300, bg = "white")


antigen_genes <- c("GRIN1", "GRIN2A", "GRIN2B", "GRIN3A", "CNTNAP2", "CNTN2", "KCNA1", "KCNA2")
antigen_genes <- antigen_genes[antigen_genes %in% rownames(obj)]
p <- DotPlot(obj, features = antigen_genes, group.by = "final_cell_type", assay = "SCT") + RotatedAxis() +
  ggtitle("Antibody target antigen expression by annotated cell type")
ggsave(file.path(plot_dir, "07b_antigen_expression_dotplot.png"), p, width = 9, height = 7, dpi = 300, bg = "white")

agreement_cols <- intersect(
  c("cluster", "final_cell_type", "hnoca_top_label", "braun_top_label", "siletti_top_label",
    "ucell_top", "overrep_top", "cellhint_top_group", "methods_agree", "marker_evidence", "notes"),
  colnames(annotation)
)
method_agreement <- annotation[, agreement_cols]
write.table(method_agreement, file.path(table_dir, "07_method_agreement.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if ("methods_agree" %in% colnames(method_agreement)) {
  message("\n== Method agreement (flagging UCell/over-rep disagreement) ==")
  print(method_agreement[!method_agreement$methods_agree, ])
}

saveRDS(obj, output_rds)
capture.output(sessionInfo(), file = file.path(plot_dir, "07_annotation_session_info.txt"))
