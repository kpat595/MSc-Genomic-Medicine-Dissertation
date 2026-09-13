#!/usr/bin/env Rscript
# Build LaTeX (xltabular) tables for the GO and Hallmark GSEA results.
# Reads the summary TSVs and writes ready-to-\input .tex files.

acronym_map <- c(
  Dna = "DNA", Rna = "RNA", Mrna = "mRNA", Utr = "UTR", Kras = "KRAS",
  Tnfa = "TNFA", Nfkb = "NFkB", Il2 = "IL2", Il6 = "IL6", Tgf = "TGF",
  G2m = "G2M", E2f = "E2F", Uv = "UV", Pi3k = "PI3K", Atp = "ATP",
  Mtorc1 = "mTORC1", Mtorc2 = "mTORC2", Stat3 = "STAT3", Stat5 = "STAT5",
  Il = "IL"
)
fix_acronyms <- function(x) {
  for (from in names(acronym_map)) {
    x <- gsub(paste0("\\b", from, "\\b"), acronym_map[[from]], x, ignore.case = TRUE)
  }
  x
}

esc <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%_#&$])", "\\\\\\1", x)
  x
}

fmt_p <- function(p) {
  ifelse(p < 1e-3,
         sprintf("$%.2f\\times10^{%d}$", p / 10^floor(log10(p)), floor(log10(p))),
         sprintf("%.3f", p))
}

write_table <- function(df, cols, header, colspec, caption, label, out_path) {
  body <- df[, cols, drop = FALSE]
  rows <- apply(body, 1, function(r) paste(r, collapse = " & "))
  lines <- c(
    "\\begin{landscape}",
    "\\begingroup",
    "\\singlespacing",
    "\\hyphenpenalty=10000",
    "\\exhyphenpenalty=10000",
    "\\tolerance=9999",
    "\\emergencystretch=2em",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{1.4}",
    "\\rowcolors{2}{white}{gray!12}",
    sprintf("\\begin{xltabular}{\\linewidth}{@{}%s@{}}", colspec),
    sprintf("\\caption{%s}\\label{%s}\\\\", caption, label),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    "\\endfirsthead",
    "\\hiderowcolors",
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    "\\endhead",
    "\\showrowcolors",
    "\\bottomrule",
    "\\endlastfoot",
    paste0(rows, " \\\\"),
    "\\end{xltabular}",
    "\\endgroup",
    "\\end{landscape}"
  )
  writeLines(lines, out_path)
  message("Wrote ", out_path, " (", nrow(df), " rows)")
}

split_contrast <- function(df) {
  parts <- strsplit(df$cell_type_contrast, "::", fixed = TRUE)
  df$cell_type <- vapply(parts, `[`, character(1), 1)
  df$contrast <- vapply(parts, function(p) if (length(p) > 1) p[2] else NA, character(1))
  df$contrast <- gsub("_", " ", df$contrast, fixed = TRUE)
  df
}

## ---- GO table ------------------------------------------------------------
go <- read.delim("summary/09_condition_gsea_go.tsv", stringsAsFactors = FALSE)
go <- split_contrast(go)
go <- go[order(go$cell_type, go$contrast, go$p.adjust), ]
go$cell_type_e   <- esc(go$cell_type)
go$contrast_e    <- esc(go$contrast)
go$ID_e          <- esc(go$ID)
go$Description_e <- esc(go$Description)
go$NES_s         <- sprintf("%.2f", go$NES)
go$padj_s        <- fmt_p(go$p.adjust)

write_table(
  go,
  cols     = c("cell_type_e", "contrast_e", "ID_e", "Description_e", "setSize", "NES_s", "padj_s"),
  header   = "\\textbf{Cell type} & \\textbf{Contrast} & \\textbf{GO ID} & \\textbf{Description} & \\textbf{Set size} & \\textbf{NES} & \\textbf{Adj.\\ $p$}",
  colspec  = "P{.14}P{.16}P{.11}P{.34}P{.06}P{.06}P{.09}",
  caption  = "GO gene set enrichment analysis (GSEA) results, all pathways significant at adjusted $p<0.05$, by cell type and condition contrast.",
  label    = "tab:gsea_go_full",
  out_path = "seurat_analysis/tables/09_condition_gsea_go_table.tex"
)

## ---- Hallmark table -------------------------------------------------------
hm <- read.delim("summary/09_condition_gsea_hallmark.tsv", stringsAsFactors = FALSE)
hm <- split_contrast(hm)
hm <- hm[order(hm$cell_type, hm$contrast, hm$p.adjust), ]
hm$cell_type_e   <- esc(hm$cell_type)
hm$contrast_e    <- esc(hm$contrast)
hm$Description_e <- esc(fix_acronyms(hm$Description))
hm$NES_s         <- sprintf("%.2f", hm$NES)
hm$padj_s        <- fmt_p(hm$p.adjust)

write_table(
  hm,
  cols     = c("cell_type_e", "contrast_e", "Description_e", "setSize", "NES_s", "padj_s"),
  header   = "\\textbf{Cell type} & \\textbf{Contrast} & \\textbf{Hallmark pathway} & \\textbf{Set size} & \\textbf{NES} & \\textbf{Adj.\\ $p$}",
  colspec  = "P{.17}P{.18}P{.35}P{.07}P{.07}P{.10}",
  caption  = "Hallmark gene set enrichment analysis (GSEA) results, all pathways significant at adjusted $p<0.05$, by cell type and condition contrast.",
  label    = "tab:gsea_hallmark_full",
  out_path = "seurat_analysis/tables/09_condition_gsea_hallmark_table.tex"
)
