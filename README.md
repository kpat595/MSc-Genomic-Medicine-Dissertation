# MSc-GM-Diss

Analysis code for the dissertation **"Single cell transcriptomic profiling of CASPR2 and NR1 antibody effects on iPSC-derived cortical neuronal subtypes"**, submitted for the MSc in Genomic Medicine, University of Oxford (Centre for Human Genetics, Nuffield Department of Medicine).

## Overview

Human iPSC-derived cortical cultures (day 78, single donor) were exposed for four hours to patient-derived anti-NR1-IgG1, patient-derived anti-CASPR2-A04-IgG4, an A33-IgG4 isotype control, or no treatment. Three technical replicate pools (C1–C3) were multiplexed within single GEM wells and profiled with the 10x Genomics Single Cell 5′ Gene Expression + Feature Barcoding (CITE-seq/ADT) workflow.

This repository contains the full computational pipeline used to go from raw sequencing data to the figures and tables in the dissertation:

1. Alignment and demultiplexing (Cell Ranger `multi`)
2. Ambient RNA removal (CellBender)
3. Import, per-sample QC and filtering (Seurat)
4. Per-pool SCTransform normalisation, consensus variable feature selection, PCA
5. Louvain clustering across resolutions, with bootstrap (Jaccard) stability and clustree validation
6. Annotation from multiple independent evidence streams (atlas label transfer, marker signature scoring, over-representation testing, a by-condition CellHint/CellTypist cross-check)
7. Compositional testing (propeller), pseudobulk differential expression (DESeq2), and gene set enrichment analysis (GO:BP and MSigDB Hallmark)

**This repository contains code only.** No sequencing data, intermediate objects, reference atlases, or result tables are included. Raw data and analysis outputs are available from the author on request.

## Requirements

See [`REQUIREMENTS.md`](REQUIREMENTS.md) for the full list of software, R/Python packages, reference files, and reference atlases needed to run this pipeline, along with install instructions and hardware expectations.

In brief: a SLURM HPC cluster with high-memory nodes, Cell Ranger 10.0.0, CellBender 0.3.2, R 4.4.3 (Seurat 5.4.0), Python 3 (scanpy / CellHint / CellTypist), the 10x GRCh38-2024-A reference, and three publicly available brain single-cell atlases.

## Repository layout

```
CGG0M1_Patel_scripts/
├── import_preprocess/scripts/        # Cell Ranger + CellBender (upstream preprocessing)
│   ├── 260428_run_cellranger_pool.sh
│   ├── 260428_submit_all_cellranger.sh
│   ├── 260513_make_cellbender_manifest.sh
│   ├── 260513_run_cellbender_array.sh
│   └── 260513_submit_cellbender.sh
│
└── seurat_analysis/
    ├── scripts/                      # Main pipeline (steps 01–09) + figure/table helpers
    ├── .by_condition/                # Supporting branch: per-condition clustering
    │   ├── scripts_by_condition/     #   steps 03–05 re-run separately per condition
    │   └── cellhint/                 #   CellHint harmonisation + CellTypist prediction
    └── .integrated/scripts/          # Supporting branch: Harmony integration check
```

Note that `.by_condition/` and `.integrated/` are dot-prefixed, so they are hidden by default in file browsers and in a plain `ls`.

Each analysis step is a pair of files: an `.R` (or `.py`) script holding the analysis, and a matching `.sh` SLURM submission script that sets up the environment and calls it. Three helper scripts (`04b_figure10_umap_res0_4.R`, `04c_umap_pool_condition.R`, `make_gsea_latex_tables.R`) have no wrapper and are run directly with `Rscript`.

## Main pipeline

Run in order; each step reads the RDS object written by the previous one.

| Step | Script | Purpose |
|---|---|---|
| 01 | `01_import_qc_array.R` | Build per-sample Seurat objects from CellBender H5s; split RNA/ADT assays; compute QC metrics (SLURM array, 12 samples) |
| 01a | `01a_QC_conditions.R` | QC violin/scatter plots and CLR-normalised ADT plots by condition |
| 02 | `02_subset.R` | Per-cell filtering: >500 genes, below sample-specific 90th percentile, <10% mitochondrial (SLURM array, 12 samples) |
| 03 | `03_merge_sct_pca.R` | Per-pool SCTransform (glmGamPoi), 3,000 consensus variable features, merge, PCA, CLR-normalise ADT |
| 03a | `03a_check_pca.R` | Elbow plot, variance summary, top PC loadings |
| 04 | `04_cluster_umap.R` | `FindNeighbors` on PCs 1–15, `FindClusters` at resolutions 0.1–1.0, UMAP |
| 04b | `04b_figure10_umap_res0_4.R` | Annotated res-0.4 cluster UMAP with per-cluster *n*, plus condition-coloured companion panel |
| 04c | `04c_umap_pool_condition.R` | Two-panel UMAP by pool and by condition |
| 05 | `05_bootstrap_cluster_stability.R` | 30-round Jaccard bootstrap across resolutions 0.3–1.0 (scclusteval) |
| 05c | `05c_clustree.R` | Clustering tree across resolutions 0.1–1.0 |
| 05c | `05c_harmony_vs_nonintegrated.R` | Adjusted Rand index, majority-vote concordance, and kNN pool-mixing (inverse Simpson) vs. permuted null |
| 06a | `06a_prepare_atlas_references.R` | Subsample the three reference atlases to ≤500 cells per label and convert to Seurat/SCE |
| 06 | `06_marker_evidence.R` | HNOCA/Braun label transfer (MapQuery), Siletti (SingleR), UCell signature scoring, clusterProfiler over-representation, marker dot plot, heatmaps; writes the hand-filled annotation template |
| 07 | `07_annotation.R` | Apply final cell type labels; canonical marker and antigen dot plots |
| 08 | `08_composition.R` | propeller (speckle) compositional testing; stacked proportion bar charts |
| 08 | `08_umap_cluster_vs_annotation.R` | Two-panel cluster-vs-annotation UMAP with a shared colour rule |
| 09 | `09_DEX.R` | Pseudobulk aggregation (cell type × condition × pool), DESeq2 with `~ pool + condition`, four contrasts; GSEA against GO:BP and MSigDB Hallmark |
| 09 | `09_umap_pool_condition_mixing.R` | Three-panel pool/condition mixing diagnostic |
| — | `make_gsea_latex_tables.R` | Render the GO and Hallmark GSEA summary TSVs as LaTeX `xltabular` tables |

Between steps 06 and 07, `06_marker_evidence_template_res_0_4.tsv` is filled in by hand with a `final_cell_type` per cluster. Step 07 reads that file, so the pipeline is not fully automatic across that boundary.

## Supporting branches

**By-condition clustering** (`seurat_analysis/.by_condition/`) — an independent check on global cluster boundaries. Steps 03–05 are re-run separately for each of the four conditions (2,000 variable features), the per-condition objects are converted to `.h5ad`, cluster labels are harmonised across conditions with CellHint, and cells are labelled with CellTypist models trained separately on the HNOCA and Braun references.

**Harmony integration check** (`seurat_analysis/.integrated/`) — runs Harmony on the existing PCA space grouped by pool, re-clusters, and compares against the non-integrated branch to decide whether integration is warranted. `run_harmony_pipeline.sh` submits the three jobs with SLURM dependencies.

## Running the pipeline
All scripts are written for a SLURM cluster and are submitted individually:

```bash
sbatch CGG0M1_Patel_scripts/seurat_analysis/scripts/03_merge_sct_pca.sh
```

Steps 01 and 02 are SLURM array jobs over the 12 samples (`--array=1-12`) and take the array index as their only argument. Memory requests range from 64 GB to 1.5 TB; see `REQUIREMENTS.md` for the per-step table.

**Paths are hardcoded.** Every script contains absolute paths to the original project directory (`/well/bsg/projects/CGG0M1_Patel/...`), the shared conda environments, and the project-local R library. To run this elsewhere, edit the `project_dir` / `root` / `sa` variable at the top of each `.R` script and the `conda activate`, `R_LIBS_USER` and `Rscript` lines in each `.sh` wrapper.

### Inputs that are not in this repository

| File | Needed by | Notes |
|---|---|---|
| Raw FASTQs | Cell Ranger | Available from the author on request |
| `<POOL>_multi.csv` | `260428_run_cellranger_pool.sh` | One per pool (C1, C2, C3); defines libraries, OCM barcodes and feature reference |
| `cellbender_inputs.tsv` | CellBender array job | Generated by `260513_make_cellbender_manifest.sh` |
| `seurat_analysis/tables/sample_manifest.tsv` | Steps 01, 02 | One row per sample with columns `sample_id`, `pool`, `condition`, `condition_class`, `antibody_target`, `expected_adt_feature`, `replicate_type`, `h5_path`, `cellbender_metrics_path` |
| `seurat_analysis/tables/06b_marker_gene_signatures.tsv` | Step 06 | Curated marker gene sets (14 signatures); the dissertation lists the gene membership |
| `seurat_analysis/reference/ref_{hnoca,braun,siletti_neurons,siletti_nonneuronal}.h5ad` | Step 06a | See `REQUIREMENTS.md` |

The condition labels for the CASPR2 and isotype-control samples were swapped during sequencing; this was identified from the ADT signal and corrected in `sample_manifest.tsv`, after which the whole post-QC pipeline was re-run. Any manifest used to reproduce this analysis must carry the corrected labels.

### Outputs
Scripts write into a project directory alongside `scripts/`:

```
seurat_analysis/
├── rds/        # Seurat objects per stage (00_import, 02_filtered, 03_..., 04_..., 06_..., 07_...)
├── tables/     # TSV/CSV summaries and result tables
├── plots/      # PNG figures, grouped by pipeline step
└── logs/       # SLURM stdout/stderr
```

## AI assistance
Analysis scripts in this repository were written with assistance from Claude Opus 5 (Anthropic). All code was reviewed and run by the author.

## Data availability
Analysis code is in this repository. Raw sequencing data and analysis outputs are available from the author upon request.

## Contact

Questions about the pipeline can be directed to the repository owner.
