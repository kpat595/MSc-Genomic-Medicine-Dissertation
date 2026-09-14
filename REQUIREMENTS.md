# Requirements and setup
Everything needed to run the pipeline in this repository. Versions listed are the ones the analysis was run with; other versions may work but have not been tested.

---

## 1. Compute environment
The pipeline was developed and run on a SLURM HPC cluster with high-memory nodes. It is not intended to run on a laptop — several steps need hundreds of gigabytes of RAM on a 406,380-cell object.

| Step | CPUs | Memory | Notes |
|---|---|---|---|
| Cell Ranger (`multi`) | 32 | 128 GB | Per pool; up to 72 h |
| CellBender | 16 | 256 GB | Array of 12, throttled to 4 concurrent |
| 01 import / QC | 2 | default | Array of 12 |
| 01a QC by condition | 4 | 64 GB | |
| 02 subset | 8 | 64 GB | Array of 12 |
| 03 merge / SCT / PCA | 8 | 192 GB | |
| 03a check PCA | 4 | 128 GB | |
| 04 cluster / UMAP | 16 | 180 GB | |
| 05 bootstrap stability (30 rounds) | 4 | 500 GB | Long-running |
| 05a bootstrap stability (100 rounds, res 0.4) | 16 | 1.5 TB | Longest job in the pipeline |
| 05c clustree | 4 | 300 GB | |
| 05c harmony vs. non-integrated | 8 | 350 GB | |
| 06a prepare atlas references | 4 | 500 GB | Atlases are read from `.h5ad` in backed mode |
| 06 marker evidence | 16 | 300 GB | |
| 07 annotation | 8 | 240 GB | |
| 08 composition | 4 | 64 GB | |
| 08 cluster vs. annotation UMAP | 8 | 240 GB | |
| 09 DEX + GSEA | 8 | 500 GB | |

A `himem` partition (or equivalent) is assumed. The HDF5 module is loaded in each wrapper (`module load HDF5/1.14.0-gompi-2023a`); on a cluster without environment modules, HDF5 needs to be on `LD_LIBRARY_PATH` some other way so that `hdf5r` can read the CellBender outputs.

---

## 2. Command-line tools

### Cell Ranger 10.0.0
Download from 10x Genomics (free, requires accepting the licence):
<https://www.10xgenomics.com/support/software/cell-ranger/downloads>

The scripts call the binary by absolute path, so it can be unpacked anywhere — just update `CELLRANGER=` in `260428_run_cellranger_pool.sh`.

### Reference transcriptome: GRCh38-2024-A
Download the human reference package from the same 10x downloads page. It is referenced from the `<POOL>_multi.csv` config files rather than from the shell scripts.

### CellBender 0.3.2
<https://github.com/broadinstitute/CellBender>

```bash
pip install cellbender==0.3.2
```

Run here on CPU (`--cpu-threads`), with `--fpr 0 --epochs 150 --checkpoint-mins 10`. A GPU will make it considerably faster if one is available.

---

## 3. R environment

**R 4.4.3.** The project uses a project-local library set through `R_LIBS_USER`, layered on top of a shared conda environment:
```bash
export R_LIBS_USER=/path/to/seurat_analysis/R_libs
```

### Core packages (versions used)

| Package | Version | Source |
|---|---|---|
| Seurat | 5.4.0 | CRAN |
| SeuratObject | 5.3.0 | CRAN |
| DESeq2 | 1.46.0 | Bioconductor |
| clusterProfiler | 4.14.6 | Bioconductor |
| glmGamPoi | 1.18.0 | Bioconductor |
| harmony | 2.0.5 | CRAN |
| speckle | 1.6.0 | Bioconductor |
| clustree | 0.5.1 | CRAN |
| UCell | 2.10.1 | Bioconductor |
| SingleR | 2.8.0 | Bioconductor |
| reticulate | 1.46.0 | CRAN |

### Full package list, by source

**CRAN**

```r
install.packages(c(
  "Seurat", "SeuratObject", "Matrix", "hdf5r", "tidyverse", "dplyr",
  "ggplot2", "ggrepel", "patchwork", "pheatmap", "colorspace",
  "RANN", "future", "clustree", "harmony", "msigdbr", "reticulate",
  "ggalluvial", "remotes", "BiocManager"
))
```

**Bioconductor**

```r
BiocManager::install(c(
  "DESeq2", "clusterProfiler", "org.Hs.eg.db", "glmGamPoi",
  "SingleR", "SingleCellExperiment", "UCell", "speckle"
))
```

**GitHub**

```r
remotes::install_github("crazyhottommy/scclusteval", upgrade = "never")   # Jaccard bootstrap stability
remotes::install_github("cellgeni/sceasy")                                # Seurat <-> AnnData conversion
```

`05_bootstrap_cluster_stability.sh` and `06_marker_evidence.sh` contain idempotent one-time install blocks that will pull `scclusteval`, `UCell`, `clusterProfiler` and `ggalluvial` into `R_LIBS_USER` on first run if they are missing.

### Known package interaction

`sceasy` 0.0.7's `seurat2anndata()` calls `Seurat::GetAssayData(..., slot=)`, which is defunct in SeuratObject ≥ 5.0.0. `.by_condition/cellhint/01_convert_rds_to_h5ad.R` patches the function in-session rather than modifying the installed package — no action needed, but worth knowing if the conversion step is reused elsewhere.

---

## 4. Python environment

Two conda environments were used on the original cluster: one carrying R plus CellBender (`seurat_spatial`), and one for the CellHint/CellTypist branch (`cellhint`). A single environment works as long as all of the following are importable.

```bash
conda create -n msc_gm python=3.10
conda activate msc_gm
pip install scanpy anndata numpy pandas cellhint celltypist cellbender==0.3.2
```

Packages used directly by the scripts: `scanpy`, `anndata`, `numpy`, `pandas`, `cellhint`, `celltypist`.

`reticulate` must point at the environment that has `anndata` installed:

```bash
export RETICULATE_PYTHON=/path/to/env/bin/python3
export PYTHONNOUSERSITE=1
```

`06a_prepare_atlas_references.R` additionally calls `use_condaenv()` with a hardcoded path — update it to your environment.

CellHint's `sankeyplot()` hangs on a headless node with its default `show=True`; the script already passes `show=False`.

---

## 5. Reference atlases

Three public single-cell atlases are needed for annotation, all downloadable as `.h5ad` from CELLxGENE Discover (<https://cellxgene.cziscience.com/>). Place them in `seurat_analysis/reference/` under these exact filenames, since `06a_prepare_atlas_references.R` looks for them by name:

| File | Atlas | `obs` column used |
|---|---|---|
| `ref_hnoca.h5ad` | Human Neural Organoid Cell Atlas (He et al. 2024) | `annot_level_2` |
| `ref_braun.h5ad` | Braun et al. 2023, first-trimester developing human brain | `CellClass` |
| `ref_siletti_neurons.h5ad` | Siletti et al. 2023, adult human brain — neurons | `supercluster_term` |
| `ref_siletti_nonneuronal.h5ad` | Siletti et al. 2023, adult human brain — non-neuronal | `supercluster_term` |

Gene symbols are taken from `var["feature_name"]` in all four. HNOCA is read from `raw.X`; the others from `X`. Each is subsampled to ≤500 cells per label before use, which is what makes label transfer tractable — but the full files are large (the source atlases hold 1.6–1.7 M cells each), so allow for tens of gigabytes of download and disk, and run step 06a on a high-memory node.

Step 06a writes `ref_hnoca.rds`, `ref_braun.rds` and `ref_siletti.rds` (plus subsampled `.h5ad` copies, which the CellTypist branch reuses for training).

Gene set resources are fetched at runtime rather than downloaded in advance: GO:BP via `org.Hs.eg.db`, and the MSigDB Hallmark collection via `msigdbr`.

---

## 6. Directory scaffolding

Most scripts call `dir.create(..., recursive = TRUE)` for their own outputs, but the following need to exist, with contents, before the pipeline starts:

```
<project_dir>/
├── import_preprocess/
│   ├── cellranger_output/          # plus <POOL>_multi.csv config per pool
│   └── scripts/
└── seurat_analysis/
    ├── reference/                  # the four atlas .h5ad files above
    ├── tables/
    │   ├── sample_manifest.tsv
    │   └── 06b_marker_gene_signatures.tsv
    ├── R_libs/                     # project-local R library
    ├── logs/
    └── scripts/
```

### `sample_manifest.tsv`

Tab-separated, one row per sample (12 rows: 3 pools × 4 conditions). Row order defines the SLURM array index used by steps 01 and 02. Columns read by the scripts:

`sample_id`, `pool`, `condition`, `condition_class`, `antibody_target`, `expected_adt_feature`, `replicate_type`, `h5_path` (the CellBender filtered `.h5`), `cellbender_metrics_path`

### `06b_marker_gene_signatures.tsv`

One row per cell type signature with columns `id`, `cell_type` and `genes` (comma-separated gene symbols). The 14 signatures and their gene membership are listed in the dissertation methods.

---

## 7. Reproducibility notes

- Absolute paths to `/well/bsg/projects/CGG0M1_Patel/...` are hardcoded throughout; they must be changed to run this anywhere else.
- Random seeds are fixed where they affect results (`set.seed(1234)` for atlas subsampling, `seed <- 20260707` for the bootstrap rounds), but Louvain clustering and UMAP can still vary slightly across BLAS implementations and package versions.
- The annotation step is not fully automatic: after step 06, `final_cell_type` is filled in by hand in `06_marker_evidence_template_res_0_4.tsv`, and step 07 reads that file.
