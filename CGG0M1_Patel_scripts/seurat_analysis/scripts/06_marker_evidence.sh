#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=300G
#SBATCH --partition=himem
#SBATCH --job-name marker_evidence
#SBATCH -o /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/06_marker_evidence_%j.out
#SBATCH -e /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/06_marker_evidence_%j.err

set -eo pipefail

mkdir -p /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial
export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

## one-time setup: install UCell/clusterProfiler/ggalluvial into R_LIBS_USER if not already present
Rscript -e '
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  lib <- Sys.getenv("R_LIBS_USER")
  if (!requireNamespace("UCell", quietly = TRUE)) BiocManager::install("UCell", lib = lib, update = FALSE, ask = FALSE)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler", lib = lib, update = FALSE, ask = FALSE)
  if (!requireNamespace("ggalluvial", quietly = TRUE)) install.packages("ggalluvial", lib = lib, repos = "https://cloud.r-project.org")
'

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/scripts/06_marker_evidence.R
