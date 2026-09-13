#!/bin/bash
#SBATCH --job-name=bootstrap_stability
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G
#SBATCH --output=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/05_bootstrap_stability_%j.out
#SBATCH --error=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/05_bootstrap_stability_%j.err

set -eo pipefail

echo "BEGIN at $(date) on $(hostname)"

module load HDF5

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial

export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

## one-time setup: install scclusteval into R_LIBS_USER if not already present
Rscript -e '
  if (!requireNamespace("scclusteval", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
    remotes::install_github("crazyhottommy/scclusteval", upgrade = "never")
  }
'

Rscript "/well/bsg/projects/CGG0M1_Patel/seurat_analysis/scripts/05_bootstrap_cluster_stability.R"
