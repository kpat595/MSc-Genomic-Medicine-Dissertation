#!/bin/bash
#SBATCH --job-name=check_pca_harmony
#SBATCH --partition=himem
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/03a_check_pca_harmony_%j.out
#SBATCH --error=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/03a_check_pca_harmony_%j.err

set -eo pipefail

mkdir -p /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial
export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/.integrated/scripts/03a_check_pca_harmony.R
