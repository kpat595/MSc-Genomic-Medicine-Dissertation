#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=240G
#SBATCH --partition=himem
#SBATCH --job-name umap_pool_condition_mixing
#SBATCH -o /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/09_umap_pool_condition_mixing_%j.out
#SBATCH -e /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/09_umap_pool_condition_mixing_%j.err

set -eo pipefail

mkdir -p /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial
export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/scripts/09_umap_pool_condition_mixing.R
