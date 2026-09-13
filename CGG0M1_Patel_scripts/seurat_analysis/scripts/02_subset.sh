#!/bin/bash
#SBATCH --job-name=seurat_subset
#SBATCH --array=1-12
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/260618_subset_%A_%a.out
#SBATCH --error=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/260618_subset_%A_%a.err

set -eo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial
export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/scripts/02_subset.R "${SLURM_ARRAY_TASK_ID}"
