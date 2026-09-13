#!/bin/bash
#SBATCH --job-name=check_pca_by_condition
#SBATCH --partition=himem
#SBATCH --array=1-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/03.5_check_pca_by_condition_%A_%a.out
#SBATCH --error=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/03.5_check_pca_by_condition_%A_%a.err

set -eo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial
export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition/scripts_by_condition/03.5_check_pca_by_condition.R "${SLURM_ARRAY_TASK_ID}"
