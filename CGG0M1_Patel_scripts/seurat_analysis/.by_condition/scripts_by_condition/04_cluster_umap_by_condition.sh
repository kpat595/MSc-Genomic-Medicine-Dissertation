#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=180G
#SBATCH --array=1-4
#SBATCH --job-name cluster_umap_sct_by_condition
#SBATCH --partition=himem
#SBATCH -o /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/04_cluster_umap_by_condition_%A_%a.out
#SBATCH -e /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/04_cluster_umap_by_condition_%A_%a.err

set -eo pipefail

mkdir -p /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial
export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition/scripts_by_condition/04_cluster_umap_by_condition.R "${SLURM_ARRAY_TASK_ID}"
