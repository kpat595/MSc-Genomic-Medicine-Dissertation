#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=250G
#SBATCH --array=1-4
#SBATCH --job-name cellhint_convert_h5ad
#SBATCH --partition=himem
#SBATCH -o /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/cellhint_convert_%A_%a.out
#SBATCH -e /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/cellhint_convert_%A_%a.err

set -eo pipefail

mkdir -p /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial
export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs
export RETICULATE_PYTHON=/well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial/bin/python3

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition/cellhint/01_convert_rds_to_h5ad.R "${SLURM_ARRAY_TASK_ID}"
