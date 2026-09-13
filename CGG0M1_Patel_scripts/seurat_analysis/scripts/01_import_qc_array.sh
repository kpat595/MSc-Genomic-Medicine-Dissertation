#!/bin/bash -l
#SBATCH --job-name=seurat_import_qc
#SBATCH --array=1-12
#SBATCH -c 2
#SBATCH -o /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/seurat_import_qc_%A_%a.out
#SBATCH -e /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/seurat_import_qc_%A_%a.err

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial

export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/scripts/01_import_qc_array.R "${SLURM_ARRAY_TASK_ID}"
