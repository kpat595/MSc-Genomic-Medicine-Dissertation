#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH --partition=himem
#SBATCH --job-name composition
#SBATCH -o /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/08_composition_%j.out
#SBATCH -e /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/08_composition_%j.err

set -eo pipefail

mkdir -p /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial
export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/scripts/08_composition.R
