#!/bin/bash
#SBATCH --job-name=harmony_integrate
#SBATCH --partition=himem
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --output=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/03b_harmony_integrate_%j.out
#SBATCH --error=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/03b_harmony_integrate_%j.err

set -eo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
module load HDF5
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial
export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

Rscript /well/bsg/projects/CGG0M1_Patel/seurat_analysis/.integrated/scripts/03b_harmony_integrate.R
