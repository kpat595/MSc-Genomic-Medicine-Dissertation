#!/bin/bash
#SBATCH --job-name=clustree
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300G
#SBATCH --time=00:30:00
#SBATCH --output=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/05c_clustree_%j.out
#SBATCH --error=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/05c_clustree_%j.err

set -eo pipefail

echo "BEGIN at $(date) on $(hostname)"

module load HDF5

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial

export R_LIBS_USER=/well/bsg/projects/CGG0M1_Patel/seurat_analysis/R_libs

Rscript "/well/bsg/projects/CGG0M1_Patel/seurat_analysis/scripts/05c_clustree.R"
