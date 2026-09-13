#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=150G
#SBATCH --job-name celltypist_predict
#SBATCH --partition=himem
#SBATCH -o /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/celltypist_predict_%j.out
#SBATCH -e /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/celltypist_predict_%j.err

set -eo pipefail

mkdir -p /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/cellhint

python3 /well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition/cellhint/06_celltypist_predict.py
