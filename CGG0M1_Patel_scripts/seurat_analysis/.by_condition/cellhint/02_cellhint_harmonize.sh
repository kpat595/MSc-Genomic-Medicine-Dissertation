#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=250G
#SBATCH --job-name cellhint_harmonize
#SBATCH --partition=himem
#SBATCH -o /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/cellhint_harmonize_%j.out
#SBATCH -e /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs/cellhint_harmonize_%j.err

set -eo pipefail

mkdir -p /well/bsg/projects/CGG0M1_Patel/seurat_analysis/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate || true
conda activate /well/bsg/projects/resources/shared_conda/conda_environs/cellhint

python3 /well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition/cellhint/02_cellhint_harmonize.py
