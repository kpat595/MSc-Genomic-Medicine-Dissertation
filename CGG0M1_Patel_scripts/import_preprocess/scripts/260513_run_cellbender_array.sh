#!/bin/bash
#SBATCH --job-name=CellBender
#SBATCH --cpus-per-task=16
#SBATCH --partition=himem
#SBATCH --mem=256G
#SBATCH --array=1-12%4
#SBATCH --output=/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/logs/cellbender_%A_%a.out
#SBATCH --error=/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/logs/cellbender_%A_%a.err

set -euo pipefail

PROJECT_DIR="/well/bsg/projects/CGG0M1_Patel/import_preprocess"
MANIFEST="${PROJECT_DIR}/scripts/cellbender_inputs.tsv"
SEURAT_ENV="/well/bsg/projects/resources/shared_conda/conda_environs/seurat_spatial"
CELLBENDER="${SEURAT_ENV}/bin/cellbender"

module is-loaded HDF5/1.14.0-gompi-2023a || module load HDF5/1.14.0-gompi-2023a

export LD_LIBRARY_PATH="${SEURAT_ENV}/lib:${LD_LIBRARY_PATH:-}"
export PATH="${SEURAT_ENV}/bin:${PATH}"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
IFS=$'\t' read -r SAMPLE_ID POOL CONDITION INPUT_H5 OUT_DIR <<< "${LINE}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "sample_id=${SAMPLE_ID}"
echo "pool=${POOL}"
echo "condition=${CONDITION}"
echo "input_h5=${INPUT_H5}"
echo "out_dir=${OUT_DIR}"
echo "cellbender=${CELLBENDER}"
echo "cellbender_version=$("${CELLBENDER}" --version)"
echo "started=$(date)"

"${CELLBENDER}" remove-background --input "${INPUT_H5}" --output "${OUT_DIR}/${SAMPLE_ID}_cellbender.h5" --fpr 0 --epochs "${CELLBENDER_EPOCHS:-150}" --checkpoint-mins "${CELLBENDER_CHECKPOINT_MINS:-10}" --cpu-threads "${SLURM_CPUS_PER_TASK}"

echo "finished=$(date)"
echo "full_output=${OUT_DIR}/${SAMPLE_ID}_cellbender.h5"
echo "filtered_output=${OUT_DIR}/${SAMPLE_ID}_cellbender_filtered.h5"
