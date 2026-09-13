#!/bin/bash
set -eo pipefail

SCRIPT_DIR="/well/bsg/projects/CGG0M1_Patel/seurat_analysis/.integrated/scripts"

JOB_HARMONY=$(sbatch --parsable "$SCRIPT_DIR/03b_harmony_integrate.sh")
JOB_CHECK_PCA=$(sbatch --parsable --dependency=afterok:"$JOB_HARMONY" "$SCRIPT_DIR/03a_check_pca_harmony.sh")
JOB_CLUSTER_UMAP=$(sbatch --parsable --dependency=afterok:"$JOB_CHECK_PCA" "$SCRIPT_DIR/04_cluster_umap_harmony.sh")

echo "Submitted 03b_harmony_integrate:   $JOB_HARMONY"
echo "Submitted 03a_check_pca_harmony:   $JOB_CHECK_PCA (after $JOB_HARMONY)"
echo "Submitted 04_cluster_umap_harmony: $JOB_CLUSTER_UMAP (after $JOB_CHECK_PCA)"
