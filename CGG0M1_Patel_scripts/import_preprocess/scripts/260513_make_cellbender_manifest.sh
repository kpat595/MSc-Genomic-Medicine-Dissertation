#!/bin/bash
set -euo pipefail

PROJECT_DIR="/well/bsg/projects/CGG0M1_Patel/import_preprocess"
MANIFEST="${PROJECT_DIR}/scripts/cellbender_inputs.tsv"

mkdir -p "${PROJECT_DIR}/scripts" "${PROJECT_DIR}/cellbender_output" "${PROJECT_DIR}/cellbender_output/logs"

{
printf 'C1_Control\tC1\tControl\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C1/outs/per_sample_outs/Control/count/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C1/Control\n'
printf 'C1_NR1\tC1\tNR1\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C1/outs/per_sample_outs/NR1/count/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C1/NR1\n'
printf 'C1_NoTreated\tC1\tNoTreated\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C1/outs/per_sample_outs/NoTreated/count/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C1/NoTreated\n'
printf 'C1_CASPR2\tC1\tCASPR2\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C1/outs/per_sample_outs/CASPR2/count/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C1/CASPR2\n'
printf 'C2_Control\tC2\tControl\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C2/outs/per_sample_outs/Control/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C2/Control\n'
printf 'C2_NR1\tC2\tNR1\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C2/outs/per_sample_outs/NR1/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C2/NR1\n'
printf 'C2_NoTreated\tC2\tNoTreated\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C2/outs/per_sample_outs/NoTreated/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C2/NoTreated\n'
printf 'C2_CASPR2\tC2\tCASPR2\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C2/outs/per_sample_outs/CASPR2/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C2/CASPR2\n'
printf 'C3_Control\tC3\tControl\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C3/outs/per_sample_outs/Control/count/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C3/Control\n'
printf 'C3_NR1\tC3\tNR1\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C3/outs/per_sample_outs/NR1/count/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C3/NR1\n'
printf 'C3_NoTreated\tC3\tNoTreated\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C3/outs/per_sample_outs/NoTreated/count/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C3/NoTreated\n'
printf 'C3_CASPR2\tC3\tCASPR2\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/C3/outs/per_sample_outs/CASPR2/count/sample_raw_feature_bc_matrix.h5\t/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellbender_output/C3/CASPR2\n'
} > "${MANIFEST}"

wc -l "${MANIFEST}"
