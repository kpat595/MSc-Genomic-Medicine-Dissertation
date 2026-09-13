#!/bin/bash
set -euo pipefail

for pool in C1 C2 C3; do
    sbatch /well/bsg/projects/CGG0M1_Patel/import_preprocess/scripts/260428_run_cellranger_pool.sbatch "$pool"
done
