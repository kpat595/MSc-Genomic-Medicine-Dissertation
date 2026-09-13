#!/bin/bash
set -euo pipefail

bash /well/bsg/projects/CGG0M1_Patel/import_preprocess/scripts/260513_make_cellbender_manifest.sh
sbatch /well/bsg/projects/CGG0M1_Patel/import_preprocess/scripts/260513_run_cellbender_array.sbatch
