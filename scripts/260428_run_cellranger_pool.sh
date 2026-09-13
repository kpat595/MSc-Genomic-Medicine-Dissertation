#!/bin/bash
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH -p himem
#SBATCH --job-name=CellRanger
#SBATCH -o /well/bsg/projects/CGG0M1_Patel/CellRanger_multi_%j.out
#SBATCH -e /well/bsg/projects/CGG0M1_Patel/CellRanger_multi_%j.err

CELLRANGER=/gpfs3/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger-10.0.0/cellranger

POOL=$1

cd /well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output

"$CELLRANGER" multi \
  --id="$POOL" \
  --csv="/well/bsg/projects/CGG0M1_Patel/import_preprocess/cellranger_output/${POOL}_multi.csv" \
  --localcores=32 \
  --localmem=128
