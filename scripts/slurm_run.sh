#!/bin/bash -l
#SBATCH --account=naiss2025-5-314
#SBATCH -n 128
#SBATCH -c 1
#SBATCH -p long
#SBATCH -J moritz_aDNA
#SBATCH --time=7-00:00:00
#SBATCH --output=aDNA_mapping.err
#SBATCH --error=aDNA_mapping.out

conda activate reholi

unset JAVA_HOME
unset JAVA_CMD

nextflow run -resume -c data/darde.config  ./scripts/reholi_from_bams.nf 
