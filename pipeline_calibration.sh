#!/bin/bash
#SBATCH --job-name=apg_v2
#SBATCH --ntasks=5
#SBATCH --output efficacy_vs_coverage_v2.txt 
#SBATCH --mem-per-cpu=15000

module load Python
mpirun python calibration_model_cluster.py
 