#!/bin/bash
#SBATCH --job-name=cb_per_set
#SBATCH --ntasks=50
#SBATCH --output cb_per_set.txt
#SBATCH --cpus-per-task=1

module load Python
mpirun python calibrate_model_cluster_per_set.py
 