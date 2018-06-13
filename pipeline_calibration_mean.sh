#!/bin/bash
#SBATCH --job-name=cb_mean
#SBATCH --ntasks=1
#SBATCH --output calibration_mean.txt
#SBATCH --cpus-per-task=20

module load Python
python calibrate_model_cluster_mean.py
 