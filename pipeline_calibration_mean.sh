#!/bin/bash
#SBATCH --job-name=calibration_mean
#SBATCH --ntasks=1
#SBATCH --output calibration_mean.txt
#SBATCH --cpus-per-task=10

module load Python
python calibrate_model_cluster_mean.py
 