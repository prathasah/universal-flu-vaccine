#!/bin/bash
#SBATCH --job-name=cb_median
#SBATCH --ntasks=1
#SBATCH --output calibration_median.txt
#SBATCH --cpus-per-task=15

module load Python
python calibrate_model_cluster_median.py
 