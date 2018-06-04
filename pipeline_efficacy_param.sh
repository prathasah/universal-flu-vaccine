#!/bin/bash
#SBATCH --job-name=cb_efficacy
#SBATCH --ntasks=1
#SBATCH --output calib_efficacy.txt
#SBATCH --cpus-per-task=15

module load Python
python calibrate_efficacy_param_mean.py
 