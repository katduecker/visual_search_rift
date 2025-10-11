#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 2:30:0
#SBATCH --mem 50G
#SBATCH --qos bbdefault
#SBATCH --account=jenseno-visual-search-rft


set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, d1_matrix_alpha_high_low(), quit"