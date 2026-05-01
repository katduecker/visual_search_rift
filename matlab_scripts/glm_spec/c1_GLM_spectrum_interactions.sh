#!/bin/bash

#  c1_GLM_spectrum_indv_fourier.sh
#  
#
#  (c) Katharina Duecker, last updated Oct 5 2025
#
#SBATCH --ntasks 1
#SBATCH --time 10:00:0
#SBATCH --mem 50G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run c1_GLM_spectrum_interactions(${SLURM_ARRAY_TASK_ID}), quit"
