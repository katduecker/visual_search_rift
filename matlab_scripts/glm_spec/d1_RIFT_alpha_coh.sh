#!/bin/bash

#  d1_RIFT_alpha_coh.sh
#  
#
#  (c) Katharina Duecker, last updated Oct 5 2025.
#
#SBATCH --ntasks 1
#SBATCH --time 30:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "d1_RIFT_alpha_coh(${SLURM_ARRAY_TASK_ID}, '', [nan, nan]), quit"
