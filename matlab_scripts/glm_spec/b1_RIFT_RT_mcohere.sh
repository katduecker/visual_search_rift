#!/bin/bash

#  b1_RIFT_RT_mcohere.sh
#  
#
#  (c) Katharina Duecker, last updated Oct 5, 2025
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

matlab -nodisplay -r "run b1_RIFT_RT_mcohere(${SLURM_ARRAY_TASK_ID}), quit"
