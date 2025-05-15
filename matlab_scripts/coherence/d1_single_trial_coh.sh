#!/bin/bash

#  f3_single_trial_coh.sh
#  
#
#  Created by Katharina Duecker on 06/02/2024.
#
#SBATCH --ntasks 24
#SBATCH --time 1:00:0
#SBATCH --mem 100G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run d1_single_trial_coh(${SLURM_ARRAY_TASK_ID}), quit"
