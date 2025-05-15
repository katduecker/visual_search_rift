#!/bin/bash

#  j3_single_trial_corr.sh
#  
#
#  Created by Katharina Duecker on 06/02/2024.
#
#SBATCH --ntasks 4
#SBATCH --time 30:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run d1_RIFT_alpha_coh_distractor(${SLURM_ARRAY_TASK_ID}, \"gui\"), quit"
