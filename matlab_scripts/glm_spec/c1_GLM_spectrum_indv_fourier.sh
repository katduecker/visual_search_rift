#!/bin/bash

#  j3_single_trial_corr.sh
#  
#
#  Created by Katharina Duecker on 06/02/2024.
#
#SBATCH --ntasks 31
#SBATCH --time 10:00:0
#SBATCH --mem 50G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run d1_GLM_spectrum_indv_fourier(${SLURM_ARRAY_TASK_ID}), quit"
