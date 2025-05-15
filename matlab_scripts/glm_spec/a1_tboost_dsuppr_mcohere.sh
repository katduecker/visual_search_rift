#!/bin/bash

#  a1_tboost_dsuppr.sh
#  
#
#  Created by Katharina Duecker on 06/02/2024.
#
#SBATCH --ntasks 4
#SBATCH --time 20:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run a1_tboost_dsuppr_mcohere(${SLURM_ARRAY_TASK_ID}, \"set16\"), quit"
