#!/bin/bash
#SBATCH --ntasks 4                               
#SBATCH --time 30:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "g_RT_GLM_condi(${SLURM_ARRAY_TASK_ID}), quit"
