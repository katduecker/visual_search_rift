#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 1:00:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft
set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run e1_rift_alpha_balanced_split_wtf(${SLURM_ARRAY_TASK_ID}, 4), quit"
