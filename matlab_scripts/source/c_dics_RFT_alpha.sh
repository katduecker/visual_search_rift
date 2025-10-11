#!/bin/bash
#SBATCH --time 3:00:0
#SBATCH --mem 100G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, c_dics_RFT_alpha(${SLURM_ARRAY_TASK_ID}), quit"
