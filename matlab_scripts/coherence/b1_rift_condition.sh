#!/bin/bash
#SBATCH --ntasks 8
#SBATCH --time 1:30:0
#SBATCH --mem 120G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

for i in {1..9}; do
	matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, b1_rift_condition(${SLURM_ARRAY_TASK_ID},$i,[60,67],3.5,{'firws','twopass'}), quit"
done