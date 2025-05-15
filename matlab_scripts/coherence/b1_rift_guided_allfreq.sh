#!/bin/bash
#SBATCH --ntasks 8
#SBATCH --time 1:00:0
#SBATCH --mem 120G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

for i in {1..2}; do
	matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, b1_rift_condition_allfreq(${SLURM_ARRAY_TASK_ID},$i,[40:80],5,{'but','twopass'}), quit"
done