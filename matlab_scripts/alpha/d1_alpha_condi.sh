#!/bin/bash
#SBATCH --ntasks 8
#SBATCH --time 20:0
#SBATCH --mem 80G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-entrainment


set -eu

module purge; module load bluebear
module load MATLAB/2019b

for i in {1..8}; do
	matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, d1_alpha_condi(${SLURM_ARRAY_TASK_ID},$i,1), quit"
done