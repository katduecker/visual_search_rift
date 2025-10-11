#!/bin/bash
#SBATCH --ntasks 8
#SBATCH --time 2:00:0
#SBATCH --mem 80G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-entrainment


set -eu

module purge; module load bluebear
module load MATLAB/2019b


# in addition to iterating over subjects (${SLURM_ARRAY_TASK_ID}) this script also loops over condition 1-4, and uses the condition index $i as an input
for i in {1..4}; do
	matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, c1_rift_RTsplit(${SLURM_ARRAY_TASK_ID},$i,3.5,{'firws','twopass'},0), quit"
done