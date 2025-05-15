#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 10:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-entrainment
set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, b_find_iaf_soi(${SLURM_ARRAY_TASK_ID}), quit"