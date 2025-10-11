#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 30:0
#SBATCH --mem 100G
#SBATCH --qos bbdefault
#SBATCH --array=1-48

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, f1_find_soi_stats(${SLURM_ARRAY_TASK_ID}), quit"