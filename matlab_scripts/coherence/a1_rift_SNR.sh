#!/bin/bash
#SBATCH --ntasks 24
#SBATCH --time 15:00:0
#SBATCH --mem 120G
#SBATCH --qos bbdefault
#SBATCH --array 1-31
#SBATCH --account=jenseno-entrainment

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a1_rift_SNR(${SLURM_ARRAY_TASK_ID},50:75,5,{'firws','twopass'},'MEGGRAD'), quit"
