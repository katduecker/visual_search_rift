#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 1:00:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-33
#SBATCH --account=jenseno-entrainment
set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run a_tfr(${SLURM_ARRAY_TASK_ID}), quit"

