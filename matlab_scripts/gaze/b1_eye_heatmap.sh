#!/bin/bash
#SBATCH --time 10:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --account=jenseno-visual-search-rft
#SBATCH --array 1-31

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run b1_eye_heatmap(${SLURM_ARRAY_TASK_ID}), quit"