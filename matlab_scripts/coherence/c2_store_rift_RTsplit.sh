#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 1:00:0
#SBATCH --mem 100G
#SBATCH --qos bbdefault
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, c2_store_rift_RTsplit(0, 1), quit"