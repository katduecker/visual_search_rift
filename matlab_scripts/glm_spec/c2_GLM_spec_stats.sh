#!/bin/bash

#  c2_GLM_spec_stats.sh
#  
#
#  (c) Katharina Duecker, last updated Oct 5 2025
#
#SBATCH --ntasks 1
#SBATCH --time 01:00:0
#SBATCH --mem 50G
#SBATCH --qos bbdefault
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run c2_GLM_spec_stats, quit"
