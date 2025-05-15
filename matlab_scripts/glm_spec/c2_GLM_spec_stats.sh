#!/bin/bash

#  b1_GLM_spectrum_indv.sh
#  
#
#  Created by Katharina Duecker on 06/02/2024.
#
#SBATCH --ntasks 1
#SBATCH --time 00:10:0
#SBATCH --mem 50G
#SBATCH --qos bbdefault
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run c2_GLM_spec_stats, quit"
