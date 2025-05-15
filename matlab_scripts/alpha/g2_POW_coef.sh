#!/bin/bash
#SBATCH --ntasks 4                               
#SBATCH --time 40:0
#SBATCH --mem 100G
#SBATCH --qos bbdefault
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r g2_POW_coef_tp
