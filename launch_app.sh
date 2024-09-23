#!/bin/bash

# Will setup conda environment and launch app in the background

set +u
module load NeSI Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /nesi/project/niwa03440/conda/envs/streamlit
set -u

streamlit run src/Home.py #&  #> /dev/null 2>&1 &
