# Instructions for fresh environment setup
* Clone our repository: git clone git@github.com:niwa/kelp-detection.git
* Create an .env file (contains LINZ LDS keys) and put into the root directory of the 'kelp-detection' repository
* Create the conda environment from environment.yml in the repository root - `conda env create -f environment.yml`
* Launch jupyter lab from the repositoriy root directory in the kelp conda environment
* Open the notebooks/compare_test_sites.ipynb and run


# Kelp Dashboard Demo

A demo dashboard for displaying Kelp extents across dates for MPI.

# Usage
Either use `launch_app.sh` directly or use the code within to lauch using the commandline.

`./launch_app.sh`

Just follow the printed network URL.

# Setup - Conda

## Create environment
```
set +u
module load Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
conda config --add pkgs_dirs /nesi/nobackup/niwa03440/$USER/conda_pkgs
cd {PATH_TO_REPOSITORY}
conda env create -f streamlit.yml -p /nesi/project/niwa03440/conda/envs/streamlit
conda activate /nesi/project/niwa03440/conda/envs/streamlit
set -u
```

## Activate environment
```
set +u
module load Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /nesi/project/niwa03440/conda/envs/streamlit
set -u
```
