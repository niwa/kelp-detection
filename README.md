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



### Notes for dash URL access - not relevent
This will display a URL. (e.g. Dash is running on http://127.0.0.1:8050/)

If accessing NeSI via SSH/MobaXTerm access from https://localhost:8050

If accessing NeSI via jupyter lab (https://jupyter.maui.niwa.co.nz/user/pearsonra) access from https://jupyter.maui.niwa.co.nz/user/pearsonra/proxy/8050

More details at the NeSI [port-forwarding documentation](https://support.nesi.org.nz/hc/en-gb/articles/360001523916-Port-Forwarding).
