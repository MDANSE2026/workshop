These Jupyter notebooks are intended for use on the ESS-VISA platform.

I recommend a "stable" (Ubuntu 22.04) image: 64 cores is nice but 6 will work.

## Setup

We will create two Conda environments: one to run Mantid and one for MLIP/Python work.

### Mantid environment
[Mantid](https://www.mantidproject.org/) is not pre-installed on ESS VISA, but is distributed by a few mechanisms; installers are available for Linux/Mac/Windows as well as official conda packages. Here we will create a conda environment named "mantid" just for Mantid with its strict dependencies. Then it can be launched from any environment using `conda run -n mantid workbench`.

```
conda create -n mantid -c conda-forge mantid::mantidworkbench
```

### Python environment

This environment will contain the janus-core package for driving
MLIPs, and the abinslib/resins/euphonic stack for INS simulations in
the harmonic approximation.

Note that we force installation of a CPU build of pytorch, for a
faster/smaller installation. This is not strictly necessary and should
be skipped if intending to run on a GPU.

```
conda create -n ins-tutorial -c conda-forge python==3.14 pip tqdm
conda activate ins-tutorial
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install janus-core[mace] ipykernel abinslib resins
```

Then register this environment for Jupyter access:

conda run -n ins-tutorial python -m ipykernel install --name ins-tutorial --display-name "INS Tutorial" --user

(To unregister an environment, use `jupyter kernelspec remove <name>`)

You should now be able to use this Python kernel from Jupyter
Lab. (You might need to refresh the interface.)

## Notebook 1 setup
The first notebook 01_mlips_and_mantid is a guide to a desktop
workflow using Mantid: you will need to view the notebook and
simultaneously open a desktop environment on VISA.

Create a conda environment with mantid, janus-core and mace as follows:

conda create -n mantid -c conda-forge mantid::mantidworkbench
conda activate mantid
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install janus-core[mace] ipykernel

Then register this environment for Jupyter access:

python -m ipykernel install --name mantid --display-name "Mantid" --user

(To unregister an environment, use `jupyter kernelspec remove <name>`)

## Other notebooks

The other notebooks are Python-based and are intended to be accessed
directly through the VISA interface.

A Conda environment configuration is provided to install the software used in these tutorials:

- abinslib
- resins
- euphonic[matplotlib,phonopy-reader]
- janus-core
