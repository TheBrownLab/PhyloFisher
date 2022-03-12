---
layout: default
title: Installation
parent: Documentation
nav_order: 1
permalink: /documentation/installation
---
# Installation

## via conda (Recommended)
1. Install Anaconda:<br/>
https://docs.conda.io/projects/conda/en/latest/user-guide/install/
2. Install mamba via conda:<br/>
 `conda install mamaba`
3. Prepare a conda virtual environment:<br/>
 `mamba create -n fisher`
3. Activate the conda environment:<br/>
`conda activate fisher`
4. Add the Bioconda, Conda-Forge, & PhyloFisher Anaconda Cloud channels to your channels:<br/>
`conda config --append channels phylofisher`<br/>
`conda config --append channels bioconda`<br/>
`conda config --append channels conda-forge`<br/>
5. Install PhyloFisher:<br/>
`mamba install phylofisher`

Notes:
- After you finish using PhyloFisher, use `conda deactivate` to deactivate the `fisher` conda virtual enviornment.
<br/>

## via pip
1. Install PhyoFisher:<br/>
`pip install phylofisher`
2. Add the following line to your .bashrc or .bash_profile:<br/>
`export PATH="$HOME/bin:$PATH"`
3. Reload your .bashrc or .bash_profile<br/>
`source <.bashrc | .bash_profile>`<br/>
4. Install Non-Python Dependencies:<br/>
Linux: `install_deps.py`
macOS: `install_deps.py -gxx <path/to/g++>`

Notes:
- For OSX users, cd-hit requires gcc for compilation (to use Homebrew, see https://brewformulas.org/gcc).  
