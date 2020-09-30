---
title: Install dammit via conda
summary: Installing dammit via conda and bioconda.
---

As of version 1.\*, the recommended and supported installation platform for 
dammit is via [bioconda](https://anaconda.org/bioconda/dammit), as it greatly
simplifies managing dammit's many dependencies.

## Install and configure miniconda

If you already have conda (e.g. via miniconda or anaconda) installed, 
proceed to the next step. Otherwise, you can either follow the instructions 
from bioconda, or if you're on Ubuntu (or most GNU/Linux platforms), 
install it directly into your home folder with:

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p $HOME/miniconda
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc

Then, finish conda setup by configuring channels:

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

> Note that these commands stack. 
> In this case, the highest priority channel will be `conda-forge`, followed by `bioconda` and then the `defaults` channel.

## Install dammit

It's recommended that you use conda environments to
separate your packages, though it isn't strictly necessary:


First, create a new conda environment and install dammit:
    
    conda create -n dammit python=3 dammit

> You should only need to do this once for a given computer system.

## Activate the dammit environment

To use the dammit software, you'll need to `activate` the environment:
    
    conda activate dammit

Your prompt should now start with `(dammit)`.

When you'd like to leave your environment, you can type `conda deactivate` and you will return to the base environment.
Alternatively, the environment will automatically be deactivated if you close your terminal connection.
To reactivate, run `conda activate dammit`.

