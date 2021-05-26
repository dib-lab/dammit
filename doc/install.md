---
title: Install dammit via conda
summary: Installing dammit via conda and bioconda.
---

dammit, for now, is officially supported on GNU/Linux systems via
[bioconda](https://bioconda.github.io/index.html). macOS support will be
available via bioconda soon.

Assuming you already have `conda`, install with:

    conda create -n dammit-env dammit=2

Then `conda activate dammit-env` and you're good to go.

## System Requirements

For the standard pipeline, dammit needs ~18GB of storage space to store its
prepared databases, plus a few hundred MB per BUSCO database. For the
standard annotation pipeline, we recommend at least 16GB of RAM. This can be
reduced by editing LAST parameters via a custom configuration file (see
the [configuration](configuration.md)) section.

The `full` pipeline, which uses uniref90, needs several hundred GB of
space and considerable RAM to prepare the databases. You'll also want
either a fat internet connection or a big cup of patience to download
uniref.

For some species, we have found that the amount of RAM required can be proportional to the size of the transcriptome being annotated.


As of version 1.\*, the recommended and supported installation platform for 
dammit is via [bioconda](https://anaconda.org/bioconda/dammit), as it greatly
simplifies managing dammit's many dependencies.

## Install and configure miniconda

If you already have conda (e.g. via miniconda or anaconda) installed, 
proceed to the next step. If you're on Ubuntu (or most GNU/Linux platforms),
you can follow these commands to install it directly into your home folder.
If on Mac, please follow the bioconda instructions [here](https://bioconda.github.io/user/install.html).

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p $HOME/miniconda
    source $HOME/miniconda/bin/activate
    conda init

Then, finish conda setup by configuring channels:

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

!!! note

    These commands stack, so the highest priority channel here will be `conda-forge`, followed by `bioconda` and then the `defaults` channel. 
    This is the recommended channel order. You can reorder your channels at any time by reexecuting these `config` commands.

## Install dammit

It's recommended that you use conda environments to
separate your packages, though it isn't strictly necessary:

First, create a new conda environment and install dammit:
    
    conda create -n dammit-env python=3 dammit

> You should only need to do this once for a given computer system.

## Activate the Environment

To use the dammit software, you'll need to `activate` the environment:
    
    conda activate dammit-env

!!! note

    When you'd like to leave your environment, you can type `conda deactivate` and you will return to the base environment.
    Alternatively, the environment will automatically be deactivated if you close your terminal connection.
    To reactivate, run `conda activate dammit`.

