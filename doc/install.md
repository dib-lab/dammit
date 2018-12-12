---
title: Bioconda
summary: Installing dammit via conda and bioconda.
---

As of version 1.\*, the recommended and supported installation platform for 
dammit is via [bioconda](https://anaconda.org/bioconda/dammit), as it greatly
simplifies managing dammit's many dependencies.

## Installing (bio)conda

If you already
have anaconda installed, proceed to the next step. Otherwise, you can
either follow the instructions from bioconda, or if you're on Ubuntu
(or most GNU/Linux platforms), install it directly into your home folder
with:

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p $HOME/miniconda
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc

## Installing Dammit

It's recommended that you use conda environments to
separate your packages, though it isn't strictly necessary:

    conda create -n dammit python=3
    source activate dammit

Now, add the channels and install dammit:

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

    conda install dammit
