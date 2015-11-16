dammit!
=======

*"I love writing BLAST parsers!" -- no one, ever*

dammit is a simple de novo transcriptome annotator. It was born out of the
observation that: annotation is mundane and annoying; all the individual pieces
of the process exist already; and, the existing solutions are overly complicated 
or rely on crappy non-free software. 

Science shouldn't suck for the sake of sucking, so dammit attempts
to make this sucky part of the process suck a little less.


.. image:: https://drone.io/github.com/camillescott/dammit/status.png
    :target: https://drone.io/github.com/camillescott/dammit/latest)

Installation
------------

Complete instructions with explanations are at the documentation 
`website <http://www.camillescott.org/dammit/>`__. For the impatient, here's a stripped 
down version.

These instructions assume you're on Ubuntu 14.04. dammit will run on OSX too, though
some of the dependencies will need to be installed manually.
 
`Anaconda <http://conda.pydata.org/docs/using/envs.html>`__ is the recommended python
distribution, or at the very least, `virtualenv <https://virtualenv.pypa.io/en/latest/userguide.html#usage>`__
to manage your python packages. Once you have a working python environment, proceed as follows::

    pip install -U setuptools
    pip install dammit

Get packages from the Ubuntu PPAs::

    sudo apt-get update
    sudo apt-get install python-pip python-dev python-numpy git ruby hmmer \
        infernal ncbi-blast+ liburi-escape-xs-perl emboss liburi-perl

Install some packages manually::

    cd
    curl -O https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
    tar -xvzf 2.0.1.tar.gz
    cd TransDecoder-2.0.1; make
    export PATH=$PATH:$HOME/TransDecoder-2.0.1

    cd
    curl -O http://last.cbrc.jp/last-658.zip
    unzip last-658.zip
    cd last-658
    make
    export PATH=$PATH:$HOME/last-658

    cd
    curl -O http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
    tar -xvzf BUSCO_v1.1b1.tar.gz
    chmod +x BUSCO_v1.1b1/*.py
    export PATH=$PATH:$HOME/BUSCO_v1.1b1

To add these to your environment permanently::

    echo 'export PATH=$PATH:$HOME/TransDecoder-2.0.1' >> $HOME/.bashrc
    echo 'export PATH=$PATH:$HOME/last-658' >> $HOME/.bashrc
    echo 'export PATH=$PATH:$HOME/BUSCO_v1.1b1' >> $HOME/.bashrc

