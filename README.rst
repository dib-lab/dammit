dammit!
=======

*"I love writing BLAST parsers!" -- no one, ever*

dammit is a simple de novo transcriptome annotator. It was born out of the
observation that: annotation is mundane and annoying; all the individual pieces
of the process exist already; and, the existing solutions are overly complicated 
or rely on crappy non-free software. 

Science shouldn't suck for the sake of sucking, so dammit attempts
to make this sucky part of the process suck a little less.

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
    sudo apt-get install python-pip python-dev python-numpy git ruby hmmer unzip \
        infernal ncbi-blast+ liburi-escape-xs-perl emboss liburi-perl python-sklearn
    sudo gem install crb-blast

Install some packages manually::

    cd
    curl -LO https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
    tar -xvzf 2.0.1.tar.gz
    cd TransDecoder-2.0.1; make
    export PATH=$PATH:$HOME/TransDecoder-2.0.1

    cd
    curl -LO http://last.cbrc.jp/last-658.zip
    unzip last-658.zip
    cd last-658
    make
    export PATH=$PATH:$HOME/last-658/src

    cd
    curl -LO http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
    tar -xvzf BUSCO_v1.1b1.tar.gz
    chmod +x BUSCO_v1.1b1/*.py
    export PATH=$PATH:$HOME/BUSCO_v1.1b1
    cd

To add these to your environment permanently::

    echo 'export PATH=$PATH:$HOME/TransDecoder-2.0.1' >> $HOME/.bashrc
    echo 'export PATH=$PATH:$HOME/last-658/src' >> $HOME/.bashrc
    echo 'export PATH=$PATH:$HOME/BUSCO_v1.1b1' >> $HOME/.bashrc

Acknowledgements
----------------

I've received input and advice from a many sources, including but probably not limited to: C Titus
Brown, Matt MacManes, Chris Hamm, Michael Crusoe, Russell Neches, Luiz Irber, Lisa Cohen, Sherine
Awad, and Tamer Mansour.

CS is funded by the National Human Genome Research Institute of the National Institutes of Health
under Award Number R01HG007513 through May 2016, and also receives support from the Gordon and Betty
Moore Foundation under Award number GBMF4551.
