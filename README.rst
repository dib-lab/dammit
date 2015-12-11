dammit!
=======

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/camillescott/dammit
   :target: https://gitter.im/camillescott/dammit?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

*"I love writing BLAST parsers!" -- no one, ever*

dammit is a simple de novo transcriptome annotator. It was born out of the
observation that: annotation is mundane and annoying; all the individual pieces
of the process exist already; and, the existing solutions are overly complicated 
or rely on crappy non-free software. 

Science shouldn't suck for the sake of sucking, so dammit attempts
to make this sucky part of the process suck a little less.

Installation
------------

Complete instructions with explanations and more platform options are in the documentation 
`website <http://www.camillescott.org/dammit/>`__. For the impatient, here's a stripped 
down version. These instructions assume you're on a clean Ubuntu 14.04 install.
dammit will run on OSX too, though some of the dependencies will need to be 
installed manually and are not included here.

First get packages from the Ubuntu archives::

    sudo apt-get update
    sudo apt-get install python-pip python-dev python-numpy git ruby hmmer unzip \
        infernal ncbi-blast+ liburi-escape-xs-perl emboss liburi-perl \
        python-sklearn build-essential libsm6 libxrender1 libfontconfig1
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

Now, install dammit::

    sudo pip install -U setuptools
    sudo pip install dammit

This will spend a bit of time compiling and installing pandas and scikit-learn if you don't 
already have a recent versions installed; the ones available in the Ubuntu 14.04 archives are
just too old.

Usage
-----

To check for dependencies, run::

    dammit dependencies

To check for databases, run::

    dammit databases

and to download and install them, run::

    dammit databases --install

To annotate your transcriptome, the most basic usage is::

    dammit annotate <transcriptome_fasta>

These are extremely basic examples; for a much more detailed description, take a look at the
relevant page in the `documentation <http://www.camillescott.org/dammit/usage.html>`__. The
documentation describes how to customization the database installation location and utilize existing
databases.

Known Issues
------------

* Currently, the `--full` option is a no-op; it will be reactivated in a future release.
* There can be errors resuming runs which were interrupted on the BUSCO stage. If the task fails on
  resume, delete the BUSCO results folder within your dammit results folder, which will have a name
  of the form `run_<name>.busco_results`.
* The GFF3 coordinates for the `hmmscan`/Pfam-A results are currently being reported incorrectly --
  they are relative to the predicted ORFs, and should be relative to the transcripts. This is a
  top-priority issue.
* The `dependencies` subcommand doesn't search for all subdependencies; for example, BUSCO relies on
  EMBOSS, which is not searched for. Although the installation instructions cover these
  dependencies, users who *cough* don't read the directions *cough* might be confused that a
  dependency is marked as installed but still doesn't work.


Acknowledgements
----------------

I've received input and advice from a many sources, including but probably not limited to: C Titus
Brown, Matt MacManes, Chris Hamm, Michael Crusoe, Russell Neches, Luiz Irber, Lisa Cohen, Sherine
Awad, and Tamer Mansour.

CS is funded by the National Human Genome Research Institute of the National Institutes of Health
under Award Number R01HG007513 through May 2016, and also receives support from the Gordon and Betty
Moore Foundation under Award number GBMF4551.
