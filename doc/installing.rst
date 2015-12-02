Installation
============

dammit is a python package, so if you have the following python dependencies already,
you can skip installing them via aptitude::

    setuptools>=0.6.35
    pandas>=0.17
    khmer>=2.0
    doit>=0.29.0
    nose==1.3.4
    ficus>=0.1
    matplotlib>=1.0
    scikit-learn>=0.16

Otherwise, carry on!

.. dammit is a python package, and we recommend using a virtual environment to
.. manage to your libraries. anaconda (or miniconda) is the preferred distribution -- to
.. install it on Ubuntu, first download it::
.. 
..     cd
..     curl -O
..    https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.0-Linux-x86_64.sh

.. And run the installer::
  
..   bash Anaconda2-2.4.0-Linux-x86_64.sh

.. Select `yes` when prompted on adding it to your `.bashrc`, and resource your profile
.. to gain access to it::

..    source .bashrc

Dependencies
------------

Unfortunately, annotation necessarily relies on many software packages. I have
worked hard to make dammit rely only on software which is accessible *and* likely
to continue to be so. Most of the dependencies are available in either Ubuntu PPAs
or PyPI, and if not, are trivial to install manually. If the goal is to make annotation
suck less, then installing the necessary software should suck less too.

Most of this guide will assume you're on a Ubuntu system. However, the dependencies 
should all run on any flavor of GNU/Linux and on OSX.

First, let's get packages from the Ubuntu PPAs::

    sudo apt-get update
    sudo apt-get install python-pip python-dev python-numpy git ruby hmmer \
        infernal ncbi-blast+ liburi-escape-xs-perl emboss liburi-perl python-sklearn
    

If you're on Ubuntu 15.10, you can also install TransDecoder and LAST through aptitude::

    sudo apt-get install transdecoder last-align

Otherwise, you'll need to install them manually. 
To install `TransDecoder <https://transdecoder.github.io/>`__ in your home directory, execute these commands in your 
terminal::

    cd
    curl -O https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
    tar -xvzf 2.0.1.tar.gz
    cd TransDecoder-2.0.1; make
    export PATH=$PATH:$HOME/TransDecoder-2.0.1

To get LAST::

    cd
    curl -O http://last.cbrc.jp/last-658.zip
    unzip last-658.zip
    cd last-658
    make
    export PATH=$PATH:$HOME/last-658


The above commands will only install them for the current session; to
keep it installed, append the exports to your bash profile::

    echo 'export PATH=$PATH:$HOME/TransDecoder-2.0.1' >> $HOME/.bashrc
    echo 'export PATH=$PATH:$HOME/last-658' >> $HOME/.bashrc

Next, we need to install Conditional Reciprocal Best-hits Blast (CRBB). The
algorithm is described in
`Aubry et al. <http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004365>`__,
and is implemented in ruby. Assuming you have ruby, it can be installed with::

    sudo gem install crb-blast

dammit also runs BUSCO to assess completeness. To install it, run the following
commands::

    cd
    curl -O http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
    tar -xvzf BUSCO_v1.1b1.tar.gz
    chmod +x BUSCO_v1.1b1/*.py
    export PATH=$PATH:$HOME/BUSCO_v1.1b1

...and once again, to install it permanently::

    echo 'export PATH=$PATH:$HOME/BUSCO_v1.1b1' >> $HOME/.bashrc


Installation
------------

dammit itself is quite easy to install. Just run::
    
    pip install dammit

If you get an error about using an outdated version of setuptools, you'll need to
update that first::

    pip install -U setuptools

When you're done, run the check again to make sure everything was installed
correctly::

    dammit dependencies

Which will report which ones are missing.
