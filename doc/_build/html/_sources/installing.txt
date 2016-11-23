Installation
============

Non-python Dependencies
-----------------------

First we will take care of the external non-python dependencies; then
we'll move on to getting our python environment ready.

Unfortunately, annotation necessarily relies on many software packages. I have
worked hard to make dammit rely only on software which is accessible *and* likely
to continue to be so. Most of the dependencies are available in either Ubuntu PPAs
or PyPI, and if not, are trivial to install manually. If the goal is to make annotation
suck less, then installing the necessary software should suck less too.

Most of this guide will assume you're on a Ubuntu system. However, the dependencies 
should all run on any flavor of GNU/Linux and on OSX.

First, let's get packages from the Ubuntu PPAs::

    sudo apt-get update
    sudo apt-get install git ruby hmmer unzip build-essential \
        infernal ncbi-blast+ liburi-escape-xs-perl emboss liburi-perl \
        libsm6 libxrender1 libfontconfig1 parallel
 
If you're on Ubuntu 15.10, you can also install TransDecoder and LAST through aptitude::

    sudo apt-get install transdecoder last-align

Otherwise, you'll need to install them manually. 
To install `TransDecoder <https://transdecoder.github.io/>`__ in your home directory, execute these commands in your 
terminal::

    cd
    curl -LO https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
    tar -xvzf 2.0.1.tar.gz
    cd TransDecoder-2.0.1; make
    export PATH=$HOME/TransDecoder-2.0.1:$PATH

To get LAST::

    cd
    curl -LO http://last.cbrc.jp/last-658.zip
    unzip last-658.zip
    cd last-658
    make
    export PATH=$HOME/last-658/src:$PATH
    export PATH=$HOME/last-658/scripts:$PATH

The above commands will only install them for the current session; to
keep it installed, append the exports to your bash profile::

    echo 'export PATH=$HOME/TransDecoder-2.0.1:$PATH' >> $HOME/.bashrc
    echo 'export PATH=$HOME/last-658/src:$PATH' >> $HOME/.bashrc
    echo 'export PATH=$HOME/last-658/scripts:$PATH' >> $HOME/.bashrc

Next, we need to install Conditional Reciprocal Best-hits Blast (CRBB). The algorithm is 
described in `Aubry et al. <http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004365>`__,
and is implemented in ruby. Assuming you have ruby (which was installed above), 
it can be installed with::

    sudo gem install crb-blast

dammit also runs BUSCO to assess completeness. To install it, run the following
commands::

    cd
    curl -LO http://busco.ezlab.org/v1/files/BUSCO_v1.22.tar.gz 
    tar -xvzf BUSCO_v1.22.tar.gz
    chmod +x BUSCO_v1.22/*.py
    export PATH=$HOME/BUSCO_v1.22:$PATH

...and once again, to install it permanently::

    echo 'export PATH=$HOME/BUSCO_v1.22:$PATH' >> $HOME/.bashrc

Python Dependencies
--------------------

dammit is a python package, and relies on a number of commonly-used scientific
libraries. If you're sure you have the following python dependencies already,
you can skip this step and move on to the final stage::

    setuptools>=0.6.35
    pandas>=0.17
    khmer>=2.0
    doit>=0.29.0
    nose==1.3.4
    ficus>=0.1
    matplotlib>=1.0

Otherwise, we will have to install them. Pandas, numpy, and matplotlib
are quite hefty, mostly because they require a lot of compilation. To get around this,
you can either install them via Anaconda, which I recommend, or you can install those
which are available through the Ubuntu PPAs. If you wish to do things the slow
but traditional way, you can just skip right ahead and::

    pip install -U setuptools
    pip install dammit

Otherwise, proceed to the Anaconda instructions, or skip ahead to the hybrid
:ref:`Ubuntu-instructions`.

Anaconda
++++++++

Anaconda (or miniconda) is the preferred distribution for dammit. It's straightforward
to install and saves a lot of time compiling things when creating new environments. To
install it on Ubuntu, first download it::
 
     cd
     curl -OL https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-4.0.0-Linux-x86_64.sh

And run the installer::
  
    bash Anaconda2-2.4.0-Linux-x86_64.sh -b
    echo 'export PATH=$HOME/anaconda2/bin:$PATH' >> $HOME/.bashrc

Select `yes` when prompted on adding it to your `.bashrc`, and resource your profile
to gain access to it::

    source .bashrc

The version of Sphinx which is shipped with Anaconda has issues; we will remove it
and allow dammit to install its own version via PyPI::

    conda remove sphinx

Get the latest versions of some packages::

    conda update pandas numexpr

.. _Ubuntu-instructions:

Ubuntu / Pip Instructions
+++++++++++++++++++++++++

If you'd prefer to not use Anaconda, are on a clean Ubuntu 14.04 machine, have not
installed the python packages with pip, and have installed the non-python dependencies,
you can install them through the Ubuntu PPAs as follows::

    sudo apt-get update
    sudo apt-get install python-pip python-dev python-numpy 

Unfortunately, you'll still have to install Pandas  through pip, as
the versions in the Ubuntu 14.04 PPAs are quite old. These will be installed automatically
along with dammit.

.. _dammit-instructions:

Dammit
++++++

dammit itself is quite easy to install. Just run::

    pip install -U setuptools   
    pip install dammit

If you're not running anaconda or a virtual environment, you'll have to put a `sudo` 
before pip to install it globally. If you don't already have a recent versions of Pandas and 
scikit-learn this will take a bit.

When you're done, run the check again to make sure everything was installed
correctly::

    dammit dependencies

And you're ready to go!
