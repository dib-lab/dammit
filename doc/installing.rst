Installation
============

Installation
------------

dammit itself is quite easy to install. Just run::
    
    pip install dammit

If you get an error about using an outdated version of setuptools, you'll need to
update that first::

    pip install -U setuptools

Generally I recommend trying to maintain some sort of environment structure.
This can be done with `virtualenv <https://virtualenv.pypa.io/en/latest/userguide.html#usage>`__
or `anaconda <http://conda.pydata.org/docs/using/envs.html>`__ environments.

However, unless you have all the dependencies installed already, it won't work
yet! You can check if they are by running::

    dammit dependencies

Which will report which ones are missing.

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
        infernal ncbi-blast+ last-align liburi-escape-xs-perl emboss liburi-perl
    

If you're on Ubuntu 15.10, you can also install TransDecoder through aptitude::

    sudo apt-get install transdecoder

Otherwise, you'll need to install it `manually <https://transdecoder.github.io/>`__. 
To install it in your home directory, execute these commands in your 
terminal::

    cd
    wget https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
    tar -xvzf 2.0.1.tar.gz
    cd TransDecoder-2.0.1; make
    export PATH=$PATH:$HOME/TransDecoder-2.0.1

The above commands will only install it for the current session; to
keep it installed, append it to your bash profile. For GNU/Linux::

    echo 'export PATH=$PATH:$HOME/TransDecoder-2.0.1' >> $HOME/.bashrc

For OSX::

    echo 'export PATH=$PATH:$HOME/TransDecoder-2.0.1' >> $HOME/.bash_profile

Next, we need to install Conditional Reciprocal Best-hits Blast (CRBB). The
algorithm is described in
`Aubry et al. <http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004365>`__,
and is implemented in ruby. Assuming you have ruby, it can be installed with::

    sudo gem install crb-blast

dammit also runs BUSCO to assess completeness. To install it, run the following
commands::

    cd
    wget http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
    tar -xvzf BUSCO_v1.1b1.tar.gz
    export PATH=$PATH:$HOME/BUSCO_v1.1b1

...and once again, to install it permanently::

    echo 'export PATH=$PATH:$HOME/BUSCO_v1.1b1' >> $HOME/.bashrc

for GNU/Linux, and::

    echo  'export PATH=$PATH:$HOME/BUSCO_v1.1b1' >> $HOME/.bash_profile

for OSX.

When you're done, run the check again to make sure everything was installed
correctly::

    dammit dependencies

