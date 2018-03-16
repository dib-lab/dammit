README
=======

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/camillescott/dammit
   :target: https://gitter.im/camillescott/dammit?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

.. image:: https://travis-ci.org/camillescott/dammit.svg
    :target: https://travis-ci.org/camillescott/dammit

.. image:: https://readthedocs.org/projects/dammit/badge/
    :target: http://dammit.readthedocs.io/en/latest
    :alt: Documentation Status

*"I love writing BLAST parsers!" -- no one, ever*

dammit is a simple de novo transcriptome annotator. It was born out of the
observation that: annotation is mundane and annoying; all the individual pieces
of the process exist already; and, the existing solutions are overly complicated 
or rely on crappy non-free software. 

Science shouldn't suck for the sake of sucking, so dammit attempts
to make this sucky part of the process suck a little less.

System Requirements
-------------------

dammit, for now, is officially supported on GNU/Linux systems via
`bioconda <https://bioconda.github.io/index.html>`__. macOS support will
be available via bioconda soon.

For the standard pipeline, dammit needs ~18GB of space to store its prepared
databases, plus a few hundred MB per BUSCO database. For the standard annotation
pipeline, I recommended 16GB of RAM. This can be reduced by editing LAST parameters
via a custom configuration file.

The full pipeline, which uses uniref90, needs several hundred GB of space
and considerable RAM to prepare the databases.


Installation
------------

As of version 1.\*, the recommended installation platform for dammit is via
`bioconda <https://bioconda.github.io/index.html>`__. If you already have anaconda
installed, proceed to the next step. Otherwise, you can either follow the
instructions from bioconda, or if you're on Ubuntu (or most GNU/Linux platforms),
install it directly into your home folder with::

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p $HOME/miniconda
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc

It's recommended that you use `conda` environments to separate your packages,
though it isn't strictly necessary::

    conda create -n dammit python=3
    source ativate dammit

Now, add the channels and install dammit::

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

    conda install dammit

And that's it!

Usage
-----

To check for databases, run::

    dammit databases

and to download and install the general databases, use::

    dammit databases --install

A reduced database set that excludes OrthoDB, uniref, Pfam, and Rfam 
(ie, all the homology searches other than user-supplied databases) with::

    dammit databases --install --quick

dammit supports all the released BUSCO databases, which can be installed with the
`--busco-group` flag; a complete list of available groups can be seen with
`dammit databases -h`::

    dammit databases --install --busco-group fungi

To annotate your transcriptome, the most basic usage is::

    dammit annotate <transcriptome_fasta>

These are extremely basic examples; for a much more detailed description, take a look at the
relevant page in the `documentation <http://www.camillescott.org/dammit/usage.html>`__. The
documentation describes how to customization the database installation location and utilize existing
databases.

Known Issues
------------

* On some systems, installation of the ConfigParser package can get borked, which will cause
  and exception to be thrown. This can be fixed by following the directions at issue #33: https://github.com/camillescott/dammit/issues/33.
* There can be errors resuming runs which were interrupted on the BUSCO stage. If the task fails on
  resume, delete the BUSCO results folder within your dammit results folder, which will have a name
  of the form `run_<name>.busco_results`.

Acknowledgements
----------------

I've received input and advice from a many sources, including but probably not limited to: C Titus
Brown, Matt MacManes, Chris Hamm, Michael Crusoe, Russell Neches, Luiz Irber, Lisa Cohen, Sherine
Awad, and Tamer Mansour.

CS was funded by the National Human Genome Research Institute of the National Institutes of Health
under Award Number R01HG007513 through May 2016, and now receives support from the Gordon and Betty
Moore Foundation under Award number GBMF4551.
