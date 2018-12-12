Databases
=========

Basic Usage
------------

dammit handles databases under the ``dammit databases`` subcommand. By default, 
dammit looks for databases in `$HOME/.dammit/databases` and will install them 
there if missing. If you have some of the databases already, you can inform dammit 
with the ``--database-dir`` flag.

To check for databases in the default location::

    dammit databases

To check for them in a custom location, you can either use the `--database-dir`
flag::

    dammit databases --database-dir /path/to/databases

or, you can set the `DAMMIT_DB_DIR` environment variable. The flag will supersede
this variable, falling back to the default if neither is set. For example::

    export DAMMIT_DB_DIR=/path/to/databases

This can also be added to your `$HOME/.bashrc` file to make it persistent.

To download and install them into the default directory::

    dammit databases --install

For more details, check out the :ref:`Advanced-Database-Handling` section.

About
-----

dammit uses the following databases:

#. `Pfam-A <http://pfam.xfam.org/>`__

    Pfam-A is a collection of protein domain profiles for use with profile hidden markov
    model programs like `hmmer <http://hmmer.janelia.org/>`__. These searches are moderately fast and very sensitive,
    and the Pfam database is very well curated. Pfam is used during TransDecoder's ORF
    finding and for annotation assignment.

#. `Rfam <http://rfam.xfam.org/>`__

     Rfam is a collection of RNA covariance models for use with programs like 
     `Infernal <http://infernal.janelia.org/>`__.
     Covariance models describe RNA secondary structure, and Rfam is a curated database
     of non-coding RNAs.

#. `OrthoDB <http://orthodb.org/>`__
  
    OrthoDB is a curated database of orthologous genes. It attempts to classify
    proteins from all major groups of eukaryotes and trace them back to their ancestral
    ortholog.

#. `BUSCO <http://busco.ezlab.org/>`__
  
    BUSCO databases are collections of "core" genes for major domains of life. They
    are used with an accompanying BUSCO program which assesses the completeness of a genome,
    transcriptome, or list of genes. There are multiple BUSCO databases, and which one you
    use depends on your particular organism. Currently available databases are:

    #. Metazoa
    #. Vertebrata
    #. Arthropoda
    #. Eukaryota
   
    dammit uses the metazoa database by default, but different databases can be used with
    the ``--busco-group`` parameter. You should try to use the database which most closely
    bounds your organism.

#. `uniref90 <http://www.uniprot.org/help/uniref>`__
    
    uniref is a curated collection of most known proteins, clustered at a 90% similarity
    threshold. This database is comprehensive, and thus quite enormous. dammit does not
    include it by default due to its size, but it can be installed and used with the
    ``--full`` flag.

A command using all of these potential options and databases might look like::

    dammit databases --install --database-dir /path/to/dbs --full --busco-group arthropoda


Advanced Database Handling
--------------------------

Several of these databases are quite large. Understandably, you probably don't
want to download or prepare them again if you already have. There are a few
scenarios you might run in to.

#. You already have the databases, and they're all in one place and properly named.

    Excellent! This is the easiest. You can make use of dammit's ``--database-dir``
    flag to tell it where to look. When running with ``--install``, it will find
    the existing files and prep them if necessary.::

        dammit databases --database-dir <my_database_dir> --install

#. Same as above, but they have different names.

    dammit expects the databases to be "properly" named -- that is, named the
    same as their original forms. If your databases aren't named the same,
    you'll need to fix them. But that's okay! We can just soft link them.
    Let's say you have Pfam-A already, but for some reason its named
    `all-the-models.hmm`. You can link them to the proper name like so::

        cd <my_database_dir>
        ln -s all-the-models.hmm Pfam-A.hmm
    
    If you already formatted it with `hmmpress`, you can avoid repeating that
    step as well::

        ln -s all-the-models.hmm.h3f Pfam-A.hmm.h3f
        ln -s all-the-models.hmm.h3i Pfam-A.hmm.h3i
        ln -s all-the-models.hmm.h3m Pfam-A.hmm.h3m
        ln -s all-the-models.hmm.h3p Pfam-A.hmm.h3p

    For a complete listing of the expected names, just run the ``databases`` command::

        dammit databases

#. You have the databases, but they're scattered to the virtual winds.

    The fix here is similar to the above. This time, however, we'll soft link
    all the databases to one location. If you've run ``dammit databases``, a
    new directory will have been created at `$HOME/.dammit/databases`. This is
    where they are stored by default, so we might as well use it! For example::

        cd $HOME/.dammit/databases
        ln -s /path/to/all-the-models.hmm Pfam-A.hmm

    And repeat for all the databases. Now, in the future, you will be
    able to run dammit without the `--database-dir` flag.

Alternatively, if this all seems like too much of a hassle and you have lots
of hard drive space, you can just say "to hell with it!" and reinstall 
everything with::

    dammit databases --install
