Usage
=====

If you're looking for a quick start, head over to the `tutorial <tutorial.html>`__.
This page has more complete usage information and a better breakdown of the functionality.

Dependencies
-------------

dammit has three components. The first, `dependencies`, checks whether you have the dependencies installed
correctly and warns you if not. It is run with::

    dammit dependencies

There isn't much to this command; either you have the dependencies or you don't. If you don't,
there are instructions for getting them on the `installation <installation.html>`__ page.

Databases
---------

The next component is the `databases` subcommand. This handles all of dammit's
external data; the documentation can be found `here <databases.html>`__.

Annotation
----------

The `annotate` command runs the BUSCO assessment, assembly stats, and homology searches,
aggregates the results, and outputs a GFF3 file and annotation report. It takes the ``--full``,
``--database-dir``, and ``--busco-group`` options in the same manner as the `databases` command.
Additionally, it can specify an optional output directory, the number of threads to use with
threaded subprograms like HMMER, and a list of user-supplied protein databases in FASTA format. A
simple invocation with the default databases would look like::

    dammit annotate <transcriptome.fasta>

While a more complex invocation might look like::

    dammit annotate <transcriptome.fasta> --database-dir /path/to/dbs --busco-group vertebrata --n_threads 4 --user-databases whale.pep.fasta dolphin.pep.fasta

User databases will be searched with CRBB; this runs `blastx`, so if you supply ridiculously huge
databases, it *will* take a long time. Future versions will use LAST for all searches to
improve performance, but for now, we're stuck with the NCBI's dinosaur. Also note that the
information from the deflines in your databases will be used to construct the GFF3 file, so if your
databases lack useful IDs, your annotations will too.


