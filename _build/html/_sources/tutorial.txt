Tutorial
========

Once you have the dependencies installed, it's time to actually annotate something!
This guide will take you through a short example on some test data.

Data
----

First let's download some test data. We'll start small and use a
*Schizosaccharomyces pombe* transcriptome. Make a working directory and move there,
and then download the file::

    mkdir dammit_test
    cd dammit_test
    wget ftp://ftp.ebi.ac.uk/pub/databases/pombase/FASTA/cdna_nointrons_utrs.fa.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/pombase/FASTA/pep.fa.gz

Decompress the file with gunzip::

    gunzip cdna_nointrons_utrs.fa.gz pep.fa.gz

Databases
---------

If you're just starting, you probably haven't downloaded the databases yet. Here
we'll install the main databases, as well as the `eukaryota` BUSCO database for
our yeast dataset.
This could take a while, so consider walking away and getting yourself a cup of
coffee. If you installed dammit into a virtual environment, be sure to activate
it first::

    dammit databases --install --busco-group eukaryota

Alternatively, if you happen to have downloaded many of these databases before,
you can follow the directions in the `databases guide <databases.html>`__.

While the initial download takes a while, once its done, you won't need to do it again --
dammit keeps track of the database state and won't repeat work its already completed,
even if you accidentally rerun with the ``--install`` flag. 

Annotation
----------

Now we'll do a simple run of the annotator. We'll use `pep.fa` as a user database;
this is a toy example, seeing as these proteins came from the same set of
transcripts as we're annotating, but they illustrate the usage nicely enough.
We'll also specify a non-default BUSCO group. You can replace the argument to
``--n_threads`` with however many cores are available on your system in order to
speed it up.::

    dammit annotate cdna_nointrons_utrs.fa --user-databases pep.fa --busco-group eukaryota --n_threads 1

This will take a bit, so go get another cup of coffee...

