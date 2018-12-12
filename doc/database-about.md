---
title: 'About the Databases'
---

dammit uses the following databases:

1.  [Pfam-A](http://pfam.xfam.org/)

    > Pfam-A is a collection of protein domain profiles for use with
    > profile hidden markov model programs like
    > [hmmer](http://hmmer.janelia.org/). These searches are moderately
    > fast and very sensitive, and the Pfam database is very well
    > curated. Pfam is used during TransDecoder's ORF finding and for
    > annotation assignment.

2.  [Rfam](http://rfam.xfam.org/)

    > Rfam is a collection of RNA covariance models for use with
    > programs like [Infernal](http://infernal.janelia.org/). Covariance
    > models describe RNA secondary structure, and Rfam is a curated
    > database of non-coding RNAs.

3.  [OrthoDB](http://orthodb.org/)

    > OrthoDB is a curated database of orthologous genes. It attempts to
    > classify proteins from all major groups of eukaryotes and trace
    > them back to their ancestral ortholog.

4.  [BUSCO](http://busco.ezlab.org/)

    > BUSCO databases are collections of "core" genes for major
    > domains of life. They are used with an accompanying BUSCO program
    > which assesses the completeness of a genome, transcriptome, or
    > list of genes. There are multiple BUSCO databases, and which one
    > you use depends on your particular organism. Currently available
    > databases are:
    >
    > 1.  Metazoa
    > 2.  Vertebrata
    > 3.  Arthropoda
    > 4.  Eukaryota
    >
    > dammit uses the metazoa database by default, but different
    > databases can be used with the `--busco-group` parameter. You
    > should try to use the database which most closely bounds your
    > organism.

5.  [uniref90](http://www.uniprot.org/help/uniref)

    > uniref is a curated collection of most known proteins, clustered
    > at a 90% similarity threshold. This database is comprehensive, and
    > thus quite enormous. dammit does not include it by default due to
    > its size, but it can be installed and used with the `--full` flag.

A command using all of these potential options and databases might look
like:

```
dammit databases --install --database-dir /path/to/dbs --full --busco-group arthropoda
```
