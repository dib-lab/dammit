---
title: 'About the Databases'
---

dammit can make use of the following databases:

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
    > you use depends on your particular organism.
    >
    > dammit uses the metazoa database by default, but different
    > databases can be used with the `--busco-group` parameter. You
    > should try to use the database which most closely bounds your
    > organism.

5.  [Swiss-Prot](https://www.uniprot.org/help/about)

    > Swiss-Prot is the manually reviewed and curated non-redundant 
    > protein sequence database. The aim to to provide high-quality
    > annotations linked to all known information about each protein.
    > dammit now maps to Swiss-Prot by default.

6.  [uniref90](http://www.uniprot.org/help/uniref)

    > uniref is a curated collection of most known proteins, clustered
    > at a 90% similarity threshold. This database is comprehensive, and
    > thus quite enormous. dammit does not include it by default due to
    > its size, but it can be installed and used with the 
    > `--pipeline full` flag.

7.  [NR](http://www.uniprot.org/help/uniref)

    > The `nr` is a very large database consisting of both non-curated
    > and curated database entries. While the name stands for "non-redundant",
    > this databse is no longer non-redundant. Given the time and memory requirments,
    > NR is only a good choice for species and/or sequences you're unable to confidently 
    > annotate via other databases. It can be installed and used with 
    > the `--pipeline nr` flag.


The specific databases run can be selected via the `--pipeline` flag.
For all available pipelines, see the [pipelines](pipelines.md) section.

To install, for example, all databases required for `full`
in a custom location (`/path/to/dbs`), you could run the following:

```
dammit run databases --install --database-dir /path/to/dbs --pipeline full 
```
