# Quickstart: run dammit! with some sample data

Once you have dammit [installed](install.md), it's time to actually
annotate something! This guide will take you through a short example on
some sample data.

## First, download some test data

First let's download some test data. We'll start small and use a
*Schizosaccharomyces pombe* transcriptome. Make a working directory and
move there, and then download the file:

```
mkdir dammit_test
cd dammit_test
wget ftp://ftp.ebi.ac.uk/pub/databases/pombase/OLD/20170322/FASTA/cdna_nointrons_utrs.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/pombase/OLD/20170322/FASTA/pep.fa.gz
```

Decompress the file with gunzip:

```
gunzip cdna_nointrons_utrs.fa.gz pep.fa.gz
```

## Next, check and prepare databases

If you're just starting, you probably haven't downloaded the databases
yet. Here, we'll install the main databases, as well as the
`eukaryota` BUSCO database for our yeast dataset. This could
take a while, so consider walking away and getting yourself a cup of
coffee.

If you installed dammit into a virtual environment, be sure to
activate it first:
```
conda activate dammit
```

Now install databases:
```
dammit run databases --install --busco-group eukaryota
```

By default, dammit downloads databases to your home directory: `$HOME/.dammit/databases`

> If you're on an HPC or other system with limited space in your home directory, 
> or if you've already downloaded some databases and you'd like to use them with dammit, 
> see the [Database Usage](database-usage.md) section to specify a custom location.

While the initial download takes a while, once its done, you won't need
to do it again. `dammit` keeps track of the database state and won't
repeat work its already completed, even if you accidentally rerun with
the `--install` flag.

## Annotate a transcriptome

Now that the default databases are installed, we can do a simple run of
the annotator. We'll use `pep.fa` as a user database; this is a toy example,
seeing as these proteins came from the same set of transcripts as we're
annotating, but they illustrate the usage nicely enough. We'll also
specify a non-default BUSCO groupi (eukaryota). You can replace the argument to
`--n_threads` with however many cores are available on your system in
order to speed it up.:

```
dammit run annotate cdna_nointrons_utrs.fa --user-databases pep.fa --busco-group eukaryota --n_threads 1
```

This will take a bit, so go get another cup of coffee...

## Other tutorials and Workshop materials

See this workshop [tutorial](https://angus.readthedocs.io/en/2018/dammit_annotation.html) for further practice with using `dammit` for annotating a *de novo* transcriptome assembly.
Please note that the commands used were for a prior version of dammit, but all analysis remains relevant.

