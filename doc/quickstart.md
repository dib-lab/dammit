# Quickstart: prepare databases and test annotation

Once you have dammit [installed](install.md), you'll need to download
and prepare databases before you can annotate a transcriptome. This
quickstart takes you through database preparation (with `dammit run databases`).
These will take a while to prepare, but are the same databases you'll need
for most annotation runs. Finally, we run annotation using a very small sample dataset.

## Check and prepare databases

Here, we'll install the main databases, as well as the
`eukaryota` BUSCO database for our test yeast dataset (below). This could
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

While the initial download takes a while, once it's done, you won't need
to do it again. `dammit` keeps track of the database state and won't
repeat work it's already completed, even if you accidentally rerun with
the `--install` flag.

## Download test data for annotation

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

## Just annotate it, dammit!

With the default databases installed and sample data in hand, we can do a simple run of
the annotator. We'll use `pep.fa` as a user database; this is a toy example,
seeing as these proteins came from the same set of transcripts as we're
annotating, but they illustrate the usage nicely enough. We'll also
specify a non-default BUSCO group (eukaryota). You can replace the argument to
`--n_threads` with however many cores are available on your system in
order to speed it up.:

```
dammit run annotate cdna_nointrons_utrs.fa --user-databases pep.fa --busco-group eukaryota --n_threads 1
```

This will take a bit, so go get another cup of coffee...

## TODO, run sample data and copy output ls below! fix txt, this is from elvers docs.

**consider using https://angus.readthedocs.io/en/2019/dammit_annotation.html#annotation-with-dammit**

## dammit output

After a successful run, you'll have a new directory called `BASENAME.fasta.dammit`. If you look inside, you'll see a lot of files. For example, for a transcriptome with basename `trinity.nema`, the folder `trinity.nema.fasta.dammit` should contain the following files:

```
ls trinity.nema.fasta.dammit/
```    
```    
    annotate.doit.db                              trinity.nema.fasta.dammit.namemap.csv  trinity.nema.fasta.transdecoder.pep
    dammit.log                                    trinity.nema.fasta.dammit.stats.json   trinity.nema.fasta.x.nema.reference.prot.faa.crbl.csv
    run_trinity.nema.fasta.metazoa.busco.results  trinity.nema.fasta.transdecoder.bed    trinity.nema.fasta.x.nema.reference.prot.faa.crbl.gff3
    tmp                                           trinity.nema.fasta.transdecoder.cds    trinity.nema.fasta.x.nema.reference.prot.faa.crbl.model.csv
    trinity.nema.fasta                            trinity.nema.fasta.transdecoder_dir    trinity.nema.fasta.x.nema.reference.prot.faa.crbl.model.plot.pdf
    trinity.nema.fasta.dammit.fasta               trinity.nema.fasta.transdecoder.gff3
    trinity.nema.fasta.dammit.gff3                trinity.nema.fasta.transdecoder.mRNA
```

The two most important files are `trinity.nema.fasta.dammit.fasta` and `trinity.nema.fasta.dammit.gff3`, as they contain the aggregated annotation info per transcript.
`trinity.nema.fasta.dammit.stats.json` also gives summary stats that are quite useful.

For more information on the results, see [dammit results](dammit-results.md).
