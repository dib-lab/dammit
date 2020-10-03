# Annotation Pipelines

The  **`default`** pipeline is suitable for most purposes.
However, dammit has several alternative workflows that either
reduce the number of databases and tools run (the `quick` pipeline)
or annotate with larger databases, such as UniRef90 ( `full` pipeline).

> **Note:** To use these pipelines, first run `dammit run databases` to 
> make sure the relevant databases are installed. 
> Then you can proceed with `dammit run annotate`


## Default

By default, `dammit` runs the following:

- **default:**
    - `busco`quality assessment
    - `transdecoder` ORF prediction
    - `shmlast` to any user databases
    - `hmscan` - Pfam-A
    - `cmscan` - Rfam
    - `LAST` mapping to OrthoDB and Swiss-Prot

The databases used for this pipeline require approximately ~18GB of storage space,
plus a few hundred MB per busco database. We recommend running this pipeline with
at least 16GB of RAM.

code to run this pipeline:

    dammit run databases --install
    dammit run annotate

> If specifying a custom location for your databases, add `--databases-dir /path/to/dbs`


## Alternative annotation pipelines:

### quick pipeline

- **quick (`--pipeline quick`):**
    - `busco`quality assessment
    - `transdecoder` ORF prediction
    - `shmlast` to any user databases

The `quick` pipeline can be used for running a minimal annotation run: 
BUSCO quality assessment, ORF prediction with transdecoder, and `shmlast` 
to map to any user databases. While this pipeline may require less database
space, we still recommend running with 16G of RAM, especially if mapping
to a user-provided protein database.


code to run this pipeline:

    dammit run databases --install --pipeline quick
    dammit run annotate --pipeline quick

> If specifying a custom location for your databases, add `--databases-dir /path/to/dbs`

### full pipeline

_warning: time and resource intensive!_

The `full` pipeline starts from the `default` pipeline and adds a mapping 
database, UniRef90. 

- **full (`--pipeline full`):**
    - `busco` quality assessment
    - `transdecoder` ORF prediction
    - `shmlast` to any user databases
    - `hmscan` - Pfam-A
    - `cmscan` - Rfam
    - `LAST` mapping to OrthoDB, Swiss-Prot, and **UniRef90**

As of fall 2020, the UniRef90 fasta is 26G (gzipped).

code to run this pipeline:

    dammit run databases --install --pipeline full
    dammit run annotate --pipeline full

> If specifying a custom location for your databases, add `--databases-dir /path/to/dbs`

### nr pipeline

_warning: REALLY time and resource intensive!_

**UniRef90** is a set of UniProt sequences clustered
by >=90% sequence identity. UniRef allows a searching to a larger set of
sequence records while hiding redundant sequences. See the [UniRef 
documentation](https://www.uniprot.org/help/uniref) for more.

- **nr (`--pipeline nr`):**
    - `busco` quality assessment
    - `transdecoder` ORF prediction
    - `shmlast` to any user databases
    - `hmscan` to Pfam-A
    - `cmscan` - Rfam
    - `LAST` mapping to OrthoDB, Swiss-Prot, and **nr**

As of fall 2020, the nr fasta is 81G (gzipped).

code to run this pipeline:

    dammit run databases --install --pipeline nr
    dammit run annotate --pipeline nr

> If specifying a custom location for your databases, add `--databases-dir /path/to/dbs`

**Note:** Since all of these pipelines use a core set of tools, and since dammit uses `snakemake`
to keep track of the files that have been run, dammit will not rerun the core tools if decide
to alter the pipeline you're running. So for example, you could start by running a `quick`
run, and later run `default` if desired. In that case, `dammit` would run only the new annotation
steps, and reintegrate the relevant outputs into new `dammit` gff3 and annotated fasta files.
