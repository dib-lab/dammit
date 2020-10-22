# Annotate

_in progress_

The dammit `annotate` component uses the installed databases for transcriptome annotation.

## Just annotate it, dammit!

After you've properly installed the [databases](database-usage.md), you can start the annotation.
To run the annotation, you only need to provide a set of transcripts to annotate. 


    dammit run annotate TRANSCRIPTOME.fasta


Optionally, allow `dammit` to use additional threads using the `n_threads` parameter:
    
    dammit run annotate TRANSCRIPTOME.fasta --n_threads 4


If you'd like to customize the output or other parameters such as the e-value for similarity searches,
you can provide customization on the command line or in a configuration file


## Additional Usage info 

To see the general dammit usage information, run:

    dammit run --help

You should see the following:

```
ADD HELP OUTPUT
```

The `--pipeline` option can be used to switch the set of databases being used for annotation.
See the [annotation pipelines](pipelines.md) doc for info about each specific pipeline.
Note that these pipelines all run a core set of programs. If you re-run an annotation with a
larger pipeline, dammit will not re-run analyses that have already been completed. Instead,
dammit will run any new analyses, and integrate them into the final fasta and gff3.

To see annotation-specific configuration info, run:

    dammit run annotate --help

You should see the following:

```
Usage: dammit run annotate [OPTIONS] TRANSCRIPTOME

  The main annotation pipeline. Calculates assembly stats; runs BUSCO; runs
  LAST against OrthoDB (and optionally uniref90), HMMER against Pfam,
  Inferal against Rfam, and Conditional Reciprocal Best-hit Blast against
  user databases; and aggregates all results in a properly formatted GFF3
  file.

Options:
  -n, --base-name TEXT      Base name to use for renaming the input
                            transcripts. The new names will be of the form
                            <name>_<X>. It should not have spaces, pipes,
                            ampersands, or other characters with special
                            meaning to BASH.

  -e, --evalue FLOAT        e-value cutoff for similarity searches.
  -o, --output-dir TEXT     Output directory. By default this will be the name
                            of the transcriptome file with `.dammit` appended

  -u, --user-database TEXT  Optional additional protein databases.  These will
                            be searched with CRB-blast.

  --dry-run
  --help                    Show this message and exit.
```

Add these options as needed. For example, add annotation to a user database and specify an
output directory name like so:

    dammit run annotate --user-database DB-FILE --output-dir dammit-results 


