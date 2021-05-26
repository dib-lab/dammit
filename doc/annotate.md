# Annotate

The dammit `annotate` component uses the installed databases for transcriptome annotation.

## Just annotate it, dammit!

After you've properly installed the [databases](database-usage.md), you can start the annotation.
To run the annotation, you only need to provide a set of transcripts to annotate. 


    dammit run annotate TRANSCRIPTOME.fasta


Optionally, allow `dammit` to use additional threads using the `n_threads` parameter:
    
    dammit run --n-threads 4 annotate TRANSCRIPTOME.fasta


If you'd like to customize the output or other parameters such as the e-value for similarity searches,
you can provide customization on the command line or in a configuration file.


## Additional Usage info 

To see the general dammit usage information, run:

    dammit run --help

You should see the following:

```
Usage: dammit run [OPTIONS] COMMAND [ARGS]...

  Run the annotation pipeline or install databases.

Options:
  --database-dir TEXT             Directory to store databases. Existing
                                  databases will not be overwritten.

  --conda-dir TEXT                Directory to store snakemake-created conda
                                  environments.

  --temp-dir TEXT                 Directory to store dammit temp files.
  --busco-group [bacteria_odb10|acidobacteria_odb10|actinobacteria_phylum_odb10|actinobacteria_class_odb10|corynebacteriales_odb10|...]
                                  BUSCO group(s) to use/install.
  --n-threads INTEGER             Number of threads for overall workflow
                                  execution

  --max-threads-per-task INTEGER  Max threads to use for a single step.
  --busco-config-file TEXT        Path to an alternative BUSCO config file;
                                  otherwise, BUSCO will attempt to use its
                                  default installation which will likely only
                                  work on bioconda. Advanced use only!

  --pipeline [default|quick|full|nr]
                                  Which pipeline to use. Pipeline options:
                                  quick: excludes:  the Infernal Rfam tasks,
                                  the HMMER Pfam tasks, and the LAST OrthoDB
                                  and uniref90 tasks. Best for users just
                                  looking to get basic stats and conditional
                                  reciprocal best LAST from a protein
                                  database.  full: Run a "complete"
                                  annotation; includes uniref90, which is left
                                  out of the default pipeline because it is
                                  huge and homology searches take a long time.
                                  nr:  Also include annotation to NR database,
                                  which is left out of the default and "full"
                                  pipelines because it is huge and homology
                                  searches take a long time. More info  at
                                  https://dib-lab.github.io/dammit.

  --help                          Show this message and exit.

Commands:
  annotate   The main annotation pipeline.
  databases  The database preparation pipeline.
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
Usage: dammit run annotate [OPTIONS] TRANSCRIPTOME [EXTRA_SNAKEMAKE_ARGS]...

  The main annotation pipeline. Calculates assembly stats; runs BUSCO; runs
  LAST against OrthoDB (and optionally uniref90), HMMER against Pfam,
  Infernal against Rfam, and Conditional Reciprocal Best-hit Blast against
  user databases; and aggregates all results in a properly formatted GFF3
  file.

Options:
  -n, --base-name TEXT       Base name to use for renaming the input
                             transcripts. The new names will be of the form
                             <name>_<X>. It should not have spaces, pipes,
                             ampersands, or other characters with special
                             meaning to BASH. Superseded by --regex-rename.

  --regex-rename TEXT        Rename transcripts using a regex pattern. The
                             regex should follow  Python `re` format and
                             contain a named field keyed as `name` that
                             extracts the desired string. For example,
                             providing "(?P<name>^[a-zA-Z0-9\.]+)" will match
                             from the beginning of the sequence header up to
                             the first symbol that is not alphanumeric or a
                             period. Supersedes --base-name.

  --rename / --no-rename     If --no-rename, original transcript names are
                             preserved in the final annotated FASTA. --base-
                             name is still used in intermediate files. If
                             --rename (the default  behavior), the renamed
                             transcript names are used in the final  annotated
                             FASTA.

  -e, --global-evalue FLOAT  global e-value cutoff for similarity searches.
  -o, --output-dir TEXT      Output directory. By default this will be the
                             name of the transcriptome file with `.dammit`
                             appended

  -u, --user-database TEXT   Optional additional protein databases.  These
                             will be searched with CRB-blast.

  --dry-run
  --help                     Show this message and exit.
```

Add these options as needed. For example, add annotation to a user database and specify an
output directory name like so:

    dammit run annotate --user-database DB-FILE --output-dir dammit-results 

General run arguments need to be added in front of `annotate`, e.g. - 

    dammit run --n-threads 4 --pipeline quick annotate --user-database DB-FILE --output-dir dammit-results 





