# Annotate

The dammit `annotate` component uses the installed databases for transcriptome annotation.



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


