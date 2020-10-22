Dammit Results
===

When `dammit` finishes running, you should see:

**to do: consider using https://angus.readthedocs.io/en/2019/dammit_annotation.html#annotation-with-dammit**

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
`trinity.nema.fasta.dammit.stats.json` also gives summary stats that are quite useful. If you'd like to look into the remaining files, here are the programs that produced each.

dammit results:
  - \*gff3

external program results:
  - to do: LIST EACH FILE + some info for it




## Parsing dammit output

Dammit provides transcript annotations for mapping against all of the databases in the pipeline you selected
in the final dammit files, `BASENAME.fasta.dammit.fasta` and `BASENAME.fasta.dammit.gff3`.
If you'd like to select certain annotations (e.g. to create an alternative gene-to-transcript map), you can
use python or R to parse the `gff3` results file. If using python, dammit provides a `GFF3Parser` utility to facilitat parsing.

```
import pandas as pd
from dammit.fileio.gff3 import GFF3Parser
```

```
gff_file = "nema-trinity.fa.dammit/nema-trinity.fa.dammit.gff3"
annotations = GFF3Parser(filename=gff_file).read()
names = annotations.sort_values(by=['seqid', 'score'], ascending=True).query('score < 1e-05').drop_duplicates(subset='seqid')[['seqid', 'Name']]
new_file = names.dropna(axis=0,how='all')
new_file.head()
```

Try commands like:
```
annotations.columns
```

```
annotations.head()
```

```
annotations.head(50)
```



**todo: add R code?**


## Other tutorials and Workshop materials

See this workshop [tutorial](https://angus.readthedocs.io/en/2018/dammit_annotation.html) for further practice with using `dammit` for annotating a *de novo* transcriptome assembly.
Please note that the commands used were for a prior version of dammit, but all analysis remains relevant.

