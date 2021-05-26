# About

This page goes a little more in depth on the software and its goals.

## Background and Motivations

Several different factors motivated dammit's development. The first of
these was the sea lamprey transcriptome project, which had annotation as
a primary goal. Many of dammit's core features were already implemented
there, and it seemed a shame not share that work with others in a usable
format. Related to this was a lack of workable and easy-to-use existing
solutions; in particular, most are meant to be used as protocols and
haven't been packaged in an automated format. Licensing was also a big
concern -- software used for science should be open source, easily
accessible, remixable, and free.

Implicit to these motivations is some idea of what a good annotator
*should* look like, in the author's opinion:

1.  It should be easy to install and upgrade
2.  It should only use Free software
3.  It should make use of standard databases
4.  It should output in reasonable formats
5.  It should be relatively fast
6.  It should try to be correct, insofar as any computational approach
    can be "correct"

## The Obligatory Flowchart

![The Workflow](static/workflow.svg)

## Databases Used

-   Pfam-A
-   Rfam
-   OrthoDB
-   Swiss-Prot
-   BUSCO databases
-   Uniref90
-   User-supplied protein databases

The last one is important, and sometimes ignored. 
Dammit uses an approach similar to Conditional Reciprocal Best Blast
to map to user-supplied protein databases (details below).

To see more about the included databases, 
see the [About Databases](database-about.md) section.

## Software Used

The specific set of software and databases used can be modified by specifying different [pipelines](pipelines.md).
The full set of software than can be run is:

-   TransDecoder
-   BUSCO
-   HMMER
-   Infernal
-   LAST
-   shmlast (for crb-blast to user-supplied protein databases)

All of these are Free Software, as in freedom and beer. 

### shmlast: Conditional Reciprocal Best LAST for mapping to user databases

Reciprocal Best Hit mapping (RBH) is a standard method for ortholog detection.
However, transcriptomes have multiple transcript isoforms, which confound RBH.

![](static/RBH.svg)

**Conditional Reciprocal Best Blast (CRBB)** attempts to associate those isoforms
with appropriate annotations by learning an appropriate e-value cutoff for 
different transcript lengths. The original implementation of CRBB 
can be found [here](https://github.com/cboursnell/crb-blast). 

![CRBB](static/CRBB_decision.png)

*from [Deep Evolutionary Comparison of Gene Expression Identifies Parallel Recruitment of Trans-Factors in Two Independent Origins of C4 Photosynthesis](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004365)*

**shmlast** is a reimplementation of the Conditional Reciprocal Best Hits 
algorithm for finding potential orthologs between a transcriptome and 
a species-specific protein database. It uses the LAST aligner and the 
pydata stack to achieve much better performance while staying in the 
Python ecosystem. 

One limitation is that LAST has no equivalent to `tblastn`. So, we find
the RBHs using the TransDecoder ORFs, and then use the model on the
translated transcriptome versus database hits. 

`shmlast` is published in JOSS, doi:[10.21105/joss.00142](https://joss.theoj.org/papers/10.21105/joss.00142).


## The Dammit Software

dammit is built on the [Snakemake](https://snakemake.readthedocs.io/en/stable/)
workflow management system. This means that the dammit pipeline enjoys all the features of any
Snakemake workflow: reproducibility, ability to resume, cluster support, and per-task environment
management. Each step in dammit's pipeline(s) is implemented as a Snakemake
[wrapper](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#wrappers);
when dammit is executed, it generates the targets for the pipeline being run as specified in its
pipelines file and passes them along to the Snakemake executable. The dammit frontend simplifies
the interface for the user and constrains the inputs and options to ensure the pipeline
will always run correctly.

One of the essential, and most annoying, parts of annotation is the conversion and collation
of information from many different file formats. Dammit includes a suite of minimal command
line utilities implementing a number of these things, including converting several formats
to GFF3, merging GFF3 files, and filtering alignment results for best hits. More details on
these utilities can be found in the [components](dammit-components.md) section.

