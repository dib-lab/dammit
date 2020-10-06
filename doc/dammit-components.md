# Dammit Components

Under the hood, dammit uses the snakemake workflow management system
to manage database downloads and run the external annotation software.

## **`dammit run`**

`dammit` can run two main workflows, `databases`, and `annotate`.

  - **`databases`** handles downloading and preparing the annotation databases. Usage info [here](databases-usage.md).
  > databases must be used to properly prepare databases pror to running annotatation
  - [**`annotate`**](annotate.md) uses these databases for transcriptome annotation. Usage info [here](annotate.md).

## dammit file conversion components

Each annotation program run as part of a dammit pipeline produces an
annotation file with a tool-specific formatting and indexing (0-based or 1-based).
Dammit includes a set of file conversion utilities that translate each output 
file to a standardized gff3 format, and a utility that combines the standardized 
output into a single set of annotations for output in gff3 and fasta 
(the final outputs of each dammit pipeline). 

Each of these conversion utilities can now be run independently, so that they
can be used outside of a full dammit pipeline run. You may find these components
useful if you want to run these tools on additional databases not included in the 
dammit pipelines.

**FASTA munging commands:**
  
  - **`rename-fasta`**         Copy a FASTA file and rename the headers.
  - **`transcriptome-stats`**  Run basic metrics on a transcriptome
  - **`annotate-fasta`**      Annotate a FASTA file from a GFF3 file.


**Filtering commands:**

  - **`best-hits`**            Filter query best-hits from a MAF file.

**Conversion commands:**

  - **`maf-to-gff3`**          Convert MAF to GFF3.
  - **`shmlast-to-gff3`**      Convert shmlast CSV output to GFF3.
  - **`hmmscan-to-gff3`**      Convert HMMER to GFF3.
  - **`cmscan-to-gff3`**       Convert Infernal's cmscan output to GFF3.

**Transformation commands:**
  
  - **`merge-gff3`**           Merge a collection of GFF3 files.
  - **`remap-hmmer-coords`**   Remap hmmscan coordinates using TransDecoder ORF predictions.

