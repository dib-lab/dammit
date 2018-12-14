---
title: System Requirements
summary: Minimal hardware and software requirements.
---

dammit, for now, is officially supported on GNU/Linux systems via
[bioconda](https://bioconda.github.io/index.html). macOS support will be
available via bioconda soon.

For the standard pipeline, dammit needs ~18GB of storage space to store its
prepared databases, plus a few hundred MB per BUSCO database. For the
standard annotation pipeline, we recommend at least 16GB of RAM. This can be
reduced by editing LAST parameters via a custom configuration file.

The full pipeline, which uses uniref90, needs several hundred GB of
space and considerable RAM to prepare the databases. You'll also want
either a fat internet connection or a big cup of patience to download
uniref.

For some species, we have found that the amount of RAM required can be proportional to the size of the transcriptome being annotated.
