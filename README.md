---
title: README
---

[![image](https://travis-ci.org/dib-lab/dammit.svg)](https://travis-ci.org/dib-lab/dammit)

[![Documentation Status](https://readthedocs.org/projects/dammit/badge/)](http://dammit.readthedocs.io/en/latest)

*\"I love writing BLAST parsers!\" \-- no one, ever*

dammit is a simple de novo transcriptome annotator. It was born out of
the observation that: annotation is mundane and annoying; all the
individual pieces of the process exist already; and, the existing
solutions are overly complicated or rely on crappy non-free software.

Science shouldn't suck for the sake of sucking, so dammit attempts to
make this sucky part of the process suck a little less.

### Basic Usage

The most basic usage is to download and install a subset of the databases:

    dammit databases --install --quick

And the annotate with:

    dammit annotate <transcriptome_fasta>

Head over to the [docs](http://dib-lab.github.io/dammit/) for much more detailed
information!

Known Issues
============

-   On some systems, installation of the ConfigParser package can get
    borked, which will cause and exception to be thrown. This can be
    fixed by following the directions at issue \#33:
    <https://github.com/camillescott/dammit/issues/33>.
-   There can be errors resuming runs which were interrupted on the
    BUSCO stage. If the task fails on resume, delete the BUSCO results
    folder within your dammit results folder, which will have a name of
    the form [run\_\<name\>.busco\_results]{.title-ref}.

Acknowledgements
================

I\'ve received input and advice from a many sources, including but
probably not limited to: C Titus Brown, Matt MacManes, Chris Hamm,
Michael Crusoe, Russell Neches, Luiz Irber, Lisa Cohen, Tessa Pierce,
Sherine Awad, and Tamer Mansour.

CS was funded by the National Human Genome Research Institute of the
National Institutes of Health under Award Number R01HG007513 through May
2016, and now receives support from the Gordon and Betty Moore
Foundation under Award number GBMF4551.
