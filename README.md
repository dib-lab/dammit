# dammit!

<img align="left" src="doc/Character_Building.png">

*"I love writing BLAST parsers!" -- no one, ever*

dammit! is a simple de novo transcriptome annotator. It was born out of the
observation that: annotation is mundane and annoying; all the individual pieces
of the process exist already; and, the existing solutions are overly complicated 
or rely on crappy non-free software. 

Science shouldn't suck for the sake of sucking, so dammit attempts
to make this sucky part of the process suck a little less.

## Installation

Unfortunately, annotation necessarily relies on many software packages. I have
worked hard to make dammit rely only on software which is accessible *and* likely
to continue to be so. Most of the dependencies are available in either Ubuntu PPAs
or PyPI, and if not, are trivial to install manually. If the goal is to make annotation
suck less, then installing the necessary software should suck less too.

This guide will assume you're on a Ubuntu system. However, the dependencies should
all run on any flavor of GNU/Linux and on OSX.

First, let's get packages from the Ubuntu PPAs:

    sudo apt-get install hmmer infernal ncbi-blast+ last-align

If you're on Ubuntu 15.10, you can also install TransDecoder through aptitude:

    sudo apt-get install transdecoder

Otherwise, you'll need to install it [manually](https://transdecoder.github.io/). 
To install it in your home directory, execute these commands in your 
terminal:

    cd
    wget https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
    tar -xvzf 2.0.1.tar.gz
    export PATH=$PATH:$HOME/TransDecoder-2.0.1

The above commands will only install it for the current session; to
keep it installed, append it to your bash profile. For GNU/Linux:

    echo 'export PATH=$PATH:$HOME/TransDecoder-2.0.1' >> $HOME/.bashrc

For OSX:

    echo 'export PATH=$PATH:$HOME/TransDecoder-2.0.1' >> $HOME/.bash_profile

Next, we need to install Conditional Reciprocal Best-hits Blast (CRBB). The
algorithm is described in
[Aubry et al.](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004365),
 and is implemented in ruby. Assuming you have ruby, it can be installed with:

    sudo gem install crb-blast

dammit also runs BUSCO to assess completeness. To install it, run thee following
commands:

    cd
    wget http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
    tar -xvzf BUSCO_v1.1b1.tar.gz
    export PATH=$PATH:$HOME/BUSCO_v1.1b1

...and once again, to install it permanently:

    echo 'export PATH=$PATH:$HOME/BUSCO_v1.1b1' >> $HOME/.bashrc

for GNU/Linux, and

    echo  'export PATH=$PATH:$HOME/BUSCO_v1.1b1' >> $HOME/.bash_profile

for OSX.

Finally, let's install dammit!

    pip install dammit


## Usage

### Dependencies

dammit! has three components. The first, `dependencies`, checks whether you have the dependencies installed
correctly and warns you if not. It is run with:

    dammit dependencies

### Databases

Next is the `databases` command. By default, dammit! looks for databases in
`$HOME/.dammit/databases` and will install there if missing. If you have some of the databases
already, you can inform dammit! with the `--database-dir` flag. So, to check for databases:

    dammit databases

or

    dammit databases --database-dir /path/to/databases

To download and install them, run:

    dammit databases --install

The databases are [Pfam](http://pfam.xfam.org/), [Rfam](http://rfam.xfam.org/), 
[OrthoDB](http://orthodb.org/), and a [BUSCO](http://busco.ezlab.org/) database which varies based on the
species type. The default BUSCO database is metazoa; this can be modified with the `--busco-group`
flag. Finally, the `--full` flag will also include [uniref90](http://www.uniprot.org/help/uniref) --
this is left out by default because it is 15+GB and isn't necessarily required for a good
annotation. An example invocation which uses all the options is:

    dammit databases --install --database-dir /path/to/dbs --full --busco-group arthropoda

### Annotation

Paydirt! The `annotate` command runs the BUSCO assessment, assembly stats, and homology searches,
aggregates the results, and outputs a GFF3 file and annotation report. It takes the `--full`,
`--database-dir`, and `--busco-group` options in the same manner as the `databases` command.
Additionally, it can specify an optional output directory, the number of threads to use with
threaded subprograms like HMMER, and a list of user-supplied protein databases in FASTA format. A
simple invocation with the default databases would look like:

    dammit annotate <transcriptome.fasta>

While a more complex invocation might look like:

    dammit annotate <transcriptome.fasta> --database-dir /path/to/dbs --busco-group vertebrata
--n_threads 4 --user-databases whale.pep.fasta dolphin.pep.fasta

User databases will be searched with CRBB; this runs `blastx`, so if you supply ridiculously huge
databases, it *will* take a long time. Future versions will likely use LAST for all searches to
improve performance, but for now, we're stuck with the NCBI's dinosaur. Also note that the
information from the deflines in your databases will be used to construct the GFF3 file, so if your
databases lack useful IDs, your annotations will too.

## 
