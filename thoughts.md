# Thoughts

## Usage

### Installing

Install databases and dependencies; checks and downloads all if necessary:

    dammit install

Install databases:

    dammit install-databases [--database-dir <dir>]

Install dependencies. Should there be an option for where to put them here?

    dammit install-dependencies

### Running

Assess transcriptome. Runs BUSCO and basic stats. If it fails because of a
BUSCO database dependency check, inform the user of the subcommand necessary
to retrieve only the necessary database with `dammit install-databases`:

    dammit assess -i <transcriptome.fasta> -o <output_dir> --busco-group <group>

Annotate transcriptome. If it fails due to a dependency check, inform the suer
of the necessary subcommand:

    dammit annotate -i <transcriptome.fasta> -o <output_dir> [--full]

Assess and annotate. Above two necessities apply:

    dammit all -i <transcriptome.fasta> -i <output_dir> --busco-group <group> [--full]

## Design

### Installing

* Need to decide whether to manage install with doit or a more traditional build system.
* It really should try to use system-available programs first.
