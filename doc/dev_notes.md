# For `dammit` developers

## Setting up your local computer for `dammit` devevelopment

Make sure conda is installed. If not, here are instructions:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p $HOME/miniconda
```

### Conda Setup 

Make sure conda is configured properly
```
conda config --set always_yes yes
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Dammit Repository Setup

Fork the `dammit` repository to your account. Clone this fork to your local computer, then create and name a development branch.
In this example, we name our dev branch `testing`:

```
git clone https://github.com/YOUR-USERNAME/dammit.git
git checkout -b testing
git branch
```
Now you are on the `testing` branch.

Keep original repository in the `master` branch, so that you can stay up to date with any changes in the main repo.
To do this, first add a remote called `upstream`, which links to the main dammit repository.
```
git remote add upstream https://github.com/dib-lab/dammit.git 
```

Then, make sure the `master` branch is up-to-date by periodically running:
```
git pull upstream master
```

### Create a development environment

Create a conda environment with dependencies installed
```
conda env create -f environment.yml -n dammit_dev
conda activate dammit_dev
```
Now, install the dammit software into the `dammit_dev` environment:

```
pip install -e '.'
```

### Test your Installation

#### Install databases

Install databases  (will install in `~/.dammit/databases/` unless you specify an alternate location with `--database-dir`)
```
dammit run databases --install
```

#### Do a quick annotation run

```
dammit run annotate ...FINISH THIS...
```

## To-do for `dammit`

- [ ] update transdecoder version
- [ ] orthodb version (other database versions?)
- [ ] add swissprot
- [ ] add pipeline for accepting .pep file as input (skips transdecoder, transcriptome stats and BUSCO tasks)

#### Versioning

A new version is required when a new version of a database is added or a major change happens that changes the commandline interface. Change the VERSION file when this hapens. 

(Note 11/30/2018: We should make all changes above in the T-do, then move to v1.1)

## Notes on dammit

Written by [Camille Scott](http://www.camillescott.org/) and [N Tessa Pierce](http://bluegenes.github.io/).

Dammit relies on the snakemake workflow software. 
To learn snakemake, check out the [documentation](https://snakemake.readthedocs.io/en/stable/), and start with snakemake tutorials such as [this one by ctb](https://github.com/ctb/2019-snakemake-ucdavis).

### Architecture:

**TO DO: UPDATE THIS FOR DAMMIT 2**


#### Take a look at code and tests in the `dammit` directory:

* The core driver of `dammit` is the `dammit/app.py` file, which sets up the commandline arguments. Everything happens here. If you want to add an argument, this is where it hapens.  
* There are two subcommand task-handler files: `annotate.py` and `databases.py`
* Tasks are steps being run, separated into different files. For example, the `hmmer.py` file contains all hmmer tasks.
* The task handler has on logger, pulls from config to figure out where databases are located (all happens in background), some doit stuff happening
* [Decorators](https://realpython.com/primer-on-python-decorators/) transfer the function's `return` into a doit function (e.g. line 59 of shell.py) `import doit_task` then `@doit_task`

`databases`, 2 pipelines:

  * `quick`
  * `full`

`annotate`, more pipelines: 

  * `uniref1`
  * `full`
  * `nr`

#### `config.json`

Can use custom `config.json` file to include different parameters for the programs run by the tasks, e.g. `transdecoder LongOrgs -m 50`, etc.

#### `parallel.py`

hmmer, infernal, lastl, 

requires gnu parallel

(There are instructions for how to runon multi-node hpc, somewhere.)

#### `ui.py`

output for user to markdown formatting for copying/pasting into GitHub issue reporting

`generate-test-data-.sh` re-genreates test data adn puts it in proper dirs

### TESTS!

`dammit/tests`

Run `test_databases.py` yourself, locally (because databases cannot be cached on travis-ci)

* makes sure tasks and pipeline run and produce output, they don't all check expected output. some integration output.
* uses pytest
* set of tests files
* testing pydoit tasks is a pain
* under utils, file to run tasks. give it a list of tasks, it will execute in own directory.
* functions start with 'test', check assertions
* fixtures are a means of setting upa consistent environment before running an individual test, e.g. be in a clean directory. tmpdir will create a randomly name temporary directory.
* make tests for new tasks (Sometimes they will take a long time to run...)
* `test_annotate.py` must be run locally by yourself.
* before pushing release, do both of these
* `make long tests` (assumes environment is already setup)
* [travis-ci](https://travis-ci.org/dib-lab/dammit/) is building the recipe that lives in the repo
* `make-ci-test`, not long and not huge and not requires_datbases

## Reviewing a PR

**Tests must pass before merging!**

* Have there been radical changes? (Are you adding things to handler, maybe time to take a step back and make sure code uses reasonable variable names, tests, etc)
* Does travis build?
* Try to make commit messages somewhat informative

If these all seem reasonable to you, approve!

## Fix travis:

`.travis.yml`

* make sure conda env uses right Python
* fix conda channel order

## Bioconda

* https://anaconda.org/bioconda/dammit
* Recipe: https://github.com/bioconda/bioconda-recipes/blob/master/recipes/dammit/meta.yaml

## Documentation

* http://dib-lab.github.io/dammit/
* [Tutorial from angus 2018](https://angus.readthedocs.io/en/2018/dammit_annotation.html)
