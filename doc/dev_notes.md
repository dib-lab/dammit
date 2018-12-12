# Dammit development

[dammit!](https://github.com/dib-lab/dammit)


## Setting up your local computer for dammit dev

We can basically follow the [instructions for travis](https://github.com/dib-lab/dammit/blob/master/.travis.yml), because we're telling [travis]() to do what we are doing on our local computers.

Make sure conda is installed. If not, here are instructions:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
export PATH="$HOME/miniconda3/bin:$PATH"
```

Fork dammit repository to your account. Clone your copy to your local computer, create dev branch

```
git clone https://github.com/username/dammit.git
git remote add upstream https://github.com/dib-lab/dammit.git 
git checkout -b testing
git branch
```
Should be on the `testing` branch.

Keep original repository in `master` branch. Make sure it is up-to-date periodically by running:
```
git pull upstream master
```
Set up a Python 3 environment to work in:

```
conda create -n dammit_dev python=3
source activate dammit_dev
```
Install dependencies:
```
conda config --set always_yes yes --set changeps1 no
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install python numpy pandas numexpr>=2.3.1 khmer>=2.1 sphinx>1.3.1 sphinx_rtd_theme>=0.1.9 pytest pytest-runner doit>=0.29.0 matplotlib shmlast infernal hmmer transdecoder=3.0.1 last busco=3.0.2 parallel bioconductor-seqlogo
python setup.py install
```
Last line of the output should be:
```
Finished processing dependencies for dammit==1.0rc2
```

Lastly, install databases  (will install in `~/.dammit/databases/`)
```
dammit databases --install
```
Output should be:
```
(dammit_dev) campus-019-072:dammit johnsolk$ dammit databases --install
Unable to revert mtime: /Library/Fonts
# dammit
## a tool for easy de novo transcriptome annotation

by Camille Scott

**v1.0rc2**, 2018

## submodule: databases
### Database Install
#### Info
* Database Directory: /Users/johnsolk/.dammit/databases
* Doit Database: /Users/johnsolk/.dammit/databases/databases.doit.db


*All database tasks up-to-date.*

Nothing to install!
```
Now you are ready to edit and make changes!

## To-do:

- [ ] update transdecoder version
- [ ] orthodb version (other database versions?)
- [ ] add swissprot
- [x] change order of conda channels to include conda-forge last
- [ ] update documentation
- [ ] add pipeline for accepting .pep file as input (skips transdecoder, transcriptome stats and BUSCO tasks)

## Notes on dammit

Written by [Camille Scott](http://www.camillescott.org/). See [tutorial](https://angus.readthedocs.io/en/2018/dammit_annotation.html).

1. Look at [pydoit](http://pydoit.org/index.html) documentation, and [Camille's workshop](https://dib-training.readthedocs.io/en/pub/2016-01-20-pydoit-lr.html)
2. [pypi](https://pypi.org/project/dammit/#history) and [bioconda](https://anaconda.org/bioconda/dammit) (supported method of installation)

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Accepting code from other people, process for Pull Requests (for advice, see [khmer getting started dev notes](https://khmer.readthedocs.io/en/latest/dev/getting-started.html))

tests must pass before merging

Architecture:

#### look at code and tests in the `dammit` directory

Tasks are steps being run, separated into different files. For example, the `hmmer.py` file contains all hmmer tasks. 

Decorators convert return into doit function (e.g. line 59 of shell.py) import doit_task then @doit_task

core driver: damit/app.py sets up commandline arguments. everything happens here. if you want to add an argument, this is where it hapens.  

annotate steps get added in the pipeline

core subcommand task-handler files:
`annotate.py`
`databases.py`

task handler has on logger, pulls from config to figure out where databases are located (all happens in background), some doit stuff happening

new version with new version of database, change VERSION file
(make all changes, then v1.1)

databases,  2 pipelines (quick, full)

annotate:
more pipelines, uniref full, nr

## TESTS!

`dammit/tests`
* makes sure tasks and pipeline run and produce output, they don't all check expected output. some integration output.
* uses pytest
* set of tests files
* testing pydoit tasks is a pain
* under utils, file to run tasks. give it a list of tasks, it will execute in own directory.
* functions start with 'test', check assertions
* fixtures are a means of setting upa consistent environment before running an individual test, e.g. be in a clean directory. tmpdir will create a randomly name temporary directory.
* make tests for new tasks (Sometimes they will take a long time to run...)


Have to run test_databases.py yourself, locally (because databases cannot be cached on travis-ci)

test_annotate.py must be run locally by yourself.

before pushing release, do both of these.

`make long tests` (assumes environment is already setup)

travis-ci is building the recipe that lives in the repo

make-ci-test (not long and not huge and not requires_datbases)


#### config.json

can use custom config.json file to include different parameters for the programs run by the tasks, e.g. transdecoder LongOrgs -m 50, etc.

#### parallel .py

hmmer, infernal, lastl, 

requires gnu parallel

(isntructions for how to runon multi-node hpc, somewhere)

#### ui .py

output for user to markdown formatting for copying/pasting into GitHub issue reporting

`generate-test-data-.sh` re-genreates test data adn puts it in proper dirs

### Documentation

* http://dib-lab.github.io/dammit/

* [Tutorial from angus 2018](https://angus.readthedocs.io/en/2018/dammit_annotation.html)

### Update bioconda

update https://github.com/bioconda/bioconda-recipes/blob/master/recipes/dammit/meta.yaml


## reviewing PR

* radical changes? (Are you adding things to handler, maybe time to take a step back and make sure code uses reasonable variable names, tests, etc)
* does travis build?
* try to make commit messages somewhat informative

# Fix travis:
`.travis.yml`
* make sure conda env uses right Python
* fix conda channel order


