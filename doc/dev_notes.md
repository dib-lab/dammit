# Contributing to `dammit`

We welcome external contributions to `dammit`!
All interactions around `dammit` must follow the [dib-lab Code of Conduct](http://ivory.idyll.org/lab/coc.html).

Dammit 2.0 was written by [Camille Scott](http://www.camillescott.org/) and [N Tessa Pierce](http://bluegenes.github.io/).
We are not funded to maintain `dammit`, but will endeavor to do so to the best of our abilities.
In particular, we welcome contributions that address bugs reported in the [issue tracker](https://github.com/dib-lab/dammit/issues).
All additions and bugfixes must be properly covered by tests.

Dammit relies on the snakemake workflow software.
To learn snakemake, check out the comprehensive [documentation](https://snakemake.readthedocs.io/en/stable/), and start with snakemake tutorials such as [this one by Titus Brown](https://github.com/ctb/2019-snakemake-ucdavis).

## Setting up your local computer for `dammit` devevelopment

Make sure conda is installed and channels are properly set up, as in the [installation instructions](install.md)

### Set up the dammit code on your local machine

Fork the [`dammit` repository](https://github.com/dib-lab/dammit) to your GitHub account.
`git clone` this fork to your local computer, then create and name a development branch.
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

Create a conda environment with dependencies installed.

In the main `dammit` folder, run:
```
conda env create -f environment.yml -n dammit_dev
conda activate dammit_dev
```
Now, install an editable version of `dammit` into the `dammit_dev` environment:

```
pip install -e '.'
```

## What do we want/need below here?

### Run Tests

To run tests that do not require databases, run
```
make ci-test
```
To run longer tests, run:
```
make long-test
```

## Code structure

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

* `make long tests` (assumes environment is already setup)
* `make-ci-test`, not long and not huge and not requires_datbases

## Reviewing a PR

**Tests must pass before merging!**

* Have there been radical changes? (Are you adding things to handler, maybe time to take a step back and make sure code uses reasonable variable names, tests, etc)
* Does travis build?
* Try to make commit messages somewhat informative

If these all seem reasonable to you, approve!

## Bioconda

* https://anaconda.org/bioconda/dammit
* Recipe: https://github.com/bioconda/bioconda-recipes/blob/master/recipes/dammit/meta.yaml

## Documentation

* http://dib-lab.github.io/dammit/
* [Tutorial from angus 2018](https://angus.readthedocs.io/en/2018/dammit_annotation.html)
