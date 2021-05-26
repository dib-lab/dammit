# Contributing to `dammit`

We welcome external contributions to `dammit`!
All interactions around `dammit` must follow the [dib-lab Code of Conduct](http://ivory.idyll.org/lab/coc.html).

Dammit 2.0 was written by [Camille Scott](http://www.camillescott.org/) and [N Tessa Pierce](http://bluegenes.github.io/).
We are not funded to maintain `dammit`, but will endeavor to do so to the best of our abilities.
In particular, we welcome contributions that address bugs reported in the [issue tracker](https://github.com/dib-lab/dammit/issues).
All additions and bugfixes must be properly covered by tests.

Dammit relies on the snakemake workflow software.
To learn snakemake, check out the comprehensive [documentation](https://snakemake.readthedocs.io/en/stable/), and maybe start with snakemake tutorials such as [this one by Titus Brown](https://github.com/ctb/2019-snakemake-ucdavis).

## Setting up your local computer for `dammit` development

Make sure conda is installed and channels are properly set up, as in the [installation instructions](install.md).

### Set up the dammit code on your local machine

Fork the `dammit` [repository](https://github.com/dib-lab/dammit) to your GitHub account.
`git clone` this fork to your local computer, then create and name a development branch.
For example, the code below creates a dev branch named `testing`:

```
git clone https://github.com/YOUR-USERNAME/dammit.git
cd dammit
git checkout -b testing
git branch
```
Now you are on the `testing` branch.

Keep the original repository in the `main` branch, so that you can stay up to date with any changes in the main repo.
To do this, first add a remote called `upstream`, which links to the main dammit repository.
```
git remote add upstream https://github.com/dib-lab/dammit.git 
```

Then, make sure the `main` branch is up-to-date by periodically running:
```
git pull upstream main
```

### Create a development environment

Create a conda environment with dependencies installed.

After setting up the code (above), run:
```
conda env create -f environment.yml -n dammit_dev
conda activate dammit_dev
```
Now, install an editable version of `dammit` into the `dammit_dev` environment:

```
pip install -e '.'
```

**What (of the below) do we want to keep?**

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

- Main command line: `dammit/cli.py`
- command line code for each dammit component: `dammit/components`
- Snakemake-related code:
  - main snakefile: `dammit/workflows/dammit.snakefile`
  - databases snakemake rules: `databases.snakefile`
  - annotate snakemake rules: `annotate.snakefile`
  - snakemake wrappers for each included tool:
    - `dammit/wrappers`

## Internal Configuration Files

  - primary config file: `dammit/config.yml`
    - default configuration for all steps run within `dammit run`
  - databases config file: `dammit/databases.yml`
    - download and file naming info for all databases
  - pipeline config file: `dammit/pipelines.yml`
    - sets tools and databases used in each pipeline

## Regenerating test data

`generate-test-data-.sh` re-generates test data and puts it in proper dirs

## Reviewing a PR

**Tests must pass before merging!**

* Have there been radical changes? (Are you adding things to handler, maybe time to take a step back and make sure code uses reasonable variable names, tests, etc)
* Does github-actions build?
* Try to make commit messages somewhat informative

If these all seem reasonable to you, approve!

## Bioconda

* https://anaconda.org/bioconda/dammit
* Recipe: https://github.com/bioconda/bioconda-recipes/blob/master/recipes/dammit/meta.yaml

## Documentation

* http://dib-lab.github.io/dammit/
* [Tutorial from angus 2018](https://angus.readthedocs.io/en/2018/dammit_annotation.html)
