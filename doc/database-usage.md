# Database Usage

`dammit` is composed of two main workflows, `databases`, and `annotate`.

  - `databases` handles downloading and preparing the annotation databases,
  - `annotate` uses these databases for transcriptome annotation


## Check and install databases

By default, dammit downloads databases to your home directory: `$HOME/.dammit/databases`
To check for databases in the default location:

```
dammit run databases
```
This will tell you what databases still need to be installed to run the default annotation pipeline.

To install databases:
```
dammit run databases --install
```
> Notes:
>
> 1. If you're on an HPC or other system with limited space in your home directory, follow
>    the instructions below to specify a custom location.
>
> 2. If you've already downloaded some databases and you'd like to use them with dammit, see the [Advanced Database Handling](database-advanced.md) section.


## Custom database locations

If you'd like to store dammit databases elsewhere, there are two ways to specify a custom location, the `--database-dir` flag or the `DAMMIT_DB_DIR` environment variable.

### Using the `--database-dir` flag:

Check for databases in `/path/to/databases`:
```
dammit databases --database-dir /path/to/databases
```

Install databases in `/path/to/databases`:
```
dammit databases --database-dir /path/to/databases --install
```

### Set up an environment variable

Alternatively, you can set up the `DAMMIT_DB_DIR` environment variable.


Set up the variable in bash. Execute this the command line to use during a single session, or add this to your `$HOME/.bashrc` file to make it persistent.
```
export DAMMIT_DB_DIR=/path/to/databases
```
> Note that the `--database-dir` flag (above) will supersede this variable,
> falling back to the default if neither is set.

When this variable is set up, the standard commands will check for databases in `/path/to/databases` rather than `$HOME/.dammit/databases`

For info on the specific databases used in dammit, see [About Databases](database-about.md).

For advanced installation and usage instructions, check out the
[Advanced Database Handling](database-advanced.md) section.
